#!/usr/bin/env python3
"""
Design principles:
  1. reading this script should not require a high level of Python experience
  2. pipeline steps should correspond directly to functions
  3. there should be a single function that represents the pipeline
  4. for development the pipeline should be 'restartable'

  To run:
  python pipeline/pipeline.py -c ./data/configs/config.txt
"""
import argparse
import glob
import gzip
import itertools
import logging
import os
import re
import shutil
import sys
import yaml
import time
from Bio import SeqIO, Seq

#TODO CHANGE BACK TO cluster_16S.
from utils import create_output_dir, gzip_files, ungzip_files, run_cmd, \
     PipelineException, get_sorted_file_list, get_forward_fastq_files, \
     get_associated_reverse_fastq_fp, get_trimmed_forward_fastq_files, \
     fastq_to_fasta


def main():
    logging.basicConfig(level=logging.INFO)
    args = get_args()
    Pipeline(**args.__dict__).run()
    return 0


def get_args():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument(
        '-c',
        '--config',
        required=True,
        help='path to config file containing parameters/arguments')
    arg_parser.add_argument(
        '-i',
        '--input_dir',
        required=True,
        help='path to input dir (step 02 output directory of the QC pipeline)')
    arg_parser.add_argument(
        '-o',
        '--out_dir',
        required=True,
        help='path to output directory')
    arg_parser.add_argument(
        '-t',
        '--threads',
        required=True,
        help='number of threads to use')

    return arg_parser.parse_args()


class Pipeline:
    def __init__(self, config, input_dir, out_dir, threads):
        self.vsearch_executable_fp = shutil.which('vsearch')
        self.frag_executable_fp = shutil.which('FragGeneScan')
        #self.interproscan_executable_fp = glob.glob('/groups/bhurwitz/tools/interproscan-*/interproscan.sh')[0]
        self.interproscan_executable_fp = glob.glob('./tools/interproscan-*/interproscan.sh')[0]
        #self.interproscan_executable_fp = shutil.which('interproscan.sh')
        print(f"FragGeneScan = {self.frag_executable_fp}\nInterProScan = {self.interproscan_executable_fp}")
        self.config_fp = config
        self.read_config()
        self.input_dir = input_dir
        self.out_dir = out_dir
        self.threads = threads


    def read_config(self):
        """
        Takes a path to a config file containing parameters for the pipeline
        :param fp: string path to config file
        :return:
        """ 
        config_in = open(self.config_fp, "r")
        config = yaml.load(config_in)
        self.debug = int(config["pipeline"]["debug"])
        self.trim_adapter_fasta = config["pipeline"]["adapter_fasta"]
        self.trim_seed_mismatches = config["pipeline"]["seed_mismatches"]
        self.trim_palindrome_clip_thresh = config["pipeline"][
            "palindrome_clip_thresh"]
        self.trim_simple_clip_thresh = config["pipeline"]["simple_clip_thresh"]
        self.trim_min_adapter_length = config["pipeline"]["min_adapter_length"]
        self.trim_keep_both_reads = config["pipeline"]["keep_both_reads"]
        self.trim_min_quality = config["pipeline"]["min_quality"]
        self.trim_min_len = config["pipeline"]["trim_min_length"]
        self.vsearch_filter_maxee = config["pipeline"]["vsearch_filter_maxee"]
        self.vsearch_filter_minlen = config["pipeline"]["vsearch_filter_minlen"]
        self.frag_train_file = config["pipeline"]["frag_train_file"]
        self.delete_intermediates = int(config["pipeline"]["delete_intermediates"])
        config_in.close()
        return


    def run(self):
        """
        Runs the entire pipeline
        :param input_dir: string path to input directory
        :return:
        """
        output_dir_list = []
        output_dir_list.append(
            self.step_04_get_gene_reads(input_dir=self.input_dir))
        
        output_dir_list.append(
            self.step_05_chunk_reads(input_dir=output_dir_list[-1]))
        output_dir_list.append(
            self.step_06_get_orfs(input_dir=output_dir_list[-1]))
        output_dir_list.append(
            self.step_07_combine_tsv(input_dir=output_dir_list[-1]))
        
        print(f"delete intermediates = {self.delete_intermediates}")
        if self.delete_intermediates == 1:
            print(f"deleting intermediates")
            for i in range(len(output_dir_list) - 1):
                for fn in os.listdir(output_dir_list[i]):
                    if fn == "log" or fn == "trim_log":
                        continue
                    fp = os.path.join(output_dir_list[i], fn)
                    os.remove(fp)
        return output_dir_list

    def initialize_step(self):
        """
        Sets up the logger and creates the output directory
        :return:
        """
        function_name = sys._getframe(1).f_code.co_name
        log = logging.getLogger(name=function_name)
        if self.debug:
            log.setLevel(logging.DEBUG)
        else:
            log.setLevel(logging.WARNING)
        output_dir = create_output_dir(output_dir_name=function_name,
                                       parent_dir=self.out_dir,
                                       debug=self.debug)
        return log, output_dir

    def complete_step(self, log, output_dir):
        """
        Checks that the output directory contains files
        :param log:
        :param output_dir:
        :return:
        """
        output_dir_list = sorted(os.listdir(output_dir))
        if len(output_dir_list) == 0:
            raise PipelineException(
                'ERROR: no output files in directory "{}"'.format(output_dir))
        return


    def step_04_get_gene_reads(self, input_dir):
        """
        Uses FragGeneScan to find reads containing fragments of genes
        :param input_dir: string path to input files
        :return: string path to output directory
        """
        log, output_dir = self.initialize_step()
        start_time = time.time()
        if len(os.listdir(output_dir)) > 0:
            log.warning(
                'output directory "%s" is not empty, this step will be skipped',
                output_dir)
        else:
            #input_fps = glob.glob(f"{input_dir}/*.fastq.gz")
            #TODO CHANGE BACK
            input_fp_list = glob.glob(f"{input_dir}/*.fasta")
            #input_fps = glob.glob(f"{input_dir}/*.fasta")
            if len(input_fp_list) == 0:
                raise PipelineException(
                    f'found no fasta files in directory "{input_dir}"')
            log.info(f"input files = {input_fp_list}")
            for fp in input_fp_list:
                out_fp = os.path.join(
                    output_dir,
                    re.sub(string=os.path.basename(fp),
                           pattern='\.fasta',
                           repl='_frags'))
                log.info(f"writing output of {fp} to {out_fp}")
                run_cmd(
                    [
                        self.frag_executable_fp, f"-genome={fp}",
                        f"-out={out_fp}", "-complete=0",
                        f"-train={self.frag_train_file}",
                        f"-thread={str(self.threads)}"
                        # INCLUDE MORE ARGS
                    ],
                    debug=self.debug)
        end_time = time.time()
        log.info(f"Time taken for this step: {int((end_time - start_time))}s")
        self.complete_step(log, output_dir)
        return output_dir

    def step_05_chunk_reads(self, input_dir):
        log, output_dir = self.initialize_step()
        start_time = time.time()
        if len(os.listdir(output_dir)) > 0:
            log.info(
                'output directory "%s" is not empty, this step will be skipped',
                output_dir)
        else:
            log.info('input directory listing:\n\t%s',
                     '\n\t'.join(os.listdir(input_dir)))
            input_files_glob = os.path.join(
                input_dir,
                f'*_trimmed_qcd_frags.faa'
            )
            #log.info('input file glob: "%s"', input_files_glob)
            input_fp_list = sorted(glob.glob(input_files_glob))
            if len(input_fp_list) == 0:
                raise PipelineException(
                    f'found no _trimmed_qcd_frags.faa files in directory "{input_dir}"'
                )
            log.info(f"input file list: {input_fp_list}")
            chunk_size = 40000
            for input_fp in input_fp_list:
                i = 0
                log.info(f"reading input file {input_fp}")
                fname, ext = input_fp.rsplit('.', 1)
                _, fname = os.path.split(fname)
                written = False
                with open(input_fp, "r") as infile:
                    outfilepath = f"{output_dir}/{fname}_{i}.{ext}"
                    log.info(f"writing chunk to {outfilepath}")
                    tmp_count = 0
                    outfile = open(outfilepath, 'w')
                    for record in SeqIO.parse(infile, "fasta"):
                        if tmp_count == chunk_size:
                            outfile.close()
                            i += 1
                            outfilepath = f"{output_dir}/{fname}_{i}.{ext}"
                            log.info(f"writing chunk to {outfilepath}")
                            outfile = open(outfilepath, 'w')
                            tmp_count = 0
                        tmp_seq = str(record.seq)
                        tmp_id = record.id
                        tmp_seq = re.sub("\*", "X", tmp_seq)
                        outfile.write(f">{tmp_id}\n{tmp_seq}\n")
                        tmp_count += 1
                try:
                    outfile.close()
                except:
                    pass
            #gzip_files(glob.glob(os.path.join(output_dir, '*.fasta')))
        end_time = time.time()
        log.info(f"Time taken for this step: {int((end_time - start_time))}s")
        self.complete_step(log, output_dir)
        return output_dir

    def step_06_get_orfs(self, input_dir):
        """
        Uses Interproscan to get ORFS and connect them to GO terms
        :param input_dir: string path to input files
        :return: string path to output directory
        """
        log, output_dir = self.initialize_step()
        start_time = time.time()
        num_chunks = len(os.listdir(input_dir))
        if len(os.listdir(output_dir)) == num_chunks * 4:
            log.warning(
                'output directory "%s" has finished every chunk, this step will be skipped',
                output_dir)
        else:
            input_fp_list = sorted(glob.glob(f"{input_dir}/*.faa"))
            if len(input_fp_list) == 0:
                raise PipelineException(
                    f'found no faa files in directory "{input_dir}"')
            log.info(f"input files = {input_fp_list}")
            for fp in input_fp_list:
                out_basename = os.path.join(
                    output_dir,
                    re.sub(string=os.path.basename(fp),
                           pattern='\.faa',
                           repl='_interpro'))
                if os.path.isfile(f"{out_basename}.tsv") and os.path.getsize(f"{out_basename}.tsv") > 0:
                    log.info(f"{out_basename} already ran to completion, skipping")
                    continue
                log.info(f"writing output of {fp} to {out_basename}")
                run_cmd(
                    [
                        self.interproscan_executable_fp, "-appl", "Pfam", "-i",
                        fp, "-b", out_basename, "-goterms", "-iprlookup",
                        "-dra", "-cpu",
                        str(self.threads)
                        # INCLUDE MORE ARGS
                    ],
                    debug=self.debug)
                if not os.path.isfile(f"{out_basename}.tsv"):
                    log.error("ERROR: interproscan failed for file '{out_basename}.tsv', exiting...")
                    exit(1)
        end_time = time.time()
        log.info(f"Time taken for this step: {int((end_time - start_time))}s")
        self.complete_step(log, output_dir)
        return output_dir

    def step_07_combine_tsv(self, input_dir):
        log, output_dir = self.initialize_step()
        start_time = time.time()
        if len(os.listdir(output_dir)) > 0:
            log.warning(
                'output directory "%s" is not empty, this step will be skipped',
                output_dir)
        else:
            input_fp_list = glob.glob(f"{input_dir}/*.tsv")
            if len(input_fp_list) == 0:
                raise PipelineException(
                    f'found no tsv files in directory "{input_dir}"')
                exit()
            log.info(f"input files = {input_fp_list}")
            out_basename = "_".join(os.path.basename(input_fp_list[0]).split("_")[0:-2])
            out_fp = os.path.join(output_dir, f"{out_basename}_interpro_combined.tsv")
            log.info(f"out_fp = {out_fp}")
            with open(out_fp, "w") as out:
                for fp in input_fp_list:
                    with open(fp, "r") as f:
                        for l in f:
                            out.write(l)
        end_time = time.time()
        log.info(f"Time taken for this step: {int((end_time - start_time))}s")
        self.complete_step(log, output_dir)
        return output_dir

if __name__ == '__main__':
    main()
