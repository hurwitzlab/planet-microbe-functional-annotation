"""
Design principles:
  1. reading this script should not require a high level of Python experience
  2. pipeline steps should correspond directly to functions
  3. there should be a single function that represents the pipeline
  4. for development the pipeline should be 'restartable'
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
import configparser

#TODO CHANGE BACK TO cluster_16S.
from utils import create_output_dir, gzip_files, ungzip_files, run_cmd, PipelineException, get_sorted_file_list, \
    get_forward_fastq_files, get_associated_reverse_fastq_fp


def main():
    logging.basicConfig(level=logging.INFO)
    args = get_args()
    Pipeline(**args.__dict__).run()
    return 0

def get_args():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-c', '--config', required=True,
                            help='path to config file containing parameters/arguments')
    return arg_parser.parse_args()


class Pipeline:
    def __init__(self,
                 config
                 ):
        self.java_executable_fp = os.environ.get('JAVA', default='java')
        self.fastqc_executable_fp = os.environ.get('FASTQC', default='fastqc')
        self.trim_executable_fp = os.environ.get('TRIMMOMATIC-0.39.JAR', default='/home/matt/Trimmomatic-0.39/trimmomatic-0.39.jar')
        print(f"trim = {self.trim_executable_fp}")
        self.frag_executable_fp = os.environ.get('FRAGGENESCAN', default='fraggenescan')
        self.interproscan_executable_fp = os.environ.get('INTERPROSCAN', default='interproscan')
        self.config_fp = config
        self.read_config()


    def read_config(self):
        """
        Takes a path to a config file containing parameters for the pipeline
        :param fp: string path to config file
        :return:
        """
        config = configparser.ConfigParser()
        config.read(self.config_fp)
        self.paired_ends = config["DEFAULT"]["paired_ends"]
        self.debug = config["DEFAULT"]["debug"]
        self.in_dir = config["DEFAULT"]["in_dir"]
        self.out_dir = config["DEFAULT"]["out_dir"]
        self.threads = config["DEFAULT"]["threads"]
        self.trim_adapter_fasta = config["DEFAULT"]["adapter_fasta"]
        self.trim_seed_mismatches = config["DEFAULT"]["seed_mismatches"]
        self.trim_palindrom_clip_thresh = config["DEFAULT"]["palindrome_clip_thresh"]
        self.trim_simple_clip_thresh = config["DEFAULT"]["simple_clip_thresh"]
        self.trim_min_adapter_length = config["DEFAULT"]["min_adapter_length"]
        self.trim_keep_both_reads = config["DEFAULT"]["keep_both_reads"]
        self.trim_min_quality = config["DEFAULT"]["min_quality"]
        self.trim_min_len = config["DEFAULT"]["min_length"]
        return


    def run(self):
        """
        Runs the entire pipeline
        :param input_dir: string path to input directory
        :return:
        """
        output_dir_list = []
        output_dir_list.append(self.step_01_trimming(input_dir=self.in_dir))
        if self.paired_ends:
            output_dir_list.append(self.step_01_1_merge_paired_end_reads(input_dir=self.in_dir))
        output_dir_list.append(self.step_02_fastqc(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_03_get_gene_reads(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_04_get_orfs(input_dir=output_dir_list[-1]))
        return output_dir_list


    def initialize_step(self):
        """
        Sets up the logger and creates the output directory
        :return:
        """
        function_name = sys._getframe(1).f_code.co_name
        log = logging.getLogger(name=function_name)
        if self.debug is True:
            log.setLevel(logging.DEBUG)
        else:
            log.setLevel(logging.WARNING)
        output_dir = create_output_dir(output_dir_name=function_name, parent_dir=self.out_dir, debug=self.debug)
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
            raise PipelineException('ERROR: no output files in directory "{}"'.format(output_dir))
        return


    def step_01_trimming(self, input_dir):
        """
        Uses Trimmomatic to trim reads
        :param input_dir: path to input files
        :return: path to output directory
        """
        log, output_dir = self.initialize_step()
        input_fps = []
        types = ['*.fastq*', '*.fq*']
        if self.paired_ends:
            input_fps = get_forward_fastq_files(input_dir, self.debug)
        else:
            for t in types:
                input_fps.extend(glob.glob(os.path.join(input_dir, t)))
        log.info(f"input files = {input_fps}")
        if len(input_fps) == 0:
            raise PipelineException(f'found no fastq files in directory "{input_dir}"')
        for fp in input_fps:
            out_fp = os.path.join(output_dir, re.sub(
                                                    string=os.path.basename(fp),
                                                    pattern='\.fq',
                                                    repl='.fastq'))
            run_arr = [self.java_executable_fp, "-jar", self.trim_executable_fp]
            if self.paired_ends:
                run_arr.append("PE")
            else:
                run_arr.append("SE")
            out_base = re.sub(
                string=os.path.basename(fp),
                pattern=r'_([0R])1',
                repl="")
            trim_log = f"{output_dir}/trim_log"
            cmd_log = f"{output_dir}/cmd_log"
            run_arr.extend(["-threads", self.threads, "-trimlog", trim_log, "-basein", fp, "-baseout",
                            os.path.join(output_dir, out_base)])
            illuminaclip_str = f"ILLUMINACLIP:{self.trim_adapter_fasta}:{self.trim_seed_mismatches}:" \
                              f"{self.trim_palindrom_clip_thresh}:{self.trim_simple_clip_thresh}:" \
                              f"{self.trim_min_adapter_length}:{self.trim_keep_both_reads}"
            leading_str = f"LEADING:{self.trim_min_quality}"
            trailing_str = f"TRAILING:{self.trim_min_quality}"
            minlen_str = f"MINLEN:{self.trim_min_len}"
            run_arr.append(illuminaclip_str)
            run_arr.append(leading_str)
            run_arr.append(trailing_str)
            run_arr.append(minlen_str)
            log.info(f"writing output of {fp} to {output_dir}/{out_base}")
            #run_arr = [self.java_executable_fp, "-jar", self.trim_executable_fp]
            run_cmd(
                run_arr,
                log_file=cmd_log,
                debug=self.debug
            )
        self.complete_step(log, output_dir)
        return output_dir


    def step_01_1_merge_paired_end_reads(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.warning('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.info('PEAR executable: "%s"', self.pear_executable_fp)

            for compressed_forward_fastq_fp in get_forward_fastq_files(input_dir=input_dir, debug=self.debug):
                compressed_reverse_fastq_fp = get_associated_reverse_fastq_fp(forward_fp=compressed_forward_fastq_fp)

                forward_fastq_fp, reverse_fastq_fp = ungzip_files(
                    compressed_forward_fastq_fp,
                    compressed_reverse_fastq_fp,
                    target_dir=output_dir
                )

                joined_fastq_basename = re.sub(
                    string=os.path.basename(forward_fastq_fp),
                    pattern=r'_([0R]1)',
                    repl=lambda m: '_merged'.format(m.group(1)))[:-6]

                joined_fastq_fp_prefix = os.path.join(output_dir, joined_fastq_basename)
                log.info('joining paired ends from "%s" and "%s"', forward_fastq_fp, reverse_fastq_fp)
                log.info('writing joined paired-end reads to "%s"', joined_fastq_fp_prefix)
                run_cmd([
                        self.pear_executable_fp,
                        '-f', forward_fastq_fp,
                        '-r', reverse_fastq_fp,
                        '-o', joined_fastq_fp_prefix,
                        '--min-overlap', str(self.pear_min_overlap),
                        '--max-assembly-length', str(self.pear_max_assembly_length),
                        '--min-assembly-length', str(self.pear_min_assembly_length),
                        '-j', str(self.core_count)
                    ],
                    log_file = os.path.join(output_dir, 'log'),
                    debug=self.debug
                )

                # delete the uncompressed input files
                os.remove(forward_fastq_fp)
                os.remove(reverse_fastq_fp)
                gzip_files(glob.glob(joined_fastq_fp_prefix + '.*.fastq'), debug=self.debug)
                with open(os.path.join(output_dir, 'log'), 'r') as logcheck:
                    num_assembled = 0
                    num_discarded = 0
                    forward_fp = ""
                    reverse_fp = ""
                    for l in logcheck:
                        if 'Forward reads file' in l:
                            forward_fp = l.split(' ')[-1]
                        elif 'Reverse reads file' in l:
                            reverse_fp = l.split(' ')[-1]
                        elif 'Assembled reads' in l and 'file' not in l:
                            num_assembled = int(l.split(' ')[3].replace(',', ''))
                        elif 'Discarded reads' in l and 'file' not in l:
                            num_discarded = int(l.split(' ')[3].replace(',', ''))
                            log.info("num_assembled = {}, num_discarded = {}".format(num_assembled, num_discarded))
                            if num_discarded > num_assembled:
                                log.warning("More sequences discarded than kept by PEAR for files '{}' and '{}'".format(forward_fp, reverse_fp))

        self.complete_step(log, output_dir)
        return output_dir


    def step_02_fastqc(self, input_dir):
        """
        Uses FastQC to filter out bad reads
        :param input_dir: string path to input files
        :return: string path to output directory
        """
        log, output_dir = self.initialize_step()
        if self.paired_ends:
            input_fps = glob.glob(f"{input_dir}/*.assembled*.fastq.gz")
        else:
            input_fps = glob.glob(f"{input_dir}/*.fastq.gz")
        log.info(f"input files = {input_fps}")
        if len(input_fps) == 0:
            raise PipelineException(f'found no fastq files in directory "{input_dir}"')
        for fp in input_fps:
            log.info(f"qc on {fp}")
            run_cmd([
                self.fastqc_executable_fp,
                fp,
                f"--outdir={output_dir}"
                # INCLUDE MORE ARGS
            ],
                log_file=os.path.join(output_dir, 'log'),
                debug=self.debug
            )
        self.complete_step(log, output_dir)
        return output_dir


    def step_03_get_gene_reads(self, input_dir):
        """
        Uses FragGeneScan to find reads containing fragments of genes
        :param input_dir: string path to input files
        :return: string path to output directory
        """
        log, output_dir = self.initialize_step()
        input_fps = glob.glob(f"{input_dir}/*.fastq")
        log.info(f"input files = {input_fps}")
        if len(input_fps) == 0:
            raise PipelineException(f'found no fastq files in directory "{input_dir}"')
        for fp in input_fps:
            out_fp = os.path.join(output_dir, re.sub(
                                                    string=os.path.basename(fp),
                                                    pattern='\.fastq\.gz',
                                                    repl='frags.fastq.gz'))
            log.info(f"writing output of {fp} to {out_fp}")
            run_cmd([
                self.frag_executable_fp,
                f"-genome={fp}",
                f"-out={out_fp}",
                "-complete=0",
                f"thread={self.threads}"
                # INCLUDE MORE ARGS
            ],
                log_file=os.path.join(output_dir, 'log'),
                debug=self.debug
            )
        self.complete_step(log, output_dir)
        return output_dir


    def step_04_get_orfs(self, input_dir):
        """
        Uses Interproscan to get ORFS and connect them to GO terms
        :param input_dir: string path to input files
        :return: string path to output directory
        """
        log, output_dir = self.initialize_step()
        input_fps = glob.glob(f"{input_dir}/*.fastq")
        log.info(f"input files = {input_fps}")
        if len(input_fps) == 0:
            raise PipelineException(f'found no fastq files in directory "{input_dir}"')
        for fp in input_fps:
            out_fp = os.path.join(output_dir, re.sub(
                                                    string=os.path.basename(fp),
                                                    pattern='\.fq',
                                                    repl='.fastq'))
            log.info(f"writing output of {fp} to {out_fp}")
            run_cmd([
                self.interproscan_executable_fp_executable_fp
                # INCLUDE MORE ARGS
            ],
                log_file=os.path.join(output_dir, 'log'),
                debug=self.debug
            )
        self.complete_step(log, output_dir)
        return output_dir


if __name__ == '__main__':
    main()

