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
        self.java_executable_fp = os.environ.get('JAVA', default='/groups/bhurwitz/tools/jdk-11.0.8/bin')
        self.vsearch_executable_fp = os.environ.get('vsearch',
                                                    default='vsearch')
        # self.trim_executable_fp = os.environ.get('TRIMMOMATIC-0.39.JAR',
        # default='/home/matt/Trimmomatic-0.39/trimmomatic-0.39.jar')
        self.trim_executable_fp = os.environ.get(
            'TRIMMOMATIC-0.39.JAR',
            default=
            '/home/u29/mattmiller899/Trimmomatic-0.39/trimmomatic-0.39.jar')

        #print(f"trim = {self.trim_executable_fp}")
        self.frag_executable_fp = os.environ.get('run_FragGeneScan.pl',
                                                 default='/groups/bhurwitz/tools/FragGeneScan1.31/run_FragGeneScan.pl')
        self.interproscan_executable_fp = os.environ.get(
            'INTERPROSCAN.SH', default='/groups/bhurwitz/tools/interproscan-5.46-81.0/interproscan.sh')
        self.pear_executable_fp = os.environ.get('PEAR', default='pear')
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
        #self.paired_ends = int(config["DEFAULT"]["paired_ends"])
        #self.in_dir = config["DEFAULT"]["in_dir"]
        #self.out_dir = config["DEFAULT"]["out_dir"]
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
        self.pear_min_overlap = config["pipeline"]["pear_min_length"]
        self.pear_max_assembly_length = config["pipeline"]["pear_max_assembly"]
        self.pear_min_assembly_length = config["pipeline"]["pear_min_assembly"]
        self.vsearch_filter_maxee = config["pipeline"]["vsearch_filter_maxee"]
        self.vsearch_filter_minlen = config["pipeline"]["vsearch_filter_minlen"]
        #self.vsearch_filter_minlen = config["DEFAULT"]["vsearch_filter_minlen"]
        #self.vsearch_filter_trunclen = config["DEFAULT"]["vsearch_filter_trunclen"]
        #self.vsearch_derep_minuniquesize = config["DEFAULT"]["vsearch_derep_minuniquesize"]
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
        #output_dir_list.append(self.step_01_copy_and_compress(input_dir=self.in_dir))
        #output_dir_list.append(self.step_02_trimming(input_dir=output_dir_list[-1]))
        #_ = self.step_03_centrifuge_taxonomy(input_dir=output_dir_list[-1])
        output_dir_list.append(
            self.step_04_get_gene_reads(input_dir=self.input_dir))
        """
        output_dir_list.append(
            self.step_05_chunk_reads(input_dir=output_dir_list[-1]))
        output_dir_list.append(
            self.step_06_get_orfs(input_dir=output_dir_list[-1]))
        output_dir_list.append(
            self.step_07_combine_tsv(input_dir=output_dir_list[-1]))
        """
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

    """
    def step_01_copy_and_compress(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.warning('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.debug('output_dir: %s', output_dir)
            # Check for fasta and qual files
            fasta_files = []
            types = ['*.fasta*', '*.fa*']
            for t in types:
                fasta_files.extend(glob.glob(os.path.join(input_dir, t)))
            fasta_files = sorted(fasta_files)
            log.info('possible fasta files: "%s"', fasta_files)
            qual_file_glob = os.path.join(input_dir, '*.qual*')
            qual_files = sorted(glob.glob(qual_file_glob))
            log.info('qual files: "%s"', qual_files)
            fastas_no_slashes = [fasta.replace('\\', ' ').replace('/', ' ').split()[-1] for fasta in fasta_files]
            quals_no_slashes = [qual.replace('\\', ' ').replace('/', ' ').split()[-1] for qual in qual_files]
            fastas = [fasta.split(".")[0] for fasta in fastas_no_slashes]
            quals = [qual.split('.')[0] for qual in quals_no_slashes]
            for fasta in fastas:
                for qual in quals:
                    if fasta == qual:
                        tmpfasta = os.path.join(input_dir, fasta + ".fasta")
                        tmpqual = os.path.join(input_dir, qual + ".qual")
                        tmpfastq = os.path.join(input_dir, fasta + ".fastq")
                        log.info("Creating %s", tmpfastq)
                        fasta_qual_to_fastq(tmpfasta, tmpqual, tmpfastq)
                        log.info('%s created', tmpfastq)
                        quals.remove(qual)
                        break
            input_fp_list = []
            types = ['*.fastq*', '*.fq*']
            for t in types:
                input_fp_list.extend(glob.glob(os.path.join(input_dir, t)))
            log.info('input files: %s', input_fp_list)

            if len(input_fp_list) == 0:
                raise PipelineException('found no fastq files in directory "{}"'.format(input_dir))

            for input_fp in input_fp_list:
                destination_name = re.sub(
                                        string=os.path.basename(input_fp),
                                        pattern='\.fq',
                                        repl='.fastq')
                destination_fp = os.path.join(output_dir, destination_name)
                if input_fp.endswith('.gz'):
                    log.info(f"already compressed, copying {input_fp} to {destination_fp}")
                    with open(input_fp, 'rb') as f, open(destination_fp, 'wb') as g:
                        shutil.copyfileobj(fsrc=f, fdst=g)
                else:
                    destination_fp = destination_fp + '.gz'
                    log.info(f"compressing {input_fp} to {destination_fp}")
                    with open(input_fp, 'rt') as f, gzip.open(destination_fp, 'wt') as g:
                        shutil.copyfileobj(fsrc=f, fdst=g)

        self.complete_step(log, output_dir)
        return output_dir
    """

    def step_01_trimming(self, input_file):
        """
        Uses Trimmomatic to trim reads
        :param input_dir: path to input files
        :return: path to output directory
        """
        log, output_dir = self.initialize_step()
        start_time = time.time()
        if len(os.listdir(output_dir)) > 0:
            log.warning(
                'output directory "%s" is not empty, this step will be skipped',
                output_dir)
        else:
            input_fp_list = []
            """
            types = ['*.fastq*', '*.fq*']
            if self.paired_ends:
                input_fp_list = get_forward_fastq_files(input_dir, self.debug)
            else:
                for t in types:
                    input_fp_list.extend(glob.glob(os.path.join(input_dir, t)))
            if len(input_fp_list) == 0:
                raise PipelineException(
                    f'found no fastq files in directory "{input_dir}"')
            log.info(f"input files = {input_fp_list}")
            """
            if not os.path.isfile(input_file):
                raise PipelineException(f'input file {input_file} is not a file or does not exist')
                exit()
            run_arr = [
                self.java_executable_fp, "-jar", self.trim_executable_fp
            ]
            trim_log = f"{output_dir}/trim_log"
            """
            if self.paired_ends:
                run_arr.append("PE")
                out_base = re.sub(string=os.path.basename(input_file),
                                  pattern=r'_1\.(fastq|fq)',
                                  repl=".fastq")
                run_arr.extend([
                    "-threads", self.threads, "-trimlog", trim_log,
                    "-basein", fp, "-baseout",
                    os.path.join(output_dir, out_base)
                ])
            
            else:
            """
            run_arr.append("SE")
            out_base = re.sub(string=os.path.basename(input_file),
                              pattern=r'.(fastq|fq)',
                              repl=".fastq")
            run_arr.extend([
                "-threads",
                str(self.threads), "-trimlog", trim_log, input_file,
                os.path.join(output_dir, out_base)
            ])
            illuminaclip_str = (f"ILLUMINACLIP:{self.trim_adapter_fasta}:"
                                f"{self.trim_seed_mismatches}:"
                                f"{self.trim_palindrome_clip_thresh}:"
                                f"{self.trim_simple_clip_thresh}:"
                                f"{self.trim_min_adapter_length}:"
                                f"{self.trim_keep_both_reads}")
            leading_str = f"LEADING:{self.trim_min_quality}"
            trailing_str = f"TRAILING:{self.trim_min_quality}"
            minlen_str = f"MINLEN:{self.trim_min_len}"
            run_arr.append(illuminaclip_str)
            run_arr.append(leading_str)
            run_arr.append(trailing_str)
            run_arr.append(minlen_str)
            log.info(f"writing output of {input_file} to {output_dir}/{out_base}")
            #run_arr = [self.java_executable_fp, "-jar", self.trim_executable_fp]
            run_cmd(run_arr,
                    log_file=os.path.join(output_dir, 'log'),
                    debug=self.debug)
            # Check the log to make sure most reads were trimmed properly
            with open(os.path.join(output_dir, 'log'), 'r') as logcheck:
                percent_surviving = 0.0
                for l in logcheck:
                    if "Surviving" in l:
                        larr = l.split(" ")
                        #if self.paired_ends:
                        #    percent_surviving = float(larr[7][1:-2])
                        #else:
                        percent_surviving = float(larr[5][1:-2])
                if percent_surviving < 10.0:
                    log.error(
                        f"Fewer than 10% ({percent_surviving}%) of reads "
                        f"survived after Trimmomatic trimming for {input_file}... "
                        f"Exiting")
                    exit(1)
                elif percent_surviving < 50.0:
                    log.warning(
                        f"Fewer than 50% ({percent_surviving}%) of reads "
                        f"survived after Trimmomatic trimming for {input_file}")
                else:
                    log.info(
                        f"{percent_surviving}% of reads survived after "
                        f"Trimmomatic trimming for {input_file}")

        end_time = time.time()
        log.info(f"Time taken for this step: {int((end_time - start_time))}s")
        self.complete_step(log, output_dir)
        return output_dir

    def step_01_1_merge_paired_end_reads(self, input_dir):
        log, output_dir = self.initialize_step()
        start_time = time.time()
        if len(os.listdir(output_dir)) > 0:
            log.warning(
                'output directory "%s" is not empty, this step will be skipped',
                output_dir)
        else:
            log.info('PEAR executable: "%s"', self.pear_executable_fp)

            #for compressed_forward_fastq_fp in get_trimmed_forward_fastq_files(input_dir=input_dir, debug=self.debug):
            #    compressed_reverse_fastq_fp = get_associated_reverse_fastq_fp(forward_fp=compressed_forward_fastq_fp)

            #forward_fastq_fp, reverse_fastq_fp = ungzip_files(
            #    compressed_forward_fastq_fp,
            #    compressed_reverse_fastq_fp,
            #    target_dir=output_dir,
            #    debug=self.debug
            #)
            for forward_fastq_fp in get_trimmed_forward_fastq_files(
                    input_dir=input_dir, debug=self.debug):
                reverse_fastq_fp = get_associated_reverse_fastq_fp(
                    forward_fp=forward_fastq_fp)

                joined_fastq_basename = re.sub(
                    string=os.path.basename(forward_fastq_fp),
                    pattern=r'_1P.fastq',
                    repl=lambda m: '_merged')

                joined_fastq_fp_prefix = os.path.join(output_dir,
                                                      joined_fastq_basename)
                log.info('joining paired ends from "%s" and "%s"',
                         forward_fastq_fp, reverse_fastq_fp)
                log.info('writing joined paired-end reads to "%s"',
                         joined_fastq_fp_prefix)
                run_cmd([
                    self.pear_executable_fp, '-f', forward_fastq_fp, '-r',
                    reverse_fastq_fp, '-o', joined_fastq_fp_prefix,
                    '--min-overlap',
                    str(self.pear_min_overlap), '--max-assembly-length',
                    str(self.pear_max_assembly_length),
                    '--min-assembly-length',
                    str(self.pear_min_assembly_length), '-j',
                    str(self.threads)
                ],
                        log_file=os.path.join(output_dir, 'log'),
                        debug=self.debug)

                #gzip_files(glob.glob(joined_fastq_fp_prefix + '.*.fastq'), debug=self.debug)
                # Check log file to make sure most reads were assembled
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
                            num_assembled = int(
                                l.split(' ')[3].replace(',', ''))
                        elif 'Discarded reads' in l and 'file' not in l:
                            num_discarded = int(
                                l.split(' ')[3].replace(',', ''))
                            log.info("num_assembled = {}, num_discarded = {}".
                                     format(num_assembled, num_discarded))
                            if num_discarded > num_assembled:
                                log.warning(
                                    "More sequences discarded than kept by PEAR for files '{}' and '{}'"
                                    .format(forward_fp, reverse_fp))
                            if num_assembled == 0:
                                log.error(
                                    f"All sequences discarded by PEAR for files '{forward_fp}' and '{reverse_fp}'... Exiting"
                                )
                                exit(1)

        end_time = time.time()
        log.info(f"Time taken for this step: {int((end_time - start_time))}s")
        self.complete_step(log, output_dir)
        return output_dir

    """
    def step_01_2_fastqc(self, input_dir):
       
         Uses FastQC to get the sequence length distribution for vsearch
        :param input_dir: string path to input files
        :return: string path to output directory
        
        log, output_dir = self.initialize_step()
        if self.paired_ends:
            input_fps = glob.glob(f"{input_dir}/*.assembled*.fastq")
        else:
            input_fps = glob.glob(f"{input_dir}/*.fastq")
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
    """

    def step_02_qc_reads_with_vsearch(self, input_dir):
        log, output_dir = self.initialize_step()
        start_time = time.time()
        if len(os.listdir(output_dir)) > 0:
            log.warning(
                'output directory "%s" is not empty, this step will be skipped',
                output_dir)
        else:
            #if self.paired_ends:
            #    input_files_glob = os.path.join(input_dir,
            #                                    '*.assembled*.fastq*')
            #else:
            input_files_glob = os.path.join(input_dir, '*.fastq*')
            input_fp_list = glob.glob(input_files_glob)
            if len(input_fp_list) == 0:
                raise PipelineException(
                    'found no .fastq files in directory "{}"'.format(
                        input_dir))
            #log.info('input file glob: "%s"', input_files_glob)
            log.info(f"input file list: {input_fp_list}")
            for input_fastq_fp in input_fp_list:
                # Uncompress if its still compressed
                if input_fastq_fp[-3:] == ".gz":
                    uncompressed_input_fp = ungzip_files(input_fastq_fp,
                                                         target_dir=input_dir,
                                                         debug=self.debug)[0]
                    os.remove(input_fastq_fp)
                    input_fastq_fp = uncompressed_input_fp
                input_file_basename = os.path.basename(input_fastq_fp)
                output_file_basename = re.sub(
                    string=input_file_basename,
                    pattern='\.fastq*',
                    repl='.ee{}minlen{}.fasta'.format(
                        self.vsearch_filter_maxee, self.vsearch_filter_minlen))
                output_fasta_fp = os.path.join(output_dir,
                                               output_file_basename)
                log.info('vsearch executable: "%s"',
                         self.vsearch_executable_fp)
                log.info('filtering "%s"', input_fastq_fp)
                run_cmd(
                    [
                        self.vsearch_executable_fp,
                        '-fastq_filter',
                        input_fastq_fp,
                        '-fastaout',
                        output_fasta_fp,
                        '-fastq_maxee',
                        str(self.vsearch_filter_maxee),
                        '-fastq_minlen',
                        str(self.vsearch_filter_minlen),
                        #'-fastq_trunclen', str(self.vsearch_filter_trunclen),
                        '-threads',
                        str(self.threads)
                    ],
                    log_file=os.path.join(output_dir, 'log'),
                    debug=self.debug)
            #gzip_glob = glob.glob(os.path.join(output_dir, '*.fastq'))
            #log.info(f"zipping files {gzip_glob}")
            #gzip_files(gzip_glob, debug=self.debug)
            with open(os.path.join(output_dir, 'log'), 'r') as logcheck:
                kept_num = 0
                discarded_num = 0
                for l in logcheck:
                    if 'sequences kept' in l:
                        l_arr = l.split(' ')
                        kept_num = int(l_arr[0])
                        discarded_num = int(l_arr[7])
                    if 'executing' in l:
                        l_arr = l.split(' ')
                        ran_fp = l_arr[3]
                        log.info("kept_num = {}, discarded_num = {}".format(
                            kept_num, discarded_num))
                        if kept_num == 0:
                            log.error(
                                "No sequences kept by vsearch qc for input file '{}'... Exiting"
                                .format(ran_fp))
                            exit(1)
                        if discarded_num > kept_num:
                            log.warning(
                                "More sequences discarded than kept by vsearch qc for input file '{}'"
                                .format(ran_fp))
        end_time = time.time()
        log.info(f"Time taken for this step: {int((end_time - start_time))}s")
        self.complete_step(log, output_dir)
        return output_dir

    """
    def step_03_centrifuge_taxonomy(self, input_dir):
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
                f'*.ee{self.vsearch_filter_maxee}minlen{self.vsearch_filter_minlen}*.fasta'
            )
            #log.info('input file glob: "%s"', input_files_glob)
            input_fp_list = sorted(glob.glob(input_files_glob))
            if len(input_fp_list) == 0:
                raise PipelineException(
                    f'found no .ee{self.vsearch_filter_maxee}minlen{self.vsearch_filter_minlen}.fasta files in directory "{input_dir}"'
                )
            log.info(f"input file list: {input_fp_list}")
            for input_fp in input_fp_list:
                input_basename = os.path.splitext(
                    os.path.basename(input_fp))[0]
                cent_results_fp = f"{output_dir}/{input_basename}_centrifuge_hits.tsv"
                cent_report_fp = f"{output_dir}/{input_basename}_centrifuge_report.tsv"
                fast_type = "f"
                kreport_fp = f"{output_dir}/{input_basename}_kreport.tsv"
                log.info(
                    f"running centrifuge on {input_fp}, outputting results to {cent_results_fp} and report to "
                    f"{cent_report_fp}")
                run_cmd(
                    [
                        self.centrifuge_executable_fp, f"-x",
                        self.centrifuge_db, f"-U", input_fp, f"-S",
                        cent_results_fp, f"--report-file", cent_report_fp,
                        f"-p",
                        str(self.threads), f"-{fast_type}"
                        # INCLUDE MORE ARGS
                    ],
                    log_file=os.path.join(output_dir, 'log'),
                    debug=self.debug)
                log.info(f"creating kraken report from {cent_results_fp}")
                run_cmd([
                    self.centrifuge_kraken_executable_fp, "-x",
                    self.centrifuge_db, f"{cent_results_fp}"
                ],
                        log_file=kreport_fp,
                        debug=False)
        end_time = time.time()
        log.info(f"Time taken for this step: {int((end_time - start_time))}s")
        self.complete_step(log, output_dir)
        return output_dir
    """

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
            #log.info("uncompressing input files")
            #uncompressed_input_fps = ungzip_files(*input_fps, target_dir=input_dir, debug=self.debug)
            #for fp in uncompressed_input_fps:
            for fp in input_fp_list:
                #fasta_fp = re.sub(
                #                    string=fp,
                #                    pattern='\.fastq',
                #                    repl='.fasta')
                #log.info(f"converting fastq {fp} to fasta {fasta_fp}")
                #fastq_to_fasta(fp, fasta_fp)
                #os.remove(fp)
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
                log.info(f"writing output of {fp} to {out_basename}")
                run_cmd(
                    [
                        self.interproscan_executable_fp, "-appl", "Pfam", "-i",
                        fp, "-b", out_basename, "-goterms", "-iprlookup",
                        "-dra", "-cpu",
                        str(self.threads)
                        # INCLUDE MORE ARGS
                    ],
                    log_file=os.path.join(output_dir, 'log'),
                    debug=self.debug)
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
