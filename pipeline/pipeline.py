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
import configparser

#TODO CHANGE BACK TO cluster_16S.
from utils import create_output_dir, gzip_files, ungzip_files, run_cmd, PipelineException, get_sorted_file_list, \
    get_forward_fastq_files, get_associated_reverse_fastq_fp, get_trimmed_forward_fastq_files, fastq_to_fasta


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
        self.vsearch_executable_fp = os.environ.get('vsearch', default='vsearch')
        self.trim_executable_fp = os.environ.get('TRIMMOMATIC-0.39.JAR', default='/home/matt/Trimmomatic-0.39/trimmomatic-0.39.jar')
        #print(f"trim = {self.trim_executable_fp}")
        self.frag_executable_fp = os.environ.get('run_FragGeneScan.pl', default='run_FragGeneScan.pl')
        self.interproscan_executable_fp = os.environ.get('INTERPROSCAN.SH', default='interproscan.sh')
        self.pear_executable_fp = os.environ.get('PEAR', default='pear')
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
        self.debug = int(config["DEFAULT"]["debug"])
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
        self.trim_min_len = config["DEFAULT"]["trim_min_length"]
        self.pear_min_overlap = config["DEFAULT"]["pear_min_length"]
        self.pear_max_assembly_length = config["DEFAULT"]["pear_max_assembly"]
        self.pear_min_assembly_length = config["DEFAULT"]["pear_min_assembly"]
        self.vsearch_filter_maxee = config["DEFAULT"]["vsearch_filter_maxee"]
        self.vsearch_filter_minlen = config["DEFAULT"]["vsearch_filter_minlen"]
        self.frag_train_file = config["DEFAULT"]["frag_train_file"]
        return


    def run(self):
        """
        Runs the entire pipeline
        :param input_dir: string path to input directory
        :return:
        """
        output_dir_list = []
        output_dir_list.append(self.step_01_copy_and_compress(input_dir=self.in_dir))
        output_dir_list.append(self.step_02_trimming(input_dir=output_dir_list[-1]))
        if self.paired_ends:
            output_dir_list.append(self.step_02_1_merge_paired_end_reads(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_03_qc_reads_with_vsearch(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_04_get_gene_reads(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_05_get_orfs(input_dir=output_dir_list[-1]))
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


    def step_02_trimming(self, input_dir):
        """
        Uses Trimmomatic to trim reads
        :param input_dir: path to input files
        :return: path to output directory
        """
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.warning('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
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
                run_arr = [self.java_executable_fp, "-jar", self.trim_executable_fp]
                if self.paired_ends:
                    run_arr.append("PE")
                else:
                    run_arr.append("SE")
                out_base = re.sub(
                    string=os.path.basename(fp),
                    pattern=r'_1\.(fastq|fq)',
                    repl=".fastq")
                trim_log = f"{output_dir}/trim_log"
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
                    log_file=os.path.join(output_dir, 'log'),
                    debug=self.debug
                )
        self.complete_step(log, output_dir)
        return output_dir


    def step_02_1_merge_paired_end_reads(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.warning('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.info('PEAR executable: "%s"', self.pear_executable_fp)

            for compressed_forward_fastq_fp in get_trimmed_forward_fastq_files(input_dir=input_dir, debug=self.debug):
                compressed_reverse_fastq_fp = get_associated_reverse_fastq_fp(forward_fp=compressed_forward_fastq_fp)

                forward_fastq_fp, reverse_fastq_fp = ungzip_files(
                    compressed_forward_fastq_fp,
                    compressed_reverse_fastq_fp,
                    target_dir=output_dir,
                    debug=self.debug
                )

                joined_fastq_basename = re.sub(
                    string=os.path.basename(forward_fastq_fp),
                    pattern=r'_1P.fastq',
                    repl=lambda m: '_merged')

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
                        '-j', str(self.threads)
                    ],
                    log_file = os.path.join(output_dir, 'log'),
                    debug=self.debug
                )

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


    def step_03_qc_reads_with_vsearch(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.warning('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            if self.paired_ends:
                input_files_glob = os.path.join(input_dir, '*.assembled*.fastq.gz')
            else:
                input_files_glob = os.path.join(input_dir, '*.fastq.gz')
            input_file_list = glob.glob(input_files_glob)
            if len(input_file_list) == 0:
                raise PipelineException('found no .fastq.gz files in directory "{}"'.format(input_dir))
            log.info('input file glob: "%s"', input_files_glob)
            for assembled_fastq_fp in input_file_list:
                input_file_basename = os.path.basename(assembled_fastq_fp)
                output_file_basename = re.sub(
                    string=input_file_basename,
                    pattern='\.fastq\.gz',
                    repl='.ee{}minlen{}.fastq.gz'.format(self.vsearch_filter_maxee, self.vsearch_filter_minlen)[:-3]
                )
                output_fastq_fp = os.path.join(output_dir, output_file_basename)
                uncompressed_input_fp = ungzip_files(assembled_fastq_fp, target_dir=input_dir, debug=self.debug)[0]
                
                log.info('vsearch executable: "%s"', self.vsearch_executable_fp)
                log.info('filtering "%s"', uncompressed_input_fp)
                run_cmd([
                        self.vsearch_executable_fp,
                        '-fastq_filter', uncompressed_input_fp,
                        '-fastqout', output_fastq_fp,
                        '-fastq_maxee', str(self.vsearch_filter_maxee),
                        '-fastq_minlen', str(self.vsearch_filter_minlen),
                        '-threads', str(self.threads)
                    ],
                    log_file = os.path.join(output_dir, 'log'),
                    debug=self.debug
                )
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
                        log.info("kept_num = {}, discarded_num = {}".format(kept_num, discarded_num))
                        if kept_num == 0:
                            log.error("No sequences kept by vsearch qc for input file '{}'".format(ran_fp))
                        if discarded_num > kept_num:
                            log.warning("More sequences discarded than kept by vsearch qc for input file '{}'".format(ran_fp))
        self.complete_step(log, output_dir)
        return output_dir


    """
    def step_02_fastqc(self, input_dir):
        
        Uses FastQC to filter out bad reads
        :param input_dir: string path to input files
        :return: string path to output directory
        
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
    """

    def step_04_get_gene_reads(self, input_dir):
        """
        Uses FragGeneScan to find reads containing fragments of genes
        :param input_dir: string path to input files
        :return: string path to output directory
        """
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.warning('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            #input_fps = glob.glob(f"{input_dir}/*.fastq.gz")
            input_fps = glob.glob(f"{input_dir}/*.fastq")
            log.info(f"input files = {input_fps}")
            if len(input_fps) == 0:
                raise PipelineException(f'found no fastq files in directory "{input_dir}"')
            #log.info("uncompressing input files")
            #uncompressed_input_fps = ungzip_files(*input_fps, target_dir=input_dir, debug=self.debug)
            #for fp in uncompressed_input_fps:
            for fp in input_fps:
                fasta_fp = re.sub(
                                    string=fp,
                                    pattern='\.fastq',
                                    repl='.fasta')
                log.info(f"converting fastq {fp} to fasta {fasta_fp}")
                fastq_to_fasta(fp, fasta_fp)
                out_fp = os.path.join(output_dir, re.sub(
                                                        string=os.path.basename(fp),
                                                        pattern='\.fastq',
                                                        repl='_frags'))
                log.info(f"writing output of {fasta_fp} to {out_fp}")
                run_cmd([
                    self.frag_executable_fp,
                    f"-genome={fasta_fp}",
                    f"-out={out_fp}",
                    "-complete=0",
                    f"-train={self.frag_train_file}",
                    f"thread={self.threads}"
                    # INCLUDE MORE ARGS
                ],
                    log_file=os.path.join(output_dir, 'log'),
                    debug=self.debug
                )
        self.complete_step(log, output_dir)
        return output_dir


    def step_05_get_orfs(self, input_dir):
        """
        Uses Interproscan to get ORFS and connect them to GO terms
        :param input_dir: string path to input files
        :return: string path to output directory
        """
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.warning('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            input_fps = glob.glob(f"{input_dir}/*.ffn")
            log.info(f"input files = {input_fps}")
            if len(input_fps) == 0:
                raise PipelineException(f'found no fnn files in directory "{input_dir}"')
            for fp in input_fps:
                out_basename = os.path.join(output_dir, re.sub(
                                                        string=os.path.basename(fp),
                                                        pattern='\.fnn',
                                                        repl='_interpro'))
                log.info(f"writing output of {fp} to {out_basename}")
                run_cmd([
                        self.interproscan_executable_fp,
                        "-appl Pfam",
                        f"-i {fp}",
                        f"-b {out_basename}",
                        "-goterms"
                        # INCLUDE MORE ARGS
                    ],
                    log_file=os.path.join(output_dir, 'log'),
                    debug=self.debug
                )
        self.complete_step(log, output_dir)
        return output_dir


if __name__ == '__main__':
    main()

