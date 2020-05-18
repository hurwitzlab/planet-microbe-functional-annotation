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
from utils import create_output_dir, gzip_files, ungzip_files, run_cmd, PipelineException, get_sorted_file_list


def main():
    logging.basicConfig(level=logging.INFO)
    args = get_args()
    Pipeline(**args.__dict__).run(input_dir=args.input_dir)
    return 0

def get_args():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-c', '--config', required=True,
                            help='path to config file containing parameters/arguments')
    return arg_parser.parse_args()


class Pipeline:
    def __init__(self,
                 config_fp
                 ):
        self.fastqc_executable_fp = os.environ.get('FASTQC', default='fastqc')
        self.trim_executable_fp = os.environ.get('TRIMMOMATIC', default='trimmomatic')
        self.frag_executable_fp = os.environ.get('FRAGGENESCAN', default='fraggenescan')
        self.interproscan_executable_fp = os.environ.get('INTERPROSCAN', default='interproscan')
        self.config_fp = config_fp
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
        return


    def run(self, input_dir):
        """
        Runs the entire pipeline
        :param input_dir: string path to input directory
        :return:
        """
        output_dir_list = []
        output_dir_list.append(self.step_01_trimming(input_dir=input_dir))
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
            log.info(f"writing output of {fp} to {out_fp}")
            run_cmd([
                    self.trim_executable_fp
                    # INCLUDE MORE ARGS
                ],
                log_file=os.path.join(output_dir, 'log'),
                debug=self.debug
            )
        self.complete_step(log, output_dir)
        return output_dir


    def step_02_fastqc(self, input_dir):
        """
        Uses FastQC to filter out bad reads
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
                self.fastqc_executable_fp
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
                                                    pattern='\.fq',
                                                    repl='.fastq'))
            log.info(f"writing output of {fp} to {out_fp}")
            run_cmd([
                self.frag_executable_fp
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

