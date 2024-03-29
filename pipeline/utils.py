import glob
import gzip
import logging
from operator import attrgetter
import os
import re
import shutil
import subprocess
import traceback


class PipelineException(BaseException):
    pass


def get_sorted_file_list(dir_path):
    return tuple(
        sorted(
            [
                entry
                for entry
                in os.scandir(dir_path)
                if entry.is_file()
            ],
        key=attrgetter('name')
        )
    )


def create_output_dir(output_dir_name, parent_dir=None, input_dir=None, debug=False):
    if parent_dir is not None and input_dir is None:
        pass
    elif input_dir is not None and parent_dir is None:
        parent_dir, _ = os.path.split(input_dir)
    else:
        raise ValueError('exactly one of parent_dir and input_dir must be None')

    log = logging.getLogger(name=__name__)
    if debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.WARNING)
    output_dir = os.path.join(parent_dir, output_dir_name)
    if os.path.exists(output_dir):
        log.warning('directory "%s" already exists', output_dir)
    else:
        os.makedirs(output_dir)
    return output_dir


def gzip_files(file_list, debug=False):
    log = logging.getLogger(name=__name__)
    if debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.WARNING)
    log.info("in gzip_files")
    for fp in file_list:
        with open(fp, 'rt') as src, gzip.open(fp + '.gz', 'wt') as dst:
            log.info('compressing file "%s"', fp)
            shutil.copyfileobj(fsrc=src, fdst=dst)
        os.remove(fp)


def ungzip_files(*fp_list, target_dir, debug=False):
    log = logging.getLogger(name=__name__)
    if debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.WARNING)
    uncompressed_fps = []
    for fp in fp_list:
        # remove '.gz' from the file path
        uncompressed_fp = os.path.join(target_dir, os.path.basename(fp)[:-3])
        log.info('unzipping {} to {}'.format(fp, uncompressed_fp))
        uncompressed_fps.append(uncompressed_fp)
        with gzip.open(fp, 'rt') as src, open(uncompressed_fp, 'wt') as dst:
            shutil.copyfileobj(fsrc=src, fdst=dst)
    return uncompressed_fps


def run_cmd(cmd_line_list, logfile="", debug=True, **kwargs):
    log = logging.getLogger(name=__name__)
    if debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.WARNING)
    counter = 0
    while counter < 100:
        try:
            counter += 1
            if logfile != "":
                with open(logfile, "at") as out:
                    cmd_line_str = ' '.join((str(x) for x in cmd_line_list))
                    out.write(f'executing "{cmd_line_str}"')
                    output = subprocess.run(
                        cmd_line_list,
                        stdout=out,
                        stderr=subprocess.STDOUT,
                        universal_newlines=True,
                        **kwargs)
                    log.info(output)
                    return output
            else:
                cmd_line_str = ' '.join((str(x) for x in cmd_line_list))
                log.info(f'executing "{cmd_line_str}"')
                output = subprocess.run(
                    cmd_line_list,
                    #stdout=log_file,
                    #stderr=subprocess.STDOUT,
                    universal_newlines=True,
                    **kwargs)
                log.info(output)
                return output

        except subprocess.CalledProcessError as c:
            logging.exception(c)
            print(c.message)
            print(c.cmd)
            print(c.output)
            #raise c
        except Exception as e:
            logging.exception(e)
            print('blarg!')
            print(e)
            traceback.print_exc()
            #raise e


def get_forward_fastq_files(input_dir, debug=False):
    log = logging.getLogger(name=__name__)
    if debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.WARNING)
    input_glob = os.path.join(input_dir, '*_1.f*q*')
    log.info('searching for forward read files with glob "%s"', input_glob)
    forward_fastq_files = glob.glob(input_glob)
    if len(forward_fastq_files) == 0:
        raise PipelineException('found no forward reads from glob "{}"'.format(input_glob))
    return forward_fastq_files


def get_trimmed_forward_fastq_files(input_dir, debug=False):
    log = logging.getLogger(name=__name__)
    if debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.WARNING)
    input_glob = os.path.join(input_dir, '*_1P.fastq*')
    log.info('searching for forward read files with glob "%s"', input_glob)
    forward_fastq_files = glob.glob(input_glob)
    if len(forward_fastq_files) == 0:
        raise PipelineException('found no forward reads from glob "{}"'.format(input_glob))
    return forward_fastq_files


def get_associated_reverse_fastq_fp(forward_fp):
    forward_input_dir, forward_basename = os.path.split(forward_fp)
    reverse_fastq_basename = re.sub(
        string=forward_basename,
        pattern=r'_1P.fastq',
        repl=lambda m: '_2P.fastq')
    reverse_fastq_fp = os.path.join(forward_input_dir, reverse_fastq_basename)
    return reverse_fastq_fp


def fastq_to_fasta(fastq_fp, fasta_fp):
    line_count = 0
    if fastq_fp[-3:] == ".gz":
        in_file = gzip.open(fastq_fp, "rb")
        out_file = gzip.open(fasta_fp, "w")
    else:
        in_file = open(fastq_fp, "r")
        out_file = open(fasta_fp, "w")
    for line_count, l in enumerate(in_file):
        if line_count % 4 == 0:
            out_file.write(f">{l[1:]}")
        elif line_count % 4 == 1:
            out_file.write(l)
    in_file.close()
    out_file.close()
    
    
