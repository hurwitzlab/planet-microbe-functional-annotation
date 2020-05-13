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
    if debug is True:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.WARNING)
    output_dir = os.path.join(parent_dir, output_dir_name)
    if os.path.exists(output_dir):
        log.warning('directory "%s" already exists', output_dir)
    else:
        os.mkdir(output_dir)
    return output_dir


def gzip_files(file_list, debug=False):
    log = logging.getLogger(name=__name__)
    if debug is True:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.WARNING)
    for fp in file_list:
        with open(fp, 'rt') as src, gzip.open(fp + '.gz', 'wt') as dst:
            log.info('compressing file "%s"', fp)
            shutil.copyfileobj(fsrc=src, fdst=dst)
        os.remove(fp)


def ungzip_files(*fp_list, target_dir):
    uncompressed_fps = []
    for fp in fp_list:
        # remove '.gz' from the file path
        uncompressed_fp = os.path.join(target_dir, os.path.basename(fp)[:-3])
        print('unzipping {} into {}'.format(fp, uncompressed_fp))
        uncompressed_fps.append(uncompressed_fp)
        with gzip.open(fp, 'rt') as src, open(uncompressed_fp, 'wt') as dst:
            shutil.copyfileobj(fsrc=src, fdst=dst)
    return uncompressed_fps


def run_cmd(cmd_line_list, log_file, debug=False, **kwargs):
    log = logging.getLogger(name=__name__)
    if debug is True:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.WARNING)
    try:
        with open(log_file, 'at') as log_file:
            cmd_line_str = ' '.join((str(x) for x in cmd_line_list))
            log.info('executing "%s"', cmd_line_str)
            log_file.write('executing "{}"'.format(cmd_line_str))
            output = subprocess.run(
                cmd_line_list,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                **kwargs)
            log.info(output)
        return output
    except subprocess.CalledProcessError as c:
        logging.exception(c)
        print(c.message)
        print(c.cmd)
        print(c.output)
        raise c
    except Exception as e:
        logging.exception(e)
        print('blarg!')
        print(e)
        traceback.print_exc()
        raise e
