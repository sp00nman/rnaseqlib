"""
Collections of utility functions.
"""

import logging
from os import (system, remove, mkdir)
from os.path import (split, splitext, join, exists)


def run_cmd(
        message,
        command,
        debug):
    """
    Print stdout message and log commands, return non-zero exit codes
    if execution failed.
    :param message: Stdout message to be printed.
    :param command: Command to be executed.
    :param debug: Do not execute command if debug is set to 1
    :return: integer, 1 -success, 0 - fail
    """

    print message
    print command # only for debug reasons

    logging.info(message)
    logging.debug(command)
    status = 0

    if not debug:
        status = system(command)
        if status != 0:
            logging.warning("command '%s' returned non-zero "
                            "status: %d'" % (command, status))
    return status


def load_dictionary(
        input_file,
        sep=':'):
    """
    Reads in files of the following format:
    key1: value1
    key2: value2
    :param file: plain text file with
    :return: dictionary of key-value pairs
    """

    d = {}

    try:
        file_handle = open(input_file)

        for line in file_handle:
            key_value_pair = line.split(sep)
            d[key_value_pair[0].strip(' ')] = \
                key_value_pair[1].strip(' ').rstrip('\n')

    except IOError:
        print('Key value read-in file missing.')

    finally:
        file_handle.close()

    return d


def create_output_dir(
        output_dir,
        project_name):
    """
    Create project directory.
    """
    if not exists(output_dir + "/" + project_name):
        logging.info('Create folder %s' % output_dir)
        try:
            mkdir(output_dir + "/" + project_name, 0777)
        except IOError, e:
            exit('%s\nFailed to create directory', (e, output_dir))


def load_tab_delimited(
    input_file,
    sep='\t'
):

    try:
        file_handle = open(input_file)

        vcfs = [line.rstrip('\n').split(sep) for line in file_handle]

    except IOError:
        print('Patient VCF file missing')


    finally:
        file_handle.close()

    return vcfs