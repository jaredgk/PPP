import os
import sys
import subprocess
import shutil
import argparse
import glob
import logging

sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'jared')))

from vcf_reader_func import checkFormat
from logging_module import initLogger, logArgs
from plink import convert_haps_to_vcf
#from vcftools import bgzip_decompress_vcfgz
#from bcftools import convert_to_bcf, check_for_index, create_index

def check_shapeit_for_errors (shapeit_stdout, output_prefix):
    '''
        Checks the shapeit stdout for errors

        Parameters
        ----------
        shapeit_stdout : str
            shapeit stdout
        output_prefix : str
            Output filename prefix

        Raises
        ------
        Exception
            If shapeit stdout returns an error
    '''

    # Returns True if the job completed without error
    if 'Running time:' in str(shapeit_stdout):
        pass

    # Print output if not completed and no error found. Unlikely to be used, but included.
    else:
        # Remove intermediate files before reporting the error
        remove_intermediate_files(output_prefix, error_intermediates = True)
        raise Exception(str(shapeit_stdout))

def remove_intermediate_files (output_prefix, error_intermediates = False):
    '''
        Removes shapeit intermediate files

        This function is used to remove the various intermediate files created
        by shapeit. The exact intermediate files to be removed are defined by
        the error-state of shapeit. The function will also return warnings if
        the intermediate files were not found.

        Parameters
        ----------
        output_prefix : str
            Output filename prefix
        error_intermediates : bool, optional
            Defines if shapeit encountered an error

    '''
    if error_intermediates:

        # Check that the log file was created, give a warning otherwise
        if not os.path.isfile(output_prefix + '.phase.log'):
            logging.warning('shapeit intermediate file %s.phase.log does not exist' % output_prefix)
        else:
            # Remove shapeit log file
            os.remove(output_prefix + '.phase.log')

    else:

        # Check that the phase.ind.mm file was created, give a warning otherwise
        if not os.path.isfile(output_prefix + '.phase.ind.mm'):
            logging.warning('shapeit intermediate file %s.phase.ind.mm does not exist' % output_prefix)
        else:
            # Remove shapeit phase.ind.mm file
            os.remove(output_prefix + '.phase.ind.mm')

        # Check that the phase.snp.mm file was created, give a warning otherwise
        if not os.path.isfile(output_prefix + '.phase.snp.mm'):
            logging.warning('shapeit intermediate file %s.phase.snp.mm does not exist' % output_prefix)
        else:
            # Remove shapeit phase.snp.mm file
            os.remove(output_prefix + '.phase.snp.mm')

    # Check that the haps file was created, give a warning otherwise
    if not os.path.isfile(output_prefix + '.haps'):
        logging.warning('shapeit intermediate file %s.haps does not exist' % output_prefix)
    else:
        # Remove shapeit haps file
        os.remove(output_prefix + '.haps')

    # Check that the sample file was created, give a warning otherwise
    if not os.path.isfile(output_prefix + '.sample'):
        logging.warning('shapeit intermediate file %s.sample does not exist' % output_prefix)
    else:
        # Remove shapeit sample file
        os.remove(output_prefix + '.sample')

    logging.info('shapeit-related files removed')

def standard_shapeit_call (shapeit_call_args, output_prefix):
    '''
        Calls shapeit using subprocess

        This function is used to call shapeit and passes the resulting stdout
        to check_shapeit_for_errors to check for errors. The function also
        passes output_prefix to check_shapeit_for_errors to delete shapeit
        intermediate files if shapeit results in an error.

        Parameters
        ----------
        shapeit_call_args : list
            Argument list for shapeit
        output_prefix : str
            Output filename prefix

    '''

    logging.info('shapeit phasing parameters assigned')

    # Phasing subprocess call
    phase_call = subprocess.Popen(['shapeit'] + shapeit_call_args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    phase_stdout, phase_stderr = phase_call.communicate()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        phase_stdout = phase_stdout.decode()

    # Check shapeit call for errors
    check_shapeit_for_errors(phase_stdout, output_prefix)

    logging.info('shapeit phasing complete (HAPS format)')

def call_shapeit (shapeit_call_args, output_prefix, output_format):
    '''
        Calls shapeit and automates file conversions

        The function is used to call shapeit and also automates conversion to
        VCF, VCF.GZ, and BCF using plink2

        Parameters
        ----------
        shapeit_call_args : list
            Argument list for shapeit
        output_prefix : str
            Output filename prefix
        output_format : str
            Output file format

    '''

    # Standard call to beagle
    standard_shapeit_call(shapeit_call_args, output_prefix)

    # Convert haps-format to vcf
    convert_haps_to_vcf(output_prefix, output_format)

    logging.info('HAPS conversion to VCF complete')
