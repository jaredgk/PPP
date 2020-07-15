import os
import sys
import subprocess
import shutil
import argparse
import glob
import logging
import pkg_resources

#sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'pppipe')))

from pgpipe.vcf_reader_func import checkFormat
from pgpipe.logging_module import initLogger, logArgs
from pgpipe.misc import confirm_executable

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
        remove_shapeit_intermediate_files(output_prefix, error_intermediates = True)
        raise Exception(str(shapeit_stdout))

def remove_shapeit_intermediate_files (output_prefix, error_intermediates = False):
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

def check_for_shapeit_intermediate_files (output_prefix, overwrite = False):

    # Check if intermediates files should not be overwritten
    if not overwrite:

        # List to hold intermediate files
        shapeit_intermediate_files = []
    
        # Check if a phase.ind.mm file exists
        if os.path.isfile(output_prefix + '.phase.ind.mm'):

            # Add the intermediate file to the list
            shapeit_intermediate_files.append(output_prefix + '.phase.ind.mm')

        # Check if a phase.snp.mm file exists
        if os.path.isfile(output_prefix + '.phase.snp.mm'):

            # Add the intermediate file to the list
            shapeit_intermediate_files.append(output_prefix + '.phase.snp.mm')

        # Check if a haps file exists
        if os.path.isfile(output_prefix + '.haps'):

            # Add the intermediate file to the list
            shapeit_intermediate_files.append(output_prefix + '.haps')

        # Check if a sample file exists
        if os.path.isfile(output_prefix + '.sample'):

            # Add the intermediate file to the list
            shapeit_intermediate_files.append(output_prefix + '.sample')

        # Check if intermediate_files were found, and report the error
        if shapeit_intermediate_files:
            raise Exception('shapeit intermediate files exist (%s)' % ', '.join(shapeit_intermediate_files))

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

    # Confirm where the specifed executable is located
    shapeit_path = confirm_executable('shapeit')

    # Check if the executable was found
    if not shapeit_path:
        raise IOError('shapeit not found. Please confirm the executable is installed')

    #shapeit_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'bin','shapeit')

    #sys.stderr.write(str(shapeit_path)+'\n')
    phase_call = subprocess.Popen([shapeit_path] + shapeit_call_args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
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
    