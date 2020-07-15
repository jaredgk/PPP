import os
import sys
import subprocess
import shutil
import argparse
import glob
import logging

#sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'pppipe')))

from pgpipe.vcf_reader_func import checkFormat
from pgpipe.logging_module import initLogger, logArgs
from pgpipe.vcftools import bgzip_decompress_vcfgz
from pgpipe.bcftools import convert_vcf, check_for_index, create_index
from pgpipe.misc import confirm_executable

def delete_beagle_log (output_prefix):
    '''
        Delete beagle log file

        This function is used to delete beagle's log file if an error is
        encountered. A warning is produced if the log file cannot be found.

        Parameters
        ----------
        output_prefix : str
            Output file prefix
    '''

    # Check that log file exists, if not return warning
    if not os.path.isfile(output_prefix + '.log'):
        logging.warning('beagle log file %s.log does not exist' % output_prefix)
    else:
        os.remove(output_prefix + '.log')

def check_beagle_for_errors (beagle_stderr, output_prefix):
    '''
        Checks the beagle stdout for errors

        Parameters
        ----------
        beagle_stderr : str
            beagle stderr
        output_prefix : str
            Output file prefix

        Raises
        ------
        Exception
            If beagle stdout returns an error
    '''

    # Check if beagle completed without an error
    if not beagle_stderr.strip():
        pass

    # Print missing data message if that is likely
    elif 'ERROR: genotype is missing allele separator:' in str(beagle_stderr):
        # Delete the beagle log file
        delete_beagle_log(output_prefix)

        # Store reported error
        error_reported = 'ERROR: genotype is missing allele separator'
        # Store message for user about error
        user_message = 'Please confirm the input has no missing data.'
        # Report on the error
        raise Exception(error_reported + '\n' + user_message)

    # Print output for beagle if error is detected
    elif 'ERROR:' in str(beagle_stderr):
        # Delete the beagle log file
        delete_beagle_log(output_prefix)

        # Splits log into list of lines
        beagle_stderr_lines = beagle_stderr.splitlines()
        # Prints the error(s)
        raise Exception('\n'.join((output_line for output_line in beagle_stderr_lines if output_line.startswith('ERROR:'))))

    # Print output if not completed and no error found. Unlikely to be used, but included.
    else:
        # Delete the beagle log file
        delete_beagle_log(output_prefix)

        raise Exception(beagle_stderr)

def check_for_beagle_intermediate_files (output_prefix, output_format, overwrite = False):

    # Check if the output format is not vcf.gz
    if output_format != 'vcf.gz':

        # Check if the intermediate should be removed
        if not overwrite:

            # Check if an intermediate file exists and riase error
            if os.path.isfile(output_prefix + '.vcf.gz'):
                raise Exception('Beagle intermediate (%s) found. Use --overwrite to ignore')

def standard_beagle_call (beagle_path, beagle_call_args, output_prefix):
    '''
        Calls beagle using subprocess

        This function is used to call beagle under standard conditions. The
        functions then passes the stderr to check_beagle_for_errors to check
        for errors.

        Parameters
        ----------
        beagle_path : str
            Path to beagle.jar
        beagle_call_args : list
            Argument list for beagle
        output_prefix : str
            Output file prefix
    '''

    # Assign location of picard jar file
    if beagle_path is None:

        # Create a string with the picard path
        beagle_jar = confirm_executable('beagle.jar')
        
    else:
        # Use path if specified
        beagle_jar = os.path.join(picard_path, 'beagle.jar')

    # Check if executable is installed
    if not beagle_jar:
         raise IOError('beagle.jar not found. Please confirm the executable is installed')

    logging.info('beagle phasing parameters assigned')

    # Phasing subprocess call
    phase_call = subprocess.Popen(['java', '-jar', beagle_jar] + beagle_call_args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    phase_stdout, phase_stderr = phase_call.communicate()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        phase_stderr = phase_stderr.decode()

    # Check beagle call for errors
    check_beagle_for_errors(phase_stderr, output_prefix)

    logging.info('beagle phasing complete')

def call_beagle (beagle_path, beagle_call_args, output_prefix, output_format):
    '''
        Automates beagle calls

        This function passes the argument list to standard_beagle_call. Once the
        beagle call has finished, the function will automatically convert the
        bgzip compressed output of beagle to BCF and VCF, if either format is
        specified.

        Parameters
        ----------
        beagle_path : str
            Path to beagle.jar
        beagle_call_args : list
            Argument list for beagle
        output_prefix : str
            Output file prefix
        output_format : str
            Output file format
    '''

    # Standard call to beagle
    standard_beagle_call(beagle_path, beagle_call_args, output_prefix)

    # Check if the desired format is VCF
    if output_format == 'vcf':

        # Decompress the VCF file
        bgzip_decompress_vcfgz(output_prefix + '.vcf.gz')

    # Check if the desired format is BCF
    elif output_format == 'bcf':

        # Check if there is an index file
        if check_for_index(output_prefix + '.vcf.gz') == False:
            # Create an index if not found
            create_index(output_prefix + '.vcf.gz')
        
        # Convert vcf.gz to bcf
        convert_vcf(output_prefix + '.vcf.gz', output_prefix, output_format, overwrite = True)
