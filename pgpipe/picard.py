import os
import sys
import subprocess
import shutil
import argparse
import logging
import copy

# Call PPP-based scripts
#sys.path.insert(0, os.path.abspath(os.path.join(os.pardir,'pppipe')))

from pgpipe.logging_module import initLogger, logArgs
from pgpipe.fasta import check_format
from pgpipe.misc import confirm_executable

def check_ref_file_association (ref_filename, dict_filename, dictonary_ext = 'dict', gunzip_ext = 'gz'):

    # Check the format of the reference file
    ref_format = check_format(ref_filename)

    # Get the basename of the reference
    ref_basename = os.path.basename(ref_filename)

    # Check if the file is fasta formatted
    if ref_format == 'fasta':

        # Check if the file has an extension
        if ref_basename.count(os.extsep) != 0:

            # Split the filename to seperate an ext
            ref_filename_wo_ext, ref_ext = ref_filename.split(os.extsep, 1)

            # Check if the dictonary file is correctly associated
            if dict_filename == ref_filename_wo_ext + os.extsep + dictonary_ext:
                return True

    # Check if the file is fasta formatted
    elif ref_format == 'gzip':

        # Check if the file has the gunzip file extension
        if ref_basename.endswith(gunzip_ext) and ref_basename.count(os.extsep) >= 2: 

            # Split the filename to seperate an ext
            ref_filename_wo_ext, ref_ext, gz_ext = ref_filename.split(os.extsep, 2)

            # Save the filename of the possible dictonary file
            dictonary_filename = ref_filename_wo_ext + os.extsep + dictonary_ext

            # Check if the dictonary file is correctly associated
            if dict_filename == ref_filename_wo_ext + os.extsep + dictonary_ext:
                return True

    # Return None, if association failed
    return None

def get_ref_dictonary (ref_filename, dictonary_ext = 'dict'):

    # Save the filename of the possible dictonary file
    dictonary_filename = ref_filename + os.extsep + dictonary_ext

    # Check if the dictonary file exists
    if os.path.isfile(dictonary_filename):

        # Return the filename
        return dictonary_filename

    # Save the filename before split
    previous_filename = copy.deepcopy(ref_filename)

    # Split the filename to seperate an ext
    ref_filename_wo_ext, ref_ext = ref_filename.split(os.extsep, 1)

    # Loop the reference filename until there are no more extentions
    while ref_filename_wo_ext != previous_filename:

        # Save the filename of the possible dictonary file
        dictonary_filename = ref_filename_wo_ext + os.extsep + dictonary_ext

        # Check if the dictonary file exists
        if os.path.isfile(dictonary_filename):

            # Return the filename
            return dictonary_filename

        # Split the file
        else:

            # Save the filename before split
            previous_filename = copy.deepcopy(ref_filename_wo_ext)

            # Split the filename to seperate an ext
            ref_filename_wo_ext, ref_ext = ref_filename.split(os.extsep, 1)

    # Return None, if no file is found
    return None

def produce_picard_log (output, filename):
    '''
        Creates the picard log file

        This function will create a log file from the picard stderr. Please
        run `check_picard_for_errors` prior to check that picard finished
        without error.

        Parameters
        ----------
        output : str
            picard stderr
        filename : str
            Specifies the filename for the log file

        Returns
        -------
        output : file
            picard log file

    '''

    picard_log_file = open(filename + '.log','w')
    picard_log_file.write(str(output))
    picard_log_file.close()

def remove_picard_index (out_filename, out_format):
    '''
        Removes the picard index file

        This function will remove the index file created by picard.

        Parameters
        ----------
        out_filename : str
            Specifies the filename for the vcf output file
        out_format : str 
            Specifies the file format of the vcf output

    '''
    
    # Check if the file is a gzipped VCF
    if out_format == 'vcf.gz':

        # Remove the index, if found
        if os.path.isfile(out_filename + '.tbi'):
            os.remove(out_filename + '.tbi')

    # Check if the file is a gzipped VCF
    if out_format == 'vcf':

        # Remove the index, if found
        if os.path.isfile(out_filename + '.idx'):
            os.remove(out_filename + '.idx')

    # Check if the file is a gzipped VCF
    if out_format == 'bcf':

        # Remove the index, if found
        if os.path.isfile(out_filename + '.idx'):
            os.remove(out_filename + '.idx')


def check_picard_for_errors (picard_stderr):
    '''
        Checks the picard stderr for errors

        Parameters
        ----------
        picard_stderr : str
            picard stderr

        Raises
        ------
        IOError
            If picard stderr returns an error
    '''

    # Check if there was no error
    if 'Elapsed time:' in picard_stderr:
        pass

    # Print exception output
    elif 'Exception' in picard_stderr:

        # Loop stderr lines
        for picard_stderr_line in picard_stderr.splitlines():

            # Check if the current line has the exception
            if 'Exception' in picard_stderr_line:
                raise Exception(picard_stderr_line)

    # Print error output. Build up as errors are discovered
    elif picard_stderr:
        raise Exception(picard_stderr)

def standard_picard_call (picard_path, picard_call_args, output_filename):
    '''
        Standard call of picard

        The function calls picard. Returns the stderr of picard to create a log
        file of the call.

        Parameters
        ----------
        picard_path : str
            Path to picard.jar
        picard_call_args : list
            picard arguments
        output_filename : str
            Output filename string

        Raises
        ------
        Exception
            If picard stderr returns an error
    '''

    # Assign location of picard jar file
    if picard_path is None:

        # Create a string with the picard path
        picard_jar = confirm_executable('picard.jar')

        #picard_jar = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'bin','picard.jar')

    else:

        # Use path if specified
        picard_jar = os.path.join(picard_path, 'picard.jar')

    # Check if executable is installed
    if not picard_jar:
         raise IOError('picard.jar not found. Please confirm the executable is installed')

    logging.info('picard parameters assigned')

    # Phasing subprocess call
    picard_call = subprocess.Popen(['java', '-jar', picard_jar] + picard_call_args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    picard_stdout, picard_stderr = picard_call.communicate()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        picard_stderr = picard_stderr.decode()

    # Check picard call for errors
    check_picard_for_errors(picard_stderr)

    logging.info('picard call complete')

    # Create the log file
    produce_picard_log (picard_stderr, output_filename)

    logging.info('picard log created')


def call_picard (picard_path, picard_call_args, output_filename):
    '''
        Automates picard calls

        This function passes the argument list to standard_picard_call.

        Parameters
        ----------
        picard_path : str
            Path to picard.jar
        picard_call_args : list
            Argument list for picard
        output_filename : str
            Output filename string
    ''' 

    # Standard call to picard
    standard_picard_call(picard_path, picard_call_args, output_filename)




