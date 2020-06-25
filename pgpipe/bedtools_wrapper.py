import os
import sys
import subprocess
import logging
import fileinput

from pgpipe.misc import confirm_executable

def log_bedtools_reference (out_filename, append_mode = False, ref_header = True):

    # Check if the file is to be written in append mode
    if append_mode:

        # Open the file
        log_file = open(out_filename + '.log', 'a')

        # Check if the ref header should be added
        if ref_header:
            log_file.write('\nPlease Reference alongside the PPP:\n')

    else:

        # Open the file
        log_file = open(out_filename + '.log', 'w')

        # Check if the ref header should be added
        if ref_header:
            log_file.write('Please Reference alongside the PPP:\n')

    log_file.write('')
    log_file.close()

    logging.info('Reference assigned')

def check_bedtools_for_errors (bedtools_stderr):
    '''
        Checks the bedtools stderr for errors

        Parameters
        ----------
        bedtools_stderr : str
            bedtools stderr

        Raises
        ------
        IOError
            If bedtools stderr returns an error
    '''

    # Print output if error found. Build up as errors are discovered
    if bedtools_stderr:

        # Splits log into list of lines
        bedtools_stderr_lines = bedtools_stderr.splitlines()

        # Prints the error(s)
        raise Exception('\n'.join((output_line.split(':')[1] for output_line in bedtools_stderr_lines if 'ERROR:' in output_line)))

def merge_bed_files (bed_files, bed_output_filename, optional_merge_args):
    '''
    Merge BED files

    The function calls bedtools

    Parameters
    ----------
    bedtools_call_args : list
    bedtools arguments

    Raises
    ------
    Exception
    If bedtools stderr returns an error
    '''

    # Confirm where the specifed executable is located
    bedtools_path = confirm_executable('bedtools')

    # Check if the executable was found
    if not bedtools_path:
        raise IOError('bedtools not found. Please confirm the executable is installed')

    # Assign the BEDtools input using fileinput
    bed_input = fileinput.input(files = bed_files, mode = 'rb')

    # Create the output file
    bed_output_file = open(bed_output_filename, 'w')

    # Call BEDtools to sort the input files
    bedtools_sort_call = subprocess.Popen([bedtools_path, 'sort', '-i', '-'], stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)

    # Loop the BEDtools input
    for bed_line in bed_input:

    # Pipe the BEDtools input to stdin
        bedtools_sort_call.stdin.write(bed_line)

    # Call BEDtools to merge the sorted input
    bedtools_merge_call = subprocess.Popen([bedtools_path, 'merge', '-i', '-'] + list(map(str, optional_merge_args)), stdin = bedtools_sort_call.stdout, stdout = bed_output_file, stderr = subprocess.PIPE)

    # Close the sort stdin
    bedtools_sort_call.stdin.close()

    # Close the sort stdout
    bedtools_sort_call.stdout.close()

    # Read the sort stderr
    bedtools_sort_stderr = bedtools_sort_call.stderr.read()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:

    # Convert bytes to string
        bedtools_sort_stderr = bedtools_sort_stderr.decode()

    # Check BEDtools sort for possible errors
    check_bedtools_for_errors(bedtools_sort_stderr)

    # Wait for sort call to finish
    bedtools_sort_call.wait()

    # Save the stderr from bgzip, stdout = None
    bedtools_merge_stdout, bedtools_merge_stderr = bedtools_merge_call.communicate()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:

        # Convert bytes to string
        bedtools_merge_stderr = bedtools_merge_stderr.decode()

    # Check BEDtools sort for possible errors
    check_bedtools_for_errors(bedtools_merge_stderr)

    # Close BEDtools output
    bed_output_file.close()

    # Close BEDtools output
    bed_output_file.close()

def standard_bedtools_call (bedtools_call_args, bed_output_filename):
    '''
        Standard call of bedtools

        The function calls bedtools. Returns the stderr of bedtools to create a log
        file of the call.

        Parameters
        ----------
        bedtools_call_args : list
            bedtools arguments

        Raises
        ------
        Exception
            If bedtools stderr returns an error
    '''

    # Create the output file
    bed_output_file = open(bed_output_filename, 'w')

    # Confirm where the specifed executable is located
    bedtools_path = confirm_executable('bedtools')

    # Check if the executable was found
    if not bedtools_path:
        raise IOError('bedtools not found. Please confirm the executable is installed')

    #bedtools_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'bin', 'bedtools')

    # bedtools subprocess call
    bedtools_call = subprocess.Popen([bedtools_path] + list(map(str, bedtools_call_args)), stdout = bed_output_file, stderr = subprocess.PIPE)

    # Wait for bedtools to finish
    bedtools_out, bedtools_err = bedtools_call.communicate()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        bedtools_err = bedtools_err.decode()

    logging.info('bedtools call complete')

    # Check that the log file was created correctly
    check_bedtools_for_errors(bedtools_err)
