import os
import sys
import copy
import logging
import subprocess
import string
import random

from pgpipe.vcf_reader_func import checkFormat
from pgpipe.bcftools import check_bcftools_for_errors
from pgpipe.misc import confirm_executable

def log_vcftools_reference ():

    # Write the log header
    logging.info('Please Reference alongside the PPP:\n')

    # Write the reference
    logging.info('Danecek, P. et al. The variant call format and VCFtools. '
                 'Bioinformatics 27, 2156-2158 (2011).')

def check_for_vcftools_output (vcftools_output):
    '''
        Checks for the previous vcftools output

        Confirms that neither a previous vcftools log or output file exists.

        Parameters
        ----------
        vcftools_output : str
            Specifies the output filename to be checked

        Raises
        ------
        IOError
            If the vcftools output file exists
        IOError
            If the vcftools log file exists

    '''
    # Check if output file already exists
    if os.path.isfile(vcftools_output):
        raise IOError('VCF output file already exists')

    logging.info('Output file assigned')

    # Check if log file already exists
    if os.path.isfile(vcftools_output + '.log'):
        raise IOError('Log file already exists')

    logging.info('Log file assigned')

def assign_vcftools_unique_prefix (output_prefix, output_format, random_seed = None, string_size = 6, string_limit = 10, assignment_limit = 100):

    '''
        Assigns a unique vcftools filename prefix

        Used to assign a unique filename prefix for vcftools jobs. If an
        output file with the same prefix and suffix, either from previous 
        or ongoing jobs, this function will generate a unique prefix. This 
        function is only used if no prefix has been specified by the user.

        Parameters
        ----------
        output_prefix : str
            Specifies the filename prefix
        output_format : str
            Specifies the file format suffix
        random_seed : int, str, optional
            Specifies the random seed
        string_size : int, optional
            Specifies the number of character in the random string
        string_limit  : int, optional
            Specifies the character limit of the random string
        assignment_limit : int, optional
            Specifies a limit the number of assignment attempts

        Returns
        -------
        updated_prefix: str
            Specifies the unqiue prefix


        Raises
        ------
        Exception
            If unable to assign a unique filename

    '''

    # Assign the random seed, with a string or int
    random.seed(random_seed)

    # Assign the assignment counter, used to avoid infinite loops
    assignment_counter = 0

    # Generate a random string
    random_string = ''.join(random.choice(string.ascii_uppercase + string.digits) for digit in range(string_size))

    # Holds the updated prefix that will be returned when unique
    updated_prefix = output_prefix + '.' + random_string

    # Loop until a unique filename is created
    while os.path.isfile(updated_prefix + output_format):

        # Add to the assignment counter
        assignment_counter += 1

        # Check if assignment counter has reached it's limit
        if assignment_counter > assignment_limit:

            # Set the assignment counter to zero
            assignment_counter = 0

            # Add to the string size
            string_size += 1

            # Check if string size has reached it's limit and report the error
            if string_size > string_limit:
                raise Exception('Unable to assign unique intermediate output')

        # Generate a random string
        random_string = ''.join(random.choice(string.ascii_uppercase + string.digits) for digit in range(string_size))

        # Holds the updated prefix that will be returned when unique
        updated_prefix = output_prefix + '.' + random_string

    # Return the intermediate output filename
    return updated_prefix

def assign_vcftools_filename_prefix (output_prefix, output_format, output_filename):

    '''
        Assigns a vcftools prefix using a filename

        Used to assign a unique prefix for vcftools jobs using the user
        specified filename. This is to avoid output file with the same 
        prefix, either from previous or ongoing jobs. This function is 
        only used if no prefix has been specified by the user.

        Parameters
        ----------
        output_prefix : str
            Specifies the filename prefix
        output_format : str
            Specifies the file format suffix
        output_filename : str
            Specifies the filename specified by the user

        Returns
        -------
        unique_prefix: str
            Specifies the unqiue prefix

        Raises
        ------
        Exception
            If unable to assign a unique filename

    '''

    # List to hold the complete file path
    file_path_list = []

    # Split the filename
    split_file_path = os.path.split(output_filename)

    while split_file_path[1]:

        # Add the split path section to the file path list
        file_path_list = [split_file_path[1]] + file_path_list

        # Split the filename
        split_file_path = os.path.split(split_file_path[0])

    # Save the updated prefix
    updated_prefix = ''.join(file_path_list).replace('.','')

    # Check if the file already exists 
    if os.path.isfile(updated_prefix + output_format):
        raise Exception('Unable to assign prefix output. %s already exists' % (updated_prefix + output_format))

    return updated_prefix

def check_vcftools_for_errors (vcftools_stderr):
    '''
        Checks the vcftools stderr for errors

        Parameters
        ----------
        vcftools_stderr : str
            vcftools stderr

        Raises
        ------
        IOError
            If vcftools stderr returns an error
    '''

    # Returns True if the job completed without error
    if 'Run Time' in str(vcftools_stderr):
        pass

    # Print output for vcftools if error is detected
    elif 'Error' in str(vcftools_stderr):
        # Splits log into list of lines
        vcftools_stderr_lines = vcftools_stderr.splitlines()
        # Prints the error(s)
        raise Exception('\n'.join((output_line for output_line in vcftools_stderr_lines if output_line.startswith('Error'))))

    # Print output if not completed and no error found. Unlikely to be used, but included.
    else:
        raise Exception(vcftools_stderr)

def produce_vcftools_output (output, filename, append_mode = False, strip_header = False):
    '''
        Creates the vcftools output file

        This function will create an output file from the vcftools stdout.
        Please run `check_vcftools_for_errors` prior to check that vcftools
        finished without error.

        Parameters
        ----------
        output : str
            vcftools stdout
        filename : str
            Specifies the filename for the output file
        append_mode : bool
            Used to create a single output file from multiple calls
        strip_header : bool
            Used to remove the header if not needed

        Returns
        -------
        output : file
            vcftools output file

    '''

    # Check if the header should be stripped
    if strip_header:
        output = ''.join(output.splitlines(True)[1:])

    # Check if single log file is required from multiple calls
    if append_mode:
        vcftools_log_file = open(filename,'a')
    else:
        vcftools_log_file = open(filename,'w')

    vcftools_log_file.write(str(output))
    vcftools_log_file.close()

def produce_vcftools_log (output, filename, append_mode = False):
    '''
        Creates the vcftools log file

        This function will create a log file from the vcftools stderr. Please
        run `check_vcftools_for_errors` prior to check that vcftools finished
        without error.

        Parameters
        ----------
        output : str
            vcftools stderr
        filename : str
            Specifies the filename for the log file
        append_mode : bool
            Used to create a single log file from multiple calls

        Returns
        -------
        output : file
            vcftools log file

    '''
    # Check if single log file is required from multiple calls
    if append_mode:
        vcftools_log_file = open(filename + '.log','a')
    else:
        vcftools_log_file = open(filename + '.log','w')

    vcftools_log_file.write(str(output))
    vcftools_log_file.close()

def assign_vcftools_input_arg (filename):
    '''
        Confirms file format for vcftools

        Parameters
        ----------
        filename : str
            Specifies the input filename of unknown format

        Returns
        -------
        list
            Returns vcftools input command for `filename`

        Raises
        ------
        IOError
            If filename is an unknown file format
    '''

    # True if file extensions is recognized by vcftools
    if filename.endswith('.vcf') or filename.endswith('.vcf.gz') or filename.endswith('.bcf'):
        # Assign the associated input command
        if filename.endswith('.vcf'):
            return ['--vcf', filename]
        elif filename.endswith('.vcf.gz'):
            return ['--gzvcf', filename]
        elif filename.endswith('.bcf'):
            return ['--bcf', filename]

    # True if file extension is unknown or not recognized
    else:

        # Checks if the file is unzipped, bgzipped, or gzipped
        vcfname_format = checkFormat(filename)

        # Assign the associated input command, or return an error.
        if vcfname_format == 'vcf':
            return ['--vcf', filename]
        elif vcfname_format == 'bgzip':
            return ['--gzvcf', filename]
        elif vcfname_format == 'bcf':
            return ['--bcf', filename]
        else:
            raise Exception('Unknown VCF file format')

def check_bgzip_for_errors (bgzip_stderr):
    '''
        Checks the bgzip stderr for errors

        Parameters
        ----------
        bgzip_stderr : str
            bgzip stderr

        Raises
        ------
        IOError
            If bgzip stderr returns an error
    '''

    if bgzip_stderr:
        raise IOError('Error occured while compressing the vcf file')

def bgzip_decompress_vcfgz (vcfgz_filename, out_prefix = None, keep_original = False):
        '''
            Converts a vcf.gz to vcf

            The function automates bgzip to decompress a vcf.gz file into a vcf

            Parameters
            ----------
            vcfgz_filename : str
                The file name of the vcf.gz file to be decompressed
            out_prefix : str
                Output file prefix (i.e. filename without extension)
            keep_original : bool
                Specifies if the original file should be kept

            Raises
            ------
            IOError
                Error in creating the compressed file
        '''

        # Run bgzip with stdout piped to file
        if keep_original or out_prefix:

            if out_prefix:

                # Assign the bgzip filename
                vcf_filename = out_prefix + '.vcf'

            else:

                # Seperate into path and filename
                split_path, split_filename = os.path.split(vcfgz_filename)

                # Remove any file extensions
                vcf_basename = split_filename.split(os.extsep)[0] + '.vcf'

                # Join path and filename
                vcf_filename = os.path.join(split_path, vcf_basename)

            # Create the output file
            vcf_file = open(vcf_filename, 'w')

            # bgzip subprocess call
            bgzip_call = subprocess.Popen(['bgzip', '-dc', vcfgz_filename], stdout = vcf_file, stderr = subprocess.PIPE)

        # Run bgzip normally
        else:

            # bgzip subprocess call
            bgzip_call = subprocess.Popen(['bgzip', '-d', vcfgz_filename], stdout = subprocess.PIPE, stderr = subprocess.PIPE)

        # Save the stdout and stderr from bgzip
        bgzip_out, bgzip_err = bgzip_call.communicate()

        # Check that output file was compressed correctly
        check_bgzip_for_errors(bgzip_err)

        # Delete input when also using an output prefix
        if out_prefix and not keep_original:
            os.remove(vcfgz_filename)

def bgzip_compress_vcf (vcf_filename, out_prefix = '', keep_original = False):
        '''
            Converts a vcf to vcf.gz

            The function automates bgzip to compress a vcf file into a vcf.gz

            Parameters
            ----------
            vcf_filename : str
                The file name of the vcf file to be compressed
            keep_original : bool
                Specifies if the original file should be kept

            Raises
            ------
            IOError
                Error in creating the compressed file
        '''

        # Confirm where the specifed executable is located
        bgzip_path = confirm_executable('bgzip')

        # Check if the executable was found
        if not bgzip_path:
            raise IOError('bgzip not found. Please confirm the executable is installed')

        # Compress and keep the original file
        if keep_original or out_prefix:

            if out_prefix:

                # Assign the filename
                vcfgz_filename = out_prefix + '.vcf.gz'

            else:

                # Seperate into path and filename
                split_path, split_filename = os.path.split(vcfgz_filename)

                # Remove any file extensions
                vcfgz_basename = split_filename.split(os.extsep)[0] + '.vcf.gz'

                # Join path and filename
                vcfgz_filename = os.path.join(split_path, vcfgz_basename)


            # Create the output file
            vcfgz_file = open(vcfgz_filename, 'w')

            # bgzip subprocess call
            bgzip_call = subprocess.Popen([bgzip_path, '-c', vcf_filename], stdout = vcfgz_file, stderr = subprocess.PIPE)

        # Compress and do not keep the original file
        else:

            # bgzip subprocess call
            bgzip_call = subprocess.Popen([bgzip_path, vcf_filename], stdout = subprocess.PIPE, stderr = subprocess.PIPE)

        # Save the stdout and stderr from bgzip
        bgzip_out, bgzip_err = bgzip_call.communicate()

        # Check that output file was compressed correctly
        check_bgzip_for_errors(bgzip_err)

def cvt_vcftools_site_to_bed (vcftools_out_str):
    # Check if str in the header
    if 'CHROM' not in vcftools_out_str or 'POS' not in vcftools_out_str:
        # Split the line into a list
        vcftools_out_data = vcftools_out_str.strip().split('\t')
        # Calc chromEnd
        chrom_end = copy.deepcopy(vcftools_out_data[1])
        # Convert the chromStart to chromStart0
        vcftools_out_data[1] = int(vcftools_out_data[1]) - 1
        # Add chrom_end to the list
        vcftools_out_data = vcftools_out_data + [chrom_end]
        # Return the list as a string (with newline element)
        return '\t'.join(map(str, vcftools_out_data)) + '\n'
    else:
        # Remove the header
        return ''

def pipe_vcftools_to_bed_file (vcftools_call_args, output_filename):

    '''
        Pipes site-file output of vcftools to a bed formmated file

        The purpose of this function is to avoid creating large uncompressed
        vcf files by directly piping the output of vcftools to bgzip. This
        results in creating a vcf.gz file without any intermediates.

        Parameters
        ----------
        vcftools_call_args : list
            vcftools arguments
        output_filename : str
            Filename of the bed file

    '''

    # Confirm where the specifed executable is located
    vcftools_path = confirm_executable('vcftools')

    # Check if the executable was found
    if not vcftools_path:
        raise IOError('vcftools not found. Please confirm the executable is installed')

    # Open vcftools pipe
    vcftools_call = subprocess.Popen([vcftools_path, '--stdout'] + list(map(str, vcftools_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Create the bed file
    bed_output = open(output_filename, 'w')

    try:
        # Iterate the vcftools stdout unless error occurs
        for vcftools_stdout_line in iter(vcftools_call.stdout.readline, b''):
            bed_output.write(cvt_vcftools_site_to_bed(vcftools_stdout_line))
        # Close the bed file
        bed_output.close()
    except:
        # Close the bed file
        bed_output.close()
        # Delete the file
        os.remove(output_filename)

    # Wait for vctools to finish
    vcftools_call.wait()

    # Close the vcftools stdout
    vcftools_call.stdout.close()

    # Read the vcftools stderr
    vcftools_stderr = vcftools_call.stderr.read()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        vcftools_stderr = vcftools_stderr.decode()

    # Check that the log file was created correctly
    check_vcftools_for_errors(vcftools_stderr)

    logging.info('vcftools call complete')

    return vcftools_stderr

def pipe_vcftools_bgzip (vcftools_call_args, output_filename):
    '''
        Pipes the output of vcftools to bgzip

        The purpose of this function is to avoid creating large uncompressed
        vcf files by directly piping the output of vcftools to bgzip. This
        results in creating a vcf.gz file without any intermediates.

        Parameters
        ----------
        vcftools_call_args : list
            vcftools arguments
        output_filename : str
            Filename of the compressed vcf file

    '''

    # Confirm where the specifed executable is located
    vcftools_path = confirm_executable('vcftools')

    # Check if the executable was found
    if not vcftools_path:
        raise IOError('vcftools not found. Please confirm the executable is installed')

    vcftools_call = subprocess.Popen([vcftools_path, '--stdout'] + list(map(str, vcftools_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Create bgzip output file
    bgzip_output = open(output_filename, 'wb')

    # Confirm where the specifed executable is located
    bgzip_path = confirm_executable('bgzip')

    # Check if the executable was found
    if not bgzip_path:
        raise IOError('bgzip not found. Please confirm the executable is installed')

    # bgzip subprocess call
    bgzip_call = subprocess.Popen([bgzip_path], stdin = vcftools_call.stdout, stdout = bgzip_output, stderr = subprocess.PIPE)

    # Wait for vctools to finish
    vcftools_call.wait()

    # Close the vcftools stdout
    vcftools_call.stdout.close()

    # Read the vcftools stderr
    vcftools_stderr = vcftools_call.stderr.read()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        vcftools_stderr = vcftools_stderr.decode()

    # Check that the log file was created correctly
    check_vcftools_for_errors(vcftools_stderr)

    # Wait for bgzip to finish
    bgzip_call.wait()

    # Close the compressed vcf file
    bgzip_output.close()

    # Save the stderr from bgzip, stdout = None
    bgzip_stdout, bgzip_stderr = bgzip_call.communicate()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        bgzip_stderr = bgzip_stderr.decode()

    # Check that output file was compressed correctly
    check_bgzip_for_errors(bgzip_stderr)

    logging.info('vcftools and bgzip calls complete')

    return vcftools_stderr

def pipe_vcftools_bcftools (vcftools_call_args, output_filename):
    '''
        Pipes the output of vcftools to bcftools

        The purpose of this function is to avoid the vcftools command
        --recode-bcf that may result in malformed BCF files. To avoid large
        uncompressed intermediates, this function pipes the stdout of vcftools
        to bcftools.

        Parameters
        ----------
        vcftools_call_args : list
            vcftools arguments
        output_filename : str
            Filename of the BCF file

    '''

    # Confirm where the specifed executable is located
    vcftools_path = confirm_executable('vcftools')

    # Check if the executable was found
    if not vcftools_path:
        raise IOError('vcftools not found. Please confirm the executable is installed')

    vcftools_call = subprocess.Popen([vcftools_path, '--stdout'] + list(map(str, vcftools_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Holds the arguments to convert to BCF format
    convert_args = ['view', '-O', 'b']

    # Create bgzip output file
    bcftools_output = open(output_filename, 'wb')

    # Confirm where the specifed executable is located
    bcftools_path = confirm_executable('bcftools')

    # Check if the executable was found
    if not bcftools_path:
        raise IOError('bcftools not found. Please confirm the executable is installed')

    # bcftools subprocess call
    bcftools_call = subprocess.Popen([bcftools_path] + convert_args, stdin = vcftools_call.stdout, stdout = bcftools_output, stderr = subprocess.PIPE)

    # Wait for vctools to finish
    vcftools_call.wait()

    # Close the vcftools stdout
    vcftools_call.stdout.close()

    # Read the vcftools stderr
    vcftools_stderr = vcftools_call.stderr.read()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        vcftools_stderr = vcftools_stderr.decode()

    # Check that the log file was created correctly
    check_vcftools_for_errors(vcftools_stderr)

    # Wait for bgzip to finish
    bcftools_call.wait()

    # Save the stderr from bgzip, stdout = None
    bcftools_stdout, bcftools_stderr = bcftools_call.communicate()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        bcftools_stderr = bcftools_stderr.decode()

    # Check that output file was compressed correctly
    check_bcftools_for_errors(bcftools_stderr)

    logging.info('vcftools and bcftools calls complete')

    return vcftools_stderr

def pipe_vcftools_to_file (vcftools_call_args, output_filename, append_output = False):
    '''
        Pipes file output of vcftools to a standard file

        The function calls vcftools. Returns the stderr of vcftools to
        create log file of the call. The function may be used to append multiple
        calls to vcftools to a single file

        Parameters
        ----------
        vcftools_call_args : list
            vcftools arguments
        append_output : bool
            The output file should be written in append mode

        Returns
        -------
        vcftools_err : str
            vcftools log output

        Raises
        ------
        Exception
            If vcftools stderr returns an error
    '''

    # Confirm where the specifed executable is located
    vcftools_path = confirm_executable('vcftools')

    # Check if the executable was found
    if not vcftools_path:
        raise IOError('vcftools not found. Please confirm the executable is installed')

    # Open vcftools pipe
    vcftools_call = subprocess.Popen([vcftools_path, '--stdout'] + list(map(str, vcftools_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Check if the output should be opened in append mode
    if append_output:
        # Create the output file (in append mode)
        output_file = open(output_filename, 'a')
    else:
        # Create the output file (in write mode)
        output_file = open(output_filename, 'w')

    try:
        # Create iterator of the vcftools stdout
        stdout_iter = iter(vcftools_call.stdout.readline, b'')

        # Check if the output is being appended and the file is empty
        if append_output and os.stat(output_filename).st_size != 0:
            # Skip the header if the file isn't empty and appending
            next(stdout_iter)

        # Iterate the vcftools stdout
        for vcftools_stdout_line in stdout_iter:

            # Check if code is running in python 3
            if sys.version_info[0] == 3:
                # Convert bytes to string
                vcftools_stdout_line = vcftools_stdout_line.decode()

            output_file.write(vcftools_stdout_line)

        # Close the bed file
        output_file.close()

    except:
        # Close the bed file
        output_file.close()
        # Delete the file
        os.remove(output_filename)

        raise Exception('vcftools to python pipe error')

    # Wait for vctools to finish
    vcftools_call.wait()

    # Close the vcftools stdout
    vcftools_call.stdout.close()

    # Read the vcftools stderr
    vcftools_stderr = vcftools_call.stderr.read()

    # Close the vcftools stderr
    vcftools_call.stderr.close()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        vcftools_stderr = vcftools_stderr.decode()

    # Check that the log file was created correctly
    check_vcftools_for_errors(vcftools_stderr)

    logging.info('vcftools call complete')

    return vcftools_stderr

def standard_vcftools_call (vcftools_call_args):
    '''
        Calls vcftools

        The function calls vcftools. Returns the stderr of vcftools to
        create log file of the call.

        Parameters
        ----------
        vcftools_call_args : list
            vcftools arguments

        Returns
        -------
        vcftools_out : str
            vcftools call output
        vcftools_err : str
            vcftools log output

        Raises
        ------
        Exception
            If vcftools stderr returns an error
    '''

    # Confirm where the specifed executable is located
    vcftools_path = confirm_executable('vcftools')

    # Check if the executable was found
    if not vcftools_path:
        raise IOError('vcftools not found. Please confirm the executable is installed')

    # vcftools subprocess call without stdout
    vcftools_call = subprocess.Popen([vcftools_path] + list(map(str, vcftools_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Wait for vcftools to finish
    vcftools_stdout, vcftools_stderr = vcftools_call.communicate()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        vcftools_stderr = vcftools_stderr.decode()

    logging.info('vcftools call complete')

    # Check that the log file was created correctly
    check_vcftools_for_errors(vcftools_stderr)

    return vcftools_stderr

def call_vcftools (vcftools_call_args, output_format = None, output_filename = None, append_mode = False):
    '''
        Calls vcftools

        The function calls vcftools. Returns the stderr of vcftools to
        create log file of the call.

        Parameters
        ----------
        vcftools_call_args : list
            vcftools arguments
        output_format : str, optional
            Allows for specific formats to be given which may result in
            the stdout of vcftools being piped to another script or 
            program
        output_filename : str, optional
            The output filename to be assigned if the stdout of vcftools 
            is piped
        append_mode : bool, optional
            Allows for vcftools stdout to append a file. Currently 
            incompatible with output_format

        Returns
        -------
        vcftools_err : str
            vcftools log output

        Raises
        ------
        Exception
            If vcftools stderr returns an error
    ''' 

    # Check if an output format was specified
    if output_format and output_filename: 

        # Check if the output is a bgzipped vcf
        if output_format == 'vcf.gz':

            # Pipe vcftools stdout to bgzip to create a bgzipped vcf
            vcftools_err = pipe_vcftools_bgzip(vcftools_call_args, output_filename)

        # Check if the output is a bcf
        elif output_format == 'bcf':

            # Pipe vcftools stdout to bgzip to create a bgzipped vcf
            vcftools_err = pipe_vcftools_bcftools(vcftools_call_args, output_filename)

        # Check if the output is a bed-based file
        elif output_format in ['removed_bed','kept_bed']:

            # Pipe vcftools stdout to bed file
            vcftools_err = pipe_vcftools_to_bed_file(vcftools_call_args, output_filename)

        # Check if the output is another format, that does not require a pipe
        else:
            
            # Pipe the vcftools stdout to a standard file
            vcftools_err = pipe_vcftools_to_file(vcftools_call_args, output_filename)

    # Check if a filename but no format was specified 
    elif output_filename:

        # Check if the output should be written in append mode
        if append_mode:

            # Pipe the vcftools stdout to a standard file and allow appending
            vcftools_err = pipe_vcftools_to_file(vcftools_call_args, output_filename, append_output = True)

        else:

            # Pipe the vcftools stdout to a standard file
            vcftools_err = pipe_vcftools_to_file(vcftools_call_args, output_filename)

    # If no format or filename is given, run vcftools with just the passed arguments
    else:

        # Call vcftools under standard conditions, if not 
        vcftools_err = standard_vcftools_call(vcftools_call_args)

    # Return the log
    return vcftools_err

