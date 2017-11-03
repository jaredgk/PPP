import os
import sys
import logging
import subprocess

sys.path.insert(0, os.path.abspath(os.path.join(os.pardir,'jared')))

import vcf_reader_func

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
        logging.error('Error occured while compressing the vcf file')
        raise IOError('Error occured while compressing the vcf file')

def bgzip_compress_vcf (vcf_filename):
        '''
            Converts a vcf to vcf.gz

            The function automates bgzip to compress a vcf file into a vcf.gz

            Parameters
            ----------
            vcf_filename : str
                The file name of the vcf file to be compressed

            Raises
            ------
            IOError
                Error in creating the compressed file
        '''

        # bgzip subprocess call
        bgzip_call = subprocess.Popen(['bgzip', vcf_filename], stdout = subprocess.PIPE, stderr = subprocess.PIPE)

        # Save the stdout and stderr from bgzip
        bgzip_out, bgzip_err = bgzip_call.communicate()

        # Check that output file was compressed correctly
        check_bgzip_for_errors(bgzip_err)


def call_vcftools_bgzip (vcftools_call_args, vcf_gz_filename):
    '''
        Calls vcftools and bgzip in tandem

        The function calls vcftools and pipes the output to bgzip to compress to
        create a vcf.gz file. Returns the stderr of vcftools to create log file
        of the call.

        Parameters
        ----------
        vcftools_call_args : list
            vcftools arguments
        vcf_gz_filename : str
            The output name of the compressed vcf file

        Returns
        -------
        vcftools_err : str
            vcftools log output

        Raises
        ------
        Exception
            If vcftools stderr returns an error
        Exception
            If bgzip stderr returns an error
    '''

    # vcftools subprocess call
    vcftools_call = subprocess.Popen(['vcftools', '--stdout'] + list(map(str, vcftools_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Create bgzip output file
    bgzip_output = open(vcf_gz_filename, 'wb')

    # bgzip subprocess call
    bgzip_call = subprocess.Popen(['bgzip'], stdin = vcftools_call.stdout, stdout = bgzip_output, stderr = subprocess.PIPE)

    # Close the vcftools stdout
    vcftools_call.stdout.close()

    # Save the stderr from bgzip, stdout = None
    bgzip_out, bgzip_err = bgzip_call.communicate()

    # Close the compressed vcf file
    bgzip_output.close()

    # Save the stderr from vcftools
    vcftools_err = vcftools_call.stderr.read()

    logging.info('vcftools and bgzip calls complete')

    # Check that the log file was created correctly
    check_vcftools_for_errors(vcftools_err)

    # Check that output file was compressed correctly
    check_bgzip_for_errors(bgzip_err)

    return vcftools_err

def call_vcftools (vcftools_call_args):
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
        vcftools_err : str
            vcftools log output

        Raises
        ------
        Exception
            If vcftools stderr returns an erro
    '''

    # vcftools subprocess call
    vcftools_call = subprocess.Popen(['vcftools'] + list(map(str, vcftools_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Wait for vcftools to finish
    vcftools_out, vcftools_err = vcftools_call.communicate()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        vcftools_out = vcftools_out.decode()
        vcftools_err = vcftools_err.decode()

    logging.info('vcftools call complete')

    # Check that the log file was created correctly
    check_vcftools_for_errors(vcftools_err)

    return vcftools_out, vcftools_err


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
        logging.error('VCF output file already exists')
        raise IOError('VCF output file already exists')

    logging.info('Output file assigned')

    # Check if log file already exists
    if os.path.isfile(vcftools_output + '.log'):
        logging.error('Log file already exists')
        raise IOError('Log file already exists')

    logging.info('Log file assigned')

def delete_vcftools_output (vcftools_output):
    '''
        Deletes previous vcftools output

        Confirms if previous vcftools output exists, and if so, deletes it

        Parameters
        ----------
        vcftools_output : str
            Specifies the output filename to be deleted

        Raises
        ------
        IOError
            If the vcftools output cannot be deleted
        IOError
            If the vcftools log cannot be deleted
    '''

    # Check if output file already exists
    if os.path.isfile(vcftools_output):
        try:
            # Delete the output
            os.remove(vcftools_output)
        except:
            logging.error('VCF output file cannot be deleted')
            raise IOError('VCF output file cannot be deleted')

    logging.info('Output file assigned')

    # Check if log file already exists
    if os.path.isfile(vcftools_output + '.log'):
        try:
            # Delete the output
            os.remove(vcftools_output + '.log')
        except:
            logging.error('Log file cannot be deleted')
            raise IOError('Log file cannot be deleted')

    logging.info('Log file assigned')

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
        logging.error('\n'.join((output_line for output_line in vcftools_stderr_lines if output_line.startswith('Error'))))
        raise Exception('\n'.join((output_line for output_line in vcftools_stderr_lines if output_line.startswith('Error'))))

    # Print output if not completed and no error found. Unlikely to be used, but included.
    else:
        logging.error(vcftools_stderr)
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
        vcfname_format = vcf_reader_func.checkFormat(filename)

        # Assign the associated input command, or return an error.
        if vcfname_format == 'nozip':
            return ['--vcf', filename]
        elif vcfname_format == 'bgzip':
            return ['--gzvcf', filename]
        elif vcfname_format == 'bcf':
            return ['--bcf', filename]
        else:
            logging.error('Unknown VCF file format')
            raise Exception('Unknown VCF file format')
