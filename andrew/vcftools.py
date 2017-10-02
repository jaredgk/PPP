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


def call_vcftools_bgzip (vcfname_arg, vcftools_call_args, vcf_gz_filename):
    '''
        Converts a vcf to vcf.gz

        The function pipes the output of vcftools to bgzip to compress a vcf
        file into a vcf.gz

        Parameters
        ----------
        vcf_call : Popen
            The vcftools subprocess call
        vcf_gz_filename : str
            The output name of the compressed vcf file

        Returns
        -------
        vcftools_err : str
            vcftools log output

        Raises
        ------
        IOError
            If vcftools stderr returns an error
        IOError
            If bgzip stderr returns an error
    '''

    # vcftools subprocess call
    vcftools_call = subprocess.Popen(['vcftools'] + vcfname_arg + list(map(str, vcftools_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

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

def call_vcftools (vcfname_arg, vcftools_call_args):

    # vcftools subprocess call
    vcftools_call = subprocess.Popen(['vcftools'] + vcfname_arg + list(map(str, vcftools_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Wait for vcftools to finish
    vcftools_out, vcftools_err = vcftools_call.communicate()

    logging.info('vcftools call complete')

    return vcftools_err



def check_for_vcftools_output (vcftools_output):
    '''
        Checks for the previous vcftools output

        Confirms that neither a previous vcftools log or output file exists.

        Parameters
        ----------
        output_prefix : str
            Specifies the output filename to be created

        Raises
        ------
        IOError
            If the vcftools standard output exists
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

def produce_vcftools_log (output, filename):
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

        Returns
        -------
        output : file
            vcftools log file

    '''

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
    if os.path.basename(filename).split('.', 1)[-1] in ['vcf', 'vcf.gz', 'bcf']:

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

        if vcfname_format == 'nozip':
            return ['--vcf', filename]
        elif vcfname_format == 'gzip':
            return ['--gzvcf', filename]
        elif vcfname_format == 'bgzip':
            return ['--bcf', filename]
        else:
            logging.error('Unknown VCF file format')
            raise Exception('Unknown VCF file format')
