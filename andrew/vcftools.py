import os
import sys
import logging

sys.path.insert(0, os.path.abspath(os.path.join(os.pardir,'jared')))

import vcf_reader_func

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
    if filename.split('.', 1)[-1] in ['vcf', 'vcf.gz', 'bcf']:

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
