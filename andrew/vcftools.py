import os
import sys
import logging
import subprocess

sys.path.insert(0, os.path.abspath(os.path.join(os.pardir,'jared')))

from vcf_reader_func import checkFormat

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

def bgzip_decompress_vcfgz (vcfgz_filename, out_prefix = '', keep_original = False):
        '''
            Converts a vcf.gz to vcf

            The function automates bgzip to decompress a vcf.gz file into a vcf

            Parameters
            ----------
            vcfgz_filename : str
                The file name of the vcf.gz file to be decompressed
            keep_original : bool
                Specifies if the original file should be kept

            Raises
            ------
            IOError
                Error in creating the compressed file
        '''

        # Compress and keep the original file
        if keep_original or out_prefix:

            if out_prefix:

                # Assign the bgzip filename
                vcf_filename = split_filename.split(os.extsep)[0] + '.vcf'

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

        # Compress and do not keep the original file
        else:

            # bgzip subprocess call
            bgzip_call = subprocess.Popen(['bgzip', '-d', vcfgz_filename], stdout = subprocess.PIPE, stderr = subprocess.PIPE)

        # Save the stdout and stderr from bgzip
        bgzip_out, bgzip_err = bgzip_call.communicate()

        # Check that output file was compressed correctly
        check_bgzip_for_errors(bgzip_err)

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
            bgzip_call = subprocess.Popen(['bgzip', '-c', vcf_filename], stdout = vcfgz_file, stderr = subprocess.PIPE)

        # Compress and do not keep the original file
        else:

            # bgzip subprocess call
            bgzip_call = subprocess.Popen(['bgzip', vcf_filename], stdout = subprocess.PIPE, stderr = subprocess.PIPE)

        # Save the stdout and stderr from bgzip
        bgzip_out, bgzip_err = bgzip_call.communicate()

        # Check that output file was compressed correctly
        check_bgzip_for_errors(bgzip_err)

def cvt_vcftools_site_to_bed (vcftools_out_str):
    # Check if str in the header
    if 'CHROM' not in vcftools_out_str or 'POS' not in vcftools_out_str:
        # Split the line into a list
        vcftools_out_data = vcftools_out_str.strip().split('\t')
        # Convert the chromStart to int
        vcftools_out_data[1] = int(vcftools_out_data[1])
        # Calc chromEnd
        chrom_end = vcftools_out_data[1] + 1
        # Add chrom_end to the list
        vcftools_out_data = vcftools_out_data + [chrom_end]
        # Return the list as a string (with newline element)
        return '\t'.join(map(str, vcftools_out_data)) + '\n'
    else:
        # Remove the header
        return ''

def pipe_vcftools (vcftools_call_args):
    '''
        Calls vcftools with pipe output

        The output of this function is the stdout and stderr of vcftools. This
        function should only be used if vcftools is being used as the stdin of
        another function. Please note that this function does not check the for
        errors in the vcftools call. Please check for errors after the call is
        closed using check_vcftools_for_errors.

        Parameters
        ----------
        vcftools_call_args : list
            vcftools arguments

        Returns
        -------
        vcftools_call : subprocess.Popen
            vcftools subprocess call
        vcftools_call.stdout : PIPE
            vcftools stdout PIPE (Results)
        vcftools_call.stderr : PIPE
            vcftools stderr PIPE (Log)

    '''

    # vcftools subprocess call
    vcftools_call = subprocess.Popen(['vcftools', '--stdout'] + list(map(str, vcftools_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    return vcftools_call

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
    # Open vcftools pipe
    vcftools_call = pipe_vcftools(vcftools_call_args)

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

    # Close the vcftools stdout
    vcftools_call.stdout.close()

    # Wait for vctools to finish
    vcftools_call.wait()

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

    vcftools_call = pipe_vcftools(vcftools_call_args)

    # Create bgzip output file
    bgzip_output = open(output_filename, 'wb')

    # bgzip subprocess call
    bgzip_call = subprocess.Popen(['bgzip'], stdin = vcftools_call.stdout, stdout = bgzip_output, stderr = subprocess.PIPE)

    # Close the vcftools stdout
    vcftools_call.stdout.close()

    # Wait for vctools to finish
    vcftools_call.wait()

    # Read the vcftools stderr
    vcftools_stderr = vcftools_call.stderr.read()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        vcftools_stderr = vcftools_stderr.decode()

    # Check that the log file was created correctly
    check_vcftools_for_errors(vcftools_stderr)

    # Close the compressed vcf file
    bgzip_output.close()

    # Wait for bgzip to finish
    bgzip_call.wait()

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

    # vcftools subprocess call without stdout
    vcftools_call = subprocess.Popen(['vcftools'] + list(map(str, vcftools_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

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

def call_vcftools (vcftools_call_args, output_format, vcftools_output_filename):
    '''
        Calls vcftools

        The function calls vcftools. Returns the stderr of vcftools to
        create log file of the call.

        Parameters
        ----------
        vcftools_call_args : list
            vcftools arguments
        output_format : str
            The output format
        vcftools_output_filename : str
            The output filename assigned by vcftools (for piped calls)

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

    # Check if the output is a bgzipped vcf
    if output_format == 'vcf.gz':
        # Pipe vcftools stdout to bgzip to create a bgzipped vcf
        vcftools_err = pipe_vcftools_bgzip(vcftools_call_args, vcftools_output_filename)
    elif output_format == 'removed_bed' or output_format == 'kept_bed':
        # Pipe vcftools stdout to bed file
        vcftools_err = pipe_vcftools_to_bed_file(vcftools_call_args, vcftools_output_filename)
    else:
        # Call vcftools under standard conditions
        vcftools_err = standard_vcftools_call(vcftools_call_args)

    # Return the log
    return vcftools_err

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
            raise IOError('VCF output file cannot be deleted')

    logging.info('Output file assigned')

    # Check if log file already exists
    if os.path.isfile(vcftools_output + '.log'):
        try:
            # Delete the output
            os.remove(vcftools_output + '.log')
        except:
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
