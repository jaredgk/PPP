import os
import sys
import logging
import subprocess

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

def call_vcftools (vcftools_call_args, output_format, vcftools_output_filename, user_output_filename):
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
            The output filename assigned by vcftools
        user_output_filename : str
            The output filename specified by the user (if specified)

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
    elif output_format == 'bed':
        # Pipe vcftools stdout to bed file
        vcftools_err = pipe_vcftools_to_bed_file(vcftools_call_args, vcftools_output_filename)
    else:
        # Call vcftools under standard conditions
        vcftools_err = standard_vcftools_call(vcftools_call_args)

    # Return the log
    return vcftools_err

# Test arguments
input_call = ['--gzvcf', sys.argv[1], '--max-missing',  '0.5', '--recode']

#pipe_vcftools_bgzip(input_call, 'delete.vcf.gz')

input_call2 = ['--gzvcf', sys.argv[1], '--max-missing',  '0.5', '--kept-sites']

call_vcftools(input_call2, 'bed', 'delete.bed', '')
