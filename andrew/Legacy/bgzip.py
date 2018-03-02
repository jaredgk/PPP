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
