import os
import sys
import logging
import subprocess

sys.path.insert(0, os.path.abspath(os.path.join(os.pardir,'jared')))

from vcf_reader_func import checkFormat

def call_bcftools (bcftools_call_args):

    # bcftools subprocess call
    bcftools_call = subprocess.Popen(['bcftools'] + list(map(str, bcftools_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Wait for bcftools to finish
    bcftools_out, bcftools_err = bcftools_call.communicate()

    logging.info('bcftools call complete')

    return bcftools_out, bcftools_err

def check_for_index (filename):

    # Assign the file format
    file_format = checkFormat(filename)

    # Check if the file to be indexed is a vcf.gz
    if file_format == 'bgzip':
        # Check if the index (.tbi) exists
        if os.path.isfile(filename + '.tbi'):
            return True

    # Check if the file to be indexed is a bcf
    elif file_format == 'bcf':
        # Check if the index (.csi) exists
        if os.path.isfile(filename + '.csi'):
            return True

    # Return false if no index is found
    return False

def create_index (filename):

    # Assign the file format
    file_format = checkFormat(filename)

    # Check if the file to be indexed is a vcf.gz
    if file_format == 'bgzip':
        # Create a index (.tbi)
        index_out, index_err = call_bcftools(['index', '-t', filename])

    # Check if the file to be indexed is a bcf
    elif file_format == 'bcf':
        # Create a index (.csi)
        index_out, index_err = call_bcftools(['index', '-c', filename])

    # Report if file cannot be indexed
    else:
        raise Exception('Error creating index for: %s. Only .bcf and .vcf.gz (bgzip) files are supported.' % filename)

    # Check bcftools for errors
    if index_out or index_err:
        print (index_out, index_err)

def convert_to_bcf (filename):
    pass

def convert_from_bcf (filename, format):
    pass
