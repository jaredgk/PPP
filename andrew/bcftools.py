import os
import sys
import logging
import subprocess

sys.path.insert(0, os.path.abspath(os.path.join(os.pardir,'jared')))

import vcf_reader_func

def call_bcftools (bcftools_call_args):

    # bcftools subprocess call
    bcftools_call = subprocess.Popen(['bcftools'] + list(map(str, bcftools_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Wait for bcftools to finish
    bcftools_out, bcftools_err = bcftools_call.communicate()

    logging.info('bcftools call complete')

    return bcftools_err

def convert_to_bcf (filename):
    pass

def convert_from_bcf (filename, format):
    pass
