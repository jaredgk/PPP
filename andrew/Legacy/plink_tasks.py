import os
import sys
import subprocess
import shutil
import argparse
import glob
import copy
import logging

# Call PPP-based scripts
sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'pppipe')))

from logging_module import initLogger, logArgs
from plink import *

def plink_argument_parser(passed_arguments):
    '''Phase Argument Parser - Assigns arguments for vcftools from command line.
    Depending on the argument in question, a default value may be specified'''

    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError('%s not found' % value)
                setattr(args, self.dest, value)
        return customAction

    def metavar_list (var_list):
        '''Create a formmated metavar list for the help output'''
        return '{' + ', '.join(var_list) + '}'

    plink_parser = argparse.ArgumentParser()

    plink_prefix = plink_parser.add_mutually_exclusive_group()

    # Input VCF argument
    plink_parser.add_argument("--vcf", dest = 'vcf_filename', help = "Input VCF filename", type = str, action = parser_confirm_file())

    # Sets a family ID for the samples in the VCF
    plink_parser.add_argument("--vcf-fid", dest = 'vcf_fid', help = "Specifies the family ID for all samples", type = str)

    # Input HAPS arguments
    plink_parser.add_argument("--haps", dest = 'haps_filename', help = "Input HAPS filename", type = str, action = parser_confirm_file())
    plink_parser.add_argument("--sample", dest = 'sample_filename', help = "Input SAMPLE filename. Called alongside --haps", type = str, action = parser_confirm_file())
    plink_prefix.add_argument("--haps-prefix", help = "Input HAPS filename prefix", type = str)

    # Input PED arguments
    plink_parser.add_argument("--ped", dest = 'ped_filename', help = "Input PED filename", type = str, action = parser_confirm_file())
    plink_parser.add_argument("--map", dest = 'map_filename', help = "Input MAP filename. Called alongside --ped", type = str, action = parser_confirm_file())
    plink_prefix.add_argument("--ped-prefix", help = "Input PED filename prefix", type = str)

    # Input BED arguments
    plink_parser.add_argument("--binary-ped", dest = 'bed_filename', help = "Input BED filename", type = str, action = parser_confirm_file())
    plink_parser.add_argument("--fam", dest = 'fam_filename', help = "Input FAM filename. Called alongside --bed", type = str, action = parser_confirm_file())
    plink_parser.add_argument("--bim", dest = 'bim_filename', help = "Input BIM filename. Called alongside --bed", type = str, action = parser_confirm_file())
    plink_prefix.add_argument("--bed-prefix", help = "Input BED filename prefix", type = str)

    output_formats = ['vcf', 'vcf.gz', 'bcf', 'ped', 'binary-ped']
    output_default = 'ped'
    plink_parser.add_argument('--out-format', metavar = metavar_list(output_formats), help = 'Specifies the output format', type = str, choices = output_formats, default = output_default)

    # Other basic arguments. Expand as needed
    plink_parser.add_argument('--out', help = 'Defines the output filename. Only usable with vcf-based output')
    plink_parser.add_argument('--out-prefix', help = 'Defines the output filename prefix', default = 'out')

    # General arguments.
    plink_parser.add_argument('--overwrite', help = "Overwrite previous output files", action = 'store_true')
    plink_parser.add_argument('--threads', help = "Set the number of threads", type = int, default = 1)

    if passed_arguments:
        return plink_parser.parse_args(passed_arguments)
    else:
        return plink_parser.parse_args()

def run (passed_arguments = []):
    '''
        PLINK toolset automator

        Currently a basic PLINK automator for file conversion. Will be updated
        in the future if functionality needs to be expanded.

        Parameters
        ----------
        --vcf : str
            Specifies the input VCF filename
        --vcf-fid : str
            Specifies the family ID for all samples
        --ped : str
            Specifies the input PED filename
        --map : str
            Specifies the input MAP filename. Called alongside --ped
        --ped-prefix : str
            Specifies the input PED filename prefix
        --bed : str
            Specifies the input BED filename
        --map : str
            Specifies the input FAM filename. Called alongside --bed
        --bim : str
            Specifies the input BIM filename. Called alongside --bed
        --bed-prefix : str
            Specifies the input BED filename prefix
        --out : str
            Specifies the output filename. Only usable with vcf-based output
        --out-prefix : str
            Specifies the output filename prefix
        --out-format : str
            Specifies the output format {vcf, vcf.gz, bcf, ped, binary-ped} (Default: ped)

        Returns
        -------
        output : file
            Statistic file output
        log : file
            Log file output

        Raises
        ------
        IOError
            Input file does not exist

    '''

    # Grab plink arguments from command line
    plink_args = plink_argument_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(plink_args, func_name = 'plink_tasks')

    # List to hold the general argument(s)
    plink_general_args = []

    # Assign the number of threads
    plink_general_args.extend(['--threads', plink_args.threads])

    # List to hold the input argument(s)
    plink_input_args = []

    if plink_args.ped_prefix:

        # Assign the ped based files
        plink_input_args = assign_plink_input_from_prefix(plink_args.ped_prefix, 'ped')

    elif plink_args.bed_prefix:

        # Assign the bed based files
        plink_input_args = assign_plink_input_from_prefix(plink_args.bed_prefix, 'binary-ped')

    elif plink_args.haps_prefix:

        # Assign the bed based files
        plink_input_args = assign_plink_input_from_prefix(plink_args.haps_prefix, 'haps')

    else:
        # Assign the general input
        plink_input_args = assign_plink_input_from_command_args(**vars(plink_args))

    logging.info('Input parameters assigned')

    plink_output_args = assign_plink_output_args(plink_args.out_prefix, plink_args.out_format)

    # Assign the expect output filename
    if plink_args.out_format == 'binary-ped':
        plink_expected_output = plink_args.out_prefix + '.bed'
    else:
        plink_expected_output = '%s.%s' % (plink_args.out_prefix, plink_args.out_format)

    # Check if previous output should not be overwritten
    if not plink_args.overwrite:

        # Check if previous output exists
        if os.path.exists(plink_expected_output):
            raise Exception('%s already exists. Add --overwrite to ignore'  % plink_expected_output)

    logging.info('Output parameters assigned')

    # Call plink with the selected arguments
    call_plink(plink_input_args + plink_output_args + plink_general_args, plink_args.out_prefix, plink_args.out_format)

    # Add reference to log file
    log_plink_reference(plink_args.out_prefix + '.log')

    # Rename output to plink_args.out, if specified
    if plink_args.out:
        shutil.move(plink_expected_output, plink_args.out)
        shutil.move(plink_args.out_prefix + '.log', plink_args.out + '.log')

    # Rename log using plink_args.out_prefix
    else:
        shutil.move(plink_args.out_prefix + '.log', '%s.log' % plink_expected_output)

    logging.info('plink log file created')

if __name__ == "__main__":
    initLogger()
    run()
