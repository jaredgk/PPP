import os
import sys
import subprocess
import shutil
import argparse
import glob
import copy
import logging

# Call PPP-based scripts
sys.path.insert(0, os.path.abspath(os.path.join(os.pardir,'jared')))

import bcftools
from vcf_reader_func import checkFormat
from logging_module import initLogger, logArgs

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

    # Input VCF argument
    plink_parser.add_argument("--vcf", help = "Input VCF filename", type = str, action = parser_confirm_file())

    # Input PED arguments
    plink_parser.add_argument("--ped", help = "Input PED filename", type = str, action = parser_confirm_file())
    plink_parser.add_argument("--map", help = "Input MAP filename. Called alongside --ped", type = str, action = parser_confirm_file())

    # Input BED arguments
    plink_parser.add_argument("--bed", help = "Input BED filename", type = str, action = parser_confirm_file())
    plink_parser.add_argument("--fam", help = "Input FAM filename. Called alongside --bed", type = str, action = parser_confirm_file())
    plink_parser.add_argument("--bim", help = "Input BIM filename. Called alongside --bed", type = str, action = parser_confirm_file())

    plink_prefix = plink_parser.add_mutually_exclusive_group()
    plink_prefix.add_argument("--ped-prefix", help = "Input PED filename prefix", type = str)
    plink_prefix.add_argument("--bed-prefix", help = "Input BED filename prefix", type = str)

    output_formats = ['vcf', 'vcf.gz', 'bcf', 'ped', 'binary-ped']
    output_default = 'ped'
    plink_parser.add_argument('--out-format', metavar = metavar_list(output_formats), help = 'Specifies the output format', type = str, choices = output_formats, default = output_default)

    # Other basic arguments. Expand as needed
    plink_parser.add_argument('--out', help = 'Defines the output filename. Only usable with vcf-based output')
    plink_parser.add_argument('--out-prefix', help = 'Defines the output filename prefix', default = 'out')

    if passed_arguments:
        return plink_parser.parse_args(passed_arguments)
    else:
        return plink_parser.parse_args()

def check_plink_for_errors (plink_stderr):
    '''
        Checks the plink stderr for errors

        Parameters
        ----------
        plink_stderr : str
            plink stderr

        Raises
        ------
        IOError
            If plink stderr returns an error
    '''

    # Print output if error found. Build up as errors are discovered
    if plink_stderr:
        raise Exception(plink_stderr)

def call_plink (plink_call_args):
    '''
        Calls plink

        The function calls plink. Returns the stderr of plink to create a log
        file of the call.

        Parameters
        ----------
        plink_call_args : list
            plink arguments

        Raises
        ------
        Exception
            If plink stderr returns an error
    '''

    # vcftools subprocess call
    plink_call = subprocess.Popen(['plink'] + list(map(str, plink_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Wait for vcftools to finish
    plink_out, plink_err = plink_call.communicate()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        plink_out = plink_out.decode()
        plink_err = plink_err.decode()

    logging.info('plink call complete')

    # Check that the log file was created correctly
    check_plink_for_errors(plink_err)

def convert_vcfgz_to_bcf (output_prefix):

    '''
        Converts vcf.gz files to bcf files

        Parameters
        ----------
        output_prefix : str
            Specifies the filename prefix

        Raises
        ------
        IOError
            If the vcf.gz cannot be found
    '''

    # Assign the VCFGZ filename
    vcfgz_filename = output_prefix + '.vcf.gz'

    # Confirm the file exists, if not return error
    if not os.path.isfile(vcfgz_filename):
        raise IOError('Cannot locate vcf.gz intermediate file for bcf converion')

    # Convert the file
    bcftools.convert_to_bcf(vcfgz_filename, output_prefix)

    # Delete the intermediate VCFGZ file
    os.remove(vcfgz_filename)

def assign_plink_output_args (output_prefix, output_format):
    '''
        Assigns the output arguments

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

    # Check if the output is only a single file
    if output_format in ['vcf', 'vcf.gz', 'bcf']:

        if output_format == 'vcf':
            return ['--recode', 'vcf-iid', '--out', output_prefix]

        if output_format == 'vcf.gz':
            return ['--recode', 'vcf-iid', 'bgz', '--out', output_prefix]

        if output_format == 'bcf':
            return ['--recode', 'vcf-iid', 'bgz', '--out', output_prefix]

    # Check if the output if a ped or bed file-set
    elif output_format in ['ped', 'binary-ped']:

        if output_format == 'ped':
            return ['--recode', '--out', output_prefix]

        if output_format == 'binary-ped':
            return ['--make-bed', '--out', output_prefix]

    else:
        raise IOError('Unknown file format. This error should NOT show up')

def assign_plink_input_from_prefix (input_prefix, input_format):
    '''
        Assigns input arguments based on the specified files

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

    # Determine the input files from the prefix
    prefix_files = glob.glob(input_prefix + '*')

    # Check if files were identified, if not return error message
    if not prefix_files:
        raise IOError('Unable to assign input files Please confirm the prefix is specified correctly')

    # Check for ped input
    if input_format == 'ped':

        # Plink PED and Map files
        plink_ped = None
        plink_map = None

        # Check if expected files are found
        for plink_filename in prefix_files:

            # Check for PED and MAP files
            if '.ped' in plink_filename:
                plink_ped = plink_filename
            elif '.map' in plink_filename:
                plink_map = plink_filename
        # Check that input files included a ped file
        if plink_ped:

            # Check that ped and map files are both found
            if plink_ped and plink_map:

                # Return the ped-based files
                return ['--file', plink_ped[:-4]]

            else:
                raise IOError('Unable to assign map file. Please confirm the file is named correctly')

        # Return error message if no ped/bed input was found
        else:
            raise IOError('Unable to assign ped file. Please confirm the prefix and/or command is specified correctly')

    # Check for ped input
    elif input_format == 'binary-ped':

        # Plink BED, BIM, and FAM files
        plink_bed = None
        plink_bim = None
        plink_fam = None

        # Check if expected files are found
        for plink_filename in prefix_files:

            # Check for BED and associated files
            if '.bed' in plink_filename:
                plink_bed = plink_filename
            elif '.bim' in plink_filename:
                plink_bim = plink_filename
            elif '.fam' in plink_filename:
                plink_fam = plink_filename

        # Check that input files included a bed file
        if plink_bed:

            # Check that bed, bim, and fam files were found
            if plink_bed and plink_bim and plink_fam:

                # Return the bed-based files
                return ['--bfile', plink_bed[:-4]]

            else:
                raise IOError('Unable to assign bim and/or fam files. Please confirm the files are named correctly')

        # Return error message if no ped/bed input was found
        else:
            raise IOError('Unable to assign bed file. Please confirm the prefix and/or command is specified correctly')

def assign_plink_input_from_command_args (plink_args):

    input_command_args = []

    # Check if a ped file is assigned
    if plink_args.ped:

        # Check if map file is also assigned
        if plink_args.ped and plink_args.map:

            # Add the ped associated commands
            input_command_args.extend(['--ped', plink_args.ped, '--map', plink_args.map])
            # Return the input command
            return input_command_args

        else:
            raise IOError('Unable to assign ped file. Please confirm the file is named correctly')

    # Check if a bed file is assigned
    elif plink_args.bed:

        # Check if bim and fam file are also assigned
        if plink_bed and plink_bim and plink_fam:

            # Add the bed associated commands
            input_command_args.extend(['--bed', plink_args.bed, '--bim', plink_args.bim, '--fam', plink_args.fam])
            # Return the input command
            return input_command_args

        else:
            raise IOError('Unable to assign bim and/or fam files. Please confirm the files are named correctly')

    # Check if a vcf file is assigned
    elif plink_args.vcf:

        # Assign the vcf file format
        vcfname_format =  checkFormat(plink_args.vcf)

        input_command_args.extend(['--double-id', '--allow-extra-chr'])

        # Assign the associated input command, or return an error.
        if vcfname_format == 'vcf' or vcfname_format == 'bgzip':

            # Add the bed associated commands
            input_command_args.extend(['--vcf', plink_args.vcf])
            # Return the input command
            return input_command_args

        elif vcfname_format == 'bcf':

            # Add the bed associated commands
            input_command_args.extend(['--bcf', plink_args.vcf])
            # Return the input command
            return input_command_args

        else:
            raise Exception('Unknown VCF file format')

def run (passed_arguments = []):
    '''
        PLINK toolset automator

        Currently a basic PLINK automator for file conversion. Will be updated
        in the future if functionality needs to be expanded.

        Parameters
        ----------
        --vcf : str
            Specifies the input VCF filename
        --ped : str
            Specifies the input PED filename
        --map : str
            Specifies the input MAP filename. Called alongside --ped
        --bed : str
            Specifies the input BED filename
        --map : str
            Specifies the input FAM filename. Called alongside --bed
        --bim : str
            Specifies the input BIM filename. Called alongside --bed
        --ped-prefix : str
            Specifies the input PED filename prefix
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

    # List to hold the input argument(s)
    plink_input_args = []

    # Assign ped/bed input from prefix
    if plink_args.ped_prefix or plink_args.bed_prefix:

        if plink_args.ped_prefix:
            # Assign the ped based files
            plink_input_args = assign_plink_input_from_prefix(plink_args.ped_prefix, 'ped')

        elif plink_args.bed_prefix:
            # Assign the bed based files
            plink_input_args = assign_plink_input_from_prefix(plink_args.bed_prefix, 'binary-ped')

    else:
        # Assign the general input
        plink_input_args = assign_plink_input_from_command_args(plink_args)

    logging.info('Input parameters assigned')

    plink_output_args = assign_plink_output_args(plink_args.out_prefix, plink_args.out_format)

    logging.info('Output parameters assigned')

    # Call plink with the selected arguments
    call_plink(plink_input_args + plink_output_args)

    # Convert VCFGZ to BCF, as PLINK cannot directly create a BCF
    if plink_args.out_format == 'bcf':

        # Convert to BCF
        convert_vcfgz_to_bcf(plink_args.out_prefix)

if __name__ == "__main__":
    initLogger()
    run()
