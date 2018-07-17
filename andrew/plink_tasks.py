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

from vcf_reader_func import checkFormat
from logging_module import initLogger, logArgs
from plink import call_plink

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
    plink_parser.add_argument("--vcf", help = "Input VCF filename", type = str, action = parser_confirm_file())

    # Input HAPS arguments
    plink_parser.add_argument("--haps", help = "Input HAPS filename", type = str, action = parser_confirm_file())
    plink_parser.add_argument("--sample", help = "Input SAMPLE filename. Called alongside --haps", type = str, action = parser_confirm_file())
    plink_prefix.add_argument("--haps-prefix", help = "Input HAPS filename prefix", type = str)

    # Input PED arguments
    plink_parser.add_argument("--ped", help = "Input PED filename", type = str, action = parser_confirm_file())
    plink_parser.add_argument("--map", help = "Input MAP filename. Called alongside --ped", type = str, action = parser_confirm_file())
    plink_prefix.add_argument("--ped-prefix", help = "Input PED filename prefix", type = str)
    # Input BED arguments
    plink_parser.add_argument("--bed", help = "Input BED filename", type = str, action = parser_confirm_file())
    plink_parser.add_argument("--fam", help = "Input FAM filename. Called alongside --bed", type = str, action = parser_confirm_file())
    plink_parser.add_argument("--bim", help = "Input BIM filename. Called alongside --bed", type = str, action = parser_confirm_file())
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
        if not plink_ped:
            # Return error message if no ped input was found
            raise IOError('Unable to assign ped file. Please confirm the prefix and/or command is specified correctly')

        if not plink_map:
            # Return error message if no map input was found
            raise IOError('Unable to assign map file. Please confirm the file is named correctly')

        # Return the ped-based files
        return ['--file', plink_ped[:-4]]

    # Check for binary-ped input
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
        if not plink_bed:
            # Return error message if no ped/bed input was found
            raise IOError('Unable to assign bed file. Please confirm the prefix and/or command is specified correctly')

        # Check that input files included a bed file
        if not plink_bim:
            # Return error message if no bim input was found
            raise IOError('Unable to assign bim file. Please confirm the file is named correctly')

        # Check that input files included a bed file
        if not plink_fam:
            # Return error message if no fam input was found
            raise IOError('Unable to assign fam file. Please confirm the file is named correctly')

        # Return the binary-ped files
        return ['--bfile', plink_bed[:-4]]

    # Check for haps input
    elif input_format == 'haps':

        # Plink haps and sample files
        plink_haps = None
        plink_sample = None

        # Check if expected files are found
        for plink_filename in prefix_files:

            # Check for haps and associated files
            if '.haps' in plink_filename:
                plink_haps = plink_filename
            elif '.sample' in plink_filename:
                plink_sample = plink_filename

        # Check that input files included a haps file
        if not plink_haps:
            # Return error message if no haps input was found
            raise IOError('Unable to assign haps file. Please confirm the prefix and/or command is specified correctly')

        # Check that input files included a sample file
        if not plink_sample:
            # Return error message if no sample input was found
            raise IOError('Unable to assign sample file. Please confirm the file is named correctly')

        # Return the haps-based files
        return ['--haps', plink_haps, 'ref-first', '--sample', plink_sample]

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
        if plink_args.bed and plink_args.bim and plink_args.fam:

            # Add the bed associated commands
            input_command_args.extend(['--bed', plink_args.bed, '--bim', plink_args.bim, '--fam', plink_args.fam])
            # Return the input command
            return input_command_args

        else:
            raise IOError('Unable to assign bim and/or fam files. Please confirm the files are named correctly')

    # Check if a haps file is assigned
    elif plink_args.haps:

        # Check if bim and fam file are also assigned
        if plink_args.haps and plink_args.sample:

            # Add the bed associated commands
            input_command_args.extend(['--haps', plink_args.haps, 'ref-first', '--sample', plink_args.sample])
            # Return the input command
            return input_command_args

        else:
            raise IOError('Unable to assign sample file. Please confirm the file is named correctly')

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

    raise Exception('No input specified')

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
        plink_input_args = assign_plink_input_from_command_args(plink_args)

    logging.info('Input parameters assigned')

    plink_output_args = assign_plink_output_args(plink_args.out_prefix, plink_args.out_format)

    logging.info('Output parameters assigned')

    # Call plink with the selected arguments
    call_plink(plink_input_args + plink_output_args, plink_args.out_prefix, plink_args.out_format)

    # Rename output to plink_args.out, if specified
    if plink_args.out:
        shutil.move('%s.%s' % (plink_args.out_prefix, plink_args.out_format), plink_args.out)
        shutil.move(plink_args.out_prefix + '.log', plink_args.out + '.log')
    # Rename log using plink_args.out_prefix
    else:
        shutil.move(plink_args.out_prefix + '.log', '%s.%s.log' % (plink_args.out_prefix, plink_args.out_format))

    logging.info('plink log file created')

if __name__ == "__main__":
    initLogger()
    run()
