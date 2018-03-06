import os
import sys
import subprocess
import shutil
import argparse
import glob
import copy
import logging

# Call PPP-based scripts
#sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'jared')))

#from logging_module import initLogger

def plink_argument_parser(passed_arguments):
    '''Phase Argument Parser - Assigns arguments for vcftools from command line.
    Depending on the argument in question, a default value may be specified'''

    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError # File not found
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
    plink_prefix.add_argument("--ped-prefix", help = "Input PED/BED filename prefix", type = str)
    plink_prefix.add_argument("--bed-prefix", help = "Input PED/BED filename prefix", type = str)

    output_formats = ['vcf', 'vcf.gz', 'ped', 'binary-ped']
    output_default = 'ped'
    plink_parser.add_argument('--out-format', metavar = metavar_list(output_formats), help = 'Specifies the output format', type = str, choices = output_formats, default = output_default)

    # Other basic arguments. Expand as needed
    plink_parser.add_argument('--out', help = 'Defines the output filename')
    plink_parser.add_argument('--out-prefix', help = 'Defines the output filename prefix', default = 'out')

    if passed_arguments:
        return plink_parser.parse_args(passed_arguments)
    else:
        return plink_parser.parse_args()

def logArgs(args, pipeline_function):
    '''Logs arguments from argparse system. May be replaced with a general function in logging_module'''
    logging.info('Arguments for %s:' % pipeline_function)
    for k in vars(args):
        logging.info('Arguments %s: %s' % (k, vars(args)[k]))

def assign_plink_output_args (output_format, output_prefix):
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
    if output_format in ['vcf', 'vcf.gz']:

        if output_format == 'vcf':
            return ['--recode', 'vcf', '--out', output_prefix]

        if output_format == 'vcf.gz':
            return ['--recode', 'vcf', 'bgz', '--out', output_prefix]

    elif output_format in ['ped', 'binary-ped']:

        if output_format == 'ped':
            return ['--recode', '--out', output_prefix]

        if output_format == 'binary-ped':
            return ['--make-bed', '--out', output_prefix]

    else:
        raise IOError('Unknown file format. This error should NOT show up')


def assign_plink_input_from_prefix (input_prefix):
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

    # Plink PED and Map files
    plink_ped = None
    plink_map = None

    # Plink BED, BIM, and FAM files
    plink_bed = None
    plink_bim = None
    plink_fam = None

    # Check if expected files are found
    for plink_filename in prefix_files:

        # Check for PED and MAP files
        if '.ped' in plink_filename:
            plink_ped = plink_filename
        elif '.map' in plink_filename:
            plink_map = plink_filename

        # Check for BED and associated files
        elif '.bed' in plink_filename:
            plink_bed = plink_filename
        elif '.bim' in plink_filename:
            plink_bim = plink_filename
        elif '.fam' in plink_filename:
            plink_fam = plink_filename

    # Check that input files included a ped file
    if plink_ped:

        # Check that ped and map files are both found
        if plink_ped and plink_map:

            # Return the ped-based files
            return ['--file', plink_ped[:-4]]

        else:
            raise IOError('Unable to assign map file. Please confirm the file is named correctly')

    else:
        raise IOError('Unable to assign ped file. Please confirm the file is named correctly')

    # Check that input files included a bed file
    if plink_bed:

        # Check that bed, bim, and fam files were found
        elif plink_bed and plink_bim and plink_fam:

            # Return the bed-based files
            return ['--bfile', plink_bed[:-4]]

        else:
            raise IOError('Unable to assign bim and/or fam files. Please confirm the files are named correctly')

    else:
        raise IOError('Unable to assign bed file. Please confirm the file is named correctly')


def assign_plink_input_from_command_args (plink_args):

    # Check if a ped file is assigned
    if plink_args.ped:
        # Check if map file is also assigned
        if plink_args.ped and plink_args.map:
            # Return the input command
            return ['--ped', plink_args.ped, '--map', plink_args.map]

        else:
            raise IOError('Unable to assign ped file. Please confirm the file is named correctly')


    # Check if a bed file is assigned
    elif plink_args.bed:
        # Check if bim and fam file are also assigned
        if plink_bed and plink_bim and plink_fam:
            # Return the input command
            return ['--bed', plink_args.bed, '--bim', plink_args.bim, '--fam', plink_args.fam]

        else:
            raise IOError('Unable to assign bim and/or fam files. Please confirm the files are named correctly')

    # Check if a vcf file is assigned
    elif plink_args.vcf:
        pass



def run (passed_arguments = []):
    '''
        Summary Line

        Complete Summary

        Parameters
        ----------

        Returns
        -------

        Raises
        ------

    '''

    # Grab plink arguments from command line
    plink_args = plink_argument_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(plink_args, 'plink_tasks')

    # List to hold the input argument(s)
    plink_input_args = []

    # Assign ped input from prefix
    if plink_args.ped_prefix:
        # Assign the ped-based files
        plink_input_args = assign_plink_input_from_prefix(plink_args.ped_prefix)

    # Assign bed input from prefix
    elif plink_args.bed_prefix:
        # Assign the ped-based files
        plink_input_args = assign_plink_input_from_prefix(plink_args.bed_prefix)

    else:
        # Assign the general input
        plink_input_args = assign_plink_input_from_command_args(plink_args)

    plink_output_args = assign_plink_output_args(plink_args.out_format, plink_args.out_prefix)

    print ' '.join(map(str, plink_input_args + plink_output_args))

if __name__ == "__main__":
    #initLogger()
    run()
