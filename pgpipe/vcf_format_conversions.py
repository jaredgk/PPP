#!/usr/bin/env python
'''
    Automates various simple file conversions. Currently the function is 
    capable of converting between VCF-based formats (i.e. VCF, compressed-VCF, 
    and BCF) and PLINK-based formats (i.e. PED and Binary-PED). Additional 
    formats will be added as needed. 

    ############################
    Input Command-line Arguments
    ############################
    **--vcf** *<vcf_filename>*
        Argument used to define the filename of the VCF file.
    **--vcf-fid** *<fid_str>*
        Argument used to define the family ID for all VCF samples.
    **--ped-prefix** *<ped_prefix>*
        Argument used to define the filename prefix of both PED and MAP files.
    **--ped** *<ped_filename>*
        Argument used to define the filename of the PED file. Called alongside
        **--map**.
    **--map** *<map_filename>*
        Argument used to define the filename of the MAP file. Called alongside
        **--ped**.
    **--binary-ped-prefix** *<binary_ped_prefix>*
        Argument used to define the filename prefix of the Binary-PED, FAM,
        and BIM files.
    **--binary-ped** *<ped_filename>*
        Argument used to define the filename of the Binary-PED (i.e. BED) file. 
        Called alongside **--fam** and **--bim**.
    **--fam** *<fam_filename>*
        Argument used to define the filename of the FAM file. Called alongside
        **--binary-ped** and **--bim**.
    **--bim** *<bim_filename>*
        Argument used to define the filename of the BIM file. Called alongside
        **--binary-ped** and **--fam**.
   
    #############################
    Output Command-line Arguments
    #############################
    **--out** *<output_filename>*
        Argument used to define the complete output filename, overrides **--out-prefix**.
    **--out-prefix** *<output_prefix>*
        Argument used to define the output prefix (i.e. filename without file extension).
    **--out-format** *<vcf, vcf.gz, bcf, ped, ped-12, binary-ped, eigenstrat>*
        Argument used to define the desired output format. Formats include: uncompressed 
        VCF (vcf); compressed VCF (vcf.gz); BCF (bcf); PLINK text file (ped); PLINK "12" 
        coded text file (ped-12); binary PLINK file (binary-ped); and eigenstrat file
        (eigenstrat).
    **--overwrite**
        Argument used to define if previous output should be overwritten.

    ############################
    Other Command-line Arguments
    ############################
    **--delete-original**
        Argument used to define that the original file should be deleted once converted.
    **--threads** *<thread_int>*
        Argument used to define the number of threads. This argument is currently only 
        supported by conversions to/from PED and Binary-PED.
'''

import os
import sys
import subprocess
import shutil
import argparse
import glob
import logging

from pgpipe.eigenstrat_wrapper import *
from pgpipe.bcftools import convert_vcf, log_bcftools_reference
from pgpipe.plink import convert_ped, convert_bed, convert_vcf_to_plink, log_plink_reference
from pgpipe.vcf_reader_func import checkFormat
from pgpipe.logging_module import initLogger, logArgs
from pgpipe.misc import argprase_kwargs


def convert_argument_parser(passed_arguments = []):
    '''
    Convert Argument Parser

    Assign the parameters for the convert function using argparse.

    Parameters
    ----------
    passed_arguments : list, optional
        Parameters passed by another function. sys.argv is used if
        not given. 

    Raises
    ------
    IOError
        If the input, or other specified files do not exist
    '''

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

    # Create the parser
    convert_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    # Input VCF argument
    convert_parser.add_argument("--vcf", help = "Defines the filename of the VCF file", type = str, action = parser_confirm_file())

    # Sets a family ID for the samples in the VCF
    convert_parser.add_argument("--vcf-fid", dest = 'vcf_fid', help = "Defines the family ID for all VCF samples", type = str)

    # Input FASTA argument
    #convert_parser.add_argument("--fasta", help = "Input FASTA filename", type = str, action = parser_confirm_file())

    # Input IM argument
    #convert_parser.add_argument("--im", help = "Input IM filename", type = str, action = parser_confirm_file())

    # Input PED arguments
    convert_parser.add_argument("--ped", dest = 'ped_filename', help = "Defines the filename of the PED file. Called alongside --map", type = str, action = parser_confirm_file())
    convert_parser.add_argument("--map", dest = 'map_filename', help = "Defines the filename of the MAP file. Called alongside --ped", type = str, action = parser_confirm_file())

    # Input BED arguments
    convert_parser.add_argument("--binary-ped", dest = 'bed_filename', help = "Defines the filename of the Binary-PED (i.e. BED) file. Called alongside --fam and --bim", type = str, action = parser_confirm_file())
    convert_parser.add_argument("--fam", dest = 'fam_filename', help = "Defines the filename of the FAM file. Called alongside --binary-ped and --bim", type = str, action = parser_confirm_file())
    convert_parser.add_argument("--bim", dest = 'bim_filename', help = "Defines the filename of the BIM file. Called alongside --binary-ped and --fam", type = str, action = parser_confirm_file())

    # Input PED/BED prefix arguments
    convert_prefix = convert_parser.add_mutually_exclusive_group()
    convert_prefix.add_argument("--ped-prefix", help = "Defines the filename prefix of both PED and MAP files", type = str)
    convert_prefix.add_argument("--binary-ped-prefix", dest = 'bed_prefix', help = "Defines the filename prefix of the Binary-PED, FAM, and BIM files", type = str)

    # Output arguments.
    output_formats = ['vcf', 'vcf.gz', 'bcf', 'ped', 'ped-12', 'binary-ped', 'eigenstrat']
    convert_parser.add_argument('--out-format', metavar = metavar_list(output_formats), help = 'Defines the desired output format', type = str, choices = output_formats, required = True)
    convert_parser.add_argument('--out', help = 'Defines the complete output filename, overrides --out-prefix')
    convert_parser.add_argument('--out-prefix', help = 'Defines the output prefix (i.e. filename without file extension)', default = 'out')

    # Other basic arguments. Expand as needed
    convert_parser.add_argument('--delete-original', help = 'Defines that the original file should be deleted once converted', action='store_false')
    convert_parser.add_argument('--overwrite', help = "Defines if previous output should be overwritten", action = 'store_true')
    convert_parser.add_argument('--threads', help = "Set the number of threads. Only supported with ped", type = int, default = 1)

    if passed_arguments:
        return vars(convert_parser.parse_args(passed_arguments))
    else:
        return vars(convert_parser.parse_args())

def check_if_conversion (input_format, out_format):
    '''
        Checks that the input and output formats differ

        Parameters
        ----------
        input_format : str
            File format of the input

        out_format : str
            File format of the output

        Raises
        ------
        Exception
            If the input and output formats do not differ
    '''

    # Check if the input and output are the same formats
    if input_format == out_format:
        raise Exception('Input and Output formats are identical. Cannot convert')

def check_conversion_support (out_format, supported_formats):
    '''
        Confirms it is possible to convert the input into the specified output
        format

        Parameters
        ----------
        out_format : str
            File format of the output
        supported_formats : str
            The file formats the input supports (i.e. allows conversion)

        Raises
        ------
        Exception
            If the output format is not supported
    '''

    # Check that output format is supported
    if out_format not in supported_formats:
        raise Exception('Input cannot be converted into the format specified')

def run (**kwargs):
    '''
    File conversion suite

    This function uses the argparse-based function :py:func:`convert_argument_parser`
    to parse either sys.argv or passed_arguments to obtain the parameters below. 
    The parameters are then translated to their respective equivalent. Once all the 
    parameters are assigned, the conversion function of choice is called.

    Parameters
    ----------
    --vcf : str
        Filename of the VCF
    --ped : str
        Filename of the PED
    --map : str
        Filename of the MAP
    --bed : str
        Filename of the Binary-PED
    --map : str
        Filename of the FAM
    --bim : str
        Filename of the BIM
    --ped-prefix : str
        PED file prefix
    --binary-ped-prefix : str
        Binary-PED file prefix
    --out-prefix : str
        Output filename prefix
    --out-format : str
        Output format

    Raises
    ------
    Exception
        Impossible conversion given
    Exception
        Conversion to same format given
    '''

    # Update kwargs with defaults
    if __name__ != "__main__":
        kwargs = argprase_kwargs(kwargs, convert_argument_parser)

    # Assign arguments
    convert_args = argparse.Namespace(**kwargs)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(convert_args, func_name = 'convert')

    logging.info('Conversion parameters assigned')

    # Check if the format is binary-ped, i.e. .bed
    if convert_args.out_format == 'binary-ped':

        # Assign the expected output filename
        expected_out_filename = convert_args.out_prefix + '.bed'

    # Check if the format is ped-12, i.e. .ped
    elif convert_args.out_format == 'ped-12':

        # Assign the expected output filename
        expected_out_filename = convert_args.out_prefix + '.ped'
        
    else:

        # Assign the expected output filename
        expected_out_filename = '%s.%s' % (convert_args.out_prefix, convert_args.out_format)

    # Check if input is specified as a vcf
    if convert_args.vcf:

        # Determine the format of the vcf file (e.g. vcf, vcf.gz, bcf)
        vcf_format = checkFormat(convert_args.vcf)

        # Convert bgzip flag to vcf.gz. Consider changing checkFormat function?
        if vcf_format == 'bgzip':
            
            # Convert to the format used by convert_args.out_format
            vcf_format = 'vcf.gz'

        # Check that the input is in a supported file format
        if vcf_format not in ['vcf', 'vcf.gz', 'bcf']:
            raise Exception('Unknown VCF file format')

        # Check that the input and output formats differ
        check_if_conversion(vcf_format, convert_args.out_format)

        # List of formats that are supported from vcf
        supported_out_format = ['vcf', 'vcf.gz', 'bcf', 'ped', 'ped-12', 'binary-ped']

        # Check that the conversion is supported
        check_conversion_support(convert_args.out_format, supported_out_format)

        logging.info('Convertion method assigned')

        # Check if the output is a vcf
        if convert_args.out_format in ['vcf', 'vcf.gz', 'bcf']:

            convert_vcf(convert_args.vcf, convert_args.out_prefix, convert_args.out_format, convert_args.overwrite, keep_original = True)

            # Add the reference to the log
            log_bcftools_reference()

            # Check if a specific output filename was specified
            if convert_args.out:

                # Rename the vcf file
                shutil.move(expected_out_filename, convert_args.out)

                # Rename the log file
                shutil.move(expected_out_filename + '.log', convert_args.out + '.log')

        # Check if the output is a ped or binary-ped
        elif convert_args.out_format in ['ped', 'ped-12', 'binary-ped']:

            # Call the VCF to PED/BED converter
            convert_vcf_to_plink(convert_args.vcf, convert_args.vcf_fid, convert_args.out_prefix, convert_args.out_format, threads = convert_args.threads, overwrite = convert_args.overwrite)

            # Rename the plink reference
            shutil.move(convert_args.out_prefix + '.log', expected_out_filename + '.log')

            # Add the reference to the log
            log_plink_reference()

    # Check if input is a specified as a ped
    elif convert_args.ped_prefix or (convert_args.ped_filename or convert_args.map_filename):

        # Check that the input and output formats differ
        check_if_conversion('ped', convert_args.out_format)

        # List of formats that are supported from ped
        supported_out_format = ['vcf', 'vcf.gz', 'bcf', 'ped-12', 'binary-ped', 'eigenstrat']

        # Check that the conversion is supported
        check_conversion_support(convert_args.out_format, supported_out_format)

        logging.info('Convertion method assigned')

        # Check if output format is not eigenstrat
        if convert_args.out_format != 'eigenstrat':

            # Convert the ped file using plink
            convert_ped(**vars(convert_args))

            # Rename the plink reference
            shutil.move(convert_args.out_prefix + '.log', expected_out_filename + '.log')

            # Add the reference to the log
            log_plink_reference()

        # Check if output format is eigenstrat
        else:

            # Convert the ped file using admixtools
            ped_to_eigenstrat(**vars(convert_args))

    # Check if input is a specified as a bed
    elif convert_args.bed_prefix or (convert_args.bed_filename or convert_args.fam_filename or convert_args.bim_filename):

        # Check that the input and output formats differ
        check_if_conversion('binary-ped', convert_args.out_format)

        # List of formats that are supported from bed
        supported_out_format = ['vcf', 'vcf.gz', 'bcf', 'ped', 'ped-12', 'eigenstrat']

        # Check that the conversion is supported
        check_conversion_support(convert_args.out_format, supported_out_format)

        logging.info('Convertion method assigned')

        # Check if output format is not eigenstrat
        if convert_args.out_format != 'eigenstrat':

            # Convert the ped file using plink
            convert_bed(**vars(convert_args))

            # Rename the plink reference
            shutil.move(convert_args.out_prefix + '.log', expected_out_filename + '.log')

            # Add the reference to the log
            log_plink_reference()

        # Check if output format is eigenstrat
        else:

            # Convert the ped file using admixtools
            bed_to_eigenstrat(**vars(convert_args))
            

if __name__ == "__main__":
    initLogger()
    run(**convert_argument_parser())
