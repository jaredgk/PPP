import os
import sys
import subprocess
import shutil
import argparse
import glob
import logging

sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'pppipe')))

# Import PPP modules and scripts
from bcftools import convert_vcf, log_bcftools_reference
from plink import convert_ped, convert_bed, convert_vcf_to_plink, log_plink_reference
from vcf_reader_func import checkFormat
from logging_module import initLogger, logArgs

def convert_argument_parser(passed_arguments):
    '''Convert Argument Parser - Assigns arguments from command line.
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

    convert_parser = argparse.ArgumentParser()

    # Input VCF argument
    convert_parser.add_argument("--vcf", help = "Input VCF filename", type = str, action = parser_confirm_file())

    # Sets a family ID for the samples in the VCF
    convert_parser.add_argument("--vcf-fid", dest = 'vcf_fid', help = "Specifies the family ID for all samples", type = str)

    # Input FASTA argument
    convert_parser.add_argument("--fasta", help = "Input FASTA filename", type = str, action = parser_confirm_file())

    # Input IM argument
    convert_parser.add_argument("--im", help = "Input IM filename", type = str, action = parser_confirm_file())

    # Input PED arguments
    convert_parser.add_argument("--ped", dest = 'ped_filename', help = "Input PED filename", type = str, action = parser_confirm_file())
    convert_parser.add_argument("--map", dest = 'map_filename', help = "Input MAP filename. Composite file of --ped", type = str, action = parser_confirm_file())

    # Input BED arguments
    convert_parser.add_argument("--binary-ped", dest = 'bed_filename', help = "Input Binary-PED (i.e. BED) filename", type = str, action = parser_confirm_file())
    convert_parser.add_argument("--fam", dest = 'fam_filename', help = "Input FAM filename. Composite file of --binary-ped", type = str, action = parser_confirm_file())
    convert_parser.add_argument("--bim", dest = 'bim_filename', help = "Input BIM filename. Composite file of --binary-ped", type = str, action = parser_confirm_file())

    # Input PED/BED prefix arguments
    convert_prefix = convert_parser.add_mutually_exclusive_group()
    convert_prefix.add_argument("--ped-prefix", help = "Input PED filename prefix", type = str)
    convert_prefix.add_argument("--binary-ped-prefix", dest = 'bed_prefix', help = "Input Binary-PED (i.e. BED) filename prefix", type = str)

    # Output arguments.
    output_formats = ['vcf', 'vcf.gz', 'bcf', 'fasta', 'im', 'ped', 'binary-ped']
    convert_parser.add_argument('--out-format', metavar = metavar_list(output_formats), help = 'Specifies the output format', type = str, choices = output_formats, required = True)
    convert_parser.add_argument('--out', help = 'Defines the output filename. Cannot be used with composite formats')
    convert_parser.add_argument('--out-prefix', help = 'Defines the output filename prefix', default = 'out')

    # Other basic arguments. Expand as needed
    convert_parser.add_argument('--delete-original', help = 'Defines that the original file should be deleted once converted', action='store_false')
    convert_parser.add_argument('--overwrite', help = "Overwrite previous output files", action = 'store_true')
    convert_parser.add_argument('--threads', help = "Set the number of threads. Only supported with ped", type = int, default = 1)

    if passed_arguments:
        return convert_parser.parse_args(passed_arguments)
    else:
        return convert_parser.parse_args()

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

def run (passed_arguments = []):
    '''
        File conversion suite

        Allows for various file converions, using the most efficient methods.

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
        --out-prefix : str
            Specifies the output filename prefix
        --out-format : str
            Specifies the output format {vcf, vcf.gz, bcf, ped, binary-ped}

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
    convert_args = convert_argument_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(convert_args, func_name = 'convert')

    logging.info('Conversion parameters assigned')

    # Check if the format is binary-ped, i.e. .bed
    if convert_args.out_format == 'binary-ped':

        # Assign the expected output filename
        expected_out_filename = convert_args.out_prefix + '.bed'
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
        supported_out_format = ['vcf', 'vcf.gz', 'bcf', 'ped', 'binary-ped']

        # Check that the conversion is supported
        check_conversion_support(convert_args.out_format, supported_out_format)

        logging.info('Convertion method assigned')

        # Check if the output is a vcf
        if convert_args.out_format in ['vcf', 'vcf.gz', 'bcf']:

            convert_vcf(convert_args.vcf, convert_args.out_prefix, convert_args.out_format, convert_args.overwrite, keep_original = True)

            # Add the reference to the log
            log_bcftools_reference(expected_out_filename, append_mode = False)

            # Check if a specific output filename was specified
            if convert_args.out:

                # Rename the vcf file
                shutil.move(expected_out_filename, convert_args.out)

                # Rename the log file
                shutil.move(expected_out_filename + '.log', convert_args.out + '.log')

        # Check if the output is a ped or binary-ped
        elif convert_args.out_format in ['ped', 'binary-ped']:

            # Call the VCF to PED/BED converter
            convert_vcf_to_plink(convert_args.vcf, convert_args.vcf_fid, convert_args.out_prefix, convert_args.out_format, threads = convert_args.threads, overwrite = convert_args.overwrite)

            # Rename the plink reference
            shutil.move(convert_args.out_prefix + '.log', expected_out_filename + '.log')

            # Add the reference to the log
            log_plink_reference(expected_out_filename, append_mode = True)

    # Check if input is a specified as a ped
    elif convert_args.ped_prefix or convert_args.ped_filename or convert_args.map_filename:

        # Check that the input and output formats differ
        check_if_conversion('ped', convert_args.out_format)

        # List of formats that are supported from ped
        supported_out_format = ['vcf', 'vcf.gz', 'bcf', 'binary-ped']

        # Check that the conversion is supported
        check_conversion_support(convert_args.out_format, supported_out_format)

        logging.info('Convertion method assigned')

        # Convert the ped file
        convert_ped(**vars(convert_args))

        # Rename the plink reference
        shutil.move(convert_args.out_prefix + '.log', expected_out_filename + '.log')

        # Add the reference to the log
        log_plink_reference(expected_out_filename, append_mode = True)

    # Check if input is a specified as a bed
    elif convert_args.bed_prefix or convert_args.bed_filename or convert_args.fam_filename or convert_args.bim_filename:

        # Check that the input and output formats differ
        check_if_conversion('binary-ped', convert_args.out_format)

        # List of formats that are supported from bed
        supported_out_format = ['vcf', 'vcf.gz', 'bcf', 'ped']

        # Check that the conversion is supported
        check_conversion_support(convert_args.out_format, supported_out_format)

        logging.info('Convertion method assigned')

        # Convert the ped file
        convert_bed(**vars(convert_args))

        # Rename the plink reference
        shutil.move(convert_args.out_prefix + '.log', expected_out_filename + '.log')

        # Add the reference to the log
        log_plink_reference(expected_out_filename, append_mode = True)

if __name__ == "__main__":
    initLogger()
    run()
