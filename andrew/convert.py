import os
import sys
import subprocess
import shutil
import argparse
import glob
import logging

sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'pppipe')))


# Import PPP modules and scripts
import vcftools
import bcftools
import plink_tasks
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

    # Input FASTA argument
    convert_parser.add_argument("--fasta", help = "Input FASTA filename", type = str, action = parser_confirm_file())

    # Input IM argument
    convert_parser.add_argument("--im", help = "Input IM filename", type = str, action = parser_confirm_file())

    # Input PED arguments
    convert_parser.add_argument("--ped", help = "Input PED filename", type = str, action = parser_confirm_file())
    convert_parser.add_argument("--map", help = "Input MAP filename. Called alongside --ped", type = str, action = parser_confirm_file())

    # Input BED arguments
    convert_parser.add_argument("--bed", help = "Input BED filename", type = str, action = parser_confirm_file())
    convert_parser.add_argument("--fam", help = "Input FAM filename. Called alongside --bed", type = str, action = parser_confirm_file())
    convert_parser.add_argument("--bim", help = "Input BIM filename. Called alongside --bed", type = str, action = parser_confirm_file())

    # Input PED/BED prefix arguments
    convert_prefix = convert_parser.add_mutually_exclusive_group()
    convert_prefix.add_argument("--ped-prefix", help = "Input PED or BED filename prefix", type = str)
    convert_prefix.add_argument("--bed-prefix", help = "Input PED or BED filename prefix", type = str)

    # Output arguments.
    output_formats = ['vcf', 'vcf.gz', 'bcf', 'fasta', 'im', 'ped', 'binary-ped']
    convert_parser.add_argument('--out-format', metavar = metavar_list(output_formats), help = 'Specifies the output format', type = str, choices = output_formats, required = True)
    convert_parser.add_argument('--out', help = 'Defines the output filename')
    convert_parser.add_argument('--out-prefix', help = 'Defines the output filename prefix', default = 'out')

    # Other basic arguments. Expand as needed
    convert_parser.add_argument('--delete-original', help = 'Defines that the original file should be deleted once converted', action='store_false')

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

def convert_vcf (vcf_filename, out_prefix, out_format):
    '''
        Converts vcf files

        Parameters
        ----------
        vcf_filename : str
            Filename of the VCF file to be converted
        out_prefix : str
            Specified output prefix (i.e. filename without extension)
        out_format : str
            Specified output format
    '''

    logging.info('Convertion method assigned')

    # Check if the output is a bgzipped vcf
    if out_format == 'vcf.gz':

        logging.info('Beginning conversion')

        # Call the vcf to vcf.gz converter (i.e. compress)
        vcftools.bgzip_compress_vcf(vcf_filename, out_prefix = 'out', keep_original = True)

    # Check if the output is a bcf
    elif out_format == 'bcf':

        logging.info('Beginning conversion')

        # Call the vcf to bcf converter
        bcftools.convert_to_bcf(vcf_filename, out_prefix)

    # Check if the output is a ped or binary-ped
    elif out_format == 'ped' or out_format == 'binary-ped':

        logging.info('Beginning conversion')

        # Stores arguments for ped converion
        convert_plink_args = []

        # Assign input file
        convert_plink_args.extend(['--vcf', vcf_filename])

        # Assign output format
        convert_plink_args.extend(['--out-format', out_format])

        # Assign output filename prefix
        convert_plink_args.extend(['--out-prefix', out_prefix])

        # Call the vcf to ped/binary-ped converter (i.e. plink automater)
        plink_tasks.run(convert_plink_args)

    logging.info('Finished conversion')

def convert_vcfgz (vcfgz_filename, out_prefix, out_format):
    '''
        Converts vcf.gz files

        Parameters
        ----------
        vcfgz_filename : str
            Filename of the VCF.GZ file to be converted
        out_prefix : str
            Specified output prefix (i.e. filename without extension)
        out_format : str
            Specified output format
    '''

    logging.info('Convertion method assigned')

    # Check if the output is a vcf
    if out_format == 'vcf':

        logging.info('Beginning conversion')

        # Call the vcf.gz to vcf converter (i.e. decompress)
        vcftools.bgzip_decompress_vcfgz(vcfgz_filename, out_prefix = 'out', keep_original = True)

    # Check if the output is a bcf
    elif out_format == 'bcf':

        logging.info('Beginning conversion')

        # Call the vcf.gz to bcf converter
        bcftools.convert_to_bcf(vcfgz_filename, out_prefix)

    # Check if the output is a ped or binary-ped
    elif out_format == 'ped' or out_format == 'binary-ped':

        logging.info('Beginning conversion')

        # Stores arguments for ped converion
        convert_plink_args = []

        # Assign input file
        convert_plink_args.extend(['--vcf', vcfgz_filename])

        # Assign output format
        convert_plink_args.extend(['--out-format', out_format])

        # Assign output filename prefix
        convert_plink_args.extend(['--out-prefix', out_prefix])

        # Call the vcf to ped/binary-ped converter (i.e. plink automater)
        plink_tasks.run(convert_plink_args)

    logging.info('Finished conversion')

def convert_bcf (bcf_filename, out_prefix, out_format):
    '''
        Converts bcf files

        Parameters
        ----------
        bcf_filename : str
            Filename of the BCF file to be converted
        out_prefix : str
            Specified output prefix (i.e. filename without extension)
        out_format : str
            Specified output format
    '''

    logging.info('Convertion method assigned')

    # Check if the output is a vcf
    if out_format == 'vcf':

        logging.info('Beginning conversion')

        # Call the bcf to vcf converter
        bcftools.convert_to_vcf(bcf_filename, out_prefix)

    # Check if the output is a compressed vcf
    elif out_format == 'vcf.gz':

        logging.info('Beginning conversion')

        # Call the bcf to vcfgz converter
        bcftools.convert_to_vcfgz(bcf_filename, out_prefix)

    # Check if the output is a ped or binary-ped
    elif out_format == 'ped' or out_format == 'binary-ped':

        logging.info('Beginning conversion')

        # Stores arguments for ped converion
        convert_plink_args = []

        # Assign input file
        convert_plink_args.extend(['--vcf', bcf_filename])

        # Assign output format
        convert_plink_args.extend(['--out-format', out_format])

        # Assign output filename prefix
        convert_plink_args.extend(['--out-prefix', out_prefix])

        # Call the vcf to ped/binary-ped converter (i.e. plink automater)
        plink_tasks.run(convert_plink_args)


    logging.info('Finished conversion')

def convert_ped (ped_filename, map_filename, ped_prefix, out_prefix, out_format):
    '''
        Converts ped (plink) files

        Parameters
        ----------
        ped_filename : str
            Filename of the PED file to be converted
        map_filename : str
            Filename of the PED-assoicated MAP file
        ped_prefix : str
            Specified ped prefix (i.e. prefix for the ped-related files)
        out_prefix : str
            Specified output prefix (i.e. filename without extension)
        out_format : str
            Specified output format
    '''

    logging.info('Convertion method assigned')
    logging.info('Beginning conversion')

    # Stores arguments for ped converion
    convert_plink_args = []

    # Check if the files are being assigned by a prefix
    if ped_prefix:

        # Assign input file
        convert_plink_args.extend(['--ped-prefix', ped_prefix])

    # Check if the files are being assigned seperately
    elif ped_filename or map_filename:

        # Check that the required files are specified
        if ped_filename and map_filename:

            # Assign input files
            convert_plink_args.extend(['--ped', ped_filename,
                                       '--map', map_filename])

        # Return error if both files are not present
        else:
            raise Exception('Both a PED and MAP file are required for conversion')

    # Assign output format
    convert_plink_args.extend(['--out-format', out_format])

    # Assign output filename prefix
    convert_plink_args.extend(['--out-prefix', out_prefix])

    # Call the vcf to ped/binary-ped converter (i.e. plink automater)
    plink_tasks.run(convert_plink_args)

    logging.info('Finished conversion')

def convert_bed (bed_filename, fam_filename, bim_filename, bed_prefix, out_prefix, out_format):
    '''
        Converts bed (i.e. binary-ped) files

        Parameters
        ----------
        bed_filename : str
            Filename of the BED file to be converted
        fam_filename : str
            Filename of the BED-assoicated FAM file
        bim_filename : str
            Filename of the BED-assoicated BIM file
        bed_prefix : str
            Specified bed prefix (i.e. prefix for the bed-related files)
        out_prefix : str
            Specified output prefix (i.e. filename without extension)
        out_format : str
            Specified output format
    '''

    logging.info('Convertion method assigned')
    logging.info('Beginning conversion')

    # Stores arguments for bed converion
    convert_plink_args = []

    # Check if the files are being assigned by a prefix
    if bed_prefix:

        # Assign input file
        convert_plink_args.extend(['--bed-prefix', bed_prefix])

    # Check if the files are being assigned seperately
    elif bed_filename or fam_filename or bim_filename:

        # Check that the required files are specified
        if bed_filename and fam_filename and bim_filename:

            # Assign input files
            convert_plink_args.extend(['--bed', bed_filename,
                                       '--fam', fam_filename,
                                       '--bim', bim_filename])

        # Return error if both files are not present
        else:
            raise Exception('A BED, FAM, and BIM file are required for conversion')

    # Assign output format
    convert_plink_args.extend(['--out-format', out_format])

    # Assign output filename prefix
    convert_plink_args.extend(['--out-prefix', out_prefix])

    # Call the vcf to ped/binary-ped converter (i.e. plink automater)
    plink_tasks.run(convert_plink_args)

    logging.info('Finished conversion')

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

        # Check that the input is a vcf
        if vcf_format == 'vcf':

            # Convert the vcf file
            convert_vcf(convert_args.vcf, convert_args.out_prefix, convert_args.out_format)

        # Check that the input is a vcf.gz
        elif vcf_format == 'vcf.gz':

            # Convert the vcf.gz file
            convert_vcfgz(convert_args.vcf, convert_args.out_prefix, convert_args.out_format)

        # Check that the input is a bcf
        elif vcf_format == 'bcf':

            # Convert the bcf file
            convert_bcf(convert_args.vcf, convert_args.out_prefix, convert_args.out_format)

    # Check if input is a specified as a ped
    elif convert_args.ped_prefix or convert_args.ped or convert_args.map:

        # List of formats that are supported from ped
        supported_out_format = ['vcf', 'vcf.gz', 'bcf', 'binary-ped']

        # Check that the conversion is supported
        check_conversion_support(convert_args.out_format, supported_out_format)

        # Convert the ped file
        convert_ped(convert_args.ped, convert_args.map, convert_args.ped_prefix, convert_args.out_prefix, convert_args.out_format)

    # Check if input is a specified as a bed
    elif convert_args.bed_prefix or convert_args.bed or convert_args.fam or convert_args.bim:

        # List of formats that are supported from bed
        supported_out_format = ['vcf', 'vcf.gz', 'bcf', 'ped']

        # Check that the conversion is supported
        check_conversion_support(convert_args.out_format, supported_out_format)

        # Convert the ped file
        convert_bed(convert_args.bed, convert_args.fam, convert_args.bim, convert_args.bed_prefix, convert_args.out_prefix, convert_args.out_format)


if __name__ == "__main__":
    initLogger()
    run()
