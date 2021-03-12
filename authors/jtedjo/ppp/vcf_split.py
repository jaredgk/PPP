import os
import sys
import copy
import argparse
import logging
import shutil
import pandas as pd
import numpy as np

# Import basic vcftools functions
from vcftools import *

# Model file related functions
from model import read_model_file

# Insert Jared's directory path, required for calling Jared's functions. Change when directory structure changes.
sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'jared')))

from logging_module import initLogger, logArgs

def vcf_split_parser(passed_arguments):
    '''VCF Argument Parser - Assigns arguments from command line'''

    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError('%s not found' % value)
                setattr(args, self.dest, value)
        return customAction

    def parser_add_to_list ():
        '''Custom action to add items to a list'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not getattr(args, self.dest):
                    setattr(args, self.dest, value)
                else:
                    getattr(args, self.dest).extend(value)
        return customAction

    def metavar_list (var_list):
        '''Create a formmated metavar list for the help output'''
        return '{' + ', '.join(var_list) + '}'

    vcf_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    # Input arguments.
    vcf_parser.add_argument('--vcf', help = "Input VCF filename", type = str, required = True, action = parser_confirm_file())

    # Model file arguments.
    vcf_parser.add_argument('--model-file', help = 'Defines the model file', type = str, action = parser_confirm_file())
    vcf_parser.add_argument('--model', help = 'Defines the model to analyze', type = str)

    # Split arguments
    vcf_parser.add_argument('--split-file', help='Specifies the file used for splitting', required = True, type = str, action = parser_confirm_file())

    split_list = ['statistic-file', 'bed']
    split_default = 'statistic-file'

    vcf_parser.add_argument('--split-method', metavar = metavar_list(split_list), help = 'Specifies the splitting method', type=str, choices = split_list, default = split_default)

    vcf_parser.add_argument('--statistic-window-size', help = 'Size of window calculations (use if BIN_END is absent in file)', type = int)

    vcf_parser.add_argument('--window-correction', help = 'If the window should be corrected to avoid single-position overlap (i.e. n-1)', type = bool, default = True)

    # Output file arguments
    out_format_list = ['vcf', 'vcf.gz', 'bcf']
    out_format_default = 'vcf.gz'

    vcf_parser.add_argument('--out-format', metavar = metavar_list(out_format_list), help = 'Specifies the VCF output format', type = str, choices = out_format_list, default = out_format_default)

    vcf_parser.add_argument('--out-dir', help = 'Specifies the output directory', type = str,  default = 'Split_VCFs')
    vcf_parser.add_argument('--out-prefix', help = 'Specifies the output prefix', type = str,  default = 'out')

    vcf_parser.add_argument('--log-prefix', help = 'Specifies the log file file prefix', type = str,  default = 'out.split')


    # General arguments.
    vcf_parser.add_argument('--overwrite', help = "overwrite previous output files", action = 'store_true')

    ### Filters

    # Non-model file arguments
    vcf_parser.add_argument('--filter-include-indv', help = 'Individual to include. Note: This argument may be used multiple times and cannont to be used alonside --model', nargs = '+', type = str, action = parser_add_to_list())
    vcf_parser.add_argument('--filter-exclude-indv', help = 'Individual to exclude. Note: This argument may be used multiple times and cannont to be used alonside --model', nargs = '+', type = str, action = parser_add_to_list())
    vcf_parser.add_argument('--filter-include-indv-file', help = 'Individuals to include in file. Cannont to be used alonside --model', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-indv-file', help = 'Individuals to exclude in file. Cannont to be used alonside --model', action = parser_confirm_file())

    # Position-based position filters
    vcf_parser.add_argument('--filter-include-positions', help = 'Sites to include within a chr/pos tab-sperated file', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-positions', help = 'Sites to exclude within a chr/pos tab-sperated file', action = parser_confirm_file())

    # BED-based position filters
    vcf_parser.add_argument('--filter-include-bed', help = 'Set of sites to include within a BED file', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-bed', help = 'Set of sites to exclude within a BED file', action = parser_confirm_file())


    if passed_arguments:
        return vcf_parser.parse_args(passed_arguments)
    else:
        return vcf_parser.parse_args()

def return_missing_columns (sample_headers):
    '''Reports if a required column (i.e. CHROM, BIN_START, BIN_END) is missing
    from the header of the statistic file.'''

    if not all(pos_header in sample_headers for pos_header in ['CHROM', 'BIN_START', 'BIN_END']):
        # Assign missing headers
        missing_headers = [pos_header for pos_header in ['CHROM', 'BIN_START', 'BIN_END'] if pos_header not in sample_headers]

        return missing_headers

    return None

def run (passed_arguments = []):
    '''
    Split VCF file into multiple VCFs.

    Splits a single VCF file into multiple VCF files using either a statistic or
    bed file. The statistic/bed file must contain loci-based (i.e. window-based)
    data for the function to operate. If the specified statistic file does not
    contain a BIN_END column, the --statistic-window-size argument may be used.

    Parameters
    ----------
    --vcf : str
        Specifies the input VCF filename
    --model-file
        Specifies the model filename
    --model
        Specifies the model
    --split-file : str
        Specifies the file used for splitting
    --split-method : str
        Specifies the method being used to split {statistic-file, bed} (Default: statistic-file)
    --statistic-window-size : int
        Specifies the size of window calculations (use if BIN_END is absent)
    --out-format : str
        Specifies the output format {vcf, vcf.gz, bcf} (Default: vcf.gz)
    --out-dir : str
        Specifies the output directory
    --out-prefix : str
        Specifies the output prefix
    --filter-include-indv : list or str
        Individual to include. Note: This argument may be used multiple times
        and cannont to be used alonside --model
    --filter-exclude-indv : list or str
        Individual to exclude. Note: This argument may be used multiple times
        and cannont to be used alonside --model
    --filter-include-indv-file : str
        Individuals to include in file. Cannont to be used alonside --model
    --filter-exclude-indv-file : str
        Individuals to exclude in file. Cannont to be used alonside --model
    --filter-include-positions : str
        Specifies a set of sites to include within a tsv file (chromosome and position)
    --filter-exclude-positions : str
        Specifies a set of sites to exclude within a tsv file (chromosome and position)
    --filter-include-bed : str
        Specifies a set of sites to include within a BED file
    --filter-exclude-bed : str
        Specifies a set of sites to exclude within a BED file

    Returns
    -------
    output : directory
        Directory of split VCF files
    log : file
        Log file output

    Raises
    ------
    IOError
        Input VCF file does not exist
    IOError
        Model file does not exist
    IOError
        Output file already exists and --overwrite is not specified

    '''

    # Grab VCF arguments from command line
    vcf_args = vcf_split_parser(passed_arguments)

    # Argument container for vcftools
    vcftools_call_args = []

    # Check if previous output should be overwritten
    if not vcf_args.overwrite:
        # Check if the output directory is present
        if os.path.isdir(vcf_args.out_dir):
            raise IOError('%s already exists' % vcf_args.out_dir)
    else:
        # Check if the output directory is present
        if os.path.isdir(vcf_args.out_dir):
            # Delete the directory
            shutil.rmtree(vcf_args.out_dir)

    # Check if the user has specified a model file
    if vcf_args.model_file and vcf_args.model:
        # Read in the models
        models_in_file = read_model_file(vcf_args.model_file)

        # Check that the selected model was not found in the file
        if vcf_args.model not in models_in_file:
            raise IOError('Selected model "%s" not found in: %s'
                           % (vcf_args.model, vcf_args.model_file))

        # Check that individual-based filters are not also being used
        if vcf_args.filter_include_indv or vcf_args.filter_exclude_indv:
            if vcf_args.filter_include_indv:
                raise Exception('--model and --filter-include-indv arguments '
                                'are incompatible')
            if vcf_args.filter_exclude_indv:
                raise Exception('--model and --filter-exclude-indv arguments '
                                'are incompatible')

        # Check that individuals-based filters are not also being used
        if vcf_args.filter_include_indv_file or vcf_args.filter_exclude_indv_file:
            if vcf_args.filter_include_indv_file:
                raise Exception('--model and --filter_include_indv-file '
                                'arguments are incompatible')
            if vcf_args.filter_exclude_indv_file:
                raise Exception('--model and --filter-exclude-indv-file '
                                'arguments are incompatible')

        # Select model, might change this in future versions
        selected_model = models_in_file[vcf_args.model]

        # Create individuals file
        selected_model.create_ind_file(overwrite = vcf_args.overwrite)

        # Assign the individuals file to vcftools
        vcftools_call_args.extend(['--keep', selected_model.ind_file])

        logging.info('Model file assigned')

    # Individuals-based filters
    if vcf_args.filter_include_indv_file or vcf_args.filter_exclude_indv_file:
        # Used to include a file of individuals to keep
        if vcf_args.filter_exclude_indv_file:
            vcftools_call_args.extend(['--keep', vcf_args.filter_include_indv_file])

        # Used to include a file of individuals to remove
        if vcf_args.filter_exclude_indv_file:
            vcftools_call_args.extend(['--remove', vcf_args.filter_exclude_indv_file])

    # Individual-based filters
    if vcf_args.filter_include_indv or vcf_args.filter_exclude_indv:
        if vcf_args.filter_include_indv:
            for indv_to_include in vcf_args.filter_include_indv:
                vcftools_call_args.extend(['--indv', indv_to_include])
        if vcf_args.filter_exclude_indv:
            for indv_to_exclude in vcf_args.filter_exclude_indv:
                vcftools_call_args.extend(['--remove-indv', indv_to_exclude])

    # Position (vcftools output file) filters
    if vcf_args.filter_include_positions or vcf_args.filter_exclude_positions:
        if vcf_args.filter_include_positions:
            vcftools_call_args.extend(['--positions', vcf_args.filter_include_positions])
        if vcf_args.filter_exclude_positions:
            vcftools_call_args.extend(['--exclude-positions', vcf_args.filter_exclude_positions])

    # Position (BED format file) filters
    if vcf_args.filter_include_bed or vcf_args.filter_exclude_bed:
        if vcf_args.filter_include_bed:
            vcftools_call_args.extend(['--bed', vcf_args.filter_include_bed])
        if vcf_args.filter_exclude_bed:
            vcftools_call_args.extend(['--exclude-bed', vcf_args.filter_exclude_bed])

    # Used to assign the output format to the vcftools call, and assign
    # vcftools output suffix
    if vcf_args.out_format == 'bcf':
        vcftools_call_args.append('--recode')
        vcftools_output_suffix = '.recode.bcf'
    elif vcf_args.out_format == 'vcf':
        vcftools_call_args.append('--recode')
        vcftools_output_suffix = '.recode.vcf'
    elif vcf_args.out_format == 'vcf.gz':
        vcftools_call_args.append('--recode')
        vcftools_output_suffix = '.recode.vcf.gz'

    logging.info('vcftools common parameters assigned')

    # Create the vcf/bcf output directory
    if not os.path.exists(vcf_args.out_dir):
        os.makedirs(vcf_args.out_dir)

    # Assigns the file argument for vcftools
    vcfname_arg = assign_vcftools_input_arg(vcf_args.vcf)

    logging.info('Input file assigned')

    # Read in the sample file
    split_samples = pd.read_csv(vcf_args.split_file, sep = '\t')

    logging.info('Split file assigned')

    # Assign the columns
    split_columns = list(split_samples)

    # Assign missing columns
    missing_columns = return_missing_columns(split_columns)

    # Check for missing columns
    if missing_columns:

        # Report error if too many or unexpected columns are missing
        if len(missing_columns) != 1 or 'BIN_END' not in missing_columns:
            raise ValueError('Cannot find %s column(s) in file specified by ' \
                             '--split-file.' % ', '.join(missing_headers))

        # Check if statistic_window_size has been assigned
        if not vcf_args.statistic_window_size:
            raise ValueError("'BIN_END' column absent in %s. Please use " \
                             '--statistic-window-size' % vcf_args.split_file)

        # Assign window size
        window_size = vcf_args.statistic_window_size

        # Check if the overlap correction has been specified
        if vcf_args.window_correction:
            # Correct the window size
            window_size -= 1

            logging.info('Applied window correction. To disable the correction '
                         "use '--window-correction False'")

        # Create the 'BIN_END' column
        split_samples['BIN_END'] = split_samples['BIN_START'] + window_size

        logging.info("Calculated 'BIN_END' using --statistic-window-size")

    # Loop the rows of the sample to be split
    for sample_index, split_sample in split_samples.iterrows():

        # Copy common arguments for VCF call
        sample_call_args = copy.deepcopy(vcftools_call_args)

        # Assign the sample output prefix
        sample_prefix = '%s_%s' % (vcf_args.out_prefix, sample_index)

        # Add the path to the prefix
        sample_path = os.path.join(vcf_args.out_dir, sample_prefix)

        # Store the sample output prefix
        sample_call_args.extend(['--out', sample_path])

        # Store the sample position information
        sample_call_args.extend(['--chr', split_sample['CHROM'],
                                 '--from-bp', split_sample['BIN_START'],
                                 '--to-bp', split_sample['BIN_END']])

        # Assign the expected output filename
        vcftools_sample_filename = sample_path + vcftools_output_suffix

        logging.info('vcftools sample parameters assigned')

        # Call vcftools with the specifed arguments
        vcftools_err = call_vcftools(vcfname_arg + sample_call_args, output_format = vcf_args.out_format, output_filename = vcftools_sample_filename)

        # Produce the log file (in append mode, will create a single log)
        produce_vcftools_log(vcftools_err, vcf_args.log_prefix, append_mode = True)


if __name__ == "__main__":
    initLogger()
    run()
