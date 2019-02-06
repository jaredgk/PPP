import os
import sys
import argparse
import itertools
import copy
import shutil
import logging

# Import basic vcftools functions
from vcftools import *

# Model file related functions
from model import read_model_file

# Insert Jared's directory path, required for calling Jared's functions. Change when directory structure changes.
sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'pppipe')))

from logging_module import initLogger, logArgs

def vcf_calc_parser(passed_arguments):
    '''VCF Argument Parser - Assigns arguments from command line'''

    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError('%s not found' % value)
                setattr(args, self.dest, value)
        return customAction

    def parser_confirm_file_list ():
        '''Custom action to confirm file exists in list'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                # Loop the list
                for value_item in value:
                    # Check if the file exists
                    if not os.path.isfile(value_item):
                        raise IOError('%s not found' % value_item)
                if not getattr(args, self.dest):
                    setattr(args, self.dest, value)
                else:
                    getattr(args, self.dest).extend(value)
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

    # Model file arguments
    vcf_parser.add_argument('--model-file', help = 'The model filename', type = str, action = parser_confirm_file())
    vcf_parser.add_argument('--model', help = 'Model to analyze', type = str)

    # Non-model file arguments
    vcf_parser.add_argument('--pop-file', help = 'Population file. Note: This argument may be used multiple times and cannont to be used alonside --model', nargs = '+', type = str, action = parser_confirm_file_list())

    # Other file arguments. Expand as needed
    vcf_parser.add_argument('--out', help = 'Output filename. Cannot be used if multiple output files are created', type = str)
    vcf_parser.add_argument('--out-prefix', help = 'Output prefix (vcftools naming scheme)', type = str,  default = 'out')
    vcf_parser.add_argument('--out-dir', help = "Output directory. Only used if multiple output files are created", default = 'Statistic_Files')

    # General arguments
    vcf_parser.add_argument('--overwrite', help = "Overwrite previous output files", action = 'store_true')

    # Statistic based arguments
    statistic_list = ['weir-fst', 'windowed-weir-fst', 'TajimaD', 'site-pi', 'window-pi', 'freq', 'het-fit', 'het-fis', 'hardy-weinberg']
    statistic_default = 'windowed-weir-fst'

    vcf_parser.add_argument('--calc-statistic', metavar = metavar_list(statistic_list), help = 'The statistic to calculate', type = str, choices = statistic_list, default = statistic_default)

    # Statistic window options
    vcf_parser.add_argument('--statistic-window-size', help = 'Window size of relevant statistics', type = int)
    vcf_parser.add_argument('--statistic-window-step', help = 'Step size between windows of relevant statistics', type = int)

    # Position-based position filters
    vcf_parser.add_argument('--filter-include-positions', help = 'Sites to include within a file', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-positions', help = 'Sites to exclude within a file', action = parser_confirm_file())

    # BED-based position filters
    vcf_parser.add_argument('--filter-include-bed', help = 'Set of sites to include within a BED file', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-bed', help = 'Set of sites to exclude within a BED file', action = parser_confirm_file())

    # Non-model based filters
    vcf_parser.add_argument('--filter-include-indv', help = 'Individual to include. Note: This argument may be used multiple times and cannont to be used alonside --model', nargs = '+', type = str, action = parser_add_to_list())
    vcf_parser.add_argument('--filter-exclude-indv', help = 'Individual to exclude. Note: This argument may be used multiple times and cannont to be used alonside --model', nargs = '+', type = str, action = parser_add_to_list())
    vcf_parser.add_argument('--filter-include-indv-file', help = 'Individuals to include in file. Cannont to be used alonside --model', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-indv-file', help = 'Individuals to exclude in file. Cannont to be used alonside --model', action = parser_confirm_file())

    if passed_arguments:
        return vcf_parser.parse_args(passed_arguments)
    else:
        return vcf_parser.parse_args()


def run (passed_arguments = []):
    '''
    Statistic calculation using VCFTools.

    Automates the calculation of site/windowed fixation index (Fst), Tajima's D,
    nucleotide diversity (Pi), allele frequency, and heterozygosity using
    VCFTools. If no statistic is specified, windowed Fst is used by default.

    Parameters
    ----------
    --vcf : str
        Specifies the input VCF filename
    --out : str
        Specifies the output filename
    --pop-file : str
        Defines a population file for calculating Fst-based statistics. May be
        used multple times (i.e. once per file)
    --calc-statistic : str
        Specifies the statistic to calculate. Choices: weir-fst,
        windowed-weir-fst (Default), TajimaD, pi, freq, het
    --statistic-window-size : int
        Specifies the window size for window-based statistics
    --statistic-window-step : int
        Specifies step size between windows for specific window-based statistics
    --statistic-pvalue-cutoff : float
        P-value cutoff for specific statistics
    --filter-include-positions : str
        Specifies a set of sites to include within a tsv file (chromosome and position)
    --filter-exclude-positions : str
        Specifies a set of sites to exclude within a tsv file (chromosome and position)

    Returns
    -------
    output : file
        Statistic file output
    log : file
        Log file output

    Raises
    ------
    IOError
        Input VCF file does not exist
    IOError
        Output file already exists
    '''

    def calc_exception (selected_model, exc_type, exc_value, exc_traceback):

        # Check if the selected model was correctly assigned
        if selected_model: 
            logging.info('Delete temporary files due to exception')
            selected_model.delete_pop_files()
            selected_model.delete_ind_file()

        # Report the original error    
        sys.__excepthook__(exc_type, exc_value, exc_traceback)

    # Grab VCF arguments from command line
    vcf_args = vcf_calc_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(vcf_args, func_name = 'vcf_calc')

    # Argument container for vcftools
    vcftools_call_args = []

    # Used to store population information from either model or pop file(s)
    vcftools_pop_files = []

    # Checks if the user specified both a model and population files
    if vcf_args.model_file and vcf_args.pop_file:
        raise Exception('--model and --pop-file arguments are incompatible')

    # Check that individual-based filters are not being user with a model
    if vcf_args.model_file and (vcf_args.filter_include_indv or vcf_args.filter_exclude_indv):
        if vcf_args.filter_include_indv:
            raise Exception('--model and --filter-include-indv arguments are incompatible')
        if vcf_args.filter_exclude_indv:
            raise Exception('--model and --filter-exclude-indv arguments are incompatible')

    # Check that individuals-based filters are not being user with a model
    if vcf_args.model_file and (vcf_args.filter_include_indv_file or vcf_args.filter_exclude_indv_file):
        if vcf_args.filter_include_indv_file:
            raise Exception('--model and --filter_include_indv-file arguments are incompatible')
        if vcf_args.filter_exclude_indv_file:
            raise Exception('--model and --filter-exclude-indv-file arguments are incompatible')

    # Performs actions related to --model-file argument
    if vcf_args.model_file:

        # Check if a model has been specified
        if not vcf_args.model:
            raise IOError('No selected model. Please use --model to select a model')

        # Read in the model file
        models_in_file = read_model_file(vcf_args.model_file)

        # Check that the selected model was found in the file
        if vcf_args.model not in models_in_file:
            raise IOError('Selected model "%s" not found in: %s' % (vcf_args.model, vcf_args.model_file))

        # Select model, might change this in future versions
        selected_model = models_in_file[vcf_args.model]

        logging.info('%s assigned as model' % selected_model)

        # Check if the specified statistic is fst-based or fis-based
        if vcf_args.calc_statistic in ['windowed-weir-fst', 'weir-fst', 'het-fis']:

            # Create the population files
            selected_model.create_pop_files(file_ext = '.txt', overwrite = vcf_args.overwrite)

            # Store the population files
            vcftools_pop_files = selected_model.pop_files

        # Filter individuals for other statistics
        else:

            # Create individuals file
            selected_model.create_ind_file(overwrite = vcf_args.overwrite)

            # Assign the individuals file to vcftools
            vcftools_call_args.extend(['--keep', selected_model.ind_file])

        # Create a lambda function to pass the selected model along with the exception
        lambda_excepthook = lambda exc_type, exc_value, exc_traceback: calc_exception(selected_model, exc_type, exc_value, exc_traceback)

        # Set the special escape hook
        sys.excepthook = lambda_excepthook

    # Performs actions related to --pop-file argument
    if vcf_args.pop_file:

        # Check if the specified statistic is fst-based or fis-based
        if vcf_args.calc_statistic in ['windowed-weir-fst', 'weir-fst', 'het-fis']:

            # Store the population files
            vcftools_pop_files = vcf_args.pop_file

        # Return an error if --pop-file is not supported
        else:
            raise Exception('%s does not support the --pop-file argument' % vcf_args.calc_statistic)

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

    # Check if windowed Fst is specified
    if vcf_args.calc_statistic == 'windowed-weir-fst':

        # Confirm that population files has been assigned
        if len(vcftools_pop_files) < 2:
            raise Exception('%s requires at least two populations to operate' % vcf_args.calc_statistic)

        # Confirm that the window size has been assigned
        if not vcf_args.statistic_window_size:
            raise Exception('%s requires a specifed window size to operate' % vcf_args.calc_statistic)

        # Confirm that window step size has been assigned
        if not vcf_args.statistic_window_step:
            raise Exception('%s requires a specifed window step size to operate' % vcf_args.calc_statistic)

        # Assigns the required window arguments
        vcftools_call_args.extend(['--fst-window-size', vcf_args.statistic_window_size, '--fst-window-step', vcf_args.statistic_window_step])

        # Assigns the suffix for the vcftools log file
        vcftools_out_suffix = '.windowed.weir.fst'

    elif vcf_args.calc_statistic == 'weir-fst':

        # Confirm that population files has been assigned
        if len(vcftools_pop_files) < 2:
            raise Exception('%s requires at least two populations to operate' % vcf_args.calc_statistic)

        # Assigns the suffix for the vcftools log file
        vcftools_out_suffix = '.weir.fst'

    elif vcf_args.calc_statistic == 'TajimaD':

        # Confirm that the window size has been assigned
        if not vcf_args.statistic_window_size:
            raise Exception('%s requires a specifed window size to operate' % vcf_args.calc_statistic)

        # Assigns all the vcftools arguments for calculating TajimaD
        vcftools_call_args.extend(['--TajimaD', vcf_args.statistic_window_size])

        # Assigns the suffix for the vcftools log file
        vcftools_out_suffix = '.Tajima.D'

    elif vcf_args.calc_statistic == 'site-pi':

        # Assigns all the vcftools arguments for calculating pi
        vcftools_call_args.append('--site-pi')

        # Assigns the suffix for the vcftools log file
        vcftools_out_suffix = '.sites.pi'

    elif vcf_args.calc_statistic == 'window-pi':

        # Confirm that the window size has been assigned
        if not vcf_args.statistic_window_size:
            raise Exception('%s requires a specifed window size to operate' % vcf_args.calc_statistic)

        # Confirm that window step size has been assigned
        if not vcf_args.statistic_window_step:
            raise Exception('%s requires a specifed window step size to operate' % vcf_args.calc_statistic)

        # Assigns all the vcftools arguments for calculating pi
        vcftools_call_args.extend(['--window-pi', vcf_args.statistic_window_size, '--window-pi-step', vcf_args.statistic_window_step])

        # Assigns the suffix for the vcftools log file
        vcftools_out_suffix = '.windowed.pi'

    elif vcf_args.calc_statistic == 'hardy-weinberg':

        # Assigns all the vcftools arguments for the allele frequency
        vcftools_call_args.append('--hardy')

        # Assigns the suffix for the vcftools log file
        vcftools_out_suffix = '.hwe'

    elif vcf_args.calc_statistic == 'freq':

        # Assigns all the vcftools arguments for the allele frequency
        vcftools_call_args.append('--freq')

        # Assigns the suffix for the vcftools log file
        vcftools_out_suffix = '.frq'

    elif vcf_args.calc_statistic == 'het-fit':

        # Assigns all the vcftools arguments for calculating heterozygosity
        vcftools_call_args.append('--het')

        # Assigns the suffix for the vcftools log file
        vcftools_out_suffix = '.het'

    elif vcf_args.calc_statistic == 'het-fis':

        # Confirm that at least one population file has been assigned
        if not vcftools_pop_files:
            raise Exception('%s requires at least one population to operate' % vcf_args.calc_statistic)

        # Assigns all the vcftools arguments for calculating heterozygosity
        vcftools_call_args.append('--het')

        # Assigns the suffix for the vcftools log file
        vcftools_out_suffix = '.het'

    logging.info('vcftools parameters assigned')

    # Check if the user has specified a output filename
    if vcf_args.out:

        # Assign the vcftools output filename, using the output filename
        vcftools_output_filename = vcf_args.out

    else:

        # Assign the vcftools output filename, using the prefix and suffix
        vcftools_output_filename = vcf_args.out_prefix + vcftools_out_suffix

    # Check if previous output should be overwritten
    if vcf_args.overwrite:

        # Check if the ouput will require the output directory
        if vcf_args.calc_statistic in ['windowed-weir-fst', 'weir-fst'] and len(vcftools_pop_files) > 2:
            
            # Check if the output directory is present
            if os.path.exists(vcf_args.out_dir):

                # Remove the output directory
                shutil.rmtree(vcf_args.out_dir)

    # If not to be overwritten, check if previous output exists
    else:

        # Confirm the vcftools output and log file do not exist
        check_for_vcftools_output(vcftools_output_filename)

        # Check if the ouput will require the output directory
        if vcf_args.calc_statistic in ['windowed-weir-fst', 'weir-fst'] and len(vcftools_pop_files) > 2:

            # Check if the output directory is present and report the error
            if os.path.exists(vcf_args.out_dir):
                raise IOError('Statistic Directory already exists')

    # Assigns the file argument for vcftools
    vcfname_arg = assign_vcftools_input_arg(vcf_args.vcf)

    logging.info('Input file assigned')

    # Run vcftools once if the statistic isn't het-fis
    if vcf_args.calc_statistic in ['windowed-weir-fst', 'weir-fst'] and len(vcftools_pop_files) > 2:

        def return_filename (filepath):
            return os.path.basename(filepath).split(os.extsep)[0]

        # Create the output directory
        if not os.path.exists(vcf_args.out_dir):
            os.makedirs(vcf_args.out_dir)

        # Loop each population file
        for first_pop_filepath, second_pop_filepath in itertools.combinations(vcftools_pop_files, 2):

            # Create the population-specific call
            pop_call_args = copy.deepcopy(vcftools_call_args)

            # Assign the population files
            pop_call_args.extend(['--weir-fst-pop', first_pop_filepath, '--weir-fst-pop', second_pop_filepath])

            # Extract filename from first filepath
            first_pop_filename = return_filename(first_pop_filepath)

            # Extract filename from second filepath
            second_pop_filename = return_filename(second_pop_filepath)

            # Create the population prefix, and join to the output directory
            pop_prefix = os.path.join(vcf_args.out_dir, vcf_args.out_prefix)

            # Update the population prefix with the population names
            pop_prefix += '.%s.%s' % (first_pop_filename, second_pop_filename)

            # The filename population file
            pop_filename = pop_prefix + vcftools_out_suffix

            # vcftools subprocess call, with stdout
            vcftools_err = call_vcftools(vcfname_arg + pop_call_args, output_filename = pop_filename)

            # Produce the vcftools log file, in append mode
            produce_vcftools_log(vcftools_err, vcftools_output_filename, append_mode = True)

    elif vcf_args.calc_statistic == 'het-fis':

        # Loop each population file
        for vcftools_pop_file in vcftools_pop_files:

            # Create the population-specific call
            pop_call_args = vcftools_call_args + ['--keep', str(vcftools_pop_file)]

            # vcftools subprocess call
            vcftools_err = call_vcftools(vcfname_arg + pop_call_args, output_filename = vcftools_output_filename, append_mode = True)

            # Produce the vcftools log file, in append mode
            produce_vcftools_log(vcftools_err, vcftools_output_filename, append_mode = True)

    else:

        # Check if either Fst statistic is specified
        if vcf_args.calc_statistic in ['windowed-weir-fst', 'weir-fst'] and len(vcftools_pop_files) == 2:

            # Assigns the population files to the vcftools call
            vcftools_call_args.extend([pop_args for pop_file in vcftools_pop_files for pop_args in ['--weir-fst-pop', pop_file]])

        # Check if the user has not specified an output filename
        if vcf_args.out:

            # Assign a unique output prefix 
            vcf_args.out_prefix = assign_vcftools_unique_prefix(vcf_args.out_prefix, vcftools_output_suffix)

        # Add the output argument
        vcftools_call_args.extend(['--out', vcf_args.out_prefix])

        # vcftools subprocess call
        vcftools_err = call_vcftools(vcfname_arg + vcftools_call_args)

        # Produce the vcftools log file
        produce_vcftools_log(vcftools_err, vcftools_output_filename)

    # Delete any files that were created for vcftools
    if vcf_args.model_file:
        selected_model.delete_pop_files()
        selected_model.delete_ind_file()

if __name__ == "__main__":
    initLogger()
    run()
