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
sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'jared')))

from logging_module import initLogger

def vcf_calc_parser(passed_arguments):
    '''VCF Argument Parser - Assigns arguments from command line'''

    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError('Input not found.') # File not found
                setattr(args, self.dest, value)
        return customAction

    def metavar_list (var_list):
        '''Create a formmated metavar list for the help output'''
        return '{' + ', '.join(var_list) + '}'

    vcf_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    # Input arguments.
    vcf_parser.add_argument("vcfname", metavar = 'VCF_Input', help = "Input VCF filename", type = str, action = parser_confirm_file())

    # Model file arguments.
    vcf_parser.add_argument('--model-file', help = 'Defines the model file', type = str, action = parser_confirm_file())
    vcf_parser.add_argument('--model', help = 'Defines the model to analyze', type = str)

    # Other file arguments. Expand as needed
    vcf_parser.add_argument('--pop-file', help = 'Defines the population files for calculating specific statistics', type = str, action='append')
    vcf_parser.add_argument('--out', help = 'Specifies the complete output filename. Cannot be used if multiple output files are created', type = str)
    vcf_parser.add_argument('--out-prefix', help = 'Specifies the output prefix (vcftools naming scheme)', type = str,  default = 'out')
    vcf_parser.add_argument('--out-dir', help = "Specifies the output directory. Only used if multiple output files are created", default = 'Statistic_Files')

    # General arguments.
    vcf_parser.add_argument('--overwrite', help = "Specifies if previous output files should be overwritten", action = 'store_true')

    # Statistic based arguments.
    statistic_list = ['weir-fst', 'windowed-weir-fst', 'TajimaD', 'site-pi', 'window-pi', 'freq', 'het-fit', 'het-fis']
    statistic_default = 'windowed-weir-fst'

    vcf_parser.add_argument('--calc-statistic', metavar = metavar_list(statistic_list), help = 'Specifies the statistic to calculate', type = str, choices = statistic_list, default = statistic_default)

    # Statistic window options
    vcf_parser.add_argument('--statistic-window-size', help = 'Specifies the size of window calculations', type = int, default = 10000)
    vcf_parser.add_argument('--statistic-window-step', help = 'Specifies step size between windows', type = int, default = 20000)

    # Position-based position filters
    vcf_parser.add_argument('--filter-include-positions', help = 'Specifies a set of sites to include within a file (tsv chromosome and position)', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-positions', help = 'Specifies a set of sites to exclude within a file (tsv chromosome and position)', action = parser_confirm_file())

    if passed_arguments:
        return vcf_parser.parse_args(passed_arguments)
    else:
        return vcf_parser.parse_args()

def logArgs(args, pipeline_function):
    '''Logs arguments from argparse system. May be replaced with a general function in logging_module'''
    logging.info('Arguments for %s:' % pipeline_function)
    for k in vars(args):
        logging.info('Arguments %s: %s' % (k, vars(args)[k]))

def run (passed_arguments = []):
    '''
    Statistic calculation using VCFTools.

    Automates the calculation of site/windowed fixation index (Fst), Tajima's D,
    nucleotide diversity (Pi), allele frequency, and heterozygosity using
    VCFTools. If no statistic is specified, windowed Fst is used by default.

    Parameters
    ----------
    VCF_Input : str
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
        Specifies step size between windows for spcific window-based statistics
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

    # Grab VCF arguments from command line
    vcf_args = vcf_calc_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(vcf_args, 'vcf_calc')

    # Argument container for vcftools
    vcftools_call_args = []

    # Used to store population information from either model or pop file(s)
    vcftools_pop_files = []

    # Checks population assignment for the current statistic
    if vcf_args.calc_statistic in ['windowed-weir-fst', 'weir-fst', 'het-fis']:

        # Check if population assignment is possible
        if not vcf_args.model_file:
            logging.error('%s requires a model file to operate. Please use --model-file to select a file' % vcf_args.calc_statistic)
            raise IOError('%s requires a model file to operate. Please use --model-file to select a file' % vcf_args.calc_statistic)

    # Performs actions related to --model-file argument
    if vcf_args.model_file:

        # Check if a model has been specified
        if not vcf_args.model:
            logging.error('No selected model. Please use --model to select a model')
            raise IOError('No selected model. Please use --model to select a model')

        # Read in the model file
        models_in_file = read_model_file(vcf_args.model_file)

        # Check that the selected model was found in the file
        if vcf_args.model not in models_in_file:
            logging.error('Selected model "%s" not found in: %s' % (vcf_args.model, vcf_args.model_file))
            raise IOError('Selected model "%s" not found in: %s' % (vcf_args.model, vcf_args.model_file))

        # Select model, might change this in future versions
        selected_model = models_in_file[vcf_args.model]

        # Check if the specified statistic is Fst-based
        #if vcf_args.calc_statistic in ['windowed-weir-fst', 'weir-fst']:
        if vcf_args.calc_statistic in ['windowed-weir-fst', 'weir-fst', 'het-fis']:

            # Check if there are enough populations
            if selected_model.npop < 2:
                logging.error('Two or more populations requried. Please check selected model')
                raise IOError('Two or more populations requried. Please check selected model')

            # Create the population files
            selected_model.create_pop_files(file_ext = '.txt', overwrite = vcf_args.overwrite)

            # Store the population files
            vcftools_pop_files = selected_model.pop_files

            # Check if either Fst statistic is specified
            if vcf_args.calc_statistic in ['windowed-weir-fst', 'weir-fst'] and len(vcftools_pop_files) == 2:

                # Assigns the population files to the vcftools call
                vcftools_call_args.extend([pop_args for pop_file in vcftools_pop_files for pop_args in ['--weir-fst-pop', pop_file]])

        else:

            # Create individuals file
            selected_model.create_individuals_file(overwrite = vcf_args.overwrite)

            # Assign the individuals file to vcftools
            vcftools_call_args.extend(['--keep', selected_model.individuals_file])


    # Check if windowed Fst is specified
    if vcf_args.calc_statistic == 'windowed-weir-fst':

        # Assigns the required window arguments
        vcftools_call_args.extend(['--fst-window-size', vcf_args.statistic_window_size, '--fst-window-step', vcf_args.statistic_window_step])

        # Assigns the suffix for the vcftools log file
        vcftools_out_suffix = '.windowed.weir.fst'

    elif vcf_args.calc_statistic == 'weir-fst':

        # Assigns the suffix for the vcftools log file
        vcftools_out_suffix = '.weir.fst'

    elif vcf_args.calc_statistic == 'TajimaD':

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

        # Assigns all the vcftools arguments for calculating pi
        vcftools_call_args.extend(['--window-pi', vcf_args.statistic_window_size, '--window-pi-step', vcf_args.statistic_window_step])

        # Assigns the suffix for the vcftools log file
        vcftools_out_suffix = '.windowed.pi'

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

        # Assigns all the vcftools arguments for calculating heterozygosity
        vcftools_call_args.append('--het')

        # Assigns the suffix for the vcftools log file
        vcftools_out_suffix = '.het'

    # Assing the filtered positions
    if vcf_args.filter_include_positions or vcf_args.filter_exclude_positions:

        # Assigns the sites to keep
        if vcf_args.filter_include_positions:
            vcftools_call_args.extend(['--positions', vcf_args.filter_include_positions])

        # Assigns the sites to remove
        if vcf_args.filter_exclude_positions:
            vcftools_call_args.extend(['--exclude-positions', vcf_args.filter_exclude_positions])

    logging.info('vcftools parameters assigned')

    # The filename vcftools assigns to the statistic output
    vcftools_output_filename = vcf_args.out_prefix + vcftools_out_suffix

    # Check if previous output should be overwritten
    if vcf_args.overwrite:

        if vcf_args.out:
            # Confirm the vcftools-renamed output and log file do not exist
            delete_vcftools_output(vcf_args.out)
        else:
            # Confirm the vcftools output and log file do not exist
            delete_vcftools_output(vcftools_output_filename)

        # Check if the output directory is present
        if os.path.exists(vcf_args.out_dir):
            shutil.rmtree(vcf_args.out_dir)

    # If not to be overwritten, check if previous output exists
    else:
        if vcf_args.out:
            # Confirm the vcftools-renamed output and log file do not exist
            check_for_vcftools_output(vcf_args.out)
        else:
            # Confirm the vcftools output and log file do not exist
            check_for_vcftools_output(vcftools_output_filename)

        # Check if the output directory is present
        if os.path.exists(vcf_args.out_dir):
            logging.error('Statistic Directory already exists')
            raise IOError('Statistic Directory already exists')

    # Assigns the file argument for vcftools
    vcfname_arg = assign_vcftools_input_arg(vcf_args.vcfname)

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

            # Assign the output argument
            pop_call_args.extend(['--out', pop_prefix])

            # vcftools subprocess call, with stdout
            vcftools_out, vcftools_err = call_vcftools(vcfname_arg + pop_call_args)

            # Check if the user specifed the complete output filename (only log-based)
            if vcf_args.out:
                produce_vcftools_log(vcftools_err, vcf_args.out, append_mode = True)
            else:
                produce_vcftools_log(vcftools_err, vcftools_output_filename, append_mode = True)

    elif vcf_args.calc_statistic == 'het-fis':

        # Add the stdout argument, for repeated calls
        vcftools_call_args.append('--stdout')

        # Store if the header should be removed
        strip_header = False

        # Loop each population file
        for vcftools_pop_file in vcftools_pop_files:

            # Create the population-specific call
            pop_call_args = vcftools_call_args + ['--keep', str(vcftools_pop_file)]

            # vcftools subprocess call, with stdout
            vcftools_out, vcftools_err = call_vcftools(vcfname_arg + pop_call_args)

            # Check if the user specifed the complete output filename
            if vcf_args.out:
                produce_vcftools_output(vcftools_out, vcf_args.out, append_mode = True, strip_header = strip_header)
                produce_vcftools_log(vcftools_err, vcf_args.out, append_mode = True)
            else:
                produce_vcftools_output(vcftools_out, vcftools_output_filename, append_mode = True, strip_header = strip_header)
                produce_vcftools_log(vcftools_err, vcftools_output_filename, append_mode = True)

            # Will strip the header after the first loop
            strip_header = True

    else:

        # Add the output argument
        vcftools_call_args.extend(['--out', vcf_args.out_prefix])

        # vcftools subprocess call
        vcftools_out, vcftools_err = call_vcftools(vcfname_arg + vcftools_call_args)

        # Check if the user specifed the complete output filename
        if vcf_args.out:
            os.rename(vcftools_output_filename, vcf_args.out)
            produce_vcftools_log(vcftools_err, vcf_args.out)
        else:
            produce_vcftools_log(vcftools_err, vcftools_output_filename)

    # Delete any files that were created for vcftools
    if vcf_args.model_file:
        selected_model.delete_pop_files()
        selected_model.delete_individuals_file()

if __name__ == "__main__":
    initLogger()
    run()
