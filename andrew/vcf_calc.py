import os
import sys
import subprocess
import argparse
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
    vcf_parser.add_argument('--out', help = 'Specifies the complete output filename. Renames vcftools-named intermediate files', type = str)
    vcf_parser.add_argument('--out-prefix', help = 'Specifies the output prefix (vcftools naming scheme)', type = str,  default = 'out')
    #vcf_parser.add_argument('--log', help = "Specifies if the vcftools log should be saved", action = 'store_false')

    # General arguments.
    vcf_parser.add_argument('--overwrite', help = "Specifies if previous output files should be overwritten", action = 'store_true')

    # Statistic based arguments.
    statistic_list = ['weir-fst', 'windowed-weir-fst', 'TajimaD', 'pi', 'freq', 'het']
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
    vcftools_call_args = ['--out', vcf_args.out_prefix]

    # Check if the user has specified a model file
    if vcf_args.model_file and vcf_args.model:

        # Read in the models
        models_in_file = read_model_file(vcf_args.model_file)

        # Check that the selected model was not found in the file
        if vcf_args.model not in models_in_file:
            logging.error('Selected model "%s" not found in: %s' % (vcf_args.model, vcf_args.model_file))
            raise IOError('Selected model "%s" not found in: %s' % (vcf_args.model, vcf_args.model_file))

        # Select model, might change this in future versions
        selected_model = models_in_file[vcf_args.model]

        # Both fst statistics do not require this step.
        if vcf_args.calc_statistic not in ['windowed-weir-fst', 'weir-fst']:

            # Create individuals file
            selected_model.create_individuals_file(overwrite = vcf_args.overwrite)

            # Assign the individuals file to vcftools
            vcftools_call_args.extend(['--keep', selected_model.individuals_file])


    if vcf_args.calc_statistic in ['windowed-weir-fst', 'weir-fst']:

        # Confirms that the user has specified a population assignment method
        if not vcf_args.pop_file and not vcf_args.model_file:
            logging.error('windowed-weir-fst requires either: --model/--model-file or --pop-file')
            raise IOError('windowed-weir-fst requires either: --model/--model-file or --pop-file')

        # Confirms that at least two population files have been specified
        if vcf_args.pop_file and len(vcf_args.pop_file) < 2:
            logging.error('Two or more population files requried. Please assign using --pop-file')
            raise IOError('Two or more population files requried. Please assign using --pop-file')

        # Check if the user specified a model file
        if vcf_args.model_file:

            # Create the population files
            selected_model.create_population_files(file_ext = '.txt', overwrite = vcf_args.overwrite)

            # Assigns the population files
            vcftools_pop_args = [population_args for population_file in selected_model.population_files for population_args in ['--weir-fst-pop', population_file]]

        # Check if the user specified a pop file
        elif vcf_args.pop_file:

            # Check that the population files exist
            for population_file in vcf_args.pop_file:
                if not os.path.isfile(population_file):
                    logging.error('Population file: %s. Not found.' % population_file)
                    raise IOError('Population file: %s. Not found.' % population_file)

            # Assigns the population files
            vcftools_pop_args = [population_args for population_file in vcf_args.pop_file for population_args in ['--weir-fst-pop', population_file]]


        if vcf_args.calc_statistic == 'windowed-weir-fst':

            # Assignsspecific vcftools arguments for calculating fst
            vcftools_window_args = ['--fst-window-size', vcf_args.statistic_window_size, '--fst-window-step', vcf_args.statistic_window_step]

            # Assigns all the vcftools arguments for calculating windowed fst
            vcftools_call_args.extend(vcftools_pop_args + vcftools_window_args)

            # Assigns the suffix for the vcftools log file
            vcftools_log_suffix = 'windowed.weir.fst'

        elif vcf_args.calc_statistic == 'weir-fst':
            # Assigns specific vcftools arguments for calculating site-based fst
            vcftools_call_args.extend(vcftools_pop_args)

            # Assigns the suffix for the vcftools log file
            vcftools_log_suffix = 'weir.fst'

    elif vcf_args.calc_statistic == 'TajimaD':

        # Assigns all the vcftools arguments for calculating TajimaD
        vcftools_call_args.extend(['--TajimaD', vcf_args.statistic_window_size])

        # Assigns the suffix for the vcftools log file
        vcftools_log_suffix = 'Tajima.D'

    elif vcf_args.calc_statistic == 'pi':

        # Assigns all the vcftools arguments for calculating pi
        vcftools_call_args.extend(['--window-pi', vcf_args.statistic_window_size, '--window-pi-step', vcf_args.statistic_window_step])

        # Assigns the suffix for the vcftools log file
        vcftools_log_suffix = 'windowed.pi'

    elif vcf_args.calc_statistic == 'freq':

        # Assigns all the vcftools arguments for the allele frequency
        vcftools_call_args.extend(['--freq'])

        # Assigns the suffix for the vcftools log file
        vcftools_log_suffix = 'frq'

    elif vcf_args.calc_statistic == 'het':

        # Assigns all the vcftools arguments for calculating heterozygosity
        vcftools_call_args.extend(['--het'])

        # Assigns the suffix for the vcftools log file
        vcftools_log_suffix = 'het'

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
    vcftools_output_filename = vcf_args.out_prefix + '.' + vcftools_log_suffix

    # Check if previous output should be overwritten
    if not vcf_args.overwrite:
        if vcf_args.out:
            # Confirm the vcftools-renamed output and log file do not exist
            check_for_vcftools_output(vcf_args.out)
        else:
            # Confirm the vcftools output and log file do not exist
            check_for_vcftools_output(vcftools_output_filename)

    # Assigns the file argument for vcftools
    vcfname_arg = assign_vcftools_input_arg(vcf_args.vcfname)
    logging.info('Input file assigned')

    # vcftools subprocess call
    vcftools_err = call_vcftools(vcfname_arg + vcftools_call_args)

    # Check if the user specifed the complete output filename
    if vcf_args.out:
        os.rename(vcftools_output_filename, vcf_args.out)
        produce_vcftools_log(vcftools_err, vcf_args.out)
    else:
        produce_vcftools_log(vcftools_err, vcftools_output_filename)

if __name__ == "__main__":
    initLogger()
    run()
