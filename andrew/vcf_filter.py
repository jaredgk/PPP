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

def vcf_filter_parser(passed_arguments):
    '''VCF Argument Parser - Assigns arguments from command line'''

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

    vcf_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    # Input arguments.
    vcf_parser.add_argument("vcfname", metavar = 'VCF_Input', help = "Input VCF filename", type = str, action = parser_confirm_file())

    # Model file arguments.
    vcf_parser.add_argument('--model-file', help = 'Defines the model file', type = str, action = parser_confirm_file())
    vcf_parser.add_argument('--model', help = 'Defines the model to analyze', type = str)

    # Other file arguments. Expand as needed
    vcf_parser.add_argument('--out', help = 'Specifies the complete output filename', type = str)
    vcf_parser.add_argument('--out-prefix', help = 'Specifies the output prefix (vcftools naming scheme)', type = str,  default = 'out')
    #vcf_parser.add_argument('--log', help = "Specifies if the vcftools log should be saved", action = 'store_false')

    out_format_list = ['vcf', 'vcf.gz', 'bcf', 'removed_sites', 'kept_sites']
    out_format_default = 'removed_sites'

    vcf_parser.add_argument('--out-format', metavar = metavar_list(out_format_list), help = 'Specifies the output format.', type = str, choices = out_format_list, default = out_format_default)

    # General arguments.
    vcf_parser.add_argument('--overwrite', help = "Specifies if previous output files should be overwritten", action = 'store_true')

    ### Filters

    # Allele count filters
    vcf_parser.add_argument('--filter-min-alleles', help = 'Specifies that only sites with a number of allele >= to the number given should be included', type = int,  default = 2)
    vcf_parser.add_argument('--filter-max-alleles', help = 'Specifies that only sites with a number of allele <= to the number given should be included', type = int,  default = 2)

    # Missing data filter
    vcf_parser.add_argument('--filter-max-missing', help = 'Specifies to exclude sites by the proportion of missing data (0.0: include all, 1.0: no missing data)', type = float)

    # Chromosome filters
    vcf_parser.add_argument('--filter-include-chr', help = 'Specifies the chromosome(s) to include', nargs = '+', type = str)
    vcf_parser.add_argument('--filter-exclude-chr', help = 'Specifies the chromosome(s) to exclude', nargs = '+', type = str)

    # Basic position filters
    vcf_parser.add_argument('--filter-from-bp', help = 'Specifies the lower bound of sites to include (May only be used with a single chromosome)', type = int)
    vcf_parser.add_argument('--filter-to-bp', help = 'Specifies the upper bound of sites to include (May only be used with a single chromosome)', type = int)

    # Position-based position filters
    vcf_parser.add_argument('--filter-include-positions', help = 'Specifies a set of sites to include within a file (tsv chromosome and position)', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-positions', help = 'Specifies a set of sites to exclude within a file (tsv chromosome and position)', action = parser_confirm_file())

    # BED-based position filters
    vcf_parser.add_argument('--filter-include-bed', help = 'Specifies a set of sites to include within a BED file', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-bed', help = 'Specifies a set of sites to exclude within a BED file', action = parser_confirm_file())

    # Individual filters
    vcf_parser.add_argument('--filter-keep', help = 'Specifies a file of individuals to keep', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-remove', help = 'Specifies a file of individuals to remove', action = parser_confirm_file())

    # Filter-flag filters
    vcf_parser.add_argument('--filter-include-passed', help = "Specifies that only sites with the filter flag 'PASS' should be included", action = 'store_true')
    vcf_parser.add_argument('--filter-include-filtered', help = 'Specifies that all sites with the given filter flag should be included', nargs = '+', type = str)
    vcf_parser.add_argument('--filter-exclude-filtered', help = 'Specifies that all sites with the given filter flag should be excluded', nargs = '+', type = str)

    # Info-flag filters
    vcf_parser.add_argument('--filter-include-info', help = 'Specifies that all sites with the given info flag should be included', nargs = '+', type = str)
    vcf_parser.add_argument('--filter-exclude-info', help = 'Specifies that all sites with the given info flag should be excluded', nargs = '+', type = str)

    # Additional Filters
    vcf_parser.add_argument('--filter-distance', help = 'Specifies a distance that no two sites may be within', type = int)

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
    Filter VCF files using VCFTools.

    Automates various filters to VCF files using VCFtools.

    Parameters
    ----------
    VCF_Input : str
        Specifies the input VCF filename
    --out : str
        Specifies the output filename
    --out-format : str
        Specifies the output format {vcf, bcf, removed_sites, kept_sites} (Default: removed_sites)
    --filter-include-chr : list or str
        Specifies the chromosome(s) to include
    --filter-exclude-chr : list or str
        Specifies the chromosome(s) to exclude
    --filter-from-bp : int
        Specifies the lower bound of sites to include. May only be used with a single chromosome
    --filter-to-bp : int
        Specifies the upper bound of sites to include. May only be used with a single chromosome
    --filter-include-positions : str
        Specifies a set of sites to include within a tsv file (chromosome and position)
    --filter-exclude-positions : str
        Specifies a set of sites to exclude within a tsv file (chromosome and position)
    --filter-include-bed : str
        Specifies a set of sites to include within a BED file
    --filter-exclude-bed : str
        Specifies a set of sites to exclude within a BED file
    --filter-include-passed : bool
        Specifies that only sites with the filter flag 'PASS' should be included (Default: False)
    --filter-include-filtered : list or str
        Specifies that all sites with the given filter flag should be included
    --filter-exclude-filtered : list or str
        Specifies that all sites with the given filter flag should be excluded
    --filter-include-info : list or str
        Specifies that all sites with the given info flag should be included
    --filter-exclude-info : list or str
        Specifies that all sites with the given info flag should be excluded
    --filter-min-alleles : int
        Specifies that only sites with a number of allele >= to the number given should be included (Default: 2)
    --filter-min-alleles : int
        Specifies that only sites with a number of allele <= to the number given should be included (Default: 2)

    Returns
    -------
    output : file
        Filtered file output
    log : file
        Log file output

    Raises
    ------
    IOError
        Input VCF file does not exist
    IOError
        Output file already exists and --overwrite is not specified

    '''

    # Grab VCF arguments from command line
    vcf_args = vcf_filter_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(vcf_args, 'vcf_filter')

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

        # Create individuals file
        selected_model.create_individuals_file(overwrite = vcf_args.overwrite)

        # Assign the individuals file to vcftools
        vcftools_call_args.extend(['--keep', selected_model.individuals_file])

    # Holds the filename vcftools assigns to the filtered output
    vcftools_output_filename = None

    # Used to assign the output format to the vcftools call, and assign vcftools output filename
    if vcf_args.out_format == 'removed_sites':
        vcftools_call_args.append('--removed-sites')
        vcftools_output_filename = vcf_args.out_prefix + '.removed.sites'
    elif vcf_args.out_format == 'kept_sites':
        vcftools_call_args.append('--kept-sites')
        vcftools_output_filename = vcf_args.out_prefix + '.kept.sites'
    elif vcf_args.out_format == 'bcf':
        vcftools_call_args.append('--recode-bcf')
        vcftools_output_filename = vcf_args.out_prefix + '.recode.bcf'
    elif vcf_args.out_format == 'vcf':
        vcftools_call_args.append('--recode')
        vcftools_output_filename = vcf_args.out_prefix + '.recode.vcf'
    elif vcf_args.out_format == 'vcf.gz':
        vcftools_call_args.append('--recode')
        vcftools_output_filename = vcf_args.out_prefix + '.recode.vcf.gz'

    # Individual-based filters
    if vcf_args.filter_keep or vcf_args.filter_remove:
        # Used to include a file of individuals to keep
        if vcf_args.filter_keep:
            vcftools_call_args.extend(['--keep', vcf_args.filter_keep])

        # Used to include a file of individuals to remove
        if vcf_args.filter_remove:
            vcftools_call_args.extend(['--remove', vcf_args.filter_remove])

    # Chromosome-based filters
    if vcf_args.filter_include_chr or vcf_args.filter_exclude_chr:
        if vcf_args.filter_include_chr:
            for chr_to_include in vcf_args.filter_include_chr:
                vcftools_call_args.extend(['--chr', chr_to_include])
        if vcf_args.filter_exclude_chr:
            for chr_to_exclude in vcf_args.filter_exclude_chr:
                vcftools_call_args.extend(['--not-chr', chr_to_exclude])

    # Site (i.e. basepair) filters
    if vcf_args.filter_from_bp or vcf_args.filter_to_bp:
        if vcf_args.filter_include_chr:
            vcftools_call_args.extend(['--from-bp', vcf_args.filter_from_bp])
        if vcf_args.filter_exclude_chr:
            vcftools_call_args.extend(['--to-bp', vcf_args.filter_to_bp])

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

    # Flag-based Filters
    if vcf_args.filter_include_passed or vcf_args.filter_include_filtered or vcf_args.filter_exclude_filtered:
        if vcf_args.filter_include_passed:
            vcftools_call_args.append('--remove-filtered-all')
        if vcf_args.filter_include_filtered:
            for filtered_to_include in vcf_args.filter_include_filtered:
                vcftools_call_args.extend(['--keep-filtered', filtered_to_include])
        if vcf_args.filter_exclude_filtered:
            for filtered_to_exclude in vcf_args.filter_exclude_filtered:
                vcftools_call_args.extend(['--remove-filtered', filtered_to_exclude])

    # Infor-based filters
    if vcf_args.filter_include_info or vcf_args.filter_exclude_info:
        if vcf_args.filter_include_info:
            for info_to_include in vcf_args.filter_include_info:
                vcftools_call_args.extend(['--keep-INFO', info_to_include])
        if vcf_args.filter_exclude_info:
            for info_to_exclude in vcf_args.filter_exclude_info:
                vcftools_call_args.extend(['--remove-INFO', info_to_exclude])

    # Allele-based filters
    if vcf_args.filter_min_alleles or vcf_args.filter_max_alleles:
        if vcf_args.filter_min_alleles:
            vcftools_call_args.extend(['--min-alleles', vcf_args.filter_min_alleles])
        if vcf_args.filter_max_alleles:
            vcftools_call_args.extend(['--max-alleles', vcf_args.filter_max_alleles])

    # Missing data filters
    if vcf_args.filter_max_missing:
        vcftools_call_args.extend(['--max-missing', vcf_args.filter_max_missing])

    # Distance (between sites) filters
    if vcf_args.filter_distance:
        vcftools_call_args.extend(['--thin', vcf_args.filter_distance])

    logging.info('vcftools parameters assigned')

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

    # Check if the user has requested vcf.gz, if so send the stdout to bgzip
    if vcf_args.out_format == 'vcf.gz':
        # Call both vcftools and bgzip, return stderr
        vcftools_err = call_vcftools_bgzip(vcfname_arg + vcftools_call_args, vcftools_output_filename)
    else:
        # Call only vcftools. Note: vcftools_out will be empty
        vcftools_out, vcftools_err = call_vcftools(vcfname_arg + vcftools_call_args)

    # Check if the user specifed the complete output filename
    if vcf_args.out:
        os.rename(vcftools_output_filename, vcf_args.out)
        produce_vcftools_log(vcftools_err, vcf_args.out)
    else:
        produce_vcftools_log(vcftools_err, vcftools_output_filename)

    # Delete any files that were created for vcftools
    if vcf_args.model_file:
        selected_model.delete_individuals_file()

if __name__ == "__main__":
    initLogger()
    run()
