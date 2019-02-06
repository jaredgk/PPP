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
sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'pppipe')))

from logging_module import initLogger, logArgs

def vcf_filter_parser(passed_arguments):
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

    # Other file arguments. Expand as needed
    vcf_parser.add_argument('--out', help = 'Specifies the complete output filename', type = str)
    vcf_parser.add_argument('--out-prefix', help = 'Specifies the output prefix (vcftools naming scheme)', type = str,  default = 'out')
    #vcf_parser.add_argument('--log', help = "Specifies if the vcftools log should be saved", action = 'store_false')

    out_format_list = ['vcf', 'vcf.gz', 'bcf', 'removed_sites', 'kept_sites', 'removed_bed', 'kept_bed']
    out_format_default = 'removed_sites'

    vcf_parser.add_argument('--out-format', metavar = metavar_list(out_format_list), help = 'Specifies the output format.', type = str, choices = out_format_list, default = out_format_default)

    # General arguments.
    vcf_parser.add_argument('--overwrite', help = "overwrite previous output files", action = 'store_true')

    ### Filters

    # Non-model file arguments
    vcf_parser.add_argument('--filter-include-indv', help = 'Individual to include. Note: This argument may be used multiple times and cannont to be used alonside --model', nargs = '+', type = str, action = parser_add_to_list())
    vcf_parser.add_argument('--filter-exclude-indv', help = 'Individual to exclude. Note: This argument may be used multiple times and cannont to be used alonside --model', nargs = '+', type = str, action = parser_add_to_list())
    vcf_parser.add_argument('--filter-include-indv-file', help = 'Individuals to include in file. Cannont to be used alonside --model', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-indv-file', help = 'Individuals to exclude in file. Cannont to be used alonside --model', action = parser_confirm_file())

    # Allele count filters
    vcf_parser.add_argument('--filter-min-alleles', help = 'Only sites with a number of allele >= to the number given should be included', type = int)
    vcf_parser.add_argument('--filter-max-alleles', help = 'Only sites with a number of allele <= to the number given should be included', type = int)

    # Missing data filter
    vcf_parser.add_argument('--filter-max-missing', help = 'Exclude sites by the proportion of missing data (0.0: include all, 1.0: no missing data)', type = float)

    # Indel filters
    indel_filters = vcf_parser.add_mutually_exclusive_group()
    indel_filters.add_argument('--filter-include-only-indels', help = 'Include only sites that contain an indel', action = 'store_true')
    indel_filters.add_argument('--filter-exclude-indels', help = 'Exclude sites that contain an indel', action = 'store_true')

    # Chromosome filters
    vcf_parser.add_argument('--filter-include-chr', help = 'Chromosome to include. Note: This argument may be used multiple times', nargs = '+', type = str, action = parser_add_to_list())
    vcf_parser.add_argument('--filter-exclude-chr', help = 'Chromosome to exclude. Note: This argument may be used multiple times', nargs = '+', type = str, action = parser_add_to_list())

    # Basic position filters
    vcf_parser.add_argument('--filter-from-bp', help = 'Lower bound of sites to include (Only usable with a single chromosome)', type = int)
    vcf_parser.add_argument('--filter-to-bp', help = 'Upper bound of sites to include (Only usable with a single chromosome)', type = int)

    # Position-based position filters
    vcf_parser.add_argument('--filter-include-positions', help = 'Sites to include within a chr/pos tab-sperated file', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-positions', help = 'Sites to exclude within a chr/pos tab-sperated file', action = parser_confirm_file())

    # BED-based position filters
    vcf_parser.add_argument('--filter-include-bed', help = 'Set of sites to include within a BED file', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-bed', help = 'Set of sites to exclude within a BED file', action = parser_confirm_file())

    # Filter-flag filters
    vcf_parser.add_argument('--filter-include-passed', help = "Include only sites with the filter flag 'PASS'", action = 'store_true')
    vcf_parser.add_argument('--filter-include-flag', help = 'Include all sites with the given filter flag', nargs = '+', type = str, action = parser_add_to_list())
    vcf_parser.add_argument('--filter-exclude-flag', help = 'Exclude all sites with the given filter flag', nargs = '+', type = str, action = parser_add_to_list())

    # Info-flag filters
    vcf_parser.add_argument('--filter-include-info', help = 'Include all sites with the given info flag', nargs = '+', type = str, action = parser_add_to_list())
    vcf_parser.add_argument('--filter-exclude-info', help = 'Exclude all sites with the given info flag', nargs = '+', type = str, action = parser_add_to_list())

    # SNP filters
    vcf_parser.add_argument('--filter-include-snp', help = 'Include SNP(s) with matching ID. Note: This argument may be used multiple times.', nargs = '+', type = str, action = parser_add_to_list())
    vcf_parser.add_argument('--filter-include-snps', help = 'Include all SNPs with matching ID within a file.', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-snps', help = 'Include all SNPs with matching ID within a file.', action = parser_confirm_file())

    # MAF Filter
    vcf_parser.add_argument('--filter-maf-min', help = 'Include sites with equal or greater MAF values', type = float)
    vcf_parser.add_argument('--filter-maf-max', help = 'Include sites with equal or lesser MAF values', type = float)

    # MAC Filter
    vcf_parser.add_argument('--filter-mac-min', help = 'Include sites with equal or greater MAC values', type = int)
    vcf_parser.add_argument('--filter-mac-max', help = 'Include sites with equal or lesser MAC values', type = int)

    # Additional Filters
    vcf_parser.add_argument('--filter-distance', help = 'Specifies a distance that no two sites may be within', type = int)

    if passed_arguments:
        return vcf_parser.parse_args(passed_arguments)
    else:
        return vcf_parser.parse_args()

def run (passed_arguments = []):
    '''
    Filter VCF files using VCFTools.

    Automates various filters to VCF files using VCFtools.

    Parameters
    ----------
    --vcf : str
        Specifies the input VCF filename
    --out : str
        Specifies the output filename
    --out-format : str
        Specifies the output format {vcf, bcf, removed_sites, kept_sites}
        (Default: removed_sites)
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
    logArgs(vcf_args, func_name = 'vcf_filter')

    # Argument container for vcftools
    vcftools_call_args = []

    # Check if the user has specified a model file
    if vcf_args.model_file and vcf_args.model:
        # Read in the models
        models_in_file = read_model_file(vcf_args.model_file)

        # Check that the selected model was not found in the file
        if vcf_args.model not in models_in_file:
            raise IOError('Selected model "%s" not found in: %s' % (vcf_args.model, vcf_args.model_file))

        # Check that individual-based filters are not also being used
        if vcf_args.filter_include_indv or vcf_args.filter_exclude_indv:
            if vcf_args.filter_include_indv:
                raise Exception('--model and --filter-include-indv arguments are incompatible')
            if vcf_args.filter_exclude_indv:
                raise Exception('--model and --filter-exclude-indv arguments are incompatible')

        # Check that individuals-based filters are not also being used
        if vcf_args.filter_include_indv_file or vcf_args.filter_exclude_indv_file:
            if vcf_args.filter_include_indv_file:
                raise Exception('--model and --filter_include_indv-file arguments are incompatible')
            if vcf_args.filter_exclude_indv_file:
                raise Exception('--model and --filter-exclude-indv-file arguments are incompatible')

        # Select model, might change this in future versions
        selected_model = models_in_file[vcf_args.model]

        # Create individuals file
        selected_model.create_ind_file(overwrite = vcf_args.overwrite)

        # Assign the individuals file to vcftools
        vcftools_call_args.extend(['--keep', selected_model.ind_file])

    # Holds the filename suffix vcftools assigns to the filtered output
    vcftools_out_suffix = None

    # Used to assign the output format to the vcftools call, and assign vcftools output suffix
    if vcf_args.out_format == 'removed_sites':
        vcftools_call_args.append('--removed-sites')
        vcftools_out_suffix = '.removed.sites'
    elif vcf_args.out_format == 'kept_sites':
        vcftools_call_args.append('--kept-sites')
        vcftools_out_suffix = '.kept.sites'
    elif vcf_args.out_format == 'removed_bed':
        vcftools_call_args.append('--removed-sites')
        vcftools_out_suffix = '.removed.bed'
    elif vcf_args.out_format == 'kept_bed':
        vcftools_call_args.append('--kept-sites')
        vcftools_out_suffix = '.kept.bed'
    elif vcf_args.out_format == 'bcf':
        vcftools_call_args.append('--recode')
        vcftools_out_suffix = '.recode.bcf'
    elif vcf_args.out_format == 'vcf':
        vcftools_call_args.append('--recode')
        vcftools_out_suffix = '.recode.vcf'
    elif vcf_args.out_format == 'vcf.gz':
        vcftools_call_args.append('--recode')
        vcftools_out_suffix = '.recode.vcf.gz'

    # Assign expected vcftools output filename
    vcftools_output_filename = vcf_args.out_prefix + vcftools_out_suffix

    # Check if the user has specified a output filename
    if vcf_args.out:

        # Assign the vcftools output filename, using the output filename
        vcftools_output_filename = vcf_args.out

    else:

        # Assign the vcftools output filename, using the prefix and suffix
        vcftools_output_filename = vcf_args.out_prefix + vcftools_out_suffix

    # Check if previous output should be overwritten
    if not vcf_args.overwrite:

        # Confirm the vcftools output and log file do not exist
        check_for_vcftools_output(vcftools_output_filename)

    # Individuals-based filters
    if vcf_args.filter_include_indv_file or vcf_args.filter_exclude_indv_file:
        # Used to include a file of individuals to keep
        if vcf_args.filter_include_indv_file:
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

    # Indel-based filters
    if vcf_args.filter_include_only_indels or vcf_args.filter_exclude_indels:
        if vcf_args.filter_include_only_indels:
            vcftools_call_args.append('--keep-only-indels')
        if vcf_args.filter_exclude_indels:
            vcftools_call_args.append('--remove-indels')

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
    if vcf_args.filter_include_passed or vcf_args.filter_include_flag or vcf_args.filter_exclude_flag:
        if vcf_args.filter_include_passed:
            vcftools_call_args.append('--remove-filtered-all')
        if vcf_args.filter_include_flag:
            for filtered_to_include in vcf_args.filter_include_flag:
                vcftools_call_args.extend(['--keep-filtered', filtered_to_include])
        if vcf_args.filter_exclude_flag:
            for filtered_to_exclude in vcf_args.filter_exclude_flag:
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

    # SNP-based filters
    if vcf_args.filter_include_snp or vcf_args.filter_include_snps or vcf_args.filter_exclude_snps:
        if vcf_args.filter_include_snp:
            for snp_to_include in vcf_args.filter_include_snp:
                vcftools_call_args.extend(['--snp', snp_to_include])
        if vcf_args.filter_include_snps:
            vcftools_call_args.extend(['--snps', vcf_args.filter_include_snps])
        if vcf_args.filter_exclude_snps:
            vcftools_call_args.extend(['--exclude', vcf_args.filter_exclude_snps])

    # MAF-based filters
    if vcf_args.filter_maf_min or vcf_args.filter_maf_max:
        if vcf_args.filter_maf_min:
            vcftools_call_args.extend(['--maf', vcf_args.filter_maf_min])
        if vcf_args.filter_maf_max:
            vcftools_call_args.extend(['--max-maf', vcf_args.filter_maf_max])

    # MAC-based filters
    if vcf_args.filter_mac_min or vcf_args.filter_mac_max:
        if vcf_args.filter_mac_min:
            vcftools_call_args.extend(['--mac', vcf_args.filter_mac_min])
        if vcf_args.filter_mac_max:
            vcftools_call_args.extend(['--max-mac', vcf_args.filter_mac_max])

    # Distance (between sites) filters
    if vcf_args.filter_distance:
        vcftools_call_args.extend(['--thin', vcf_args.filter_distance])

    logging.info('vcftools parameters assigned')

    # Assigns the file argument for vcftools
    vcfname_arg = assign_vcftools_input_arg(vcf_args.vcf)
    logging.info('Input file assigned')

    # Check if the output format is vcf
    if vcf_args.out_format == 'vcf':

        # Call vcftools with the specifed arguments
        vcftools_err = call_vcftools(vcfname_arg + vcftools_call_args, output_format = vcf_args.out_format, output_filename = vcftools_output_filename)

    else:

        # Call vcftools with the specifed arguments
        vcftools_err = call_vcftools(vcfname_arg + vcftools_call_args, output_format = vcf_args.out_format, output_filename = vcftools_output_filename)

    produce_vcftools_log(vcftools_err, vcftools_output_filename)

    # Delete any files that were created for vcftools
    if vcf_args.model_file:
        selected_model.delete_ind_file()

if __name__ == "__main__":
    initLogger()
    run()
