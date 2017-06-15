import os
import sys
import subprocess
import argparse
import logging

from vcftools import *

# Insert Jared's directory path, required for calling Jared's functions. Change when directory structure changes.
sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'jared')))

#from logging_module import initLogger

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

    def parser_confirm_no_file ():
        '''Custom action to confirm file does not exist'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if os.path.isfile(value):
                    raise IOError # File found
                setattr(args, self.dest, value)
        return customAction

    def parser_confirm_files ():
        '''Custom action to confirm multiple file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError # File not found
                getattr(args, self.dest).append(value)
        return customAction

    vcf_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    # Input arguments.
    vcf_parser.add_argument("vcfname", metavar='VCF_Input', help = "Input VCF filename", type = str, action = parser_confirm_file())

    # Other file arguments. Expand as needed
    vcf_parser.add_argument('--out', help = 'Specifies the output filename', type = str,  default = 'out', action = parser_confirm_no_file())

    out_format_list = ['vcf', 'bcf']
    out_format_default = 'bcf'

    vcf_parser.add_argument('--out-format', metavar = '{' + ', '.join(out_format_list) + '}', help = 'Specifies the output format', type = str, choices = out_format_list, default = out_format_default)
    ### Filters

    # Chromosome filters
    vcf_parser.add_argument('--filter-include-chr', help = 'Specifies the chromosome(s) to include', nargs = '+', type = str)
    vcf_parser.add_argument('--filter-exclude-chr', help = 'Specifies the chromosome(s) to exclude', nargs = '+', type = str)

    # Basic position filters
    vcf_parser.add_argument('--filter-from-bp', help = 'Specifies the lower bound of sites to include (May only be used with a single chromosome)', type = int)
    vcf_parser.add_argument('--filter-to-bp', help = 'Specifies the upper bound of sites to include (May only be used with a single chromosome)', type = int)

    # BED-based position filters
    vcf_parser.add_argument('--filter-include-bed', help = 'Specifies a set of sites to include within a BED file', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-bed', help = 'Specifies a set of sites to exclude within a BED file', action = parser_confirm_file())

    # Filter-flag filters
    vcf_parser.add_argument('--filter-include-passed', help = "Specifies that only sites with the filter flag 'PASS' should be included", action = 'store_true')
    vcf_parser.add_argument('--filter-include-filtered', help = 'Specifies that all sites with the given filter flag should be included', nargs = '+', type = str)
    vcf_parser.add_argument('--filter-exclude-filtered', help = 'Specifies that all sites with the given filter flag should be excluded', nargs = '+', type = str)

    # Info-flag filters
    vcf_parser.add_argument('--filter-include-info', help = 'Specifies that all sites with the given info flag should be included', nargs = '+', type = str)
    vcf_parser.add_argument('--filter-exclude-info', help = 'Specifies that all sites with the given info flag should be excluded', nargs = '+', type = str)

    # Allele count filters
    vcf_parser.add_argument('--filter-min-alleles', help = 'Specifies that only sites with a number of allele >= to the number given should be included', type = int)
    vcf_parser.add_argument('--filter-max-alleles', help = 'Specifies that only sites with a number of allele <= to the number given should be included', type = int)

    # Missing data filter
    vcf_parser.add_argument('--filter-max-missing', help = 'Specifies that only sites with more than this number of genotypes among individuals should be included', type = int)

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
    VCF_Input : str
        Specifies the input VCF filename
    --out : str
        Specifies the output filename

    Returns
    -------
    output : file
        Filtered VCF file output
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
    vcf_args = vcf_filter_parser(passed_arguments)

    # Argument container for vcftools
    vcftools_call_args = ['--out', vcf_args.out]

    if vcf_args.out_format:
        if vcf_args.out_format == 'bcf':
            vcftools_call_args.append('--recode-bcf')
        elif vcf_args.out_format == 'vcf':
            vcftools_call_args.append('--recode')

    if vcf_args.filter_include_chr or vcf_args.filter_exclude_chr:
        if vcf_args.filter_include_chr:
            for chr_to_include in vcf_args.filter_include_chr:
                vcftools_call_args.extend(['--chr', chr_to_include])
        if vcf_args.filter_exclude_chr:
            for chr_to_exclude in vcf_args.filter_exclude_chr:
                vcftools_call_args.extend(['--not-chr', chr_to_exclude])

    if vcf_args.filter_from_bp or vcf_args.filter_to_bp:
        if vcf_args.filter_include_chr:
            vcftools_call_args.extend(['--from-bp', vcf_args.filter_from_bp])
        if vcf_args.filter_exclude_chr:
            vcftools_call_args.extend(['--to-bp', vcf_args.filter_to_bp])

    if vcf_args.filter_include_bed or vcf_args.filter_exclude_bed:
        if vcf_args.filter_include_bed:
            vcftools_call_args.extend(['--bed', vcf_args.filter_include_bed])
        if vcf_args.filter_exclude_bed:
            vcftools_call_args.extend(['--exclude-bed', vcf_args.filter_exclude_bed])

    if vcf_args.filter_include_passed or vcf_args.filter_include_filtered or vcf_args.filter_exclude_filtered:
        if vcf_args.filter_include_passed:
            vcftools_call_args.append('--remove-filtered-all')
        if vcf_args.filter_include_filtered:
            for filtered_to_include in vcf_args.filter_include_filtered:
                vcftools_call_args.extend(['--keep-filtered', filtered_to_include])
        if vcf_args.filter_exclude_filtered:
            for filtered_to_exclude in vcf_args.filter_exclude_filtered:
                vcftools_call_args.extend(['--remove-filtered', filtered_to_exclude])

    if vcf_args.filter_include_info or vcf_args.filter_exclude_info:
        if vcf_args.filter_include_info:
            for info_to_include in vcf_args.filter_include_info:
                vcftools_call_args.extend(['--keep-INFO', info_to_include])
        if vcf_args.filter_exclude_info:
            for info_to_exclude in vcf_args.filter_exclude_info:
                vcftools_call_args.extend(['--remove-INFO', info_to_exclude])

    if vcf_args.filter_min_alleles or vcf_args.filter_max_alleles:
        if vcf_args.filter_min_alleles:
            vcftools_call_args.extend(['--min-alleles', vcf_args.filter_min_alleles])
        if vcf_args.filter_max_alleles:
            vcftools_call_args.extend(['--max-alleles', vcf_args.filter_max_alleles])

    if vcf_args.filter_max_missing:
        vcftools_call_args.extend(['--max-missing-count', vcf_args.filter_max_missing])

    if vcf_args.filter_distance:
        vcftools_call_args.extend(['--thin', vcf_args.filter_distance])

    # Assigns the file argument for vcftools
    vcfname_arg = assign_vcftools_input_arg(vcf_args.vcfname)

    # vcftools subprocess call
    vcftools_call = subprocess.Popen(['vcftools'] + vcfname_arg + list(map(str, vcftools_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    vcftools_out, vcftools_err = vcftools_call.communicate()

    # Check that the log file was created correctly, get the suffix for the log file, and create the file
    if check_vcftools_for_errors(vcftools_err):
        produce_vcftools_log(vcftools_err, vcf_args.out, '.filter')


if __name__ == "__main__":
    #initLogger()
    run()
