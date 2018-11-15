import os
import sys
import argparse
import itertools
import copy
import shutil
import logging
import pandas as pd

# Import basic vcftools functions
from bcftools import *

# Insert Jared's directory path, required for calling Jared's functions. Change when directory structure changes.
sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'jared')))

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
    vcf_parser.add_argument('--vcf', help = "Input VCF filename (may be used multiple times)", type = str, nargs = '+', required = True, action = parser_confirm_file_list())

    # Model file arguments - Add later if needed
    #vcf_parser.add_argument('--model-file', help = 'The model filename', type = str, action = parser_confirm_file())
    #vcf_parser.add_argument('--model', help = 'Model to analyze', type = str)

    # Utility based arguments
    utility_list = ['sample-list', 'chr-list', 'concatenate', 'merge', 'sort']
    vcf_parser.add_argument('--utility', metavar = metavar_list(utility_list), help = 'The utility to use', type = str, choices = utility_list)

    # Record merging argument
    merge_list = ['none', 'snps', 'indels', 'both', 'all', 'id']
    vcf_parser.add_argument('--record-merge-mode', metavar = metavar_list(merge_list), help = 'Types of multiallelic records to create. Only usable with the merge and concatenate utilites', type = str, choices = merge_list)

    # Missing Record argument
    vcf_parser.add_argument('--record-missing-as-ref', help = 'Convert missing records as the default. Only usable with the merge and concatenate utilites', action = 'store_true')

    # Output file arguments
    out_format_list = ['vcf', 'vcf.gz', 'bcf']
    out_format_default = 'vcf.gz'

    vcf_parser.add_argument('--out-format', metavar = metavar_list(out_format_list), help = 'Specifies the format of vcf-based output', type = str, choices = out_format_list, default = out_format_default)

    vcf_parser.add_argument('--out', help = 'Output filename. If used, overrides --out-prefix', type = str)
    vcf_parser.add_argument('--out-prefix', help = 'Defines the output prefix', default = 'out')
    vcf_parser.add_argument('--overwrite', help = "Overwrite previous output files", action = 'store_true')

    if passed_arguments:
        return vcf_parser.parse_args(passed_arguments)
    else:
        return vcf_parser.parse_args()

def utility_tsv_output (utility_results, utility_header, utility_index, utility_output_filename):

    # Check if the utility didn't produced results
    if not utility_results:
        raise Exception('Utility (%s) produced no results' % vcf_args.utility)

    # Check if index values were assigned, and if so, create file with index
    if utility_index:

        # Save the utility information into a dataframe
        utility_dataframe = pd.DataFrame(utility_results, index = utility_index, columns = utility_header)

        # Save the utility dataframe into a file
        utility_dataframe.to_csv(utility_output_filename, sep = '\t')

    # Create file without index
    else:

        # Save the utility information into a dataframe
        utility_dataframe = pd.DataFrame(utility_results, columns = utility_header)

        # Save the utility dataframe into a file
        utility_dataframe.to_csv(utility_output_filename, sep = '\t', index = False)

def run (passed_arguments = []):
    '''
    Utilites for VCF-formatted files

    Automates various utilites for VCF-formatted files. This currently includes:
    obtain list of chromosomes and obtain list of samples.

    Parameters
    ----------
    --vcf : str
        Specifies the input VCF filename
    --utility : str
        Specifies the utility to be used. Choices: chr-list,
        sample-list (Default)
    --out-prefix : str
        Specifies the output filename prefix
    --out : str
        Specifies the output filename. Overrides --out-prefix
    --overwrite
        Species if previous output should be overwriten

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
    logArgs(vcf_args, func_name = 'vcf_utilities')

    # Check if the function does not support multiple vcf files
    if vcf_args.utility in ['chr-list', 'sample-list', 'sort']:

        # Check if multiple vcf files have been assigned
        if len(vcf_args.vcf) != 1:
            raise Exception('Utility (%s) only supports a single VCF file' % vcf_args.utility)

        # Convert the VCF argument to a string if no error is found
        vcf_args.vcf = vcf_args.vcf[0]

    # Check if the function does not support a single vcf file
    if vcf_args.utility in ['merge', 'concatenate']:

        # Check if multiple vcf files have been assigned
        if len(vcf_args.vcf) <= 1:
            raise Exception('Utility (%s) only supports multiple VCF files' % vcf_args.utility)

    # Assign the standard output filename
    utility_output_filename = vcf_args.out_prefix

    # Check if the user has requested a chromosome list
    if vcf_args.utility == 'chr-list':

        # Assign the suffix
        utility_output_filename += '.chr.list'

    # Check if the user has requested a sample list
    elif vcf_args.utility == 'sample-list':

        # Assign the suffix
        utility_output_filename += '.sample.list'

    # Check if the user has requested to merge or concatenate VCFs
    elif vcf_args.utility in ['merge', 'concatenate', 'sort']:

        # Assign the suffix
        utility_output_filename += '.%s' % vcf_args.out_format

    # Check if the user specified a output filename
    if vcf_args.out:

        # Override the out_prefix output filename
        utility_output_filename = vcf_args.out

    # Check if previous output should not be overwriten
    if not vcf_args.overwrite:

        # Check if previous output exists
        if os.path.exists(utility_output_filename):
            raise IOError('Utility output already exists')

    # List of hold optional utility-based arguments
    utility_optional_args = []

    if vcf_args.utility in ['merge', 'concatenate']:

        # Check if the user has requested to convert missing data into reference alleles
        if vcf_args.record_missing_as_ref:

            # Assign the optional argument
            utility_optional_args.append('--missing-to-ref')

    if vcf_args.utility == 'merge':

        # Check if the user has selected a merge mode
        if vcf_args.record_merge_mode:

            # Assign the optional argument
            utility_optional_args.extend(['--merge', vcf_args.record_merge_mode])

    # List of hold utility results
    utility_results = []

    # List of the utility index
    utility_index = []

    # Str of the utility header
    utility_header = ''

    # Check if the user has requested a chromosome list
    if vcf_args.utility == 'chr-list':

        # Use bcftools to assign the unique chromosomes within the vcf
        utility_results = get_unique_chrs(vcf_args.vcf)

        # Assign the chromosome list header
        utility_header = ['Chromosomes']

        # Create an output file for the utility
        utility_tsv_output(utility_results, utility_header, utility_index, utility_output_filename)

    # Check if the user has requested a sample list
    elif vcf_args.utility == 'sample-list':

        # Use bcftools to assign the samples within the vcf
        utility_results = get_samples(vcf_args.vcf)

        # Assign the chromosome list header
        utility_header = ['Samples']

        # Create an output file for the utility
        utility_tsv_output(utility_results, utility_header, utility_index, utility_output_filename)

    # Check if the user has requested a concatenated VCF
    elif vcf_args.utility == 'concatenate':

        # Concatenate the VCF files
        concatenate(vcf_args.vcf, vcf_args.out_prefix, vcf_args.out_format, keep_original = True, optional_args = utility_optional_args)

    # Check if the user has requested a concatenated VCF
    elif vcf_args.utility == 'merge':

        # Concatenate the VCF files
        merge(vcf_args.vcf, vcf_args.out_prefix, vcf_args.out_format, keep_original = True, optional_args = utility_optional_args)

    # Check if the user has requested a concatenated VCF
    elif vcf_args.utility == 'sort':

        # Assign the basic arguments for sort
        sort_arguments = ['sort', '-o', utility_output_filename]

        # Assign the output format arguments
        sort_arguments += return_output_format_args(vcf_args.out_format)

        # Concatenate the VCF files
        call_bcftools(sort_argumentse)

    # Check if the output needs to renamed to a user specified filename
    if vcf_args.out and vcf_args.utility in ['merge', 'concatenate']:

        # Override the out_prefix output filename
        shutil.move('%s.%s' % (vcf_args.out_prefix, vcf_args.out_format), utility_output_filename)


if __name__ == "__main__":
    #initLogger()
    run()
