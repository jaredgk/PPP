#!/usr/bin/env python
'''
    Automates various utilites for VCF-formatted files. This currently includes:
    obtain a list of the chromosomes within a VCF-based file, obtain a list of the 
    samples within a VCF-based file, concatenate multiple VCF-based files, merge
    multiple VCF-based files, and sort a VCF-based file.

    ##################
    Command-line Usage
    ##################
    The VCF utilites function may be called using the following command:

    .. code-block:: bash
        
        python vcf_utilites.py

    *************
    Example usage
    *************
    Concatenate multiple VCF files:

    .. code-block:: bash
        
        python vcf_utilites.py --vcfs chr21.vcf.gz chr22.vcf.gz --utility concatenate

    Merge multiple VCF files:

    .. code-block:: bash
        
        python vcf_utilites.py --vcfs chr22.ceu.vcf.gz chr22.yri.vcf.gz --utility merge
    
    ############
    Dependencies 
    ############
    * `BCFtools <https://samtools.github.io/bcftools/bcftools.html>`_
    
    ############################
    Input Command-line Arguments
    ############################
    **--vcf** *<input_filename>*
        Argument used to define the filename of the VCF file.
    **--vcfs** *<input_filename>* *<input1_filename, input2_filename, etc.>*
        Argument used to define the filename of the VCF file(s). May be used multiple 
        times.
    
    #############################
    Output Command-line Arguments
    #############################
    **--out** *<output_filename>*
        Argument used to define the complete output filename, overrides **--out-prefix**.
    **--out-prefix** *<output_prefix>*
        Argument used to define the output prefix (i.e. filename without file extension)
    **--overwrite**
        Argument used to define if previous output should be overwritten.
    
    ##################################
    Utility Command-line Specification
    ##################################
    **--utility** *<sample-list, chr-list, concatenate, merge, sort>*
        Argument used to define the desired utility. Current utilities include: creation
        of a file of the samples within the VCF (sample-list); creation of a file of the 
        chromosomes within the VCF (chr-list); combine multiple VCF files with different
        variants but the same samples (concatenate); combine multiple VCF files with 
        different samples but the same variants (merge); or sort a single VCF file (sort).

    *****************************************
    Additional Utility Command-line Arguments
    *****************************************
    **--record-merge-mode** *<none, snps, indels, both, all, id>*
        Argument used to define the type of multiallelic records to create. Only usable with 
        the merge utility.
    **--record-missing-as-ref**
        Argument used to define that missing records should be converted to the reference
        allele. Only usable with the merge and concatenate utilites.
    **--out-format** *<vcf, vcf.gz, bcf>*
        Argument used to define the desired output format. Formats include: uncompressed 
        VCF (vcf); compressed VCF (vcf.gz) [default]; and BCF (bcf). Only usable with the 
        merge and concatenate utilites.
'''

import os
import sys
import argparse
import itertools
import copy
import shutil
import logging
import pandas as pd

from pgpipe.bcftools import *
from pgpipe.logging_module import initLogger, logArgs
from pgpipe.misc import argprase_kwargs

def vcf_utility_parser(passed_arguments = []):
    '''
    VCF Utility Argument Parser

    Assign the parameters for VCF Utility using argparse.

    Parameters
    ----------
    passed_arguments : list, optional
        Parameters passed by another function. sys.argv is used if
        not given. 

    Raises
    ------
    IOError
        If the input, or other specified files do not exist
    '''

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
    vcf_input = vcf_parser.add_mutually_exclusive_group(required=True)
    vcf_input.add_argument('--vcf', help = "Defines the filename of the VCF", type = str, action = parser_confirm_file())
    vcf_input.add_argument('--vcfs', help = "Defines the filenames of the VCFs (may be used multiple times)", type = str, nargs = '+', action = parser_confirm_file_list())

    # Output file arguments
    vcf_parser.add_argument('--out', help = 'Defines the complete output filename, overrides --out-prefix', type = str)
    vcf_parser.add_argument('--out-prefix', help = 'Defines the output prefix (i.e. filename without file extension)', default = 'out')
    vcf_parser.add_argument('--overwrite', help = "Defines that previous output files should be overwritten", action = 'store_true')
    
    # Utility based arguments
    utility_list = ['sample-list', 'chr-list', 'concatenate', 'merge', 'sort']
    vcf_parser.add_argument('--utility', metavar = metavar_list(utility_list), help = 'Defines the utility to use', type = str, choices = utility_list)
    
    merge_list = ['none', 'snps', 'indels', 'both', 'all', 'id']
    vcf_parser.add_argument('--record-merge-mode', metavar = metavar_list(merge_list), help = 'Defines the type of multiallelic records to create. Only usable with the merge utility', type = str, choices = merge_list)
    
    vcf_parser.add_argument('--record-missing-as-ref', help = 'Defines that missing records should be converted to the reference allele. Only usable with the merge and concatenate utilites', action = 'store_true')
    
    out_format_list = ['vcf', 'vcf.gz', 'bcf']
    out_format_default = 'vcf.gz'
    vcf_parser.add_argument('--out-format', metavar = metavar_list(out_format_list), help = 'Defines the desired output format', type = str, choices = out_format_list, default = out_format_default)
    
    if passed_arguments:
        return vars(vcf_parser.parse_args(passed_arguments))
    else:
        return vars(vcf_parser.parse_args())

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

def run (**kwargs):
    '''
    Utilites for VCF-formatted files

    This function uses the argparse-based function :py:func:`vcf_utility_parser`
    to parse either sys.argv or passed_arguments to obtain the parameters below. 
    The parameters are then translated to their BCFtools equivalent. Once all the 
    parameters are assigned, BCFtools is called.

    Parameters
    ----------
    --vcf : str
        Filename of the VCF
    --vcfs : list
        List of VCF filenames
    --utility : str
        Utility to be used
    --out : str
        Complete output filename, overrides --out-prefix
    --out-prefix : str
        Output prefix
    --out-format : str
        Desired output format
    --overwrite
        Species if previous output should be overwriten
    --record-merge-mode : str
        Type of multiallelic records to create
    --record-missing-as-ref : bool
        Convert missing records to the reference allele

    Raises
    ------
    IOError
        Output file already exists
    '''

    # Update kwargs with defaults
    if __name__ != "__main__":
        kwargs = argprase_kwargs(kwargs, vcf_utility_parser)

    # Assign arguments
    vcf_args = argparse.Namespace(**kwargs)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(vcf_args, func_name = 'vcf_utilities')

    # Check if the function does not support multiple vcf files
    if vcf_args.utility in ['chr-list', 'sample-list', 'sort']:

        # Check if multiple vcf files have been assigned
        if vcf_args.vcfs:
            raise Exception('Utility (%s) only supports a single VCF file' % vcf_args.utility)

    # Check if the function does not support a single vcf file
    if vcf_args.utility in ['merge', 'concatenate']:

        # Check if multiple vcf files have been assigned
        if vcf_args.vcf:
            raise Exception('Utility (%s) requires multiple VCF files' % vcf_args.utility)

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
        concatenate(vcf_args.vcfs, vcf_args.out_prefix, vcf_args.out_format, keep_original = True, optional_args = utility_optional_args)

    # Check if the user has requested a concatenated VCF
    elif vcf_args.utility == 'merge':

        # Concatenate the VCF files
        merge(vcf_args.vcfs, vcf_args.out_prefix, vcf_args.out_format, keep_original = True, optional_args = utility_optional_args)

    # Check if the user has requested a concatenated VCF
    elif vcf_args.utility == 'sort':

        # Assign the basic arguments for sort
        sort_arguments = ['sort', '-o', utility_output_filename]

        # Assign the output format arguments
        sort_arguments += return_output_format_args(vcf_args.out_format)

        # Append the input VCF
        sort_arguments.append(vcf_args.vcf)

        # Concatenate the VCF files
        call_bcftools(sort_argumentse)

    # Check if the output needs to renamed to a user specified filename
    if vcf_args.out and vcf_args.utility in ['merge', 'concatenate']:

        # Override the out_prefix output filename
        shutil.move('%s.%s' % (vcf_args.out_prefix, vcf_args.out_format), utility_output_filename)


if __name__ == "__main__":
    initLogger()
    run(**vcf_utility_parser())
