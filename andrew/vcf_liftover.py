import os
import sys
import argparse
import itertools
import copy
import shutil
import logging

# Insert Jared's directory path, required for calling Jared's functions. Change when directory structure changes.
sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'pppipe')))

from logging_module import initLogger, logArgs
from vcf_reader_func import checkFormat
from picard import call_picard, get_ref_dictonary, check_ref_file_association, remove_picard_index
from fasta import check_format

def vcf_liftover_parser(passed_arguments):
    '''VCF Argument Parser - Assigns arguments from command line'''

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

    vcf_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    # Input arguments.
    vcf_parser.add_argument('--vcf', help = "Input VCF filename (may be used multiple times)", type = str, required = True, action = parser_confirm_file())
    vcf_parser.add_argument('--ref', help = 'Reference fasta of target file for liftover', type = str, required = True, action = parser_confirm_file())
    vcf_parser.add_argument('--ref-dict', help = 'Reference fasta dictionary, if filename differs', type = str, action = parser_confirm_file())
    vcf_parser.add_argument('--chain', help = 'Chain file for liftover', type = str, required = True, action = parser_confirm_file())

    # Output file arguments
    out_format_list = ['vcf', 'vcf.gz', 'bcf']
    out_format_default = 'vcf.gz'
    vcf_parser.add_argument('--out-format', metavar = metavar_list(out_format_list), help = 'Format of the output', type = str, choices = out_format_list, default = out_format_default)
    vcf_parser.add_argument('--out', help = 'Output filename. If used, overrides --out-prefix', type = str)
    vcf_parser.add_argument('--out-prefix', help = 'Defines the output prefix', default = 'out')

    # Output file arguments
    reject_format_list = ['vcf', 'vcf.gz', 'bcf']
    vcf_parser.add_argument('--reject-format', metavar = metavar_list(reject_format_list), help = 'Format of the reject. Same as --out-format by default', type = str, choices = reject_format_list)
    vcf_parser.add_argument('--reject', help = 'Reject filename. If used, overrides --reject-prefix', type = str)
    vcf_parser.add_argument('--reject-prefix', help = 'Defines the reject prefix', default = 'reject')

    # Optional arguments
    vcf_parser.add_argument('--liftover-min-match', help = 'Minimum percentage required for liftover of a variant', type = float) 

    # General arguments
    vcf_parser.add_argument('--overwrite', help = "Overwrite previous output files", action = 'store_true')
    vcf_parser.add_argument('--keep-index', help = "Keep the output VCF index created by picard", action = 'store_true')
    vcf_parser.add_argument('--picard-max-records', help = 'Number of records to store in memory for picard', type = int)
    vcf_parser.add_argument('--picard-path', help = "Defines path to locate picard.jar", type = str, default = 'bin/')

    if passed_arguments:
        return vcf_parser.parse_args(passed_arguments)
    else:
        return vcf_parser.parse_args()

def assign_vcf_extension (filename):

    # Checks if the file is unzipped, bgzipped, or gzipped
    vcfname_format = checkFormat(filename)

    if vcfname_format == 'vcf':
        return '.vcf'

    elif vcfname_format == 'gzip' or vcfname_format == 'bgzip':
        return '.vcf.gz'

    elif vcfname_format == 'bcf':
        return '.bcf'

    else:
        raise Exception('%s - Unknown file format' % filename)

def assign_ref_extension (filename):

    # Checks if the file is unzipped, bgzipped, or gzipped
    ref_format = check_format(filename)

    if ref_format == 'fasta':
        return '.fa'

    elif ref_format == 'gzip':
        return '.fa.gz'

    else:
        raise Exception('%s - Unknown file format' % filename)


def run (passed_arguments = []):
    '''
    Liftover for VCF-formatted files

    Automates the liftover function from picard for VCF-formatted files.

    Parameters
    ----------
    --vcf : str
        Input VCF filename (from origin assembly)
    --ref : str
        Target assembly reference fasta
    --ref-dict
        Dictionary file of reference fasta
    --chain : str
        Chain file from origin assembly to target assembly
    --out-format : str
        Format of the output file
    --out-prefix : str
        Output filename prefix
    --out : str
        Output filename. Overrides --out-prefix
    --reject-format : str
        Format of the reject file, same as --out-format by default
    --reject-prefix : str
        Reject filename prefix
    --reject : str
        Reject filename. Overrides --reject-prefix
    --overwrite
        If previous output should be overwriten
    --keep-index
        Keep the output VCF index created by picard
    --picard-max-records : int
        Number of records to store in memory for picard
    --picard-path : str
        Defines path to locate picard.jar
        
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
    vcf_args = vcf_liftover_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(vcf_args, func_name = 'vcf_liftover')

    # Check if a reference dictionary was not assigned
    if not vcf_args.ref_dict:

        # Get filename of ref dictionary, if found
        ref_dict = get_ref_dictonary(vcf_args.ref)

        # Check if the filename was found
        if not ref_dict:
            raise IOError('No dictionary associated with reference. Please create one using fasta_utilities')

        # Assign the ref dictionary filename
        vcf_args.ref_dict = ref_dict

    # Assign file extension for VCF input file
    ref_ext = assign_ref_extension(vcf_args.ref)

    # Assign the ref filename, in case of rename
    vcf_liftover_ref = copy.deepcopy(vcf_args.ref)

    # Assign the dict filename, in case of rename
    vcf_liftover_dict = copy.deepcopy(vcf_args.ref_dict)

    # Check if the reference and reference dictionary are not associated for picard
    if not check_ref_file_association(vcf_args.ref, vcf_args.ref_dict):

        # Edit the ref argument
        vcf_liftover_dict = vcf_liftover_ref + '.dict'

        # Rename the ref file
        os.rename(vcf_args.ref_dict, vcf_liftover_dict)

        # Edit the ref argument
        vcf_liftover_ref += ref_ext

        # Rename the ref file
        os.rename(vcf_args.ref, vcf_liftover_ref)

    # Assign file extension for VCF input file
    vcfname_ext = assign_vcf_extension(vcf_args.vcf)

    # Assign the vcf filename, in case of rename
    vcf_liftover_vcf = copy.deepcopy(vcf_args.vcf)

    # Confirm input has correct file extension
    if vcfname_ext not in vcf_args.vcf:

        # Edit the VCF argument
        vcf_liftover_vcf += vcfname_ext

        # Rename the file
        os.rename(vcf_args.vcf, vcf_liftover_vcf)

    # Create list to store arguments for liftover
    vcf_liftover_args = ['LiftoverVcf']

    # Assign the input VCF file
    vcf_liftover_args.append('I=' + vcf_liftover_vcf)

    # Assign the chain file
    vcf_liftover_args.append('C=' + vcf_args.chain)

    # Assign the reference FASTA file
    vcf_liftover_args.append('R=' + vcf_liftover_ref)

    # Check if a output name has been assigned
    if vcf_args.out:

        # Assign the output name
        vcf_liftover_out = vcf_args.out + os.extsep + vcf_args.out_format

    else:

        # Assign the output name
        vcf_liftover_out = vcf_args.out_prefix + os.extsep + vcf_args.out_format

    # Assign the output VCF file
    vcf_liftover_args.append('O=' + vcf_liftover_out)

    # Check if a specific reject format has not been specified
    if not vcf_args.reject_format:

        # Assign the output format to the reject format
        vcf_args.reject_format = vcf_args.out_format

    # Check if a output name has been assigned
    if vcf_args.reject:

        # Assign the output name
        vcf_liftover_reject = vcf_args.reject + os.extsep + vcf_args.reject_format

    else:

        # Assign the output name
        vcf_liftover_reject = vcf_args.reject_prefix + os.extsep + vcf_args.reject_format

    # Assign the output VCF file
    vcf_liftover_args.append('REJECT=' + vcf_liftover_reject)

    # Check if the user has assigned a record max
    if vcf_args.picard_max_records:

        # Assign the output VCF file
        #vcf_liftover_args.extend(['--MAX_RECORDS_IN_RAM', str(vcf_args.picard_max_records)])
        vcf_liftover_args.append('MAX_RECORDS_IN_RAM=' + str(vcf_args.picard_max_records))

    # Call picard with the arguments
    call_picard(vcf_args.picard_path, vcf_liftover_args, vcf_liftover_out)

    # Check if the index file created by picard should be removed
    if not vcf_args.keep_index:
        remove_picard_index(vcf_liftover_out, vcf_args.out_format)

    # Reverts the VCF input file to orginal filename
    if vcf_liftover_vcf != vcf_args.vcf:
        os.rename(vcf_liftover_vcf, vcf_args.vcf)

    # Reverts the ref files to orginal filenames
    if vcf_liftover_ref != vcf_args.ref:
        os.rename(vcf_liftover_ref, vcf_args.ref)
        os.rename(vcf_liftover_dict, vcf_args.ref_dict)

    # Check if a output name has been assigned and rename
    if vcf_args.out:
        os.rename(vcf_liftover_out, vcf_args.out)

    # Check if a output name has been assigned and rename
    if vcf_args.reject:
        os.rename(vcf_liftover_reject, vcf_args.reject)


if __name__ == "__main__":
    #initLogger()
    run()