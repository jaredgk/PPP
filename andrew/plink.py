import os
import sys
import subprocess
import shutil
import argparse
import glob
import copy
import logging
import re

# Call PPP-based scripts
sys.path.insert(0, os.path.abspath(os.path.join(os.pardir,'jared')))

from vcf_reader_func import checkFormat
from bcftools import convert_to_bcf
from logging_module import initLogger, logArgs

def assign_delim_from_ids (filename, id_column = 0, default_delim = '_'):

    delim_symbols = ['-', '|', '/', '@', '&', '*', '^', '%', '$']

    symbols_present = []

    with open(filename, 'rU') as sample_file:
        for line_pos, sample_line in enumerate(sample_file):
            if line_pos > 1:
                # Stores the current sample ID
                sample_id = sample_line.strip().split()[id_column]
                # Symbols found in the sample ID
                id_symbols = re.sub(r'[a-zA-Z0-9]','', sample_id)
                # Iterate the symbols from the ID
                for id_symbol in iter(id_symbols):
                    # Check if that symbol has been stored
                    if id_symbol not in symbols_present:
                        # Store symbols
                        symbols_present.append(id_symbol)

    if default_delim not in symbols_present:
        return default_delim
    else:
        for delim_symbol in delim_symbols:
            if delim_symbol not in symbols_present:
                return delim_symbol

    raise IOError('Cannot assign plink delimiter')

def convert_haps_to_vcf (haps_prefix, output_format, output_prefix = ''):

    # Check that input files included a haps file
    if not os.path.isfile(haps_prefix + '.haps'):
        # Return error message if no haps input was found
        raise IOError('Unable to assign haps file. Please confirm the prefix and/or command is specified correctly')

    # Check that input files included a sample file
    if not os.path.isfile(haps_prefix + '.sample'):
        # Return error message if no sample input was found
        raise IOError('Unable to assign sample file. Please confirm the prefix and/or command is specified correctly')

    haps_input_args = ['--haps', haps_prefix + '.haps', 'ref-first', '--sample', haps_prefix + '.sample']

    # If no output_prefix is assigned, use haps_prefix
    if not output_prefix:
        output_prefix = haps_prefix

    plink_delim = assign_delim_from_ids(haps_prefix + '.sample')

    if output_format == 'vcf':
        haps_output_args = ['--recode', 'vcf-iid', 'id-delim=' + plink_delim, '--out', output_prefix]

    elif output_format == 'vcf.gz':
        haps_output_args = ['--recode', 'vcf-iid', 'bgz', 'id-delim=' + plink_delim, '--out', output_prefix]

    elif output_format == 'bcf':
        haps_output_args = ['--recode', 'vcf-iid', 'bgz', 'id-delim=' + plink_delim, '--out', output_prefix]

    else:
        raise IOError('Unknown file format. This error should NOT show up')

    # Call plink2
    standard_plink2_call(haps_input_args + haps_output_args)

    # Convert VCFGZ to BCF, as PLINK cannot directly create a BCF
    if output_format == 'bcf':
        # Convert to BCF
        convert_to_bcf(output_prefix + '.vcf.gz', output_prefix)

def check_plink_for_errors (plink_stderr):
    '''
        Checks the plink stderr for errors

        Parameters
        ----------
        plink_stderr : str
            plink stderr

        Raises
        ------
        IOError
            If plink stderr returns an error
    '''

    # Print output if error found. Build up as errors are discovered
    if plink_stderr:
        raise Exception(plink_stderr)

def standard_plink2_call (plink2_call_args):
    '''
        Standard call of plink2

        The function calls plink2. Returns the stderr of plink to create a log
        file of the call.

        Parameters
        ----------
        plink2_call_args : list
            plink2 arguments

        Raises
        ------
        Exception
            If plink2 stderr returns an error
    '''


    # plink subprocess call
    plink2_call = subprocess.Popen(['plink2'] + list(map(str, plink2_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Wait for plink to finish
    plink2_out, plink2_err = plink2_call.communicate()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        plink2_out = plink2_out.decode()
        plink2_err = plink2_err.decode()

    logging.info('plink2 call complete')

    # Check that the log file was created correctly
    check_plink_for_errors(plink2_err)

def standard_plink_call (plink_call_args):
    '''
        Standard call of plink

        The function calls plink. Returns the stderr of plink to create a log
        file of the call.

        Parameters
        ----------
        plink_call_args : list
            plink arguments

        Raises
        ------
        Exception
            If plink stderr returns an error
    '''

    # plink subprocess call
    plink_call = subprocess.Popen(['plink'] + list(map(str, plink_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Wait for plink to finish
    plink_out, plink_err = plink_call.communicate()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        plink_out = plink_out.decode()
        plink_err = plink_err.decode()

    logging.info('plink call complete')

    # Check that the log file was created correctly
    check_plink_for_errors(plink_err)

def call_plink (plink_call_args, output_prefix = '', output_format = ''):
    '''
        Determines which version of plink (plink1.9 vs. plink2a) to call

        The function is required for the plink wrapper to determine the needed
        plink version and also automates conversion to BCF, which both versions
        of plink do not support.

        Parameters
        ----------
        plink_call_args : list
            plink arguments
        output_prefix : str
            Output prefix for non-plink conversions
        output_format : str
            Output format for non-plink conversions

    '''

    # Check if the call requries plink2 to operate. This function should be
    # removed once plink2 is out of alpha
    if '--haps' in plink_call_args:
        # Convert vcf id system for plink 1.9 to plink 2
        if 'vcf-iid' in plink_call_args:
            format_index = plink_call_args.index('vcf-iid')
            plink_call_args[format_index:format_index + 1] = ['vcf', 'id-paste=fid']

        # Call the plink2 subprocess function
        standard_plink2_call(plink_call_args)

    # Call plink1.9
    else:
        # Call the plink subprocess function
        standard_plink_call(plink_call_args)

    # Convert VCFGZ to BCF, as PLINK cannot directly create a BCF
    if output_format == 'bcf':
        # Convert to BCF
        convert_to_bcf(output_prefix + '.vcf.gz', output_format)

def assign_plink_output_args (output_prefix, output_format):
    '''
        Assigns the output arguments

        Parameters
        ----------
        filename : str
            Specifies the input filename of unknown format

        Returns
        -------
        list
            Returns vcftools input command for `filename`

        Raises
        ------
        IOError
            If filename is an unknown file format
    '''

    # Check if the output is only a single file
    if output_format in ['vcf', 'vcf.gz', 'bcf']:

        if output_format == 'vcf':
            return ['--recode', 'vcf-iid', '--out', output_prefix]

        if output_format == 'vcf.gz':
            return ['--recode', 'vcf-iid', 'bgz', '--out', output_prefix]

        if output_format == 'bcf':
            return ['--recode', 'vcf-iid', 'bgz', '--out', output_prefix]

    # Check if the output if a ped or bed file-set
    elif output_format in ['ped', 'binary-ped']:

        if output_format == 'ped':
            return ['--recode', '--out', output_prefix]

        if output_format == 'binary-ped':
            return ['--make-bed', '--out', output_prefix]

    else:
        raise IOError('Unknown file format. This error should NOT show up')

def assign_plink_input_from_prefix (input_prefix, input_format):
    '''
        Assigns input arguments based on the specified files

        Parameters
        ----------
        filename : str
            Specifies the input filename of unknown format

        Returns
        -------
        list
            Returns vcftools input command for `filename`

        Raises
        ------
        IOError
            If filename is an unknown file format
    '''

    # Determine the input files from the prefix
    prefix_files = glob.glob(input_prefix + '*')

    # Check if files were identified, if not return error message
    if not prefix_files:
        raise IOError('Unable to assign input files Please confirm the prefix is specified correctly')

    # Check for ped input
    if input_format == 'ped':

        # Plink PED and Map files
        plink_ped = None
        plink_map = None

        # Check if expected files are found
        for plink_filename in prefix_files:

            # Check for PED and MAP files
            if '.ped' in plink_filename:
                plink_ped = plink_filename
            elif '.map' in plink_filename:
                plink_map = plink_filename

        # Check that input files included a ped file
        if not plink_ped:
            # Return error message if no ped input was found
            raise IOError('Unable to assign ped file. Please confirm the prefix and/or command is specified correctly')

        if not plink_map:
            # Return error message if no map input was found
            raise IOError('Unable to assign map file. Please confirm the file is named correctly')

        # Return the ped-based files
        return ['--file', plink_ped[:-4]]

    # Check for binary-ped input
    elif input_format == 'binary-ped':

        # Plink BED, BIM, and FAM files
        plink_bed = None
        plink_bim = None
        plink_fam = None

        # Check if expected files are found
        for plink_filename in prefix_files:

            # Check for BED and associated files
            if '.bed' in plink_filename:
                plink_bed = plink_filename
            elif '.bim' in plink_filename:
                plink_bim = plink_filename
            elif '.fam' in plink_filename:
                plink_fam = plink_filename

        # Check that input files included a bed file
        if not plink_bed:
            # Return error message if no ped/bed input was found
            raise IOError('Unable to assign bed file. Please confirm the prefix and/or command is specified correctly')

        # Check that input files included a bed file
        if not plink_bim:
            # Return error message if no bim input was found
            raise IOError('Unable to assign bim file. Please confirm the file is named correctly')

        # Check that input files included a bed file
        if not plink_fam:
            # Return error message if no fam input was found
            raise IOError('Unable to assign fam file. Please confirm the file is named correctly')

        # Return the binary-ped files
        return ['--bfile', plink_bed[:-4]]

    # Check for haps input
    elif input_format == 'haps':

        # Plink haps and sample files
        plink_haps = None
        plink_sample = None

        # Check if expected files are found
        for plink_filename in prefix_files:

            # Check for haps and associated files
            if '.haps' in plink_filename:
                plink_haps = plink_filename
            elif '.sample' in plink_filename:
                plink_sample = plink_filename

        # Check that input files included a haps file
        if not plink_haps:
            # Return error message if no haps input was found
            raise IOError('Unable to assign haps file. Please confirm the prefix and/or command is specified correctly')

        # Check that input files included a sample file
        if not plink_sample:
            # Return error message if no sample input was found
            raise IOError('Unable to assign sample file. Please confirm the file is named correctly')

        # Return the haps-based files
        return ['--haps', plink_haps, 'ref-first', '--sample', plink_sample]

def assign_plink_input_from_command_args (ped_filename = None, map_filename = None,
                                          bed_filename = None, bim_filename = None, fam_filename = None,
                                          haps_filename = None, sample_filename = None,
                                          vcf_filename = None, vcf_fid = None, **kwargs):

    input_command_args = []

    # Check if a ped file is assigned
    if ped_filename:

        # Check if map file is also assigned
        if ped_filename and map_filename:

            # Add the ped associated commands
            input_command_args.extend(['--ped', ped_filename, '--map', map_filename])
            # Return the input command
            return input_command_args

        else:
            raise IOError('Unable to assign ped file. Please confirm the file is named correctly')

    # Check if a bed file is assigned
    elif bed_filename:

        # Check if bim and fam file are also assigned
        if bed_filename and bim_filename and fam_filename:

            # Add the bed associated commands
            input_command_args.extend(['--bed', bed_filename, '--bim', bim_filename, '--fam', fam_filename])
            # Return the input command
            return input_command_args

        else:
            raise IOError('Unable to assign bim and/or fam files. Please confirm the files are named correctly')

    # Check if a haps file is assigned
    elif haps_filename:

        # Check if bim and fam file are also assigned
        if haps_filename and sample_filename:

            # Add the bed associated commands
            input_command_args.extend(['--haps', haps_filename, 'ref-first', '--sample', sample_filename])
            # Return the input command
            return input_command_args

        else:
            raise IOError('Unable to assign sample file. Please confirm the file is named correctly')

    # Check if a vcf file is assigned
    elif vcf_filename:

        # Assign the vcf file format
        vcfname_format = checkFormat(vcf_filename)

        # Check if a fmaily ID was assigned for the vcf file
        if vcf_fid:
            # Assign the family id
            input_command_args.extend(['--const-fid', vcf_fid])

        # If no family ID was assigned, use double IDs
        else:
            input_command_args.append('--double-id')

        input_command_args.append('--allow-extra-chr')

        # Assign the associated input command, or return an error.
        if vcfname_format == 'vcf' or vcfname_format == 'bgzip':

            # Add the bed associated commands
            input_command_args.extend(['--vcf', vcf_filename])
            # Return the input command
            return input_command_args

        elif vcfname_format == 'bcf':

            # Add the bed associated commands
            input_command_args.extend(['--bcf', vcf_filename])
            # Return the input command
            return input_command_args

        else:
            raise Exception('Unknown VCF file format')

    raise Exception('No input specified')
