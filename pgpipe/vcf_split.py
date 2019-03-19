'''
    Splits a single VCF file into multiple VCF files using either a statistic or
    bed file. The statistic/bed file must contain loci-based (i.e. window-based)
    data for the function to operate. If the specified statistic file does not
    contain a BIN_END column, the --statistic-window-size argument may be used.

    ############################
    Input Command-line Arguments
    ############################
    **--vcf** *<input_filename>*
        Argument used to define the filename of the VCF file to be split.
    **--split-file** *<split_filename>*
        Argument used to define the file to be split
    **--model-file** *<model_filename>*
        Argument used to define the model file. Please note that this argument cannot be 
        used with the **--pop-file** argument or individual-based filters.
    **--model** *<model_str>*
        Argument used to define the model (i.e. the individual(s) to include and/or the 
        populations for relevant statistics). May be used with any statistic. Please note 
        that this argument cannot be used with **--pop-file** argument or the 
        individual-based filters.

    #############################
    Output Command-line Arguments
    #############################
    **--out-prefix** *<output_prefix>*
        Argument used to define the output prefix (i.e. filename without file extension)
    **--out-format** *<vcf, vcf.gz, bcf>*
        Argument used to define the desired output format. Formats include: uncompressed 
        VCF (vcf); compressed VCF (vcf.gz) [default]; and BCF (bcf).
    **--out-dir** *<output_dir_name>*
        Argument used to define the output directory.
    **--overwrite**
        Argument used to define if previous output should be overwritten.

    ############################
    Split Command-line Arguments
    ############################
    **--split-method** *<statistic-file, bed>*
        Argument used to define the splitting method. Users may spilit using either 
        a statistic-file (statistic-file) from VCF Calc (or other methods) or a BED 
        file (bed).
    **--statistic-window-size** *<statistic_window_int>*
        Argument used to define the size of window calculations. This argument is only 
        required if the BIN_END column is absent within the file.
    **--no-window-correction**
        Argument used to define if a window should not be corrected to avoid 
        an overlap of a single position (i.e. 100-200/200-300 vs. 100-199/200-299).

    #############################
    Filter Command-line Arguments
    #############################
    If using an unfiltered VCF file (e.g. reduce the creation of unnecessary large files)
    the VCF calculator is able to use either a kept or removed sites/BED file and the
    individual-based paramemeters. 
    
    **************************
    Individual-Based Arguments
    **************************
    Please note that all individual-based arguments are not compatible with either the 
    **--model** or **--model-file** command-line arguments.

    **--filter-include-indv** *<indv_str>* *<indv1_str, indv2_str, etc.>*
        Argument used to define the individual(s) to include. This argument may be used 
        multiple times if desired.
    **--filter-exclude-indv** *<indv_str>* *<indv1_str, indv2_str, etc.>*
        Argument used to define the individual(s) to exclude. This argument may be used 
        multiple times if desired.
    **--filter-include-indv-file** *<indv_filename>*
        Argument used to define a file of individuals to include.
    **--filter-exclude-indv-file** *<indv_filename>*
        Argument used to define a file of individuals to exclude.
    
    ************************
    Position-Based Arguments
    ************************
    **--filter-include-positions** *<position_filename>*
        Argument used to define a file of positions to include within a tsv file 
        (chromosome and position).
    **--filter-exclude-positions** *<position_filename>*
        Argument used to define a file of positions to exclude within a tsv file 
        (chromosome and position).
    **--filter-include-bed** *<position_bed_filename>*
        Argument used to define a BED file of positions to include.
    **--filter-exclude-bed** *<position_bed_filename>*
        Argument used to define a BED file of positions to exclude.

    ##########################
    Example Command-line Usage
    ##########################
    Command-line to split using a statistic file:

    .. code-block:: bash
        
        python vcf_split.py --vcf input.vcf --split-file input.windowed.weir.fst --split-method statistic-file
'''

import os
import sys
import copy
import argparse
import logging
import shutil
import pandas as pd
import numpy as np

from pybedtools import BedTool

# Import basic vcftools functions
from pgpipe.vcftools import *

# Model file related functions
from pgpipe.model import read_model_file

# Import logging module
from pgpipe.logging_module import initLogger, logArgs

def vcf_split_parser(passed_arguments):
    '''
    VCF Split Argument Parser

    Assign the parameters for VCF Split using argparse.

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

    def parser_add_to_list ():
        '''Custom action to add items to a list'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):

                # Clean up any commas
                value = [item.strip(',') for item in value]
                
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
    vcf_parser.add_argument('--vcf', help = "Defines the filename of the VCF", type = str, required = True, action = parser_confirm_file())

    # Model file arguments.
    vcf_parser.add_argument('--model-file', help = 'Defines the model file', type = str, action = parser_confirm_file())
    vcf_parser.add_argument('--model', help = 'Defines the model and the individual(s) to include', type = str)

    # Split arguments
    vcf_parser.add_argument('--split-file', help='Defines the file to split', required = True, type = str, action = parser_confirm_file())

    split_list = ['statistic-file', 'bed']
    vcf_parser.add_argument('--split-method', metavar = metavar_list(split_list), help = 'Defines the splitting method', type=str, choices = split_list)
    vcf_parser.add_argument('--statistic-window-size', help = 'Defines the size of window calculations (if BIN_END is absent in file)', type = int)
    vcf_parser.add_argument('--no-window-correction', help = 'Defines if the window correction should not be used', action = 'store_true')

    # Output file arguments
    out_format_list = ['vcf', 'vcf.gz', 'bcf']
    out_format_default = 'vcf.gz'
    vcf_parser.add_argument('--out-format', metavar = metavar_list(out_format_list), help = 'Defines the desired output format', type = str, choices = out_format_list, default = out_format_default)
    vcf_parser.add_argument('--out-dir', help = 'Defines the output directory', type = str,  default = 'Split_VCFs')
    vcf_parser.add_argument('--out-prefix', help = 'Defines the output prefix (i.e. filename without file extension)', type = str,  default = 'out')

    # General arguments.
    vcf_parser.add_argument('--overwrite', help = "Overwrite previous output files", action = 'store_true')

    # Galaxy Option to pipe log to stdout
    vcf_parser.add_argument('--log-stdout', help = argparse.SUPPRESS, action = 'store_true')

    ## Filters
    # Position-based position filters
    vcf_parser.add_argument('--filter-include-positions', help = 'Defines a file of positions to include within a tsv file (chromosome and position)', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-positions', help = 'Defines a file of positions to exclude within a tsv file (chromosome and position)', action = parser_confirm_file())

    # BED-based position filters
    vcf_parser.add_argument('--filter-include-bed', help = 'Defines a BED file of positions to include', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-bed', help = 'Defines a BED file of positions to exclude', action = parser_confirm_file())

    # Non-model based filters
    vcf_parser.add_argument('--filter-include-indv', help = 'Defines the individual(s) to include. This argument may be used multiple times if desired', nargs = '+', type = str, action = parser_add_to_list())
    vcf_parser.add_argument('--filter-exclude-indv', help = 'Defines the individual(s) to exclude. This argument may be used multiple times if desired', nargs = '+', type = str, action = parser_add_to_list())
    vcf_parser.add_argument('--filter-include-indv-file', help = 'Defines a file of individuals to include', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-indv-file', help = 'Defines a file of individuals to exclude', action = parser_confirm_file())

    if passed_arguments:
        return vcf_parser.parse_args(passed_arguments)
    else:
        return vcf_parser.parse_args()

def return_missing_columns (sample_headers):
    '''Reports if a required column (i.e. CHROM, BIN_START, BIN_END) is missing
    from the header of the statistic file.'''

    if not all(pos_header in sample_headers for pos_header in ['CHROM', 'BIN_START', 'BIN_END']):
        # Assign missing headers
        missing_headers = [pos_header for pos_header in ['CHROM', 'BIN_START', 'BIN_END'] if pos_header not in sample_headers]

        return missing_headers

    return None

def split_samples_iter (split_sample_data, split_method):

    # Check if the split method is a statistic file
    if split_method == 'statistic-file':

        for sample_index, split_sample in split_sample_data.iterrows():

            yield sample_index, split_sample

    # Check if the split method is a bed file
    elif split_method == 'bed':

        for sample_index, split_sample in enumerate(split_sample_data):

            yield sample_index, split_sample

def assign_position_args (sample_data, split_method):

    # Assignment variable list - chrom, start, end
    assign_list = ['', '', '']

    # Assign the statistic file columns
    if split_method == 'statistic-file':

        # Assign the headers
        assign_list = ['CHROM', 'BIN_START', 'BIN_END']

    # Assign the BED file columns
    elif split_method == 'bed':

        # Assign the headers
        assign_list = ['chrom', 'start', 'end']

    # Return the sample position information
    return ['--chr', sample_data[assign_list[0]], 
            '--from-bp', sample_data[assign_list[1]], 
            '--to-bp', sample_data[assign_list[2]]]

def run (passed_arguments = []):
    '''
    Split VCF file into multiple VCFs.

    This function uses the argparse-based function :py:func:`vcf_split_parser`
    to parse either sys.argv or passed_arguments to obtain the parameters below. 
    The parameters are then translated to their VCFtools equivalent. Once all 
    the parameters are assigned, VCFtools is called.

    Parameters
    ----------
    --vcf : str
        Filename of the VCF
    --out-prefix : str
        Output prefix
    --out-dir : str
        Output directory
    --out-format : str
        Desired output format
    --model-file : str
        Model filename
    --model : str
        Model to use
    --split-file : str
        Filename of file to split
    --split-method : str
        Splitting method
    --statistic-window-size : int
        Size of window calculations (if BIN_END column is absent)
    --no-window-correction : bool
        If the window correction should not be used
    --overwrite : bool
        If previous output should be overwritten    
    --filter-include-indv : list or str
        Individual(s) to include. May be used multiple times. Not usable w/ --model
    --filter-exclude-indv : list or str
        Individual(s) to exclude. May be used multiple times. Not usable w/ --model
    --filter-include-indv-file : str
        File of individuals to include. Not usable w/ --model
    --filter-exclude-indv-file : str
        File of individuals to exclude.  Not usable w/ --model
    --filter-include-positions : str
        File of positions to include (tsv: chromosome and position)
    --filter-exclude-positions : str
        File of positions to exclude (tsv: chromosome and position)
    --filter-include-bed : str
        BED file of positions to include 
    --filter-exclude-bed : str
        BED file of positions to exclude

    Raises
    ------
    IOError
        Output file already exists and --overwrite is not specified
    Exception
        Incompatible arguments
    '''

    # Grab VCF arguments from command line
    vcf_args = vcf_split_parser(passed_arguments)

    # Argument container for vcftools
    vcftools_call_args = []

    # Check if previous output should be overwritten
    if not vcf_args.overwrite:
        # Check if the output directory is present
        if os.path.isdir(vcf_args.out_dir):
            raise IOError('%s already exists' % vcf_args.out_dir)
    else:
        # Check if the output directory is present
        if os.path.isdir(vcf_args.out_dir):
            # Delete the directory
            shutil.rmtree(vcf_args.out_dir)

    # Check if the user has specified a model file
    if vcf_args.model_file and vcf_args.model:
        # Read in the models
        models_in_file = read_model_file(vcf_args.model_file)

        # Check that the selected model was not found in the file
        if vcf_args.model not in models_in_file:
            raise IOError('Selected model "%s" not found in: %s'
                           % (vcf_args.model, vcf_args.model_file))

        # Check that individual-based filters are not also being used
        if vcf_args.filter_include_indv or vcf_args.filter_exclude_indv:
            if vcf_args.filter_include_indv:
                raise Exception('--model and --filter-include-indv arguments '
                                'are incompatible')
            if vcf_args.filter_exclude_indv:
                raise Exception('--model and --filter-exclude-indv arguments '
                                'are incompatible')

        # Check that individuals-based filters are not also being used
        if vcf_args.filter_include_indv_file or vcf_args.filter_exclude_indv_file:
            if vcf_args.filter_include_indv_file:
                raise Exception('--model and --filter_include_indv-file '
                                'arguments are incompatible')
            if vcf_args.filter_exclude_indv_file:
                raise Exception('--model and --filter-exclude-indv-file '
                                'arguments are incompatible')

        # Select model, might change this in future versions
        selected_model = models_in_file[vcf_args.model]

        # Create individuals file
        selected_model.create_ind_file(overwrite = vcf_args.overwrite)

        # Assign the individuals file to vcftools
        vcftools_call_args.extend(['--keep', selected_model.ind_file])

        logging.info('Model file assigned')

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

    # Used to assign the output format to the vcftools call, and assign
    # vcftools output suffix
    if vcf_args.out_format == 'bcf':
        vcftools_call_args.append('--recode')
        vcftools_output_suffix = '.recode.bcf'
    elif vcf_args.out_format == 'vcf':
        vcftools_call_args.append('--recode')
        vcftools_output_suffix = '.recode.vcf'
    elif vcf_args.out_format == 'vcf.gz':
        vcftools_call_args.append('--recode')
        vcftools_output_suffix = '.recode.vcf.gz'

    logging.info('vcftools common parameters assigned')

    # Create the vcf/bcf output directory
    if not os.path.exists(vcf_args.out_dir):
        os.makedirs(vcf_args.out_dir)

    # Assigns the file argument for vcftools
    vcfname_arg = assign_vcftools_input_arg(vcf_args.vcf)

    logging.info('Input file assigned')

    # Check if a statistic file is the split method
    if vcf_args.split_method == 'statistic-file':

        # Read in the sample file
        split_samples = pd.read_csv(vcf_args.split_file, sep = '\t')

        logging.info('Statistic file for splitting assigned')

        # Assign the columns
        split_columns = list(split_samples)

        # Assign missing columns
        missing_columns = return_missing_columns(split_columns)

        # Check for missing columns
        if missing_columns:

            # Report error if too many or unexpected columns are missing
            if len(missing_columns) != 1 or 'BIN_END' not in missing_columns:
                raise ValueError('Cannot find %s column(s) in file specified by ' \
                                 '--split-file.' % ', '.join(missing_headers))

            # Check if statistic_window_size has been assigned
            if not vcf_args.statistic_window_size:
                raise ValueError("'BIN_END' column absent in %s. Please use " \
                                 '--statistic-window-size' % vcf_args.split_file)

            # Assign window size
            window_size = vcf_args.statistic_window_size

            # Check if the overlap correction has been specified
            if not vcf_args.no_window_correction:

                # Correct the window size
                window_size -= 1

                logging.info('Applied window correction. To disable the correction '
                             "use '--no-window-correction'")

            # Create the 'BIN_END' column
            split_samples['BIN_END'] = split_samples['BIN_START'] + window_size

            logging.info("Calculated 'BIN_END' using --statistic-window-size")

    elif vcf_args.split_method == 'bed':

        # Read in the BED file
        split_samples = BedTool(vcf_args.split_file)

    # Loop the rows of the sample to be split
    for sample_index, split_sample in split_samples_iter(split_samples, vcf_args.split_method):

        # Copy common arguments for VCF call
        sample_call_args = copy.deepcopy(vcftools_call_args)

        # Assign the sample output prefix
        sample_prefix = '%s_%s' % (vcf_args.out_prefix, sample_index)

        # Add the path to the prefix
        sample_path = os.path.join(vcf_args.out_dir, sample_prefix)

        # Store the sample output prefix
        sample_call_args.extend(['--out', sample_path])

        # Assign the position args using the split method
        assigned_position_args = assign_position_args(split_sample, vcf_args.split_method)

        # Store the sample position args
        sample_call_args.extend(assigned_position_args)

        # Assign the expected output filename
        vcftools_sample_filename = sample_path + vcftools_output_suffix

        logging.info('vcftools sample parameters assigned')

        # Call vcftools with the specifed arguments
        vcftools_err = call_vcftools(vcfname_arg + sample_call_args, output_format = vcf_args.out_format, output_filename = vcftools_sample_filename)

        # Check if the log should be piped to the stdout
        if vcf_args.log_stdout:

            # Write the log to stdout
            sys.stdout.write(vcftools_err)

        # Check if log should be saved as a file
        else:

            # Produce the log file (in append mode, will create a single log)
            produce_vcftools_log(vcftools_err, vcf_args.out_prefix + '.split', append_mode = True)


if __name__ == "__main__":
    initLogger()
    run()
