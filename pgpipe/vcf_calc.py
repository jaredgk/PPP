#!/usr/bin/env python
'''
    Automates the calculation of site/windowed fixation index (Fst), Tajima's D,
    nucleotide diversity (Pi), allele frequency, and heterozygosity using
    VCFTools. If no statistic is specified, windowed Fst is used by default.

    ##################
    Command-line Usage
    ##################
    The VCF statistic calculator may be called using the following command:

    .. code-block:: bash
        
        vcf_calc.py

    *************
    Example usage
    *************
    Command-line to calculate Tajima's D:

    .. code-block:: bash
        
        vcf_calc.py --vcf examples/files/merged_chr1_10000.vcf.gz --calc-statistic TajimaD --statistic-window-size 10000

    Command-line to calculate windowed Fst on the two populations within the model *2Pop*: 

    .. code-block:: bash
        
        vcf_calc.py --vcf examples/files/merged_chr1_10000.vcf.gz --model-file examples/files/input.model --model 2Pop --calc-statistic windowed-weir-fst --statistic-window-size 10000 --statistic-window-step 10000
    
    ############
    Dependencies 
    ############
    * `VCFtools <https://vcftools.github.io/index.html>`_
     
    ############################
    Input Command-line Arguments
    ############################
    **--vcf** *<input_filename>*
        Argument used to define the filename of the VCF file for calculations.
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
    **--out** *<output_filename>*
        Argument used to define the complete output filename, overrides **--out-prefix**.
        Cannot be used if multiple output files are created.
    **--out-prefix** *<output_prefix>*
        Argument used to define the output prefix (i.e. filename without file extension)
    **--out-dir** *<output_dir_name>*
        Argument used to define the output directory. Only used if multiple output files 
        are created.
    **--overwrite**
        Argument used to define if previous output should be overwritten.
    
    ####################################
    Statistic Command-line Specification
    ####################################
    **--calc-statistic** *<weir-fst, windowed-weir-fst, TajimaD, site-pi, window-pi, freq, het-fit, het-fis, hardy-weinberg>*
        Argument used to define the statistic to be calculated. Site Fst (weir-fst), 
        windowed Fst (windowed-weir-fst), Tajima's D (TajimaD), site nucleotide 
        diversity (site-pi), windowed nucleotide diversity (window-pi), allele frequency
        (freq), Fit (het-fit), Fis (het-fis), and the hardy-weinberg equilibrium 
        (hardy-weinberg).

    ***********************************
    Statistic Command-line Requirements
    ***********************************
    It should be noted that some of the statistics in the VCF calculator require additional
    arguments (i.e. **--pop-file**, **--statistic-window-size**, **--statistic-window-step**).
    These statistics may be found below with their additional requirements. If a statistic is
    not given, only the statistic specification (i.e. **--calc-statistic**) is required.

    **--calc-statistic** *weir-fst*
        Requires: **--pop-file**/**--model**.

    **--calc-statistic** *windowed-weir-fst*
        Requires: **--pop-file**/**--model**, **--statistic-window-size**, and
        **--statistic-window-step**.

    **--calc-statistic** *TajimaD*
        Requires: **--statistic-window-size**

    **--calc-statistic** *windowed-pi*
        Requires: **--statistic-window-size** and **--statistic-window-step**.

    **--calc-statistic** *het-fis*
        Requires: **--pop-file**/**--model**.
    
    *******************************************
    Additional Statistic Command-line Arguments
    *******************************************
    **--statistic-window-size** *<size_int>*
        Defines the statistic window size. Not usable with all statistics.
    **--statistic-window-step** *<step_int>*
        Defines the statistic window step size. Not usable with all statistics.
    **--pop-file** *<pop_filename>*
        Population file. This argument may be used multiple times if desired. Please note
        the this argument is not compatible with either the **--model** or **--model-file** 
        command-line arguments.
    
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
'''

import os
import sys
import argparse
import itertools
import copy
import shutil
import logging

from pgpipe.vcftools import *
from pgpipe.model import read_model_file
from pgpipe.logging_module import initLogger, logArgs
from pgpipe.misc import argprase_kwargs

def vcf_calc_parser(passed_arguments = []):
    '''
    VCF Calc Argument Parser

    Assign the parameters for the VCF Calc using argparse.

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

    # Model file arguments
    vcf_parser.add_argument('--model-file', help = 'Defines the model file', type = str, action = parser_confirm_file())
    vcf_parser.add_argument('--model', help = 'Defines the model and the individual(s)/population(s) to include', type = str)

    # Non-model file arguments
    vcf_parser.add_argument('--pop-file', help = 'Population file. May be used multiple times', nargs = '+', type = str, action = parser_confirm_file_list())

    # Other file arguments. Expand as needed
    vcf_parser.add_argument('--out', help = 'Defines the complete output filename, overrides --out-prefix. Cannot be used if multiple output files are created', type = str)
    vcf_parser.add_argument('--out-prefix', help = 'Defines the output prefix (i.e. filename without file extension)', type = str,  default = 'out')
    vcf_parser.add_argument('--out-dir', help = "Defines the output directory. Only used if multiple output files are created", default = 'Statistic_Files')

    # General arguments
    vcf_parser.add_argument('--overwrite', help = "Overwrite previous output files", action = 'store_true')

    # Galaxy Option to pipe log to stdout
    vcf_parser.add_argument('--log-stdout', help = argparse.SUPPRESS, action = 'store_true')

    # Statistic based arguments
    statistic_list = ['weir-fst', 'windowed-weir-fst', 'TajimaD', 'site-pi', 'window-pi', 'freq', 'het-fit', 'het-fis', 'hardy-weinberg']

    vcf_parser.add_argument('--calc-statistic', metavar = metavar_list(statistic_list), help = 'The statistic to calculate', required = True, type = str, choices = statistic_list)

    # Statistic window options
    vcf_parser.add_argument('--statistic-window-size', help = 'Statistic window size', type = int)
    vcf_parser.add_argument('--statistic-window-step', help = 'Statistic window step size', type = int)

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
        return vars(vcf_parser.parse_args(passed_arguments))
    else:
        return vars(vcf_parser.parse_args())

def run (**kwargs):
    '''
    Statistic calculation using VCFTools.

    This function uses the argparse-based function :py:func:`vcf_calc_parser`
    to parse either sys.argv or passed_arguments to obtain the parameters below. 
    The parameters are then translated to their VCFtools equivalent. Once all the
    VCFtools-based parameters are assigned, VCFtools is called.

    Parameters
    ----------
    --vcf : str
        Input VCF filename
    --out : str
        Complete output filename, overrides --out-prefix
    --out-prefix : str
        Output filename prefix
    --out-dir : str
        Output directory
    --model-file : str
        Model filename
    --model : str
        Model to use
    --pop-file : str
        Defines a population file for relevant statistics. May be
        used multple times
    --calc-statistic : str
        Specifies the statistic to calculate
    --statistic-window-size : int
        Specifies the window size for window-based statistics
    --statistic-window-step : int
        Specifies step size between windows for specific window-based statistics
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

    def calc_exception (selected_model, exc_type, exc_value, exc_traceback):

        # Check if the selected model was correctly assigned
        if selected_model: 
            logging.info('Delete temporary files due to exception')
            selected_model.delete_pop_files()
            selected_model.delete_ind_file()

        # Report the original error    
        sys.__excepthook__(exc_type, exc_value, exc_traceback)

    # Update kwargs with defaults
    if __name__ != "__main__":
        kwargs = argprase_kwargs(kwargs, vcf_calc_parser)

    # Assign arguments
    vcf_args = argparse.Namespace(**kwargs)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(vcf_args, func_name = 'vcf_calc')

    # Argument container for vcftools
    vcftools_call_args = []

    # Used to store population information from either model or pop file(s)
    vcftools_pop_files = []

    # Checks if the user specified both a model and population files
    if vcf_args.model_file and vcf_args.pop_file:
        raise Exception('--model and --pop-file arguments are incompatible')

    # Check that individual-based filters are not being user with a model
    if vcf_args.model_file and (vcf_args.filter_include_indv or vcf_args.filter_exclude_indv):
        if vcf_args.filter_include_indv:
            raise Exception('--model and --filter-include-indv arguments are incompatible')
        if vcf_args.filter_exclude_indv:
            raise Exception('--model and --filter-exclude-indv arguments are incompatible')

    # Check that individuals-based filters are not being user with a model
    if vcf_args.model_file and (vcf_args.filter_include_indv_file or vcf_args.filter_exclude_indv_file):
        if vcf_args.filter_include_indv_file:
            raise Exception('--model and --filter_include_indv-file arguments are incompatible')
        if vcf_args.filter_exclude_indv_file:
            raise Exception('--model and --filter-exclude-indv-file arguments are incompatible')

    # Performs actions related to --model-file argument
    if vcf_args.model_file:

        # Check if a model has been specified
        if not vcf_args.model:
            raise IOError('No selected model. Please use --model to select a model')

        # Read in the model file
        models_in_file = read_model_file(vcf_args.model_file)

        # Check that the selected model was found in the file
        if vcf_args.model not in models_in_file:
            raise IOError('Selected model "%s" not found in: %s' % (vcf_args.model, vcf_args.model_file))

        # Select model, might change this in future versions
        selected_model = models_in_file[vcf_args.model]

        logging.info('%s assigned as model' % selected_model)

        # Check if the specified statistic is fst-based or fis-based
        if vcf_args.calc_statistic in ['windowed-weir-fst', 'weir-fst', 'het-fis']:

            # Create the population files
            selected_model.create_pop_files(file_ext = '.txt', overwrite = vcf_args.overwrite)

            # Store the population files
            vcftools_pop_files = selected_model.pop_files

        # Filter individuals for other statistics
        else:

            # Create individuals file
            selected_model.create_ind_file(overwrite = vcf_args.overwrite)

            # Assign the individuals file to vcftools
            vcftools_call_args.extend(['--keep', selected_model.ind_file])

        # Create a lambda function to pass the selected model along with the exception
        lambda_excepthook = lambda exc_type, exc_value, exc_traceback: calc_exception(selected_model, exc_type, exc_value, exc_traceback)

        # Set the special escape hook
        sys.excepthook = lambda_excepthook

    # Performs actions related to --pop-file argument
    if vcf_args.pop_file:

        # Check if the specified statistic is fst-based or fis-based
        if vcf_args.calc_statistic in ['windowed-weir-fst', 'weir-fst', 'het-fis']:

            # Store the population files
            vcftools_pop_files = vcf_args.pop_file

        # Return an error if --pop-file is not supported
        else:
            raise Exception('%s does not support the --pop-file argument' % vcf_args.calc_statistic)

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

    # Check if windowed Fst is specified
    if vcf_args.calc_statistic == 'windowed-weir-fst':

        # Confirm that population files has been assigned
        if len(vcftools_pop_files) < 2:
            raise Exception('%s requires at least two populations to operate' % vcf_args.calc_statistic)

        # Confirm that the window size has been assigned
        if not vcf_args.statistic_window_size:
            raise Exception('%s requires a specifed window size to operate' % vcf_args.calc_statistic)

        # Confirm that window step size has been assigned
        if not vcf_args.statistic_window_step:
            raise Exception('%s requires a specifed window step size to operate' % vcf_args.calc_statistic)

        # Assigns the required window arguments
        vcftools_call_args.extend(['--fst-window-size', vcf_args.statistic_window_size, '--fst-window-step', vcf_args.statistic_window_step])

        # Assigns the suffix for the vcftools log file
        vcftools_out_suffix = '.windowed.weir.fst'

    elif vcf_args.calc_statistic == 'weir-fst':

        # Confirm that population files has been assigned
        if len(vcftools_pop_files) < 2:
            raise Exception('%s requires at least two populations to operate' % vcf_args.calc_statistic)

        # Assigns the suffix for the vcftools log file
        vcftools_out_suffix = '.weir.fst'

    elif vcf_args.calc_statistic == 'TajimaD':

        # Confirm that the window size has been assigned
        if not vcf_args.statistic_window_size:
            raise Exception('%s requires a specifed window size to operate' % vcf_args.calc_statistic)

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

        # Confirm that the window size has been assigned
        if not vcf_args.statistic_window_size:
            raise Exception('%s requires a specifed window size to operate' % vcf_args.calc_statistic)

        # Confirm that window step size has been assigned
        if not vcf_args.statistic_window_step:
            raise Exception('%s requires a specifed window step size to operate' % vcf_args.calc_statistic)

        # Assigns all the vcftools arguments for calculating pi
        vcftools_call_args.extend(['--window-pi', vcf_args.statistic_window_size, '--window-pi-step', vcf_args.statistic_window_step])

        # Assigns the suffix for the vcftools log file
        vcftools_out_suffix = '.windowed.pi'

    elif vcf_args.calc_statistic == 'hardy-weinberg':

        # Assigns all the vcftools arguments for the allele frequency
        vcftools_call_args.append('--hardy')

        # Assigns the suffix for the vcftools log file
        vcftools_out_suffix = '.hwe'

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

        # Confirm that at least one population file has been assigned
        if not vcftools_pop_files:
            raise Exception('%s requires at least one population to operate' % vcf_args.calc_statistic)

        # Assigns all the vcftools arguments for calculating heterozygosity
        vcftools_call_args.append('--het')

        # Assigns the suffix for the vcftools log file
        vcftools_out_suffix = '.het'

    logging.info('vcftools parameters assigned')

    # Check if the user has specified a output filename
    if vcf_args.out:

        # Assign the vcftools output filename, using the output filename
        vcftools_output_filename = vcf_args.out

    else:

        # Assign the vcftools output filename, using the prefix and suffix
        vcftools_output_filename = vcf_args.out_prefix + vcftools_out_suffix

    # Check if previous output should be overwritten
    if vcf_args.overwrite:

        # Check if the ouput will require the output directory
        if vcf_args.calc_statistic in ['windowed-weir-fst', 'weir-fst'] and len(vcftools_pop_files) > 2:
            
            # Check if the output directory is present
            if os.path.exists(vcf_args.out_dir):

                # Remove the output directory
                shutil.rmtree(vcf_args.out_dir)

    # If not to be overwritten, check if previous output exists
    else:

        # Confirm the vcftools output and log file do not exist
        check_for_vcftools_output(vcftools_output_filename)

        # Check if the ouput will require the output directory
        if vcf_args.calc_statistic in ['windowed-weir-fst', 'weir-fst'] and len(vcftools_pop_files) > 2:

            # Check if the output directory is present and report the error
            if os.path.exists(vcf_args.out_dir):
                raise IOError('Statistic Directory already exists')

    # Assigns the file argument for vcftools
    vcfname_arg = assign_vcftools_input_arg(vcf_args.vcf)

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

            # The filename population file
            pop_filename = pop_prefix + vcftools_out_suffix

            # vcftools subprocess call, with stdout
            vcftools_err = call_vcftools(vcfname_arg + pop_call_args, output_filename = pop_filename)

            # Check if the log should be piped to the stdout
            if vcf_args.log_stdout:

                # Write the log to stdout
                sys.stdout.write(vcftools_err)

            # Check if log should be saved as a file
            else:

                # Produce the vcftools log file, in append mode
                produce_vcftools_log(vcftools_err, vcftools_output_filename, append_mode = True)

    elif vcf_args.calc_statistic == 'het-fis':

        # Loop each population file
        for vcftools_pop_file in vcftools_pop_files:

            # Create the population-specific call
            pop_call_args = vcftools_call_args + ['--keep', str(vcftools_pop_file)]

            # vcftools subprocess call
            vcftools_err = call_vcftools(vcfname_arg + pop_call_args, output_filename = vcftools_output_filename, append_mode = True)

            # Check if the log should be piped to the stdout
            if vcf_args.log_stdout:

                # Write the log to stdout
                sys.stdout.write(vcftools_err)

            # Check if log should be saved as a file
            else:

                # Produce the vcftools log file, in append mode
                produce_vcftools_log(vcftools_err, vcftools_output_filename, append_mode = True)

    else:

        # Check if either Fst statistic is specified
        if vcf_args.calc_statistic in ['windowed-weir-fst', 'weir-fst'] and len(vcftools_pop_files) == 2:

            # Assigns the population files to the vcftools call
            vcftools_call_args.extend([pop_args for pop_file in vcftools_pop_files for pop_args in ['--weir-fst-pop', pop_file]])

        # vcftools subprocess call
        vcftools_err = call_vcftools(vcfname_arg + vcftools_call_args, output_filename = vcftools_output_filename)

        # Check if the log should be piped to the stdout
        if vcf_args.log_stdout:

            # Write the log to stdout
            sys.stdout.write(vcftools_err)

        # Check if log should be saved as a file
        else:

            # Produce the vcftools log file
            produce_vcftools_log(vcftools_err, vcftools_output_filename)

    # Lof the vcftools reference
    log_vcftools_reference()

    # Delete any files that were created for vcftools
    if vcf_args.model_file:
        selected_model.delete_pop_files()
        selected_model.delete_ind_file()

if __name__ == "__main__":
    initLogger()
    run(**vcf_calc_parser())
