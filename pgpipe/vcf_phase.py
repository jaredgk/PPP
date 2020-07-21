#!/usr/bin/env python
"""
    Phasing is an essental and frequently used process in population genetic 
    analyses. Given an unphased VCF file and a selected phasing algorithm, 
    vcf_phase will produce a phased VCF. Phasing may be configured using
    various general options (e.g. specifying Ne, including a genetic map)
    or algorithm-specific options (e.g. including a compatible reference 
    panel) as needed.

    .. image:: ../../PPP_assets/PPP_Phase.png
        :width: 100 %
        :align: center
    
    In this illustration of the phasing process, unphased variants (alleles
    divided diagonally) are converted into an estimated haplotypes (alleles
    divided horizontally and on seperate strands).

    ##################
    Command-line Usage
    ##################
    The VCF file phaser may be called using the following command:

    .. code-block:: bash
        
        vcf_phase.py

    *************
    Example usage
    *************
    Command-line to phase a VCF using Beagle:

    .. code-block:: bash
        
        vcf_phase.py --vcf examples/files/merged_chr1_10000.unphased.vcf.gz --phase-algorithm beagle

    Command-line to phase a VCF using SHAPEIT:

    .. code-block:: bash
        
        vcf_phase.py --vcf examples/files/merged_chr1_10000.unphased.vcf.gz --phase-algorithm shapeit

    ############
    Dependencies 
    ############
    * `SHAPEIT <https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html>`_
    * `Beagle 5.0 <https://faculty.washington.edu/browning/beagle/beagle.html>`_
    * `BCFtools <https://samtools.github.io/bcftools/bcftools.html>`_ 
    * `plink 2.0 <https://www.cog-genomics.org/plink/2.0/>`_

    ############################
    Input Command-line Arguments
    ############################
    **--vcf** *<input_filename>*
        Argument used to define the filename of the VCF file to be phased.
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
    **--out-prefix** *<output_prefix>*
        Argument used to define the output prefix (i.e. filename without file extension)
    **--out-format** *<vcf, vcf.gz, bcf>*
        Argument used to define the desired output format. Formats include: uncompressed 
        VCF (vcf); compressed VCF (vcf.gz) [default]; and BCF (bcf).
    **--overwrite**
        Argument used to define if previous output should be overwritten.
    
    ##############################
    Phasing Command-line Arguments
    ##############################
    **--phase-algorithm** *<beagle, shapeit>*
        Argument used to define the phasing algorithm. BEAGLE 5.0 (beagle) [default] and
        SHAPEIT (shapeit). Please note: Both algorithms possess algorithm-specific 
        arguments that may be found in their respective sections.
    **--Ne** *<Ne_int>*
        Argument used to define the effective population size.
    **--genetic-map** *<genetic_map_filename>*
        Argument used to define a genetic map file.
    **--phase-chr** *<chr>*
        Argument used to define a single chromosome to phase.
    **--phase-from-bp**
        Argument used to define the lower bound of positions to include. May only be 
        used with a single chromosome.
    **--phase-to-bp**
        Argument used to define the upper bound of positions to include. May only be 
        used with a single chromosome.
    **--random-seed** *<seed_int>*
        Argument used to define the seed value for the random number generator.

    **************************************
    SHAPEIT Phasing Command-line Arguments
    **************************************
    **--shapeit-ref** *<ref_haps>* *<ref_legend>* *<ref_sample>*
        Argument used to define a reference panel. Three files are required: the reference 
        haplotypes (.haps), the snp map (.legend), and the individual information (.sample)
    **--shapeit-burn-iter** *<iteration_int>*
        Argument used to define the number of burn-in iterations.
    **--shapeit-prune-iter** *<iteration_int>*
        Argument used to define the number of pruning iterations.
    **--shapeit-main-iter** *<iteration_int>*
        Argument used to define the number of main iterations.
    **--shapeit-states** *<state_int>*
        Argument used to define the number of conditioning states for haplotype 
        estimation.
    **--shapeit-window** *<Mb_float>*
        Argument used to define the model window size in Mb.

    *************************************
    BEAGLE Phasing Command-line Arguments
    *************************************
    **--beagle-ref** *<ref_vcf, ref_bref3>*
        Argument used to define a reference panel VCF or bref3.
    **--beagle-burn-iter** *<iteration_int>*
        Argument used to define the number of burn-in iterations.
    **--beagle-iter** *<iteration_int>*
        Argument used to define the number of main iterations
    **--beagle-states** *<state_int>*
        Argument used to define the number of model states for genotype 
        estimation.
    **--beagle-error** *<probability>*
        Argument used to define the HMM allele mismatch probability.
    **--beagle-window** *<cM_float>*
        Argument used to define the sliding window size in cM.
    **--beagle-overlap** *<cM_float>*
        Argument used to define the overlap between neighboring windows in cM.
    **--beagle-step** *<cM_float>*
        Argument used to define the step length in cM used for identifying short 
        IBS segments.
    **--beagle-nsteps** *<windows_int>*
        Argument used to define the number of consecutive **--beagle-steps** used 
        for identifying long IBS segments.
    **--beagle-path** *<path>*
        Argument used to define the path to locate beagle.jar.    
"""

import os
import sys
import copy
import shutil
import argparse
import glob
import logging

from pgpipe.vcf_reader_func import checkFormat
from pgpipe.logging_module import initLogger, logArgs
from pgpipe.model import read_model_file
from pgpipe.beagle import call_beagle, check_for_beagle_intermediate_files
from pgpipe.shapeit import call_shapeit, remove_shapeit_intermediate_files, check_for_shapeit_intermediate_files
from pgpipe.bcftools import get_unique_chrs, get_samples, chr_subset_file, concatenate, check_for_index, create_index
from pgpipe.plink import convert_haps_to_vcf
from pgpipe.misc import argprase_kwargs

def phase_argument_parser(passed_arguments = []):
    '''
    VCF Phase Argument Parser

    Assign the parameters for VCF Phase using argparse.

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

    def parser_confirm_files ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):

                # Clean up any commas
                value = [item.strip(',') for item in value]
                
                # Loop the files
                for item in value:
                    if not os.path.isfile(item):
                        raise IOError('%s not found' % item)

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

    phase_parser = argparse.ArgumentParser()

    # Input arguments
    phase_parser.add_argument('--vcf', help = "Defines the filename of the VCF", type = str, required = True, action = parser_confirm_file())

    # Output arguments
    phase_parser.add_argument('--out', help = 'Defines the complete output filename, overrides --out-prefix')
    phase_parser.add_argument('--out-prefix', help = 'Defines the output prefix (i.e. filename without file extension)', default = 'out')

    out_format_list = ['vcf', 'vcf.gz', 'bcf']
    out_format_default = 'vcf.gz'

    phase_parser.add_argument('--out-format', metavar = metavar_list(out_format_list), help = 'Defines the desired output format', type = str, choices = out_format_list, default = out_format_default)

    # Model file arguments
    phase_parser.add_argument('--model-file', help = 'Defines the model file', type = str, action = parser_confirm_file())
    phase_parser.add_argument('--model', help = 'Defines the model and the individual(s) to include', type = str)

    # Non-Model individual selection
    phase_parser.add_argument('--filter-include-indv', help = 'Defines the individual(s) to include. This argument may be used multiple times if desired', nargs = '+', type = str, action = parser_add_to_list())
    phase_parser.add_argument('--filter-exclude-indv', help = 'Defines the individual(s) to exclude. This argument may be used multiple times if desired', nargs = '+', type = str, action = parser_add_to_list())
    phase_parser.add_argument('--filter-include-indv-file', help = 'Defines a file of individuals to include', action = parser_confirm_file())
    phase_parser.add_argument('--filter-exclude-indv-file', help = 'Defines a file of individuals to exclude', action = parser_confirm_file())

    # General arguments
    phase_parser.add_argument('--overwrite', help = "Overwrite previous output files", action = 'store_true')
    phase_parser.add_argument('--beagle-path', help = "Defines the path to locate beagle.jar", type = str)

    # Galaxy Option to pipe log to stdout
    phase_parser.add_argument('--log-stdout', help = argparse.SUPPRESS, action = 'store_true')

    # Phase algorithm argument
    phasing_list = ['beagle', 'shapeit']
    phasing_default = 'beagle'
    phase_parser.add_argument('--phase-algorithm', metavar = metavar_list(phasing_list), help = 'Specifies the phase algorithm to be used', type = str, choices = phasing_list, default = phasing_default)

    # Common phasing arguments
    phase_parser.add_argument('--Ne', help = 'Defines the effective population size', type = int)
    phase_parser.add_argument('--random-seed', help="Defines the seed value for the random number generator", type = int)
    phase_parser.add_argument('--genetic-map', help = 'Genetic map filename', type = str, action = parser_confirm_file())
    phase_parser.add_argument('--phase-chr', help = 'Selects a single chromosome to phase', type = str)
    phase_parser.add_argument('--phase-from-bp', help = 'Lower bound of sites to include. May only be used with a single chromosome', type = int)
    phase_parser.add_argument('--phase-to-bp', help = 'Upper bound of sites to include. May only be used with a single chromosome', type = int)

    # Shapeit-specific options
    phase_parser.add_argument('--shapeit-ref', help = 'Reference panel filenames (haps, legend, and sample)', type = str, nargs = 3, action = parser_confirm_files())
    phase_parser.add_argument('--shapeit-burn-iter', help = 'Number of the burn-in iterations (shapeit)', type = int)
    phase_parser.add_argument('--shapeit-prune-iter', help = 'Number of pruning iterations (shapeit)', type = int)
    phase_parser.add_argument('--shapeit-main-iter', help = 'Number of main iterations (shapeit)', type = int)
    phase_parser.add_argument('--shapeit-states', help = 'Number of conditioning states for haplotype estimation (shapeit)', type = int)
    phase_parser.add_argument('--shapeit-window', help = 'Model window size in Mb (shapeit)', type = float)
    phase_parser.add_argument('--shapeit-use-mt', help = 'Use chrMT rather than chrM', action = 'store_true')

    # Beagle-specific options
    phase_parser.add_argument('--beagle-ref', help = 'Reference panel filename (bref3 or VCF formats)', type = str, action = parser_confirm_file())
    phase_parser.add_argument('--beagle-burn-iter', help = 'Number of the burn-in iterations (beagle)', type = int)
    phase_parser.add_argument('--beagle-iter', help = 'Number of iterations after burn-in (beagle)', type = int)
    phase_parser.add_argument('--beagle-states', help = 'Number of model states for genotype estimation (beagle)', type = int)
    phase_parser.add_argument('--beagle-window', help = 'Sliding window size in cM (beagle)', type = float)
    phase_parser.add_argument('--beagle-overlap', help = 'Overlap between neighboring sliding windows in cM (beagle)', type = float)
    phase_parser.add_argument('--beagle-error', help = 'HMM allele mismatch probability (beagle)', type = float)
    phase_parser.add_argument('--beagle-step', help = 'Step length in cM used for identifying short IBS segments (beagle)', type = float)
    phase_parser.add_argument('--beagle-nsteps', help = 'Number of consecutive --beagle-steps used for identifying long IBS segments (beagle)', type = int)

    if passed_arguments:
        return vars(phase_parser.parse_args(passed_arguments))
    else:
        return vars(phase_parser.parse_args())

def log_to_stdout (log_filename):
    '''
        Write log file contents to stdout

        Used to allow galaxy users to see the contents of the log file. 
        Once completed, delete the log file.

        Raises
        ------
        IOError
            If a log file does not exist
    '''

    # Confirm the log file exists, raise an error if not
    if not os.path.isfile(log_filename):
        raise IOError('%s not found' % log_filename)

    # Open the log file
    log_file = open(log_filename)

    # Write the log file to stdout
    sys.stdout.write(log_file.readlines())

    # Close the log file
    log_file.close()

    # Delete the log file
    os.remove(log_filename)

def concatenate_logs (log_files, new_log_filename):
    '''
        Concatenates log files

        This function is used to concatenate log files. If a single chromosome
        is being analyzed, the function is used to concatenate the log files
        from multiple wrappers. If multiple chromosomes are being analyzed the
        function is used to concatenate the logs of each chromosome.

        Raises
        ------
        IOError
            If a log file does not exist
    '''

    # Open the new log file
    with open(new_log_filename, 'wb') as out_file:
        # Loop the logs to be combined
        for log_file in log_files:
            # Confirm the log file exists, raise an error if not
            if not os.path.isfile(log_file):
                raise IOError('%s not found' % log_file)
            # Open the log file
            with open(log_file, 'rb') as in_file:
                # Copy the contents of the file
                shutil.copyfileobj(in_file, out_file)

    # Remove old log files
    for log_file in log_files:
        os.remove(log_file)

def assign_vcf_extension (filename):
    # Checks if the file is unzipped, bgzipped, or gzipped
    vcfname_format = checkFormat(filename)

    if vcfname_format == 'vcf':
        return '.vcf'

    elif vcfname_format == 'gzip' or vcfname_format == 'bgzip':
        return '.vcf.gz'

    elif vcfname_format == 'bcf':
        raise Exception('BCF not supported in current version')

    else:
        raise Exception('Unknown file format')

def write_indv_file (filename, indv_list):

    # Open the individual file
    indv_file = open(filename, 'w')

    # Loop the list of individuals
    for indv in indv_list:

        # Write the individual to the file
        indv_file.write(indv + '\n')

    # Close the file
    indv_file.close()

def read_indv_file (filename):

    # List to hold individuals
    indv_list = []

    # Open the individual file
    with open(filename, 'r') as indv_file:

        # Read the individual lines
        for indv_line in indv_file:

            # Add the individual to the list
            indv_list.append(indv_line.strip())

    # Return the individuals
    return indv_list

def assign_filename_prefix (output_filename, output_format):

    '''
        Assigns a prefix using a filename

        Used to assign a unique prefix for jobs using the user specified
        filename. This is to avoid output files with the same prefix, 
        either from previous or ongoing jobs. This function is only used 
        if no prefix has been specified by the user.

        Parameters
        ----------
        output_filename : str
            Specifies the filename specified by the user
        output_format : str
            Specifies the file format suffix
        
        Returns
        -------
        unique_prefix: str
            Specifies the unqiue prefix

        Raises
        ------
        Exception
            If unable to assign a unique filename

    '''

    # Save the updated prefix
    updated_prefix = copy.deepcopy(output_filename)

    # Check if the file already exists 
    if os.path.isfile(updated_prefix + output_format):
        raise Exception('Unable to assign prefix output. %s already exists' % (updated_prefix + output_format))

    # Return the updated prefix
    return updated_prefix

def run (**kwargs):
    '''
    Phaser for VCF files.

    This function uses the argparse-based function :py:func:`phase_argument_parser`
    to parse either sys.argv or passed_arguments to obtain the parameters below. 
    The parameters are then translated to their BEAGLE/SHAPEIT equivalent. Once all 
    the parameters are assigned, BEAGLE or SHAPEIT is called.

    Parameters
    ----------
    --vcf : str
       Filename of the VCF
    --out : str
        Complete output filename, overrides --out-prefix
    --out-prefix : str
        Output prefix
    --out-format : str
        Desired output format
    --model-file : str
        Model filename
    --model : str
        Model to use
    --phase-algorithm : str
        The phase algorithm to use
    --Ne : int
        Effective population size
    --random-seed : int
        Seed value for the random number generator
    --genetic-map : str
        Genetic map filename
    --phase-chr : str
        Chromosome to phase
    --phase-from-bp : int
        Lower bound of sites to include. Only usable with a single 
        chromosome
    --phase-to-bp : int
        Upper bound of sites to include. Only usable with a single 
        chromosome
    --shapeit-ref : list
        Reference panel filenames (haps, legend, and sample)
    --shapeit-burn-iter : int
        Number of burn-in iterations
    --shapeit-prune-iter : int
        Number of pruning iterations
    --shapeit-main-iter : int
        Number of main iterations
    --shapeit-states : int
        Number of conditioning states for haplotype estimation
    --shapeit-window : float
        Model window size in Mb
    --beagle-ref : str
        Reference panel filename (bref3 or VCF formats)
    --beagle-burn-iter : int
        Number of the burn-in iterations
    --beagle-iter : int
        Number of main iterations
    --beagle-states : int
        Number of model states for genotype estimation
    --beagle-window : float
        Sliding window size in cM
    --beagle-overlap : float
        Overlap between neighboring sliding windows in cM
    --beagle-error : float
        HMM allele mismatch probability
    --beagle-step : float
        Step length in cM used for identifying short IBS segments
    --beagle-nsteps : int
        Number of consecutive --beagle-steps used for identifying 
        long IBS segments

    Raises
    ------
    IOError
        Output file already exists and --overwrite is not specified
    Exception
        Incompatible arguments
    '''

    # Update kwargs with defaults
    if __name__ != "__main__":
        kwargs = argprase_kwargs(kwargs, phase_argument_parser)

    # Assign arguments
    phase_args = argparse.Namespace(**kwargs)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(phase_args, func_name = 'vcf_phase')

    # Assign file extension for VCF input file
    vcfname_ext = assign_vcf_extension(phase_args.vcf)

    logging.info('Input file assigned')

    # Used to confirm if the VCF input was renamed
    vcfname_renamed = False

    # Confirm input has correct file extension
    if vcfname_ext not in phase_args.vcf:

        # Confirm the file has been renamed
        vcfname_renamed = True

        # Rename the file
        os.rename(phase_args.vcf, phase_args.vcf + vcfname_ext)
        
        # Update the vcf argument with the updated name
        phase_args.vcf += vcfname_ext

    # List to hold beagle/shapeit arguments
    phase_call_args = []

    # bool to check if an eclude file was created
    exclude_file_created = False

    # bool to check if an include file was created
    include_file_created = False

    # String to hold filename of filter individuals, if created
    filter_indv_filename = ''

    # Check if the has specified the output filename, without a prefix
    if phase_args.out:
        
        # Assign a prefix based on the output filename
        phase_args.out_prefix = assign_filename_prefix(phase_args.out, phase_args.out_format)

    # Assign expected phased output filename
    phased_output = '%s.%s' % (phase_args.out_prefix, phase_args.out_format)

    # Check if previous output should not be overwriten
    if not phase_args.overwrite:

        # Check if the user has specified an output filename
        if phase_args.out:

            # Check if previous output exists
            if os.path.exists(phase_args.out):
                raise IOError('Phased output already exists')

        else:

            # Check if previous output exists
            if os.path.exists(phased_output):
                raise IOError('Phased output already exists')

    # Performs actions related to --model-file and argument assignment
    if phase_args.model_file:

        # Check if a model has been specified
        if not phase_args.model:
            raise IOError('No selected model. Please use --model to select a model')

        # Read in the model file
        models_file = read_model_file(phase_args.model_file)

        # Check that the selected model was found in the file
        if phase_args.model not in models_file:
            raise IOError('Selected model "%s" not found in: %s' % (phase_args.model, phase_args.model_file))

        logging.info('Model assigned')

        if phase_args.phase_algorithm == 'beagle':

            # Get list of individuals to include
            inds_in_model = models_file[phase_args.model].inds

            # Create file with individuals to exclude (beagle has no include option)
            models_file.create_exclude_ind_file(inds_to_include = inds_in_model, overwrite = phase_args.overwrite)

            # Assign exclude file
            phase_call_args.append('excludesamples=' + models_file.exclude_file)

            # Confirm that the exclude file was created
            exclude_file_created = True

        elif phase_args.phase_algorithm == 'shapeit':

            # Assign the model
            selected_model = models_file[phase_args.model]

            # Create file with individuals to exclude (beagle has no include option)
            selected_model.create_ind_file(overwrite = phase_args.overwrite)

            # Assign exclude file
            phase_call_args.extend(['--include-ind', selected_model.ind_file])

            # Confirm that the include file was created
            include_file_created = True

    # Check if individuals to include or exclude were assigned
    elif phase_args.filter_include_indv or phase_args.filter_exclude_indv or phase_args.filter_include_indv_file or phase_args.filter_exclude_indv_file:

        # List of individuals to include
        include_indv_list = []

        # List of individuals to exclude
        exclude_indv_list = []

        # Check for individuals to include from the command-line
        if phase_args.filter_include_indv:

            # Add the individuals to the list
            include_indv_list.extend(phase_args.filter_include_indv)

        # Check for individuals to exclude from the command-line
        if phase_args.filter_exclude_indv:

            # Add the individuals to the list
            exclude_indv_list.extend(phase_args.filter_exclude_indv)

        # Check for individuals to include from a file
        if phase_args.filter_include_indv_file:

            # Read the individuals from the file
            include_file_indvs = read_indv_file(phase_args.filter_include_indv_file)

            # Add the individuals to the list
            include_indv_list.extend(include_file_indvs)

        # Check for individuals to exclude from a file
        if phase_args.filter_exclude_indv_file:

            # Read the individuals from the file
            exclude_file_indvs = read_indv_file(phase_args.filter_exclude_indv_file)

            # Add the individuals to the list
            exclude_indv_list.extend(exclude_file_indvs)

        # Assign the complete list of individuals within the file
        indvs_in_vcf = get_samples(phase_args.vcf)

        # Assign the individuals to include, checking both include/exclude lists
        indvs_to_include = [include_indv for include_indv in include_indv_list if include_indv not in exclude_indv_list]

        if phase_args.phase_algorithm == 'beagle':

            # Assign a filename for the indv file
            filter_indv_filename = 'filter_exclude_individuals'

            # Assign the individuals to exclude, using indvs_in_vcf and indvs_to_include 
            indvs_to_exclude = [indv_in_vcf for indv_in_vcf in indvs_in_vcf if indv_in_vcf not in indvs_to_include]

            # Create the indv file
            write_indv_file(filter_indv_filename, indvs_to_exclude)

            # Assign exclude file
            phase_call_args.append('excludesamples=' + filter_indv_filename)

        elif phase_args.phase_algorithm == 'shapeit':

            # Assign a filename for the indv file
            filter_indv_filename = 'filter_include_individuals'

            # Create the indv file
            write_indv_file(filter_indv_filename, indvs_to_include)

            # Assign exclude file
            phase_call_args.extend(['--include-ind', filter_indv_filename])

    # Get the list of chromosomes within the VCF
    chrs_in_vcf = get_unique_chrs(phase_args.vcf)

    # Check if the user specified chromosome is within the file
    if phase_args.phase_chr and (phase_args.phase_chr not in chrs_in_vcf):
        raise Exception('--phase-chr %s not found in %s' % (phase_args.phase_chr, phase_args.vcf))

    # Check if bp-based arguments are possible
    if (phase_args.phase_from_bp or phase_args.phase_to_bp) and (not phase_args.phase_chr and len(chrs_in_vcf) != 1):
        raise Exception('The --phase-from-bp and --phase-to-bp arguments require the --phase-chr '
                        'argument to function if multiple chrs are within %s' % phase_args.vcf)

    # Assign general arguments and call beagle
    if phase_args.phase_algorithm == 'beagle':

        # Check for the presence of intermediate files from beagle
        check_for_beagle_intermediate_files(phase_args.out_prefix, phase_args.out_format, overwrite = phase_args.overwrite)

        # Assigns the input and output arguments
        phase_call_args.extend(['gt=' + phase_args.vcf,
                                'out=' + phase_args.out_prefix])

        # Assign the beagle reference panel filename, if specified
        if phase_args.beagle_ref:
            phase_call_args.append('ref=' + phase_args.beagle_ref)

        # Assign the burn-in iter, if specified
        if phase_args.beagle_burn_iter:
            phase_call_args.append('burnin=' + str(phase_args.beagle_burn_iter))

        # Assign the iter, if specified
        if phase_args.beagle_iter:
            phase_call_args.append('iterations=' + str(phase_args.beagle_iter))

        # Assign the state count, if specified
        if phase_args.beagle_states:
            phase_call_args.append('phase-states=' + str(phase_args.beagle_states))

        # Assign the window length, if specified
        if phase_args.beagle_window:
            phase_call_args.append('window=' + str(phase_args.beagle_window))

        # Assign the window overlap, if specified
        if phase_args.beagle_overlap:
            phase_call_args.append('overlap=' + str(phase_args.beagle_overlap))

        # Assign the HMM error, if specified
        if phase_args.beagle_error:
            phase_call_args.append('err=' + str(phase_args.beagle_error))

        # Assign the step length, if specified
        if phase_args.beagle_step:
            phase_call_args.append('step=' + str(phase_args.beagle_step))

        # Assign the number of steps, if specified
        if phase_args.beagle_nsteps:
            phase_call_args.append('nsteps=' + str(phase_args.beagle_nsteps))

        # Assign the genetic map, if specified
        if phase_args.genetic_map:
            phase_call_args.append('map=' + phase_args.genetic_map)

        # Assign the effective pop size, if specified
        if phase_args.Ne:
            phase_call_args.append('ne=' + str(phase_args.Ne))

        # Assign the random seed, if specified
        if phase_args.random_seed:
            phase_call_args.append('seed=' + str(phase_args.random_seed))

        # Check if the file has a single chromosome
        if not phase_args.phase_chr and len(chrs_in_vcf) == 1:
                    
            # Assign the chromosome, as beagle requires it to be assigned
            phase_args.phase_chr = chrs_in_vcf[0]

            logging.info('Assigned chr %s for beagle' % chrs_in_vcf[0])

        # Assign the chromosome to phase, if specified
        if phase_args.phase_chr:

            # Store the chromosome argument
            chr_arg = 'chrom=%s' % phase_args.phase_chr

            # Check if either bp position arguments were specified
            if phase_args.phase_from_bp or phase_args.phase_to_bp:

                # List of the position arguments, in their required order
                position_args = [':', phase_args.phase_from_bp, '-', phase_args.phase_to_bp]

                # Filter the position arguments to remove empty values
                filttered_position_args = filter(None, position_args)

                # Map the arguments to str and add them to the chromosome argument
                chr_arg += ''.join(map(str, filttered_position_args))

            # Add the chromosome argument the argument list
            phase_call_args.append(chr_arg)

        # Call beagle wrapper
        call_beagle(phase_args.beagle_path, list(map(str, phase_call_args)), phase_args.out_prefix, phase_args.out_format)

        # Rename log using phased_output
        shutil.move(phase_args.out_prefix + '.log', phased_output + '.log')

    # Assign general arguments and call shapeit
    elif phase_args.phase_algorithm == 'shapeit':

        # Assign the shapeit reference panel filenames, if specified
        if phase_args.shapeit_ref:
            phase_call_args.append('--input-ref')
            phase_call_args.extend(phase_args.shapeit_ref)

        # Assign the shapeit burn in iter, if specified
        if phase_args.shapeit_burn_iter:
            phase_call_args.extend(['--burn', str(phase_args.shapeit_burn_iter)])

        # Assign the shapeit prune iter, if specified
        if phase_args.shapeit_prune_iter:
            phase_call_args.extend(['--prune', str(phase_args.shapeit_prune_iter)])

        # Assign the shapeit main iter, if specified
        if phase_args.shapeit_main_iter:
            phase_call_args.extend(['--main', str(phase_args.shapeit_main_iter)])

        # Assign the number of shapeit states if specified
        if phase_args.shapeit_states:
            phase_call_args.extend(['--states', str(phase_args.shapeit_states)])

        # Assign the shapeit window size, if specified
        if phase_args.shapeit_window:
            phase_call_args.extend(['--window', str(phase_args.shapeit_window)])

        # Assign the genetic map, if specified
        if phase_args.genetic_map:
            phase_call_args.extend(['--input-map', phase_args.genetic_map])

        # Assign the effective pop size, if specified
        if phase_args.Ne:
            phase_call_args.extend(['--effective-size', str(phase_args.Ne)])

        # Assign the random seed, if specified
        if phase_args.random_seed:
            phase_call_args.extend(['--seed', str(phase_args.random_seed)])

        # Check if only a single shapeit run is required
        if phase_args.phase_chr or len(chrs_in_vcf) == 1:

            # Check for the presence of intermediate files from shapeit
            check_for_shapeit_intermediate_files(phase_args.out_prefix, overwrite = phase_args.overwrite)

            # Assign the default shapeit input filename
            shapeit_input_vcf = phase_args.vcf

            # Check if a single chromosome is selected
            if phase_args.phase_chr and len(chrs_in_vcf) > 1:

                # Assign the chromosome-specific input prefix
                chr_input_prefix = phase_args.vcf + '.' + phase_args.phase_chr

                # Assign the chr-specific file as the shapeit input filename
                shapeit_input_vcf = chr_input_prefix + '.vcf.gz'

                # Create the chromosome-specific input
                chr_subset_file(phase_args.vcf, phase_args.phase_chr, chr_input_prefix, 'vcf.gz', 
                                from_bp = phase_args.phase_from_bp, to_bp = phase_args.phase_to_bp, 
                                overwrite = phase_args.overwrite)

                logging.info('Chr %s subset created' % phase_args.phase_chr)

            # Assigns the input and output arguments for shapeit
            phase_call_args.extend(['--input-vcf', shapeit_input_vcf,
                                    '--output-max', phase_args.out_prefix,
                                    '--output-log', phase_args.out_prefix + '.phase.log'])

            # Assign the from bp position argument, if specified
            if phase_args.phase_from_bp:
                phase_call_args.extend(['--input-from', phase_args.phase_from_bp])

            # Assign the to bp position argument, if specified
            if phase_args.phase_to_bp:
                phase_call_args.extend(['--input-to', phase_args.phase_to_bp])

            # Call shapeit wrapper
            call_shapeit(list(map(str, phase_call_args)), phase_args.out_prefix, phase_args.out_format)

            # Convert haps-format to vcf
            convert_haps_to_vcf(phase_args.out_prefix, phase_args.out_format, phase_args.shapeit_use_mt)

            logging.info('HAPS conversion to VCF complete')

            # Combine the log files
            concatenate_logs([phase_args.out_prefix + '.phase.log', phase_args.out_prefix + '.log'],  phased_output + '.log')

            # Check if a chr subset file was created
            if phase_args.phase_chr and len(chrs_in_vcf) > 1:
                # Delete the chromosome-specific input
                os.remove(shapeit_input_vcf)

                logging.info('Chr %s subset deleted' % phase_args.phase_chr)

            # Remove intermediate files created by shapeit
            remove_shapeit_intermediate_files(phase_args.out_prefix)

        # Check if multiple shapeit runs are required
        else:

            # List to store the phased filenames
            phased_filename_list = []

            # List to store the phased logs
            phased_log_list = []

            logging.info('Multi-chr shapeit phasing assigned')

            for chr_in_vcf in chrs_in_vcf:

                logging.info('Chr %s assigned' % chr_in_vcf)

                # Copy the arguments for this run
                chr_call_args = copy.deepcopy(phase_call_args)

                # Assign the chromosome-specific output prefix
                chr_out_prefix = phase_args.out_prefix + '.' + chr_in_vcf

                # Check for the presence of intermediate files from shapeit
                check_for_shapeit_intermediate_files(chr_out_prefix)

                # Assign the expected chromosome-specific output filename
                chr_out_filename = '%s.%s' % (chr_out_prefix, phase_args.out_format)

                # Store the output filename
                phased_filename_list.append(chr_out_filename)

                # Assign the chromosome-specific input prefix
                chr_input_prefix = phase_args.vcf + '.' + chr_in_vcf

                # Assign the expected chromosome-specific input filename
                chr_input_filename = chr_input_prefix + '.vcf.gz'

                # Create the chromosome-specific input
                chr_subset_file(phase_args.vcf, chr_in_vcf, chr_input_prefix, 'vcf.gz', 
                                overwrite = phase_args.overwrite)

                # Assigns the input and output arguments for shapeit
                chr_call_args.extend(['--input-vcf', chr_input_filename,
                                      '--output-max', chr_out_prefix,
                                      '--output-log', chr_out_prefix + '.phase.log'])

                # Call shapeit wrapper
                call_shapeit(list(map(str, chr_call_args)), chr_out_prefix, phase_args.out_format)

                # Combine log files
                concatenate_logs([chr_out_prefix + '.phase.log', chr_out_prefix + '.log'],  chr_out_filename + '.log')

                # Store the filename of the combined logs
                phased_log_list.append(chr_out_filename + '.log')

                # Delete the chromosome-specific input
                os.remove(chr_input_filename)

                # Remove intermediate files created by shapeit
                remove_shapeit_intermediate_files(chr_out_prefix)

            # Concatenate the vcf files
            concatenate(phased_filename_list, phase_args.out_prefix, phase_args.out_format)

            logging.info('Concatenated chromosomes')

            # Combine the log files
            concatenate_logs(phased_log_list, phased_output + '.log')

    # Check if the log should be piped to the stdout
    if phase_args.log_stdout:

        # Write the log to stdout
        log_to_stdout(phased_output + '.log')

        # Rename output to phase_args.out, if specified
        if phase_args.out:
            shutil.move(phased_output, phase_args.out)

    # Check if log should be saved as a file
    else:

        # Rename output to phase_args.out, if specified
        if phase_args.out:
            shutil.move(phased_output, phase_args.out)
            shutil.move(phased_output + '.log', phase_args.out + '.log')

        logging.info('Log file created')

    # Reverts the VCF input file
    if vcfname_renamed:
        os.rename(phase_args.vcf, phase_args.vcf[:-len(vcfname_ext)])

    # Check if a filter indv file was created
    if filter_indv_filename:

        # Delete the file
        os.remove(filter_indv_filename)

if __name__ == "__main__":
    initLogger()
    run(**phase_argument_parser())
