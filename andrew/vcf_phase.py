import os
import sys
import copy
import shutil
import argparse
import glob
import logging

sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'jared')))

from vcf_reader_func import checkFormat
from logging_module import initLogger, logArgs
from model import read_model_file
from beagle import call_beagle
from shapeit import call_shapeit, remove_intermediate_files
from bcftools import pipe_bcftools_to_chr, chr_subset_file, concatenate

def phase_argument_parser(passed_arguments):
    '''Phase Argument Parser - Assigns arguments for vcftools from command line.
    Depending on the argument in question, a default value may be specified'''

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

    phase_parser = argparse.ArgumentParser()

    # Input arguments
    phase_parser.add_argument('--vcf', help = "Input VCF filename", type = str, required = True, action = parser_confirm_file())

    # Model file arguments
    phase_parser.add_argument('--model-file', help = 'Defines the model file', type = str, action = parser_confirm_file())
    phase_parser.add_argument('--model', help = 'Defines the model to analyze', type = str)

    # General arguments
    phase_parser.add_argument('--overwrite', help = "Overwrite previous output files", action = 'store_true')
    phase_parser.add_argument('--beagle-path', help = "Defines path to locate beagle.jar", type = str, default = 'bin/')

    # Phase algorithm argument
    phasing_list = ['beagle', 'shapeit']
    phasing_default = 'beagle'
    phase_parser.add_argument('--phase-algorithm', metavar = metavar_list(phasing_list), help = 'Specifies the phase algorithm to be used', type = str, choices = phasing_list, default = phasing_default)

    # Common phasing arguments
    phase_parser.add_argument('--Ne', help = 'Defines the effective population size', type = int)
    phase_parser.add_argument('--random-seed', help="Defines the random seed value for the random number generator", type = int)
    phase_parser.add_argument('--genetic-map', help = 'Genetic map filename', type = str, action = parser_confirm_file())
    phase_parser.add_argument('--phase-chr', help = 'Selects a single chromosome to phase', type = str)
    phase_parser.add_argument('--phase-from-bp', help = 'Lower bound of sites to include (Only usable with a single chromosome)', type = int)
    phase_parser.add_argument('--phase-to-bp', help = 'Upper bound of sites to include (Only usable with a single chromosome)', type = int)

    # Shapeit-specific options
    phase_parser.add_argument('--shapeit-burn-iter', help = 'Number of the burn-in iterations (shapeit)', type = int)
    phase_parser.add_argument('--shapeit-prune-iter', help = 'Number of pruning iterations (shapeit)', type = int)
    phase_parser.add_argument('--shapeit-main-iter', help = 'Number of main iterations (shapeit)', type = int)
    phase_parser.add_argument('--shapeit-states', help = 'Number of conditioning states for haplotype estimation (shapeit)', type = int)
    phase_parser.add_argument('--shapeit-window', help = 'Model window size in Mb (shapeit)', type = float)

    # Beagle-specific options
    phase_parser.add_argument('--beagle-burn-iter', help = 'Number of the burn-in iterations (beagle)', type = int)
    phase_parser.add_argument('--beagle-iter', help = 'Number of iterations after burn-in (beagle)', type = int)
    phase_parser.add_argument('--beagle-states', help = 'Number of model states for genotype estimation (beagle)', type = int)
    phase_parser.add_argument('--beagle-window', help = 'Sliding window size in cM (beagle)', type = float)
    phase_parser.add_argument('--beagle-overlap', help = 'Overlap between neighboring sliding windows in cM (beagle)', type = float)
    phase_parser.add_argument('--beagle-error', help = 'HMM allele mismatch probability (beagle)', type = float)
    phase_parser.add_argument('--beagle-step', help = 'Step length in cM used for identifying short IBS segments (beagle)', type = float)
    phase_parser.add_argument('--beagle-nsteps', help = 'Number of consecutive --beagle-steps used for identifying long IBS segments (beagle)', type = int)

    # Output arguments
    phase_parser.add_argument('--out', help = 'Defines the output filename')
    phase_parser.add_argument('--out-prefix', help = 'Defines the output prefix (used by phasing algorithms)', default = 'out')

    out_format_list = ['vcf', 'vcf.gz', 'bcf']
    out_format_default = 'vcf.gz'

    phase_parser.add_argument('--out-format', metavar = metavar_list(out_format_list), help = 'Specifies the output format.', type = str, choices = out_format_list, default = out_format_default)

    if passed_arguments:
        return phase_parser.parse_args(passed_arguments)
    else:
        return phase_parser.parse_args()

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

def run (passed_arguments = []):
    '''
        Phaser for VCF files.

        Automates the phasing process for a specified VCF file. The function
        allows users to select between multiple phasing algorithms: beagle
        (default) and shapit.

        Parameters
        ----------
        --vcf : str
            Specifies the input VCF filename
        --phase-algorithm : str
            Specifies the algorithm to be used. Choices: beagle (default) and
            shapit
        --out : str
            Specifies the output filename

        Returns
        -------
        output : file
            Phased VCF file

        Raises
        ------
        IOError
            Input VCF file does not exist
        IOError
            Output file already exists

    '''

    # Grab VCF arguments from command line
    phase_args = phase_argument_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(phase_args, func_name = 'vcf_phase')

    # Assign file extension for VCF input file
    vcfname_ext = assign_vcf_extension(phase_args.vcf)

    logging.info('Input file assigned')

    # Used to confirm if the VCF input was renamed
    vcfname_renamed = False

    # Confirm input has correct file extension
    if vcfname_ext not in phase_args.vcf:
        vcfname_renamed = True
        os.rename(phase_args.vcf, phase_args.vcf + vcfname_ext)
        phase_args.vcf += vcfname_ext

    # List to hold beagle/shapeit arguments
    phase_call_args = []

    # bool to check if an eclude file was created
    exclude_file_created = False

    # bool to check if an include file was created
    include_file_created = False

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

    # Get the list of chromosomes within the VCF
    chrs_in_vcf = pipe_bcftools_to_chr(phase_args.vcf)

    # Check if the user specified a specific chromosome
    if phase_args.phase_chr:
        # Check if the chromosome is within the file
        if phase_args.phase_chr not in chrs_in_vcf:
            raise Exception('--phase-chr %s not found in %s' % (phase_args.phase_chr, phase_args.vcf))

    # Check that a chr was specified if a bp-based argument was specified
    if (phase_args.phase_from_bp or phase_args.phase_to_bp) and not phase_args.phase_chr:
        # Check if the file has a single chromosome
        if len(chrs_in_vcf) == 1:
            # Check if the beagle algorithm is being used
            if phase_args.phase_algorithm == 'beagle':
                # Assign the chromosome, as beagle requires it to be assigned
                phase_args.phase_chr = chrs_in_vcf[0]

                logging.info('Assigned chr %s for beagle' % chrs_in_vcf[0])
        # Report error if multiple chromosomes are within the file
        else:
            raise Exception('The --phase-from-bp and --phase-to-bp arguments '
                            'require the --phase-chr argument to function if '
                            'multiple chrs are within %s' % phase_args.vcf)

    # Assign general arguments and call beagle
    if phase_args.phase_algorithm == 'beagle':

        # Assigns the input and output arguments
        phase_call_args.extend(['gt=' + phase_args.vcf,
                                'out=' + phase_args.out_prefix])

        # Assign the burn-in iter, if specified
        if phase_args.beagle_burn_iter:
            phase_call_args.append('burnin=' + phase_args.beagle_burn_iter)

        # Assign the iter, if specified
        if phase_args.beagle_iter:
            phase_call_args.append('iterations=' + phase_args.beagle_iter)

        # Assign the state count, if specified
        if phase_args.beagle_states:
            phase_call_args.append('phase-states=' + phase_args.beagle_states)

        # Assign the window length, if specified
        if phase_args.beagle_window:
            phase_call_args.append('window=' + phase_args.beagle_window)

        # Assign the window overlap, if specified
        if phase_args.beagle_overlap:
            phase_call_args.append('overlap=' + phase_args.beagle_overlap)

        # Assign the HMM error, if specified
        if phase_args.beagle_error:
            phase_call_args.append('err=' + phase_args.beagle_error)

        # Assign the step length, if specified
        if phase_args.beagle_step:
            phase_call_args.append('step=' + phase_args.beagle_step)

        # Assign the number of steps, if specified
        if phase_args.beagle_nsteps:
            phase_call_args.append('nsteps=' + phase_args.beagle_nsteps)

        # Assign the genetic map, if specified
        if phase_args.genetic_map:
            phase_call_args.append('map=' + phase_args.genetic_map)

        # Assign the effective pop size, if specified
        if phase_args.Ne:
            phase_call_args.append('ne=' + phase_args.Ne)

        # Assign the random seed, if specified
        if phase_args.random_seed:
            phase_call_args.append('seed=' + phase_args.random_seed)

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

        # Assign expected phased output filename
        phased_output = '%s.%s' % (phase_args.out_prefix, phase_args.out_format)

        # Rename output to phase_args.out, if specified
        if phase_args.out:
            shutil.move(phased_output, phase_args.out)
            shutil.move(phase_args.out_prefix + '.log', phase_args.out + '.log')
        # Rename log using phased_output
        else:
            shutil.move(phase_args.out_prefix + '.log', phased_output + '.log')

        logging.info('beagle log file created')

    # Assign general arguments and call shapeit
    elif phase_args.phase_algorithm == 'shapeit':

        # Assign the shapeit burn in iter, if specified
        if phase_args.shapeit_burn_iter:
            phase_call_args.extend(['--burn', phase_args.shapeit_burn_iter])

        # Assign the shapeit prune iter, if specified
        if phase_args.shapeit_prune_iter:
            phase_call_args.extend(['--prune', phase_args.shapeit_prune_iter])

        # Assign the shapeit main iter, if specified
        if phase_args.shapeit_main_iter:
            phase_call_args.extend(['--main', phase_args.shapeit_main_iter])

        # Assign the number of shapeit states if specified
        if phase_args.shapeit_states:
            phase_call_args.extend(['--states', phase_args.shapeit_states])

        # Assign the shapeit window size, if specified
        if phase_args.shapeit_window:
            phase_call_args.extend(['--window', phase_args.shapeit_window])

        # Assign the genetic map, if specified
        if phase_args.genetic_map:
            phase_call_args.extend(['--input-map', phase_args.genetic_map])

        # Assign the effective pop size, if specified
        if phase_args.Ne:
            phase_call_args.extend(['--effective-size', phase_args.Ne])

        # Assign the random seed, if specified
        if phase_args.random_seed:
            phase_call_args.extend(['--seed', phase_args.random_seed])

        # Check if only a single shapeit run is required
        if phase_args.phase_chr or len(chrs_in_vcf) == 1:

            # Holds shapeit input filename, as a temp file may be required
            shapeit_input_vcf = None

            # Check if a single chromosome is selected
            if phase_args.phase_chr and len(chrs_in_vcf) > 1:
                # Assign the chromosome-specific input prefix
                chr_input_prefix = phase_args.vcf + '.' + phase_args.phase_chr

                # Assign the chr-specific file as the shapeit input filename
                shapeit_input_vcf = chr_input_prefix + '.vcf.gz'

                # Create the chromosome-specific input
                chr_subset_file(phase_args.vcf,
                                phase_args.phase_chr,
                                chr_input_prefix,
                                'vcf.gz',
                                from_bp = phase_args.phase_from_bp,
                                to_bp = phase_args.phase_to_bp)

                logging.info('Chr %s subset created' % phase_args.phase_chr)

            else:
                # Assign the shapeit input filename
                shapeit_input_vcf = phase_args.vcf

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

            # Assign expected phased output filename
            phased_output = '%s.%s' % (phase_args.out_prefix, phase_args.out_format)

            # Combine the log files
            concatenate_logs([phase_args.out_prefix + '.phase.log', phase_args.out_prefix + '.log'],  phased_output + '.log')

            # Rename output to phase_args.out, if specified
            if phase_args.out:
                shutil.move(phased_output, phase_args.out)
                shutil.move(phased_output + '.log',  phase_args.out + '.log')

            logging.info('shapeit log file created')

            # Check if a chr subset file was created
            if phase_args.phase_chr and len(chrs_in_vcf) > 1:
                # Delete the chromosome-specific input
                os.remove(shapeit_input_vcf)

                logging.info('Chr %s subset deleted' % phase_args.phase_chr)

            # Remove intermediate files created by shapeit
            remove_intermediate_files(phase_args.out_prefix)

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

                # Assign the expected chromosome-specific output filename
                chr_out_filename = '%s.%s' % (chr_out_prefix, phase_args.out_format)

                # Store the output filename
                phased_filename_list.append(chr_out_filename)

                # Assign the chromosome-specific input prefix
                chr_input_prefix = phase_args.vcf + '.' + chr_in_vcf

                # Assign the expected chromosome-specific input filename
                chr_in_filename = chr_input_prefix + '.vcf.gz'

                # Create the chromosome-specific input
                chr_subset_file(phase_args.vcf,
                                chr_in_vcf,
                                chr_input_prefix,
                                'vcf.gz')

                # Assigns the input and output arguments for shapeit
                chr_call_args.extend(['--input-vcf', chr_in_filename,
                                      '--output-max', chr_out_prefix,
                                      '--output-log', chr_out_prefix + '.phase.log'])

                # Call shapeit wrapper
                call_shapeit(list(map(str, chr_call_args)), chr_out_prefix, phase_args.out_format)

                # Combine log files
                concatenate_logs([chr_out_prefix + '.phase.log', chr_out_prefix + '.log'],  chr_out_filename + '.log')

                # Store the filename of the combined logs
                phased_log_list.append(chr_out_filename + '.log')

                # Delete the chromosome-specific input
                os.remove(chr_in_filename)

                # Remove intermediate files created by shapeit
                remove_intermediate_files(chr_out_prefix)

            # Concatenate the vcf files
            concatenate(phased_filename_list, phase_args.out_prefix, phase_args.out_format)

            logging.info('Concatenated chromosomes')

            # Assign expected concatenated output filename
            phased_output = '%s.%s' % (phase_args.out_prefix, phase_args.out_format)

            # Combine the log files
            concatenate_logs(phased_log_list, phased_output + '.log')

            # Rename output to phase_args.out, if specified
            if phase_args.out:
                shutil.move(phased_output, phase_args.out)
                shutil.move(phased_output + '.log',  phase_args.out + '.log')

            logging.info('Multi-chr shapeit log created')

    # Reverts the VCF input file
    if vcfname_renamed:
        os.rename(phase_args.vcf, phase_args.vcf[:-len(vcfname_ext)])

if __name__ == "__main__":
    initLogger()
    run()
