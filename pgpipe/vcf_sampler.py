import os
import sys
import random
import argparse
import pysam
import itertools
import logging
import shutil
import copy
import numpy as np
import pandas as pd
from collections import defaultdict

# Insert Jared's directory path, required for calling Jared's functions. Change when directory structure changes.
#sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'pppipe')))

# Import log initializer
from pgpipe.logging_module import initLogger, logArgs

# Import basic vcftools functions
from pgpipe.vcftools import *

# Import basic vcf
from pgpipe.bcftools import check_for_index, create_index

# Model file related functions
from pgpipe.model import read_model_file

def sampler_parser(passed_arguments):
    '''Sampler Argument Parser - Assigns arguments from command line.'''

    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string = None):
                if not os.path.isfile(value):
                    raise IOError('%s not found' % value)
                setattr(args, self.dest, value)
        return customAction

    def metavar_list (var_list):
        '''Create a formmated metavar list for the help output'''
        return '{' + ', '.join(var_list) + '}'

    sampler_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    # Input arguments
    sampler_parser.add_argument('--vcf', help = "Input VCF filename", type = str, required = True, action = parser_confirm_file())

    # Model file arguments.
    sampler_parser.add_argument('--model-file', help = 'Defines the model file', type = str, action = parser_confirm_file())
    sampler_parser.add_argument('--model', help = 'Defines the model to analyze', type = str)

    # Other input arguments.
    sampler_parser.add_argument('--statistic-file', help='Specifies the statistic file for filtering', required = True, type = str, action = parser_confirm_file())
    sampler_parser.add_argument("--vcf-index", help = "VCF Index filename", type = str, action = parser_confirm_file())

    # Output arguents
    sampler_parser.add_argument('--sample-file', help = 'Specifies the sampled (statistic file) tsv output filename', type = str, default = 'sampled_data.tsv')

    sampler_parser.add_argument('--out-dir', help = 'Specifies the VCF output directory', type = str, default = 'Sample_Files')
    sampler_parser.add_argument('--out-prefix', help = 'Specifies the VCF output filename prefix', type = str, default = 'Sample')

    out_format_list = ['vcf', 'vcf.gz', 'bcf']
    out_format_default = 'vcf.gz'

    sampler_parser.add_argument('--out-format', metavar = metavar_list(out_format_list), help = 'Specifies the output format.', type = str, choices = out_format_list, default = out_format_default)

    # General arguments
    sampler_parser.add_argument('--overwrite', help = "Specifies if previous output files should be overwritten", action = 'store_true')

    # Statistic based arguments.
    statistic_list = ['windowed-weir-fst', 'TajimaD', 'window-pi']
    statistic_default = 'windowed-weir-fst'
    sampler_parser.add_argument('--calc-statistic', metavar = metavar_list(statistic_list), help = 'Specifies the statistic calculated ', type=str, choices = statistic_list, default = statistic_default)

    sampler_parser.add_argument('--statistic-window-size', help = 'Specifies the size of window calculations (if BIN_END is absent)', type = int, default = 10000)

    # Sampling methods. Currently mutually exclusive to only allow a single sampling method
    sampling_list = ['uniform', 'random']
    sampling_default = 'random'
    sampler_parser.add_argument('--sampling-scheme', metavar = metavar_list(sampling_list), help = 'Specifies the sampling scheme ', type=str, choices = sampling_list, default = sampling_default)

    # Sampling options
    sampler_parser.add_argument('--uniform-bins', help="Number of bins in uniform sampling", type = int, default = 10)
    sampler_parser.add_argument('--sample-size', help="Total sample size. If using the uniform sampling scheme, this number must be divisible by the bins", type = int, default = 200)

    # Random options
    sampler_parser.add_argument('--random-seed', help="Defines the random seed value for the random number generator", type = int, default = random.randint(1, 1000000000))

    # Individual filters
    sampler_parser.add_argument('--filter-keep', help = 'Specifies a file of individuals to keep', action = parser_confirm_file())
    sampler_parser.add_argument('--filter-remove', help = 'Specifies a file of individuals to remove', action = parser_confirm_file())

    # Position-based position filters
    sampler_parser.add_argument('--filter-include-positions', help = 'Specifies a set of sites to include within a file (tsv chromosome and position)', action = parser_confirm_file())
    sampler_parser.add_argument('--filter-exclude-positions', help = 'Specifies a set of sites to exclude within a file (tsv chromosome and position)', action = parser_confirm_file())

    if passed_arguments:
        return sampler_parser.parse_args(passed_arguments)
    else:
        return sampler_parser.parse_args()

def random_vcftools_sampler (vcftools_samples, sample_size):
    '''Random Sampler. Returns a list of randomly selected elements from
    vcftools_samples. The length of this list is defined by sample_size.'''

    try:
        return list(np.random.choice(vcftools_samples, int(sample_size), replace = False))
    except ValueError:
        return vcftools_samples

def uniform_vcftools_sampler (vcftools_samples, uniform_bins, sample_size):
    '''Uniform Sampler. Seperates the elements of vcftools_samples into
    equal sized bins (number of bins specified by .uniform_bins arg). The
    samples in the bins are then randomly selected (number of selected
    samples is defined by: .sample_size / .uniform_bins).'''

    binned_sample_database = defaultdict(list)

    binned_samples = pd.cut(vcftools_samples, uniform_bins, labels = False)

    for vcftools_sample, sample_bin in zip(range(len(vcftools_samples)), binned_samples):
        binned_sample_database[sample_bin].append(vcftools_sample)

    uniform_samples = []

    for current_samples in binned_sample_database.values():
        uniform_samples.extend(random_vcftools_sampler(current_samples, sample_size / uniform_bins))

    return uniform_samples

def assign_position_columns (sample_headers):
    '''Column assignment. Returns the columns that contain the following
    headers: CHROM, BIN_START, BIN_END. Will also produce an error if not
    all the columns were found.'''

    if not all(pos_header in sample_headers for pos_header in ['CHROM', 'BIN_START', 'BIN_END']):
        # Assign missing headers
        missing_headers = [pos_header for pos_header in ['CHROM', 'BIN_START', 'BIN_END'] if pos_header not in sample_headers]
        # Print error message with missing headers
        raise ValueError('Cannot find %s column(s) in file specified by --statistic-file.' % ', '.join(missing_headers))

    return sample_headers.index('CHROM'), sample_headers.index('BIN_START'), sample_headers.index('BIN_END')

def assign_statistic_column (sample_headers, statistic):
    statistic_converter = {'windowed-weir-fst':'MEAN_FST', 'TajimaD':'TajimaD', 'window-pi':'PI'}

    if statistic not in statistic_converter:
        raise Exception('Statistic not found. Statistic list needs to be updated. Please contact the PPP Team.')

    if statistic_converter[statistic] not in sample_headers:
        raise ValueError('Statistic selected not found in file specified by --statistic-file.')

    return statistic_converter[statistic]

def run (passed_arguments = []):
    '''
        Statistics sampler for VCFTools output.

        Automates the sampling process of a specified statistic from VCFTools
        output files. The function allows the user to select both the statistic
        in question (windowed Fst [default] and Tajima's D) and the sampling
        scheme(random [default] and uniform).

        Parameters
        ----------
        --vcf : str
            Specifies the input VCF filename
        --out : str
            Specifies the VCF output filename
        --sample-file : str
            Specifies the sampled (statistic file) tsv output filename
        --statistic-file : str
            Specifies the statistic file for filtering
        --calc-statistic : str
            Specifies the statistic to calculate. Choices: windowed-weir-fst
            (default) and TajimaD
        --statistic-window-size : int
            Specifies the window size of the statistic. Required for Tajima's D
        --sampling-scheme : str
            Specifies the sampling scheme to use. Choices: random (default) and
            uniform
        --uniform-bins : int
            Specifies the number of bins for the uniform sampler
        --sample-size : int
            Specifies the total sample size. Note: If using the uniform sampling
            scheme, this number must be divisible by number of uniform bins
        --random-seed : int
            Specifies the random seed value for the random number generator

        Returns
        -------
        output : file
            Sampled statistic file
        samples : file
            Sampled VCF output

        Raises
        ------
        IOError
            Input VCF file does not exist
        IOError
            Output file already exists
        ValueError
            Sample size larger than the number of samples
        ValueError
            Sample size not divisible by the bin count
        ValueError
            Window size argument conflicts with values
        TypeError
            Window size argument not defined (if necessary)
    '''

    # Get arguments from command line
    sampler_args = sampler_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(sampler_args, func_name = 'vcf_sampler')

    # Set the random seed
    np.random.seed(sampler_args.random_seed)

    # Check if user has selected to overwrite previous output
    if not sampler_args.overwrite:
        # Check if output dir already exists
        if os.path.exists(sampler_args.out_dir):
            raise IOError('Sample Directory already exists')

        # Check if output file already exists
        if os.path.isfile(sampler_args.sample_file):
            raise IOError('Sample file already exists')
    else:
        # If previous output is to be overwritten, the previous output directory needs to be removed
        if os.path.exists(sampler_args.out_dir):
            shutil.rmtree(sampler_args.out_dir)

    # Read in the sample file
    vcftools_samples = pd.read_csv(sampler_args.statistic_file, sep = '\t')

    # Confirm statistic is correct, and assign the column ID
    assigned_statistic = assign_statistic_column(list(vcftools_samples), sampler_args.calc_statistic)

    # Remove NaN datapoints for the statistic
    vcftools_samples = vcftools_samples.dropna(subset=[assigned_statistic])

    logging.info('Sample (i.e. statistic) file assigned')

    # Confirm there are enough samples
    if len(vcftools_samples) < sampler_args.sample_size:
        raise ValueError('Sample size larger than the number of samples within sample file')

    # UPDATE Planned
    # Improve from equal to too few samples (i.e. samples = 100, sample_size = 95)
    elif len(vcftools_samples) == sampler_args.sample_size:
        # Warns the user of poor sampling
        logging.warning('Sample size equal to number of samples within sample file')

    # Run the uniform sampler
    if sampler_args.sampling_scheme == 'uniform':
        if sampler_args.sample_size % sampler_args.uniform_bins != 0:
            raise ValueError('Sample size not divisible by the bin count')

        selected_samples = sorted(uniform_vcftools_sampler(list(vcftools_samples[assigned_statistic]), sampler_args.uniform_bins, sampler_args.sample_size))

        if len(selected_samples) < sampler_args.sample_size:
            logging.warning('Uniform sampler was unable to reach sample size. Likely too few samples per bin, consider reducing the number of bins')

        logging.info('Uniform sampling complete')

    # Run the random sampler
    if sampler_args.sampling_scheme == 'random':

        selected_samples = sorted(random_vcftools_sampler(list(vcftools_samples.index), sampler_args.sample_size))

        logging.info('Random sampling complete')

    # Reduce to only selected samples
    sampled_samples = vcftools_samples[vcftools_samples.index.isin(selected_samples)].copy()

    # Create selected samples TSV file, with either the default filename or a user-defined filename
    sampled_samples.to_csv(sampler_args.sample_file, sep = '\t', float_format = '%g')

    logging.info('Created selected samples file')

    # Create the vcf/bcf output directory
    if not os.path.exists(sampler_args.out_dir):
        os.makedirs(sampler_args.out_dir)

    # Check if the data was sampled using Tajima's D. Adds BIN_END columns
    # that is required to fetch intervals from vcf files. Will need to be
    # edited when we add more statistics.
    if sampler_args.calc_statistic == 'TajimaD':
        # Check that the window size used in vcf_calc has been defined
        if not sampler_args.statistic_window_size:
            raise TypeError("--statistic-window-size argument required for the Tajima's D statistic")

        # Create iterator of all unique bin combinations
        bin_start_iter = itertools.combinations(list(sampled_samples['BIN_START']), 2)

        # Check if the bin combinations are divisible by the window size (to check if the window size is likely correct)
        if not all([abs(bin_1 - bin_2) % sampler_args.statistic_window_size == 0  for bin_1, bin_2 in bin_start_iter]):
            raise ValueError("--statistic-window-size argument conflicts with values in sample file")

        # Create list of bin start positions
        bin_start_list = list(sampled_samples['BIN_START'])

        # Create list of bin end positions
        bin_end_list = [start_pos + (sampler_args.statistic_window_size - 1) for start_pos in bin_start_list]

        # Add bin end positions to sampled_samples
        sampled_samples['BIN_END'] = pd.Series(bin_end_list, index = sampled_samples.index)

    # Open the VCF file
    if sampler_args.vcf_index:
        # Read in the file using pysam
        vcf_input = pysam.VariantFile(sampler_args.vcf, index_filename = sampler_args.vcf_index)
    else:
        # Check if there is an index file
        if check_for_index(sampler_args.vcf) == False:
            # Create an index if not found
            create_index(sampler_args.vcf)
        # Read in the file using pysam
        vcf_input = pysam.VariantFile(sampler_args.vcf)

    # Check if the user has specified a model file
    if sampler_args.model_file and sampler_args.model:

        # Read in the models
        models_in_file = read_model_file(sampler_args.model_file)

        # Check that the selected model was not found in the file
        if sampler_args.model not in models_in_file:
            raise IOError('Selected model "%s" not found in: %s' % (sampler_args.model, sampler_args.model_file))

        # Select model, might change this in future versions
        selected_model = models_in_file[sampler_args.model]

        # List of the individuals within the vcf
        individuals_in_file = list(vcf_input.header.samples)

        if not set(selected_model.inds).issubset(set(individuals_in_file)):
            raise IOError('Individuals differ between model "%s" and input.' % sampler_args.model)

        # Reduce individuals to those in model
        vcf_input.subset_samples(selected_model.inds)

    # If no model was specifed, check if a keep or remove file was specifed
    elif sampler_args.filter_keep or sampler_args.filter_remove:

        # List of the individuals within the vcf
        individuals_in_file = list(vcf_input.header.samples)

        # Check if a keep file was specifed
        if sampler_args.filter_keep:

            # List of the individuals to keep
            selected_individuals = []

            # Read in the individuals within the keep file
            with open(sampler_args.filter_keep, 'rU') as keep_file:
                for keep_line in keep_file:

                    if keep_line.strip() not in individuals_in_file:
                        raise IOError('Individuals differ between --keep-file and input.')

                    selected_individuals.append(keep_line.strip())

        # Check if a remove file was specifed
        elif sampler_args.filter_remove:

            # List of the individuals to keep
            selected_individuals = copy.deepcopy(individuals_in_file)

            # Read in the individuals within the keep file
            with open(sampler_args.filter_remove, 'rU') as remove_file:
                for remove_line in remove_file:

                    if remove_line.strip() not in individuals_in_file:
                        raise IOError('Individuals differ between --remove-file and input.')

                    selected_individuals.remove(remove_line.strip())

        # Reduce individuals to those in model
        vcf_input.subset_samples(selected_individuals)

    # Sites to be included
    include_sites =  pd.DataFrame()

    # Check if the user has specified an include sites file
    if sampler_args.filter_include_positions:
        include_sites = pd.read_csv(sampler_args.filter_include_positions, sep = '\t')

    # Sites to be excluded
    exclude_sites = pd.DataFrame()

    # Check if the user has specified an exclude sites file
    if sampler_args.filter_exclude_positions:
        exclude_sites = pd.read_csv(sampler_args.filter_exclude_positions, sep = '\t')

    # Get the chromosome, start, and end columns
    chr_col, start_col, end_col = assign_position_columns(list(sampled_samples))

    # iterate the selected samples
    for sampled_count, sampled_row in enumerate(sampled_samples.values):

        # Assign filename for sample.
        if sampler_args.out_format == 'vcf.gz':
            sample_filename =  sampler_args.out_prefix + '_%s.' %sampled_count + 'vcf'
        else:
            sample_filename =  sampler_args.out_prefix + '_%s.' %sampled_count + sampler_args.out_format

        # Join the output directory and the sample filename paths
        sample_path = os.path.join(sampler_args.out_dir, sample_filename)

        # Create the VCF output file, with either the default filename or a user-defined filename
        vcf_output = pysam.VariantFile(sample_path, 'w', header = vcf_input.header)

        # Fetch positions specified from the vcf input file
        for vcf_record in vcf_input.fetch(sampled_row[chr_col], int(sampled_row[start_col]), int(sampled_row[end_col])):

            # Bool to determine if the record should be saved
            save_record = True

            # Check if the user has specified an include sites file
            if sampler_args.filter_include_positions:
                # Check if the current record is not an included site
                if not ((include_sites['CHROM'] == vcf_record.chrom) & (include_sites['POS'] == vcf_record.pos)).any():
                    save_record = False

            # Check if the user has specified an exclude sites file
            if sampler_args.filter_exclude_positions:
                # Check if the current record is an excluded site
                if ((exclude_sites['CHROM'] == vcf_record.chrom) & (exclude_sites['POS'] == vcf_record.pos)).any():
                    save_record = False

            # Check if the record is to be saved
            if save_record:
                vcf_output.write(vcf_record)

        vcf_output.close()

        # Compress vcf to vcf.gz using bgzip
        if sampler_args.out_format == 'vcf.gz':
            bgzip_compress_vcf(sample_path)

    vcf_input.close()


if __name__ == "__main__":
    initLogger()
    run()
