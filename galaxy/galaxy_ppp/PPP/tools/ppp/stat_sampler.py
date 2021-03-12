import os
import sys
import random
import argparse
import itertools
import logging
import shutil
import copy
import numpy as np
import pandas as pd
from collections import defaultdict

# Insert Jared's directory path, required for calling Jared's functions. Change when directory structure changes.
sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'jared')))

# Import log initializer
from logging_module import initLogger, logArgs

# Import basic vcftools functions
from vcftools import *

# Import basic vcf
from bcftools import check_for_index, create_index

# Model file related functions
from model import read_model_file

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

    # Other input arguments.
    sampler_parser.add_argument('--statistic-file', help='Specifies the statistic file for filtering', required = True, type = str, action = parser_confirm_file())

    # Statistic based arguments.
    statistic_list = ['windowed-weir-fst', 'TajimaD', 'window-pi']
    statistic_default = 'windowed-weir-fst'
    sampler_parser.add_argument('--calc-statistic', metavar = metavar_list(statistic_list), help = 'Specifies the statistic calculated ', type=str, choices = statistic_list, default = statistic_default)

    # Sampling methods. Currently mutually exclusive to only allow a single sampling method
    sampling_list = ['uniform', 'random']
    sampling_default = 'random'
    sampler_parser.add_argument('--sampling-scheme', metavar = metavar_list(sampling_list), help = 'Specifies the sampling scheme ', type=str, choices = sampling_list, default = sampling_default)

    # Sampling options
    sampler_parser.add_argument('--uniform-bins', help="Number of bins in uniform sampling", type = int)
    sampler_parser.add_argument('--sample-size', help="Total sample size. If using the uniform sampling scheme, this number must be divisible by the bins", type = int)

    # Random options
    sampler_parser.add_argument('--random-seed', help="Defines the random seed value for the random number generator", type = int, default = random.randint(1, 1000000000))

    # Output arguents
    sampler_parser.add_argument('--out', help = 'Specifies the filename of the sampled statistic file', type = str)

    sampler_parser.add_argument('--out-prefix', help = 'Specifies the filename prefix of the sampled statistic file', type = str, default = 'out')

    # General arguments
    sampler_parser.add_argument('--overwrite', help = "Specifies if previous output files should be overwritten", action = 'store_true')

    if passed_arguments:
        return sampler_parser.parse_args(passed_arguments)
    else:
        return sampler_parser.parse_args()

def random_vcftools_sampler (stat_file_data, sample_size, with_replacements = False):
    '''
        Random Sampler.

        Returns a list of randomly selected elements from stat_file_data. The
        length of this list is defined by sample_size.

        Parameters
        ----------
        stat_file_data : list
            List of statistic values
        sample_size : int
            Size of the total sample
        with_replacements : bool, optional
            If the function should sample with replacements

        Returns
        -------
        random_samples : list
            List of the randomly sampled statistics

    '''

    # Use numpy choice to randomly sample the stat_file_data
    random_samples = np.random.choice(stat_file_data, int(sample_size), replace = with_replacements)

    # Convert the samples to a list
    random_samples = list(random_samples)

    # Return the sample list
    return random_samples

def uniform_vcftools_sampler (stat_file_data, uniform_bins, sample_size):
    '''
        Uniform Sampler.

        Seperates the elements of stat_file_data into equal sized bins (number
        of bins specified by uniform_bins arg). The samples in the bins are
        then randomly selected. The number of selected samples is defined by
        sample_size / uniform_bins.

        Parameters
        ----------
        stat_file_data : list
            List of statistic values
        uniform_bins : int
            Number of uniform bins
        sample_size : int
            Size of the total sample

        Returns
        -------
        uniform_samples : list
            List of the uniformly sampled statistics
    '''

    # Create defaultdict to store the binned data
    binned_sample_database = defaultdict(list)

    # Bin the data and get the sizes of the uniform bin catagories
    binned_samples, bin_catagories = pd.cut(stat_file_data, uniform_bins, labels = False, retbins = True)

    # Loop the bins for logging purposes
    for bin_number, bin_edge_pos in enumerate(range(uniform_bins)):
        # Store the right edge of the bin
        right_edge = bin_catagories[bin_edge_pos]
        # Store the left edge of the bin
        left_edge = bin_catagories[bin_edge_pos + 1]
        # Log the bins for the user
        logging.info('Bin %s: %s < statistic <= %s' % (bin_number + 1, right_edge, left_edge))

    # Get the number of samples to sample per bin
    samples_per_bin = sample_size / uniform_bins

    # Loop the binned data
    for vcftools_sample, sample_bin in zip(range(len(stat_file_data)), binned_samples):
        # Populate the defaultdict with the bins
        binned_sample_database[sample_bin].append(vcftools_sample)

    # Create a list to hold the sampled data
    uniform_samples = []

    # Loop the defaultdict
    for sample_bin, current_samples in binned_sample_database.items():

        # Check if the current bin needs to be sampled
        if len(current_samples) >= samples_per_bin:

            # Randomly sample the bin
            bin_samples = random_vcftools_sampler(current_samples, samples_per_bin)

        # Check if the bin is too small for sampling
        else:
            # Return the current_samples
            bin_samples = current_samples

        uniform_samples.extend(bin_samples)

    return uniform_samples

def assign_statistic_column (sample_headers, statistic):
    '''
        Assigns the statistic column.

        Determines the column with the associated values for the specified
        statistic. This function is required only for the uniform sampler, as
        the random sampler does not need to know the values to function. Due to
        the nature of the statistic files, a conversion system is required to
        determine the associated column header for each specified statistic.
        This is done using the statistic_converter variable, which will need to
        be updated along with --calc-statistic to add more statistics.

        Parameters
        ----------
        sample_headers : list
            List of headers from the statistic file
        statistic : str
            The specified statistic

        Returns
        -------
        uniform_samples : list
            List of the uniformly sampled statistics

        Raises
        ------
        Exception
            Statistic not found. Converter likely needs to be updated
        Exception
            Specified statistic not found in headers. Likely wrong statistic
            was selected
    '''

    # Defaultdict to hold the statistic conversion data.
    # Conversion: --calc-statistic value : header value
    statistic_converter = {'windowed-weir-fst':'MEAN_FST',
                           'TajimaD':'TajimaD',
                           'window-pi':'PI'}

    # Return error if the statistic is not found
    if statistic not in statistic_converter:
        raise Exception('Statistic not found. Statistic list needs to be '
                        'updated. Please contact the PPP Team.')

    # Return error if the wrong statistic was selected
    if statistic_converter[statistic] not in sample_headers:
        raise ValueError('Statistic selected not found in file specified by '
                         '--statistic-file.')

    # Save the converted statistic
    converted_statistic = statistic_converter[statistic]

    # Return the converted statistic
    return converted_statistic

def run (passed_arguments = []):
    '''
        Statistics sampler.

        Automates the sampling process of a specified statistic output file. The
        function allows the user to select both the statistic in question and
        the sampling scheme. Please note that all sampling is done without
        replacement.

        Parameters
        ----------
        --statistic-file : str
            Specifies the statistic file for filtering
        --calc-statistic : str
            Specifies the statistic to calculate. Choices: windowed-weir-fst
            (default) and TajimaD
        --statistic-window-size : int
            Specifies the window size of the statistic if not specified in the
            file
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
        --out : str
            Specifies the sampled statistic output filename
        --out-prefix : str
            Specifies the sampled statistic output filename prefix

        Returns
        -------
        output : file
            Sampled statistic file
        log : file
            Log file output

        Raises
        ------
        IOError
            Statistic file does not exist
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
    logArgs(sampler_args, func_name = 'stat_sampler')

    # Set the random seed
    np.random.seed(sampler_args.random_seed)

    # Set the output filename using the out-prefix
    sampled_out_filename = sampler_args.out_prefix + '.sampled'

    # Check if the user has defined a specifed output filename
    if sampler_args.out:
        # Rename the output
        sampled_out_filename = sampler_args.out

    # Check if user has selected to overwrite previous output
    if not sampler_args.overwrite:

        # Check if output file already exists
        if os.path.isfile(sampled_out_filename):
            raise IOError('Sample file already exists')

    # Read in the sample file
    stat_file_data = pd.read_csv(sampler_args.statistic_file, sep = '\t')

    # Confirm statistic is correct, and assign the column ID
    assigned_statistic = assign_statistic_column(list(stat_file_data), sampler_args.calc_statistic)

    # Remove NaN datapoints for the statistic
    stat_file_data = stat_file_data.dropna(subset=[assigned_statistic])

    logging.info('Sample (i.e. statistic) file assigned')

    # Confirm there are enough samples
    if len(stat_file_data) < sampler_args.sample_size:
        raise ValueError('Sample size larger than the number of datapoints '
                         'within statistic file')

    # UPDATE Planned
    # Improve from equal to too few samples (i.e. samples = 100, sample_size = 95)
    elif len(stat_file_data) == sampler_args.sample_size:
        # Warns the user of poor sampling
        logging.warning('Sample size equal to number of datapoints within '
                        'statistic file')

    # Run the uniform sampler
    if sampler_args.sampling_scheme == 'uniform':
        if sampler_args.sample_size % sampler_args.uniform_bins != 0:
            raise ValueError('Sample size not divisible by the bin count')

        selected_samples = sorted(uniform_vcftools_sampler(list(stat_file_data[assigned_statistic]), sampler_args.uniform_bins, sampler_args.sample_size))

        if len(selected_samples) < sampler_args.sample_size:
            logging.warning('Uniform sampler was unable to reach sample size. '
                            'Likely too few samples per bin, consider reducing '
                            'the number of bins or the sample size.')

        logging.info('Uniform sampling complete')

    # Run the random sampler
    if sampler_args.sampling_scheme == 'random':

        selected_samples = sorted(random_vcftools_sampler(list(stat_file_data.index), sampler_args.sample_size))

        logging.info('Random sampling complete')

    # Reduce to only selected samples
    sampled_samples = stat_file_data[stat_file_data.index.isin(selected_samples)].copy()

    # Create selected samples TSV file with the defined output filename
    sampled_samples.to_csv(sampled_out_filename, sep = '\t', float_format = '%g')

    logging.info('Created selected samples file')

if __name__ == "__main__":
    initLogger()
    run()
