#!/usr/bin/env python
'''
    As a single statistic file may include far more loci/windows than
    a technique is capable of analyzing, it is often necessary to 
    sample the loci/windows from the file. Given a statistic file and 
    a sampling scheme, stat_sampler will generate a pseudorandomly
    sampled file.

    .. image:: ../../PPP_assets/PPP_STAT_Sample.png
        :width: 100 %
        :align: center

    In this illustration of the sampling process, the loci found within 
    Data.VCF are pseudorandomly sampled using the corrdinates found within 
    the given statistic file.

    Two pseudorandomly sampling schemes are provided: i) a random sampler 
    that will randomly select loci/windows and ii) a uniform sampler that 
    will evenly sample across equal-sized bins of the given statistic. 
    Please note that all sampling is done without replacement.

    For BED-based sampling, please see :doc:`../Utilities/bed_utilities.rst`.

    ##################
    Command-line Usage
    ##################
    The statistic sampler may be called using the following command:

    .. code-block:: bash
        
        stat_sampler.py

    *************
    Example usage
    *************
    Randomly sampling 20 windows from a windowed Fst statistic file 
    **merged_chr1_10000.windowed.weir.fst**.

    .. code-block:: bash
        
        stat_sampler.py --statistic-file examples/files/merged_chr1_10000.windowed.weir.fst --calc-statistic windowed-weir-fst --sampling-scheme random --sample-size 20

    Uniform sampling 20 windows from four bins from a windowed pi statistic file 
    **merged_chr1_10000.windowed.pi**.

    .. code-block:: bash
        
        stat_sampler.py --statistic-file examples/files/merged_chr1_10000.windowed.pi --calc-statistic window-pi --sampling-scheme uniform --uniform-bins 4 --sample-size 20

    ############################
    Input Command-line Arguments
    ############################
    **--statistic-file** *<statistic_filename>*
        Argument used to define the filename of the statistic file for sampling.

    #############################
    Output Command-line Arguments
    #############################
    **--out** *<output_filename>*
        Argument used to define the complete output filename, overrides **--out-prefix**.
        Cannot be used if multiple output files are created.
    **--out-prefix** *<output_prefix>*
        Argument used to define the output prefix (i.e. filename without file extension)
    **--overwrite**
        Argument used to define if previous output should be overwritten.

    ###############################
    Sampling Command-line Arguments
    ###############################
    **--calc-statistic** *<windowed-weir-fst, TajimaD, window-pi>*
        Argument used to define the statistic to be sampled. Windowed Fst 
        (windowed-weir-fst), Tajima's D (TajimaD), and windowed nucleotide 
        diversity (window-pi).
    **--sampling-scheme** *<random, uniform>*
        Argument used to define the sampling scheme. Random [Default]
        sampling or uniform sampling across of number of equal-sized
        bins.
    **--uniform-bins** *<bin_int>*
        Argument used to define the number of bins in uniform sampling.
    **--sample-size** *<sample_size_int>*
        Argument used to define the total sample size. If using the uniform 
        sampling scheme, this number must be divisible by the number of bins.
    **--random-seed** *<seed_int>*
        Argument used to define the seed value for the random number generator.
'''

import os
import sys
import argparse
import itertools
import logging
import shutil
import copy
import numpy as np
import pandas as pd
from collections import defaultdict

from pgpipe.logging_module import initLogger, logArgs
from pgpipe.vcftools import *
from pgpipe.bcftools import check_for_index, create_index
from pgpipe.model import read_model_file
from pgpipe.misc import argprase_kwargs

def sampler_parser(passed_arguments = []):
    '''
    Stat Sampler Phase Argument Parser

    Assign the parameters for stat sampler using argparse.

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
    sampler_parser.add_argument('--statistic-file', help='Defines the filename of the statistic file for sampling', required = True, type = str, action = parser_confirm_file())

    # Statistic based arguments.
    statistic_list = ['windowed-weir-fst', 'TajimaD', 'window-pi']
    sampler_parser.add_argument('--calc-statistic', metavar = metavar_list(statistic_list), help = 'Defines the statistic to be sampled', type=str, choices = statistic_list)

    # Sampling methods. Currently mutually exclusive to only allow a single sampling method
    sampling_list = ['uniform', 'random']
    sampling_default = 'random'
    sampler_parser.add_argument('--sampling-scheme', metavar = metavar_list(sampling_list), help = 'Defines the sampling scheme', type=str, choices = sampling_list, default = sampling_default)

    # Sampling options
    sampler_parser.add_argument('--uniform-bins', help="Defines the number of bins in uniform sampling", type = int)
    sampler_parser.add_argument('--sample-size', help="Defines the total sample size. If using the uniform sampling, this number must be divisible by the bins.", type = int)

    # Random options
    sampler_parser.add_argument('--random-seed', help="Defines the seed value for the random number generator", type = int, default = np.random.randint(1, high = 1000000000))

    # Output arguents
    sampler_parser.add_argument('--out', help = 'Defines the complete output filename, overrides **--out-prefix**', type = str)

    sampler_parser.add_argument('--out-prefix', help = 'Defines the output prefix (i.e. filename without file extension)', type = str, default = 'out')

    # General arguments
    sampler_parser.add_argument('--overwrite', help = "Defines if previous output should be overwritten", action = 'store_true')


    if passed_arguments:
        return vars(sampler_parser.parse_args(passed_arguments))
    else:
        return vars(sampler_parser.parse_args())

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

def run (**kwargs):
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

    Raises
    ------
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

    # Update kwargs with defaults
    if __name__ != "__main__":
        kwargs = argprase_kwargs(kwargs, sampler_parser)

    # Assign arguments
    sampler_args = argparse.Namespace(**kwargs)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(sampler_args, func_name = 'stat_sampler')

    # Set the random seed
    np.random.seed(sampler_args.random_seed)

    # Check if a sample size has been specified
    if not sampler_args.sample_size:
        raise Exception('No sample size specified. Please use --sample-size')

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

        # Check that uniform bins have been specified
        if not sampler_args.uniform_bins:
            raise Exception('No uniform bin count specified. Please use --uniform-bin')

        # Check that the sample size is divisible by the bins
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
    sampled_samples.to_csv(sampled_out_filename, sep = '\t', float_format = '%g', index = False)

    logging.info('Created selected samples file')

if __name__ == "__main__":
    initLogger()
    run(**sampler_parser())
