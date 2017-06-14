import os, sys, csv, argparse, random, pysam
import numpy as np
import pandas as pd
from collections import defaultdict

def sampler_parser():
    '''Sampler Argument Parser - Assigns arguments from command line.'''

    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError # File not found
                setattr(args, self.dest, value)
        return customAction

    sampler_parser = argparse.ArgumentParser()

    # Input arguments
    sampler_parser.add_argument("vcfname", metavar='VCF_Input', help = "Input VCF filename", type = str, action = parser_confirm_file())
    sampler_parser.add_argument('--statistic-file', help='Specifies the statistic file for filtering', required = True, type = str, action = parser_confirm_file())

    # Statistic based arguments.
    statistic_list = ['windowed-weir-fst', 'TajimaD']
    statistic_default = 'windowed-weir-fst'
    sampler_parser.add_argument('--calc-statistic', metavar = '{' + ', '.join(statistic_list) + '}', help = 'Specifies the statistic calculated ', type=str, choices = statistic_list, default = statistic_default)

    sampler_parser.add_argument('--statistic-window-size', help = 'Specifies the size of window calculations', type = int, default = 10000)

    # Sampling methods. Currently mutually exclusive to only allow a single sampling method
    sampling_list = ['uniform', 'random']
    sampling_default = 'random'
    sampler_parser.add_argument('--sampling-scheme', metavar = '{' + ', '.join(sampling_list) + '}', help = 'Specifies the sampling scheme ', type=str, choices = sampling_list, default = sampling_default)

    # Sampling options
    sampler_parser.add_argument('--uniform-bins', help="Number of bins in uniform sampling", type = int, default = 10)
    sampler_parser.add_argument('--sample-size', help="Total sample size. If using the uniform sampling scheme, this number must be divisible by the bins", type = int, default = 200)

    # Other options
    sampler_parser.add_argument('--random-seed', help="Defines the random seed value for the random number generator", type = int, default = random.randint(1, 1000000000))

    return sampler_parser.parse_args()

def random_vcftools_sampler (vcftools_samples, sample_size):
    '''Random Sampler. Returns a list of randomly selected elements from
    vcftools_samples. The length of this list is defined by sample_size.'''

    try:
        return random.sample(vcftools_samples, sample_size)
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
    enough columns were found.'''

    if not all(pos_header in sample_headers for pos_header in ['CHROM', 'BIN_START', 'BIN_END']):
        if all(pos_header in sample_headers for pos_header in ['CHROM', 'BIN_START']):
            return sample_headers.index('CHROM'), sample_headers.index('BIN_START'), None
        else:
            if 'CHROM' in sample_headers:
                sys.exit('Cannot find BIN_START in file')
            else:
                sys.exit('Cannot find CHROM and BIN_START in file')
    return sample_headers.index('CHROM'), sample_headers.index('BIN_START'), sample_headers.index('BIN_END')

def assign_statistic_column (sample_headers, statistic):
    statistic_converter = {'windowed-weir-fst':'MEAN_FST', 'TajimaD':'TajimaD'}

    if not statistic_converter.has_key(statistic):
        sys.exit('Statistic not found. assign_statistic_column update needed')

    if statistic_converter[statistic] not in sample_headers:
        sys.exit('Statistic selected not found in file specified by --statistic-file')

    return statistic_converter[statistic]

def run ():
    '''
        Statistics sampler for VCFTools output.

        Automates the sampling process of a specified statistic from VCFTools
        output files. The function allows the user to select both the statistic
        in question (windowed Fst [default] and Tajima's D) and the sampling
        scheme(random [default] and uniform).

        Parameters
        ----------
        VCF_Input : str
            Specifies the input VCF filename
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

    '''

    # Get arguments from command line
    sampler_args = sampler_parser()
    # Set the random seed
    random.seed(sampler_args.random_seed)
    # Read in the sample file
    vcftools_samples = pd.read_csv(sampler_args.statistic_file, sep = '\t')

    # Confirm there are enough samples
    if len(vcftools_samples) < sampler_args.sample_size:
        sys.exit('Sample size larger than the number of samples within input')

    elif len(vcftools_samples) == sampler_args.sample_size:
        # Update with warning call
        print 'Sample size equal to number of samples within input'

    # Run the uniform sampler
    if sampler_args.sampling_scheme == 'uniform':
        if sampler_args.sample_size % sampler_args.uniform_bins != 0:
            sys.exit('Sample size not divisible by the bin count')

        assigned_statistic = assign_statistic_column (list(vcftools_samples), sampler_args.calc_statistic)
        selected_samples = sorted(uniform_vcftools_sampler(list(vcftools_samples[assigned_statistic]), sampler_args.uniform_bins, sampler_args.sample_size))

    # Run the random sampler
    if sampler_args.sampling_scheme == 'random':
        selected_samples = sorted(random_vcftools_sampler(range(len(vcftools_samples)), sampler_args.sample_size))

    # Reduce to only selected samples
    sampled_samples = vcftools_samples[vcftools_samples.index.isin(selected_samples)]

    # Create TSV file of the reduced samples
    sampled_samples.to_csv(sampler_args.statistic_file + '.sampled', sep = '\t')

    # If a VCF file is specifed, create a reduced VCF file(s) from the samples
    if sampler_args.vcfname:
        # Get the chromosome, start, and end columns
        chr_col, start_col, end_col = assign_position_columns(list(vcftools_samples))
        for sampled_row in vcftools_samples.values:
            if not end_col:
                if not sampler_args.statistic_window_size:
                    sys.argv("--statistic-window-size argument required for the Tajima's D statistic")

                print sampled_row[chr_col], sampled_row[start_col], sampled_row[start_col] + (sampler_args.statistic_window_size - 1)
            else:
                print sampled_row[chr_col], sampled_row[start_col], sampled_row[end_col]

        '''
        vcf_input = pysam.VariantFile(sampler_args.vcf_file)
        for sampled_row in vcftools_samples.values:
            if end_col:
                vcf_output = pysam.VariantFile('Sampled_{0}_{1}_{2}.vcf.gz'.format(sampled_row[chr_col], sampled_row[start_col], sampled_row[end_col]), 'w', header = vcf_input.header)
                for vcf_record in vcf_input.fetch(sampled_row[chr_col], int(sampled_row[start_col]), int(sampled_row[end_col])):
                    vcf_output.write(vcf_record)
        '''

if __name__ == "__main__":
    run()
