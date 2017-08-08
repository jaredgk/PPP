import os, sys, csv, argparse, random, pysam, itertools, logging
import numpy as np
import pandas as pd
from collections import defaultdict

# Insert Jared's directory path, required for calling Jared's functions. Change when directory structure changes.
sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'jared')))

from logging_module import initLogger

def sampler_parser(passed_arguments):
    '''Sampler Argument Parser - Assigns arguments from command line.'''

    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string = None):
                if not os.path.isfile(value):
                    raise IOError('Input not found.') # File not found
                setattr(args, self.dest, value)
        return customAction

    def parser_confirm_no_file (file_details):
        '''Custom action to confirm file does not exist'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string = None):
                if os.path.isfile(value):
                    raise IOError(file_details + ' already exists.') # File found
                setattr(args, self.dest, value)
        return customAction

    def metavar_list (var_list):
        '''Create a formmated metavar list for the help output'''
        return '{' + ', '.join(var_list) + '}'

    sampler_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    # Input arguments
    sampler_parser.add_argument("vcfname", metavar='VCF_Input', help = "Input VCF filename", type = str, action = parser_confirm_file())
    sampler_parser.add_argument('--statistic-file', help='Specifies the statistic file for filtering', required = True, type = str, action = parser_confirm_file())

    # Output arguents
    sampler_parser.add_argument('--out', help = 'Specifies the VCF output filename', type = str, default = 'out.vcf.gz',action = parser_confirm_no_file('Output'))
    sampler_parser.add_argument('--sample-file', help = 'Specifies the sampled (statistic file) tsv output filename', type = str, default = 'sampled_data.tsv', action = parser_confirm_no_file('Sample file'))

    # Statistic based arguments.
    statistic_list = ['windowed-weir-fst', 'TajimaD']
    statistic_default = 'windowed-weir-fst'
    sampler_parser.add_argument('--calc-statistic', metavar = metavar_list(statistic_list), help = 'Specifies the statistic calculated ', type=str, choices = statistic_list, default = statistic_default)

    sampler_parser.add_argument('--statistic-window-size', help = 'Specifies the size of window calculations', type = int, default = 10000)

    # Sampling methods. Currently mutually exclusive to only allow a single sampling method
    sampling_list = ['uniform', 'random']
    sampling_default = 'random'
    sampler_parser.add_argument('--sampling-scheme', metavar = metavar_list(sampling_list), help = 'Specifies the sampling scheme ', type=str, choices = sampling_list, default = sampling_default)

    # Sampling options
    sampler_parser.add_argument('--uniform-bins', help="Number of bins in uniform sampling", type = int, default = 10)
    sampler_parser.add_argument('--sample-size', help="Total sample size. If using the uniform sampling scheme, this number must be divisible by the bins", type = int, default = 200)

    # Other options
    sampler_parser.add_argument('--random-seed', help="Defines the random seed value for the random number generator", type = int, default = random.randint(1, 1000000000))

    if passed_arguments:
        return sampler_parser.parse_args(passed_arguments)
    else:
        return sampler_parser.parse_args()

def logArgs(args, pipeline_function):
    '''Logs arguments from argparse system. May be replaced with a general function in logging_module'''
    logging.info('Arguments for %s:' % pipeline_function)
    for k in vars(args):
        logging.info('Arguments %s: %s' % (k, vars(args)[k]))

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
    enough columns were found.'''

    if not all(pos_header in sample_headers for pos_header in ['CHROM', 'BIN_START', 'BIN_END']):
        if all(pos_header in sample_headers for pos_header in ['CHROM', 'BIN_START']):
            return sample_headers.index('CHROM'), sample_headers.index('BIN_START'), None
        else:
            if 'CHROM' in sample_headers:
                logging.error('Cannot find BIN_START column in file specified by --statistic-file.')
                raise ValueError('Cannot find BIN_START column in file specified by --statistic-file.')
            elif 'BIN_START' in sample_headers:
                logging.error('Cannot find CHROM column in file specified by --statistic-file.')
                raise ValueError('Cannot find CHROM column in file specified by --statistic-file.')
            else:
                logging.error('Cannot find CHROM and BIN_START columns in file specified by --statistic-file.')
                raise ValueError('Cannot find CHROM and BIN_START columns in file specified by --statistic-file.')

    return sample_headers.index('CHROM'), sample_headers.index('BIN_START'), sample_headers.index('BIN_END')

def assign_statistic_column (sample_headers, statistic):
    statistic_converter = {'windowed-weir-fst':'MEAN_FST', 'TajimaD':'TajimaD'}

    if not statistic_converter.has_key(statistic):
        logging.critical('Statistic not found. Statistic list needs to be updated. Please contact the PPP Team.')
        raise Exception('Statistic not found. Statistic list needs to be updated. Please contact the PPP Team.')

    if statistic_converter[statistic] not in sample_headers:
        logging.error('Statistic selected not found in file specified by --statistic-file.')
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
        VCF_Input : str
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
    logArgs(sampler_args, 'vcf_sampler')

    # Set the random seed
    np.random.seed(sampler_args.random_seed)

    # Read in the sample file
    vcftools_samples = pd.read_csv(sampler_args.statistic_file, sep = '\t')

    logging.info('Sample (i.e. statistic) file assigned')

    # Confirm there are enough samples
    if len(vcftools_samples) < sampler_args.sample_size:
        logging.error('Sample size larger than the number of samples within sample file')
        raise ValueError('Sample size larger than the number of samples within sample file')

    # UPDATE Planned
    # Improve from equal to too few samples (i.e. samples = 100, sample_size = 95)
    elif len(vcftools_samples) == sampler_args.sample_size:
        # Warns the user of poor sampling
        logging.warning('Sample size equal to number of samples within sample file')

    # Run the uniform sampler
    if sampler_args.sampling_scheme == 'uniform':
        if sampler_args.sample_size % sampler_args.uniform_bins != 0:
            logging.error('Sample size not divisible by the bin count')
            raise ValueError('Sample size not divisible by the bin count')

        assigned_statistic = assign_statistic_column (list(vcftools_samples), sampler_args.calc_statistic)
        selected_samples = sorted(uniform_vcftools_sampler(list(vcftools_samples[assigned_statistic]), sampler_args.uniform_bins, sampler_args.sample_size))

        logging.info('Uniform sampling complete')

    # Run the random sampler
    if sampler_args.sampling_scheme == 'random':

        selected_samples = sorted(random_vcftools_sampler(range(len(vcftools_samples)), sampler_args.sample_size))

        logging.info('Random sampling complete')

    # Reduce to only selected samples
    sampled_samples = vcftools_samples[vcftools_samples.index.isin(selected_samples)]

    # Create selected samples TSV file, with either the default filename or a user-defined filename
    sampled_samples.to_csv(sampler_args.sample_file, sep = '\t')

    logging.info('Created selected samples file')

    # If a VCF file is specifed, create a reduced VCF file(s) from the samples
    if sampler_args.vcfname:

        split_vcfname = os.path.basename(sampler_args.vcfname).split(os.extsep, 1)

        # Open the VCF file
        vcf_input = pysam.VariantFile(sampler_args.vcfname)

        # Create the VCF output file, with either the default filename or a user-defined filename
        vcf_output = pysam.VariantFile(sampler_args.out, 'w', header = vcf_input.header)

        # Get the chromosome, start, and end columns
        chr_col, start_col, end_col = assign_position_columns(list(vcftools_samples))
        for sampled_row in vcftools_samples.values:
            if not end_col:
                if sampler_args.calc_statistic == 'TajimaD':

                    if not sampler_args.statistic_window_size:
                        logging.error("--statistic-window-size argument required for the Tajima's D statistic")
                        raise TypeError("--statistic-window-size argument required for the Tajima's D statistic")

                    # Create iterator of all unique bin combinations
                    bin_start_iter = itertools.combinations(list(vcftools_samples['BIN_START']), 2)

                    # Check if the bin combinations are divisible by the window size
                    if not all([abs(bin_1 - bin_2) % sampler_args.statistic_window_size == 0  for bin_1, bin_2 in bin_start_iter]):
                        logging.error("--statistic-window-size argument conflicts with values in sample file")
                        raise ValueError("--statistic-window-size argument conflicts with values in sample file")

                    for vcf_record in vcf_input.fetch(sampled_row[chr_col], int(sampled_row[start_col]), int(sampled_row[start_col]) + (sampler_args.statistic_window_size - 1)):
                        vcf_output.write(vcf_record)
            else:
                for vcf_record in vcf_input.fetch(sampled_row[chr_col], int(sampled_row[start_col]), int(sampled_row[end_col])):
                    vcf_output.write(vcf_record)

        vcf_input.close()
        vcf_output.close()

        logging.info('Created selected VCF file')


if __name__ == "__main__":
    initLogger()
    run()
