import os, sys, csv, argparse, random, pysam
import numpy as np
import pandas as pd
from collections import defaultdict

def sampler_parser():
    '''Sampler Argument Parser - Assigns arguments from command line.'''
    
    def sampler_argument_file ():
        '''Custom action to confrim the file exists prior to assignment.'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError # File not found
                setattr(args, self.dest, value)
        return customAction
    
    sampler_parser = argparse.ArgumentParser()
    
    # Input arguments
    sampler_parser.add_argument('--input', help='Defines the vcftools output file to be processed', required = True, type = str, action = sampler_argument_file())
    
    # Sampling methods. Currently mutually exclusive to only allow a single sampling method
    vcftools_sample_type = sampler_parser.add_mutually_exclusive_group()
    vcftools_sample_type.add_argument('--fst', dest = 'sample_type', help = "Defines Fst as the sample", action = 'store_const', const = 'MEAN_FST')
    vcftools_sample_type.add_argument('--TajimaD', dest = 'sample_type', help = "Defines TajimaD as the sample", action = 'store_const', const = 'TajimaD')
    vcftools_sample_type.add_argument('--user-specified', dest = 'sample_type', help = "Used to specify a specific column", type = str)
    
    # Sampling methods. Currently mutually exclusive to only allow a single sampling method
    vcftools_sampling = sampler_parser.add_mutually_exclusive_group()
    vcftools_sampling.add_argument('--uniform', dest = 'sampling', help="Samples the data uniformly", action = 'store_const', const = 'uniform')
    vcftools_sampling.add_argument('--random', dest = 'sampling', help="Samples the data randomly", action = 'store_const', const = 'random')
    
    # Sampling options
    sampler_parser.add_argument('--uniform-bins', help="Number of bins in uniform sampling", type = int, default = 10)
    sampler_parser.add_argument('--sample-size', help="Total sample size. If using --uniform, this number must be divisible by the bins", type = int, default = 200)
    sampler_parser.add_argument('--random-seed', help="Defines the random seed value for the random number generator", type = int, default = random.randint(1, 1000000000))
    
    # VCF file arguments
    vcf_file = sampler_parser.add_mutually_exclusive_group()
    vcf_file.add_argument('--vcf', dest = 'vcf_file', help = 'Defines the VCF file to be processed', type = str, action = sampler_argument_file())
    vcf_file.add_argument('--gzvcf', dest = 'vcf_file', help = 'Defines the compressed (gzipped) VCF file to be processed', type = str, action = sampler_argument_file())
    vcf_file.add_argument('--bcf', dest = 'vcf_file', help = 'Defines the BCF file to be processed', type = str, action = sampler_argument_file())
    
    # VCF Options
    sampler_parser.add_argument('--window-size', help = "Specifies the window size used for vcftools. Only required for sampling: TajimaD", type = str)
       
    sampler_parser.set_defaults(sample_type = 'MEAN_FST', sampling = 'uniform')
    
    return sampler_parser.parse_args()

def random_vcftools_sampler (vcftools_samples, sample_size):
    '''Random Sampler. Returns a list of randomly selected elements from
    vcftools_samples. The length of this list is defined by sample_size.'''
    
    try:
        return random.sample(vcftools_samples, sample_size)
    except ValueError:
        return vcftools_samples
        
def uniform_vcftools_sampler (vcftools_samples, uniform_sampling_args):
    '''Uniform Sampler. Seperates the elements of vcftools_samples into
    equal sized bins (number of bins specified by .uniform_bins arg). The
    samples in the bins are then randomly selected (number of selected
    samples is defined by: .sample_size / .uniform_bins).'''
        
    binned_sample_database = defaultdict(list)
       
    binned_samples = pd.cut(vcftools_samples, uniform_sampling_args.uniform_bins, labels = False)
    
    for vcftools_sample, sample_bin in zip(range(len(vcftools_samples)), binned_samples):
        binned_sample_database[sample_bin].append(vcftools_sample)
    
    uniform_samples = []
    
    for current_samples in binned_sample_database.values():
        uniform_samples.extend(random_vcftools_sampler(current_samples, uniform_sampling_args.sample_size / uniform_sampling_args.uniform_bins))
    
    return uniform_samples

def vcftools_sampler_logfile (sampler_args):
    '''Logger. Will be removed once logging is fully intergrated.'''
    
    sample_type_cvt = {'MEAN_FST':'fst', 'TajimaD':'TajimaD'}
       
    log_output = 'VCFtools Sampler Log\nParameters as interpreted:\n'
    for arg in ['input', 'random_seed', 'sample_size']:
        log_output += '\t--{0} {1}\n'.format(arg, getattr(sampler_args, arg))
    log_output += '\t--' + sample_type_cvt[getattr(sampler_args, 'sample_type')] + '\n'
    log_output += '\t--' + getattr(sampler_args, 'sampling')
    if sampler_args.sampling == 'uniform':
        log_output += '\n\t--uniform_bins ' + str(getattr(sampler_args, 'uniform_bins'))
    
    log_file = open(sampler_args.input + '.filtered.log', 'w')
    log_file.write(log_output)
    log_file.close()
    
def assign_sample_columns (sample_headers):
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
             
def run ():
    '''Sampler. Automates the sampling processs. Begins with checks to confirm the
    sampling is possible, then runs the specified sampling method, and produces the
    output files.'''
    
    # Get arguments from command line
    sampler_args = sampler_parser()
    # Set the random seed
    random.seed(sampler_args.random_seed)
    # Read in the sample file
    vcftools_samples = pd.read_csv(sampler_args.input, sep = '\t')
    
    # Confirm there are enough samples    
    if len(vcftools_samples) < sampler_args.sample_size:
        sys.exit('Sample size larger than the number of samples within input')
    
    elif len(vcftools_samples) == sampler_args.sample_size:
        # Update with warning call
        print 'Sample size equal to number of samples within input'
    
    # Run the uniform sampler       
    if sampler_args.sampling == 'uniform':
        if sampler_args.sample_size % sampler_args.uniform_bins != 0:
            sys.exit('Sample size not divisible by the bin count')
        selected_samples = sorted(uniform_vcftools_sampler(list(vcftools_samples[sampler_args.sample_type]), sampler_args))
    # Run the random sampler 
    if sampler_args.sampling == 'random':
        selected_samples = sorted(random_vcftools_sampler(range(len(vcftools_samples)), sampler_args.sample_size))
    
    # Reduce to only selected samples        
    vcftools_samples = vcftools_samples[vcftools_samples.index.isin(selected_samples)]
    
    # Create TSV file of the reduced samples
    vcftools_samples.to_csv(sampler_args.input + '.sampled', sep = '\t')
   
    # If a VCF file is specifed, create a reduced VCF file(s) from the samples
    if sampler_args.vcf_file:
        # Get the chromosome, start, and end columns
        chr_col, start_col, end_col = assign_sample_columns(list(vcftools_samples))
              
        vcf_input = pysam.VariantFile(sampler_args.vcf_file)
        for sampled_row in vcftools_samples.values:
            if end_col:
                vcf_output = pysam.VariantFile('Sampled_{0}_{1}_{2}.vcf.gz'.format(sampled_row[chr_col], sampled_row[start_col], sampled_row[end_col]), 'w', header = vcf_input.header)
                for vcf_record in vcf_input.fetch(sampled_row[chr_col], int(sampled_row[start_col]), int(sampled_row[end_col])):
                    vcf_output.write(vcf_record)
                   
    vcftools_sampler_logfile(sampler_args)
      
if __name__ == "__main__":
    run()
