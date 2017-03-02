import os, sys, csv, argparse, random
import numpy as np
import pandas as pd
from collections import defaultdict

def vcftools_output_parser():
    vcftools_parser = argparse.ArgumentParser()
    
    # Input arguments. 
    vcftools_parser.add_argument('--input', help='Defines the vcftools output file to be processed', required = True)
    
    # Sampling methods. Currently mutually exclusive to only allow a single sampling method
    vcftools_sample_type = vcftools_parser.add_mutually_exclusive_group()
    vcftools_sample_type.add_argument('--fst', dest = 'sample_type', help="Defines Fst as the sample", action = 'store_const', const = 'MEAN_FST')
    vcftools_sample_type.add_argument('--TajimaD', dest = 'sample_type', help="Defines TajimaD as the sample", action = 'store_const', const = 'TajimaD')
    
    # Sampling methods. Currently mutually exclusive to only allow a single sampling method
    vcftools_sampling = vcftools_parser.add_mutually_exclusive_group()
    vcftools_sampling.add_argument('--uniform', dest = 'sampling', help="Samples the data uniformly", action = 'store_const', const = 'uniform')
    vcftools_sampling.add_argument('--random', dest = 'sampling', help="Samples the data randomly", action = 'store_const', const = 'random')
    
    vcftools_parser.add_argument('--uniform-bins', help="Number of bins in uniform sampling", type = int, default = 10)
    vcftools_parser.add_argument('--sample-size', help="Total sample size. If using --uniform, this number must be divisible by the bins", type = int, default = 200)
    vcftools_parser.add_argument('--random-seed', help="Defines the random seed value for the random number generator", type = int, default = random.randint(1, 1000000000))
    
    vcftools_parser.set_defaults(sample_type = 'MEAN_FST', sampling = 'uniform')
    
    
    return vcftools_parser.parse_args()

def random_vcftools_sampler (vcftools_samples, sample_size):
    try:
        return random.sample(vcftools_samples, sample_size)
    except ValueError:
        return vcftools_samples
        
def uniform_vcftools_sampler (vcftools_samples, uniform_sampling_args):
        
    binned_sample_database = defaultdict(list)
       
    binned_samples = pd.cut(vcftools_samples, uniform_sampling_args.uniform_bins, labels = False)

    for vcftools_sample, sample_bin in zip(range(len(vcftools_samples)), binned_samples):
        binned_sample_database[sample_bin].append(vcftools_sample)
    
    uniform_samples = []
    
    for current_samples in binned_sample_database.values():
        uniform_samples.extend(random_vcftools_sampler(current_samples, uniform_sampling_args.sample_size / uniform_sampling_args.uniform_bins))
    
    return uniform_samples

def vcftools_sampler_logfile (sampler_args):
    sample_type_cvt = {'MEAN_FST':'fst', 'TajimaD':'TajimaD'}
       
    log_output = 'VCFtools Sampler Log\nParameters as interpreted:\n'
    for arg in ['input', 'random_seed', 'sample_size']:
        log_output += '\t--{0} {1}\n'.format(arg, getattr(sampler_args, arg))
    log_output += '\t--' + sample_type_cvt[getattr(sampler_args, 'sample_type')] + '\n'
    log_output += '\t--' + getattr(sampler_args, 'sampling')
    if sampler_args.sampling == 'uniform':
        log_output += '\n\t--uniform_bins' + str(getattr(sampler_args, 'uniform_bins'))
    
    log_file = open(sampler_args.input + '.filtered.log', 'w')
    log_file.write(log_output)
    log_file.close()
          
def return_parsed_output ():
    
    vcftools_args = vcftools_output_parser()
    
    random.seed(vcftools_args.random_seed)
        
    with open(vcftools_args.input, 'rb') as vcftools_file:
        vcftools_file_reader = csv.DictReader(vcftools_file, delimiter='\t')
        vcftools_file_fieldnames = vcftools_file_reader.fieldnames
        vcftools_samples = [float(vcftools_output[vcftools_args.sample_type]) for vcftools_output in vcftools_file_reader]
    
    with open(vcftools_args.input + '.filtered', 'w') as output_file:
        output_file_writer = csv.DictWriter(output_file, fieldnames=vcftools_file_fieldnames)
        output_file_writer.writeheader()
    
        if vcftools_args.sampling == 'uniform':
            if vcftools_args.sample_size % vcftools_args.uniform_bins == 0:
                selected_samples = uniform_vcftools_sampler(vcftools_samples, vcftools_args)
            else:
                print 'Error. Sample size not divisible by the bins'
                        
        if vcftools_args.sampling == 'random':
            selected_samples = random_vcftools_sampler(range(len(vcftools_samples)), vcftools_args.sample_size)
            
        with open(vcftools_args.input, 'rb') as vcftools_output_file:
            for current_line, current_data in enumerate(csv.DictReader(vcftools_output_file, delimiter='\t')):
                if current_line in selected_samples:
                    output_file_writer.writerow(current_data)
                    selected_samples.remove(current_line)
    
    vcftools_sampler_logfile(vcftools_args)
           
if __name__ == "__main__":
    return_parsed_output()
