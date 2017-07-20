import os, sys, csv, argparse, logging, random, pysam
import numpy as np
import pandas as pd
from collections import defaultdict

# Load path from Jared's directory, remove when not needed
sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'jared')))

from parse_functions import defaultsDictForFunction, getConfigFilename, makeRequiredList, getArgsWithConfig
from logging_module import initLogger

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

def logArgs(args, pipeline_function):
    logging.info('Arguments for %s:' % pipeline_function)
    for k in vars(args):
        logging.info('Arguments %s: %s' % (k, vars(args)[k]))

def run ():

    # Get arguments from command line
    sampler_args = sampler_parser()

    #required_args = ['vcfname','statistic_file','genename']
    #args = getArgsWithConfig(parser,sys_args,required_args,'vcf_ref_to_seq')

    # Set the random seed
    random.seed(sampler_args.random_seed)

    logArgs(sampler_args, 'Test')

    print sampler_args

if __name__ == "__main__":
    initLogger()
    run()
