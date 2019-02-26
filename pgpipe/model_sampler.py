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

# Insert Jared's directory path, required for calling Jared's functions. Change when directory structure changes.
#sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'pppipe')))

# Import log initializer
from pgpipe.logging_module import initLogger, logArgs

# Model file related functions
from pgpipe.model import read_model_file, write_model_file

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

    def parser_dict_list ():
        '''Custom action to house data as defaultdict int'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, values, option_string=None):

                argument_call = self.dest.replace('_','-')

                if not isinstance(values, (list, tuple)):
                    raise Exception('--%s argument assignment error. Expected iterable' % argument_call)

                if len(values) != 2:
                    raise Exception('--%s argument assignment error. Expected iterable with a length 2' % argument_call)

                # Assign the passed value
                if getattr(args, self.dest):
                    # Append the argument with the file
                    getattr(args, self.dest)[values[0]].append(values[1])

                else:
                    # Set the argument with the file (as list)
                    arg_dict = defaultdict(list)
                    arg_dict[values[0]].append(values[1])
                    setattr(args, self.dest, arg_dict)
        return customAction

    def parser_dict_int ():
        '''Custom action to house data as defaultdict int'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, values, option_string=None):

                argument_call = self.dest.replace('_','-')

                if not isinstance(values, (list, tuple)):
                    raise Exception('--%s argument assignment error. Expected iterable' % argument_call)

                if len(values) != 2:
                    raise Exception('--%s argument assignment error. Expected iterable with a length 2' % argument_call)

                try:
                    int(values[1])
                except:
                    raise Exception('--%s argument assignment error. Expected integer as second value' % argument_call)

                # Assign the passed value
                if getattr(args, self.dest):

                    if values[0] in getattr(args, self.dest):
                        raise Exception('--%s argument assignment error. Models may only be sampled once' % argument_call)

                    # Append the argument with the file
                    getattr(args, self.dest)[values[0]] = int(values[1])

                else:
                    # Set the argument with the file (as list)
                    arg_dict = defaultdict(int)
                    arg_dict[values[0]] = int(values[1])
                    setattr(args, self.dest, arg_dict)
        return customAction

    def parser_dict_dict ():
        '''Custom action to house data as defaultdict dict'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, values, option_string=None, split_char = '='):

                argument_call = self.dest.replace('_','-')

                if not isinstance(values, (list, tuple)):
                    raise Exception('--%s argument assignment error. Expected iterable' % argument_call)

                if len(values) < 2:
                    raise Exception('--%s argument assignment error. Expected iterable with length >= 2' % argument_call)

                # Create dict to hold the non-key values
                value_dict = {}

                # Loop the non-key values
                for sample_value in values[1:]:

                    # Check that the sample_value may be split using split_char
                    if split_char not in sample_value:
                        # Report an error if the data cannot be split
                        raise Exception('--%s argument assignment error. Populations could not be split' % argument_call)

                    # Split the value using the split_char argument
                    split_value = sample_value.split(split_char)
                    # Add the split value to value_dict
                    value_dict[split_value[0]] = split_value[1]

                # Assign the passed value
                if getattr(args, self.dest):

                    if values[0] in getattr(args, self.dest):
                        raise Exception('--%s argument assignment error. Models may only be sampled once' % argument_call)

                    # Append the argument with the file
                    getattr(args, self.dest)[values[0]] = value_dict

                else:
                    # Set the argument with the file (as list)
                    arg_dict = defaultdict(dict)
                    arg_dict[values[0]] = value_dict
                    setattr(args, self.dest, arg_dict)
        return customAction

    def metavar_list (var_list):
        '''Create a formmated metavar list for the help output'''
        return '{' + ', '.join(var_list) + '}'

    sampler_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    # Input arguments
    sampler_parser.add_argument('--model-file', help = 'Model file to sample from', required = True, type = str, action = parser_confirm_file())

    # Sampler arguments
    sampler_parser.add_argument('--model-copy', help = 'Copy the specified model', type = str, nargs = 2, action = parser_dict_list())

    sampler_parser.add_argument('--model-sample-pops', help = 'Sample all populations within the selected model', type = str, nargs = 2, action = parser_dict_int())

    sampler_parser.add_argument('--model-sample-pop', help = 'Sample all populations within the selected model', type = str, nargs = '+', action = parser_dict_dict())

    # New-file options
    sampler_parser.add_argument('--out', help = 'Specifies the filename of the new model file', type = str, default = 'out.model')
    sampler_parser.add_argument('--keep-unsampled', help = 'Specifies if an output file should keep unsampled models', action = 'store_true')

    # Random options
    sampler_parser.add_argument('--random-seed', help="Defines the random seed value for the random number generator", type = int, default = np.random.randint(1, high = 1000000000))

    if passed_arguments:
        return sampler_parser.parse_args(passed_arguments)
    else:
        return sampler_parser.parse_args()

def sample_model_w_dict (model_file_data, model_to_sample, sample_pop_dict):

    # Check if the selected model was not found in the file
    if model_to_sample not in model_file_data:
        raise IOError('Selected model "%s" not found in %s' % (model_to_sample, model_file_data))

    # Loop the population sampleing dict
    for sample_pop_name, sample_pop_int in sample_pop_dict.items():

        # Sample the selected pop within the specified model
        model_file_data[model_to_sample].sample_pop(sample_pop_name, sample_pop_int)

def sample_model_w_int (model_file_data, model_to_sample, sample_pops_int):

    # Check if the selected model was not found in the file
    if model_to_sample not in model_file_data:
        raise IOError('Selected model "%s" not found in %s' % (model_to_sample, model_file_data))

    # Sample the selected model using the sample pops int
    model_file_data[model_to_sample].sample_pops(sample_pops_int)

def copy_model (model_file_data, model_to_copy, model_from_copy):

    # Check if the source model was not found in the file
    if model_to_copy not in model_file_data:
        raise IOError('Source model for copy "%s" not found in %s' % (model_to_copy, model_file_data))

    # Check if the copy is found in the file
    if model_from_copy in model_file_data:
        raise IOError('Copy "%s" already in %s' % (model_from_copy, model_file_data))

    # Copy the model
    model_file_data.copy_model(model_to_copy, model_from_copy)

def run (passed_arguments = []):
    '''
        Model sampler.

        Automates the sampling process of a specified model file. The function
        allows the user to sample populations to contain a specific number of
        individuals.

        Parameters
        ----------
        --model-file : str
            Model file to sample from
        --model-copy : str str
            Copy the specified model (i.e. source) as a new model (i.e. target).
            Example: --model-copy source_model_name target_model_name
        --model-sample-pops : str int
            Sample the population within the specified model using the given
            integer. Example: --model-copy model_name k
        --model-sample-pop : str str=int
            Sample the specified population within the specified model using the
            given integer. Example: --model-copy model_name pop=k
        --random-seed : int
            Specifies the random seed value for the random number generator
        --out : str
            Filename of the sampled model output
        --keep-unsampled :
            Specifies if unsampled models should be saved in the output

        Returns
        -------
        output : file
            Sampled statistic file
        log : file
            Log file output

        Raises
        ------
        IOError
            Model file does not exist
        Exception
            If user attepts to sample a model more than once
    '''

    # Get arguments from command line
    sampler_args = sampler_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(sampler_args, func_name = 'model_sampler')

    # Check that no model is found in both --model-sample-pop/--model-sample-pops
    if set(list(sampler_args.model_sample_pop)) & set(list(sampler_args.model_sample_pops)):
        raise Exception('Argument error. Sampling processes found in both --model-sample-pop and --model-sample-pops')

    # Set the random seed
    np.random.seed(sampler_args.random_seed)

    logging.info('Random seed: %s' % sampler_args.random_seed)

    # Read in the model file
    model_file_data = read_model_file(sampler_args.model_file)

    # Check if any models will be copied, should copy prior to sampling
    if sampler_args.model_copy:

        # Loop the copy entries
        for model_to_copy, models_from_copy in sampler_args.model_copy.items():

            # Loop the models created from the copy
            for model_from_copy in models_from_copy:

                # Copy the model
                copy_model(model_file_data, model_to_copy, model_from_copy)

            logging.info('%s successfully copied to %s' % (model_to_copy, models_from_copy))

    # Check if any models will sampled using the pops method
    if sampler_args.model_sample_pops:

        # Loop the entries
        for model_to_sample, sample_pops_int in sampler_args.model_sample_pops.items():

            sample_model_w_int(model_file_data, model_to_sample, sample_pops_int)

            logging.info('%s successfully sampled by %s' % (model_to_sample, sample_pops_int))

    # Check if any models will sampled using the pop method
    if sampler_args.model_sample_pop:

        # Loop the entries
        for model_to_sample, sample_pops_dict in sampler_args.model_sample_pop.items():

            sample_model_w_dict(model_file_data, model_to_sample, sample_pops_dict)

        logging.info('%s successfully sampled by each population' % model_to_sample )

    # Check if the sampled models should be deleted
    if not sampler_args.keep_unsampled:

        # Create a complete list of the sampled models
        sampled_models = list(sampler_args.model_sample_pop) + list(sampler_args.model_sample_pops)

        # Create lite of unsampled models
        unsampled_models =  list(set(model_file_data) - set(sampled_models))

        # Loop the unsampled models
        for unsampled_model in unsampled_models:

            # Delete the unsampled model
            del model_file_data[unsampled_model]

        logging.info('Unsampled model removed')

    # Create the model output file
    write_model_file(model_file_data, sampler_args.out)

if __name__ == "__main__":
    initLogger()
    run()
