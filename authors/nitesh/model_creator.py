import os
import sys
import json
import argparse
import logging
import itertools

from collections import defaultdict
from collections import OrderedDict

# Insert Jared's directory path, required for calling Jared's functions. Change when directory structure changes.
sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'jared')))

from logging_module import initLogger

def model_creator_parser (passed_arguments):
    '''VCF Argument Parser - Assigns arguments from command line'''

    def parser_dict_list_append ():
        '''Custom action to house data as defaultdict list'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if getattr(args, self.dest):
                    # Append the argument with the file
                    getattr(args, self.dest)[value[0]].append(value[1])

                else:
                    # Set the argument with the file (as list)
                    arg_dict = defaultdict(list)
                    arg_dict[value[0]].append(value[1])
                    setattr(args, self.dest, arg_dict)
        return customAction

    def parser_dict_list_extend ():
        '''Custom action to house data as defaultdict list'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if getattr(args, self.dest):
                    # Append the argument with the file
                    getattr(args, self.dest)[value[0]].extend(value[1:])

                else:
                    # Set the argument with the file (as list)
                    arg_dict = defaultdict(list)
                    arg_dict[value[0]].extend(value[1:])
                    setattr(args, self.dest, arg_dict)
        return customAction

    def parser_dict_str ():
        '''Custom action to house data as defaultdict int'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                # Assign the passed value
                if getattr(args, self.dest):
                    # Append the argument with the file
                    getattr(args, self.dest)[value[0]] = value[1]

                else:
                    # Set the argument with the file (as list)
                    arg_dict = defaultdict(str)
                    arg_dict[value[0]] = value[1]
                    setattr(args, self.dest, arg_dict)
        return customAction

    def parser_dict_int ():
        '''Custom action to house data as defaultdict int'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):

                # Check the passed value is an int
                try:
                    int(value[1])
                except:
                    print('--%s only accepts integers' % self.dest)
                    sys.exit()

                # Assign the passed value
                if getattr(args, self.dest):
                    # Append the argument with the file
                    getattr(args, self.dest)[value[0]] = int(value[1])

                else:
                    # Set the argument with the file (as list)
                    arg_dict = defaultdict(int)
                    arg_dict[value[0]] = int(value[1])
                    setattr(args, self.dest, arg_dict)
        return customAction

    def parser_dict_file ():
        '''Custom action to house data as defaultdict int'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):

                if not os.path.isfile(value[1]):
                    raise IOError('Input not found.') # File not found

                # Assign the passed value
                if getattr(args, self.dest):
                    # Append the argument with the file
                    getattr(args, self.dest)[value[0]] = value[1]

                else:
                    # Set the argument with the file (as list)
                    arg_dict = defaultdict(str)
                    arg_dict[value[0]] = value[1]
                    setattr(args, self.dest, arg_dict)
        return customAction

    def metavar_list (var_list):
        '''Create a formmated metavar list for the help output'''
        return '{' + ', '.join(var_list) + '}'

    model_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    # Input arguments.
    #model_parser.add_argument("--vcf", metavar = 'VCF_Input', help = "Input VCF filename", type = str, action = parser_confirm_file())

    # Other file arguments. Expand as needed
    model_parser.add_argument('--out', help = 'Specifies the complete output filename.', type = str, default = 'out.model')

    # General arguments.
    model_parser.add_argument('--overwrite', help = "Specifies if previous output files should be overwritten", action = 'store_true')

    # General model arguments
    model_parser.add_argument('--model', help = 'Defines the name of a model', type = str, action = 'append', required = True)

    model_parser.add_argument('--model-tree', help = 'Defines the population tree of the model, in Newick format', type = str, nargs = 2, action = parser_dict_str(), required = True)

    # Population assignment
    model_parser.add_argument('--assign-npop', help = 'Assigns the number of populations to a model.', type = str, nargs = 2, action = parser_dict_int())

    # Pop assignment methods
    pop_group = model_parser.add_mutually_exclusive_group(required=True)
    pop_group.add_argument('--assign-pop', dest = 'pops', help = 'Assigns a population name to a model.', type = str, nargs = 2, action = parser_dict_list_append())
    pop_group.add_argument('--assign-pops', dest = 'pops', help = 'Assigns multiple population names to a model.', type = str, nargs = '+', action = parser_dict_list_extend())

    # Individual assignment
    model_parser.add_argument('--assign-nind', help = 'Assigns number of individuals to a model.', type = str, nargs = 2, action = parser_dict_int())

    # Ind Assignment
    ind_group = model_parser.add_mutually_exclusive_group(required=True)
    ind_group.add_argument('--assign-ind', dest = 'inds', help = 'Assigns an individual name to a population.', type = str, nargs = 2, action = parser_dict_list_append())
    ind_group.add_argument('--assign-inds', dest = 'inds', help = 'Assigns multiple individual names to a population.', type = str, nargs = '+', action = parser_dict_list_extend())
    ind_group.add_argument('--assign-ind-file', help = 'Assigns multiple individual names to a population using a file.', type = str, nargs = 2, action = parser_dict_file())

    if passed_arguments:
        return model_parser.parse_args(passed_arguments)
    else:
        return model_parser.parse_args()

def logArgs(args, pipeline_function):
    '''Logs arguments from argparse system. May be replaced with a general function in logging_module'''
    logging.info('Arguments for %s:' % pipeline_function)
    for k in vars(args):
        logging.info('Arguments %s: %s' % (k, vars(args)[k]))

def read_population_file (filename):
    # List of store individuals
    inds = []

    # Read the individuals from the file
    with open(filename, 'rU') as pop_file:
        for pop_line in pop_file:
            # Store the individuals
            inds.append(pop_line.strip())

    # Return individuals
    return inds

def run (passed_arguments = []):
    '''

    Parameters
    ----------

    Returns
    -------

    Raises
    ------

    '''

    # Grab VCF arguments from command line
    model_args = model_creator_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(model_args, 'model_creator')

    # Dict to store lists of individuals
    ind_dict = defaultdict(list)

    # Check if ind_files were specified by the user
    if model_args.assign_ind_file:
        for pop, ind_file in model_args.assign_ind_file.items():
            ind_dict[pop] = read_population_file(ind_file)

    else:
        ind_dict = model_args.inds

    # Loop each model specified to confirm parameters are valid
    for current_model in model_args.model:

        # Check that a tree has been assigned to the model
        if current_model not in model_args.model_tree:
            print('No tree assigned to: %s' % current_model)
            sys.exit()

        # Check that pops have been assigned to the model
        if current_model not in model_args.pops:
            print('No populations assigned to: %s' % current_model)
            sys.exit()

        # Check if the number of populations has been assigned
        if model_args.assign_npop:
            # Check if the number of populations has been assigned for the current model
            if model_args.assign_npop[current_model]:
                if len(model_args.pops[current_model]) != model_args.assign_npop[current_model]:
                    print( len(model_args.pops[current_model]),  model_args.assign_npop[current_model])
                    print('Number of populations (--npop) assigned to %s does not match the named popultions' % current_model)
                    sys.exit()

        # Loop each population in the model
        for current_pop in model_args.pops[current_model]:

            # Check that the population has been assigned to at least one model
            if current_pop not in set(itertools.chain.from_iterable(model_args.pops.values())):
                print('Population (%s) not assigned to any model' % current_pop)
                sys.exit()

            # Check that inds have been assigned to the population
            if current_pop not in ind_dict:
                print('No individuals assigned to: %s' % current_pop)
                sys.exit()

                # Check if the number of inds has been assigned
                if model_args.assign_nind:
                    # Check if the number of inds has been assigned for the current pop
                    if model_args.assign_nind[current_pop]:
                        if len(ind_dict[current_pop]) != model_args.assign_nind[current_pop]:
                            print('Number of individuals (--nind) assigned to %s does not match the named individuals' % current_pop)
                            sys.exit()

    # List to store the models
    models_dict = []

    # Loop each model specified to generate the output
    for current_model in model_args.model:

        # Create an ordered container for model information
        model_dict = OrderedDict()

        # Save the name of the model
        model_dict['name'] = current_model

        # Check if the model isn't a single population
        if len(model_args.pops[current_model]) > 1:

            # Assign the tree to the population
            model_dict['tree'] = model_args.model_tree[current_model]

        # Create container for population information
        pop_dict = defaultdict(lambda: defaultdict(list))

        # Loop each population in the model
        for current_pop in model_args.pops[current_model]:
            # Assign the individuals to the population
            pop_dict[current_pop]['inds'] = ind_dict[current_pop]

        # Append the pop_dict
        model_dict['pops'] = pop_dict

        # Append the model_dict
        models_dict.append(model_dict)

    # Open the output file
    output_file = open(model_args.out, 'w')

    # Write the json-formmated data to the output file
    output_file.write(json.dumps(models_dict, indent = 4))

    output_file.close()


if __name__ == "__main__":
    run()
