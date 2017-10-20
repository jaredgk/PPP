import os
import sys
import subprocess
import argparse
import logging
import itertools

from collections import defaultdict

from vcftools import *

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
                    print '--%s only accepts integers' % self.dest
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

    #model_parser.add_argument('--ninds', help = 'Defines the number pf individuals within a population.', type = str, nargs='+', action = 'append')
    #model_parser.add_argument('--inds-names', help = 'Defines the individuals names within a population. ', type = str, nargs='+', action = 'append')

    # Statistic based arguments.
    #statistic_list = ['weir-fst', 'windowed-weir-fst', 'TajimaD', 'pi', 'freq', 'het']
    #statistic_default = 'windowed-weir-fst'

    #model_parser.add_argument('--calc-statistic', metavar = metavar_list(statistic_list), help = 'Specifies the statistic to calculate', type = str, choices = statistic_list, default = statistic_default)

    if passed_arguments:
        return model_parser.parse_args(passed_arguments)
    else:
        return model_parser.parse_args()

def logArgs(args, pipeline_function):
    '''Logs arguments from argparse system. May be replaced with a general function in logging_module'''
    logging.info('Arguments for %s:' % pipeline_function)
    for k in vars(args):
        logging.info('Arguments %s: %s' % (k, vars(args)[k]))

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

    # Loop each model specified to confirm parameters are valid
    for current_model in model_args.model:

        # Check that a tree has been assigned to the model
        if current_model not in model_args.model_tree:
            print 'No tree assigned to: %s' % current_model
            sys.exit()

        # Check that pops have been assigned to the model
        if current_model not in model_args.pops:
            print 'No populations assigned to: %s' % current_model
            sys.exit()

        # Check if the number of populations has been assigned
        if model_args.assign_npop:
            # Check if the number of populations has been assigned for the current model
            if model_args.assign_npop[current_model]:
                if len(model_args.pops[current_model]) != model_args.assign_npop[current_model]:
                    print 'Number of populations (--npop) assigned to %s does not match the named popultions' % current_model
                    sys.exit()

        # Loop each population in the model
        for current_pop in model_args.pops[current_model]:

            # Check that the population has been assigned to at least one model
            if current_pop not in set(itertools.chain.from_iterable(model_args.pops.values())):
                print 'Population (%s) not assigned to any model' % current_pop
                sys.exit()

            # Check that inds have been assigned to the population
            if current_pop not in model_args.inds:
                print 'No individuals assigned to: %s' % current_pop
                sys.exit()

            # Check that inds have been assigned to the population
            if current_pop not in model_args.inds:
                print 'No individuals assigned to: %s' % current_pop
                sys.exit()

            # Check if the number of inds has been assigned
            if model_args.assign_nind:
                # Check if the number of inds has been assigned for the current pop
                if model_args.assign_nind[current_pop]:
                    if len(model_args.inds[current_pop]) != model_args.assign_nind[current_pop]:
                        print 'Number of individuals (--nind) assigned to %s does not match the named individuals' % current_pop
                        sys.exit()

    model_output = open(model_args.out, 'w')
    # Loop each model specified to generate the output
    for current_model in model_args.model:
        model_output.write('MODEL:\t%s\n' % current_model)
        model_output.write('  TREE:\t%s\n' % model_args.model_tree[current_model])
        model_output.write('  NPOP:\t%s\n' % len(model_args.pops[current_model]))

        # Loop each population in the model
        for current_pop in model_args.pops[current_model]:
            model_output.write('  POP:\t%s\n' % current_pop)
            model_output.write('    NIND:\t%s\n' % len(model_args.inds[current_pop]))
            model_output.write('    IND:\t%s\n' % ', '.join(model_args.inds[current_pop]))

    model_output.close()

if __name__ == "__main__":
    run()
