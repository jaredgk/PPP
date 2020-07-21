#!/usr/bin/env python
'''
    Many PPP functions were designed to automatically assign relevant populations 
    and/or individuals using a *Model* file. To enable this functionality, the 
    model_creator.py function may be used to produce *Model* files by either: i) 
    manually entering the necessary information or ii) by using files with the 
    relevant information. It is also all possible to create multiple models 
    simultaneously and assign populations to more than a single model. 

    .. image:: ../../PPP_assets/PPP_Model.png
        :width: 100 %
        :align: center

    A simple way to visualize models are as a hierarchy. Each *Model* may contain 
    one or more *Populations* and each *Population* may contain one or more 
    *Individuals*.

    ##################
    Command-line Usage
    ##################
    The model creator may be called using the following command:

    .. code-block:: bash
        
        model_creator.py

    ***********************
    Example 1: Simple Model
    ***********************
    A basic model only require a single population (pop) with a single in individual (ind).
    Only three commands are required:

    #. Create and name a model: *--model 1Pop*
    #. Assign a pop to a model: *--model-pop 1Pop Paniscus*
    #. Assign an ind to a pop:  *--pop-ind Paniscus Pan_paniscus-9731_LB502*

    .. code-block:: bash

        model_creator.py --model 1Pop --model-pop 1Pop Paniscus --pop-ind Paniscus Pan_paniscus-9731_LB502
        
    ****************************
    Example 2: Model Using Files
    ****************************
    A model may also be created using two file options:

    #. Assign multiple pops to model: *--model-pop-file 2Pop 2Pops.txt*
    #. Assign multiple inds to pop:   *--pop-ind-file Paniscus Paniscus.txt*
                                      
    .. code-block:: bash
        
        model_creator.py --model 2Pop --model-pop-file 2Pop examples/files/2Pops.txt --pop-ind-file Paniscus examples/files/Paniscus.txt --pop-ind-file Troglodytes examples/files/Troglodytes.txt

    .. code-block:: bash
       :caption: examples/files/2Pops.txt

       Paniscus
       Troglodytes

    .. code-block:: bash
       :caption: examples/files/Paniscus.txt

       Pan_paniscus-9731_LB502
       Pan_paniscus-A915_Kosana

    .. code-block:: bash
       :caption: examples/files/Troglodytes.txt

       Pan_troglodytes_troglodytes-A957_Vaillant
       Pan_troglodytes_troglodytes-A958_Doris
    
    *************************************
    Example 3: Update Model in Model File
    *************************************
    A model may be updated if desired using the following options:
    
    #. Assign the model file:                             *--model-file input.model*
    #. Assign the model to update:                        *--update-model 2Pop*
    #. Assign a pop to remove from a model:               *--model-rm-pop 2Pop Troglodytes*
    #. Assign a new or previously created pop to a model: *--model-pop 2Pop Schweinfurthii*
    #. Assign multiple inds to pop (if new):              *--pop-ind-file Schweinfurthii Schweinfurthii.txt*
                                      
    .. code-block:: bash
        
        model_creator.py --model-file examples/files/input.model --update-model 2Pop --model-pop 2Pop Schweinfurthii --model-rm-pop 2Pop Troglodytes

    ******************************************
    Example 4: Update Population in Model File
    ******************************************
    A population may be updated if desired using the following options:
    
    #. Assign the model file:             *--model-file input.model*
    #. Assign the pop to update:          *--update-pop Paniscus*
    #. Assign a ind to add to a pop:      *--pop-ind Paniscus Pan_paniscus-Unknown*
    #. Assign a ind to remove from a pop: *--pop-rm-ind Paniscus Pan_paniscus-9731_LB502*
                                      
    .. code-block:: bash
        
        model_creator.py --model-file examples/files/input.model --update-pop Paniscus --pop-ind Paniscus Pan_paniscus-Unknown --pop-rm-ind Paniscus Pan_paniscus-9731_LB502

    #####################################
    Standard Model Command-line Arguments
    #####################################
    Except for **--model-file** all other model-based arguments may be used 
    multiple times.
    
    **--model-file** *<str>*
        Argument used to define the name of a model file.
    **--model** *<model_str>*
        Argument used to define the name of a model to create.
    **--update-model** *<str>*
        Argument used to define the name of a model to update. Allows 
        for: i) tree update and ii) populations to be added.
    **--update-pop** *<str>*
        Argument used to define the name of a population to update. 
        Allows for individuals to be added.
    **--model-tree** *<model_str>* *<newick_str>*
        Argument used to assign a population tree to a model, in Newick format. 
    **--model-tree-file** *<model_str>* *<newick_file>*
        Argument used to assign a population tree file to a model, in Newick 
        format.
    **--model-pop** *<model_str>* *<pop_str>*
        Argument used to assign a population to a model.
    **--model-pops** *<model_str>* *<pop1_str>* *<pop2_str>* ..
        Argument used to assign a multiple populations to a model.
    **--model-pop-file** *<model_str>* *<pop_file>*
        Argument used to assign a multiple populations to a model using a file.
    **--pop-ind** *<pop_str>* *<ind_str>*
        Argument used to assign a individual to a population.
    **--pop-inds** *<pop_str>* *<ind1_str>* *<ind2_str>* ..
        Argument used to assign a multiple individuals to a population.
    **--pop-ind-file** *<pop_str>* *<ind_file>*
        Argument used to assign a multiple individuals to a population using a file.

    ###############################################
    Model Update: Compatible Command-line Arguments
    ###############################################
    Please note: **--update-model** is required to update a model.
    
    **--update-model** *<str>*
        Argument used to define the name of a model to update. Allows for: 
        i) tree update and ii) populations to be added and/or removed.
    **--model-file** *<str>*
        Argument used to define the name of a model file.
    **--model-tree** *<model_str>* *<newick_str>*
        Argument may be used to replace a population tree, in Newick format.
    **--model-tree-file** *<model_str>* *<newick_file>*
        Argument may be used to replace a population usng a tree file, in 
        Newick format.
    **--model-pop** *<model_str>* *<pop_str>*
        Argument used to assign a population to the model begin updated.
    **--model-pops** *<model_str>* *<pop1_str>* *<pop2_str>* ..
        Argument used to assign multiple populations to the model begin updated.
    **--model-pop-file** *<model_str>* *<pop_file>*
        Argument used to assign multiple populations to the model begin updated 
        using a file.
    **--model-rm-pop** *<model_str>* *<pop_str>*
        Argument used to remove a population to the model begin updated.
    **--model-rm-pops** *<model_str>* *<pop1_str>* *<pop2_str>* ..
        Argument used to remove multiple populations to the model begin updated.
    **--model-rm-pop-file** *<model_str>* *<pop_file>*
        Argument used to remove multiple populations to the model begin updated 
        using a file.


    ####################################################
    Population Update: Compatible Command-line Arguments
    ####################################################
    Please note: **--update-pop** is required to update a population.

    **--update-pop** *<str>*
        Argument used to define the name of a population to update. Allows for 
        individuals to be added.
    **--model-file** *<str>*
        Argument used to define the name of a model file.
    **--pop-ind** *<pop_str>* *<ind_str>*
        Argument used to assign a individual to the population begin updated.
    **--pop-inds** *<pop_str>* *<ind1_str>* *<ind2_str>* ..
        Argument used to assign multiple individuals to the population begin 
        updated.
    **--pop-ind-file** *<pop_str>* *<ind_file>*
        Argument used to assign multiple individuals to thr population begin 
        updated using a file.
    **--pop-rm-ind** *<pop_str>* *<ind_str>*
        Argument used to remove a individual to the population begin updated.
    **--pop-rm-inds** *<pop_str>* *<ind1_str>* *<ind2_str>* ..
        Argument used to remove multiple individuals to the population begin 
        updated.
    **--pop-rm-ind-file** *<pop_str>* *<ind_file>*
        Argument used to remove multiple individuals to thr population begin 
        updated using a file.

    #############################
    Output Command-line Arguments
    #############################
    **--out** *<output_filename>*
        Argument used to define the complete output filename.
    **--overwrite**
        Argument used to define if previous output should be overwritten.
'''

import os
import sys
import json
import argparse
import logging
import itertools

from collections import defaultdict, OrderedDict

from pgpipe.model import Model, ModelFile, write_model_file, read_model_file
from pgpipe.logging_module import initLogger, logArgs
from pgpipe.misc import argprase_kwargs

def model_creator_parser (passed_arguments = []):
    '''VCF Argument Parser - Assigns arguments from command line'''

    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError('%s not found' % value)
                setattr(args, self.dest, value)
        return customAction

    def parser_dict_list_append ():
        '''Custom action to house data as defaultdict list'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, values, option_string=None):
                if getattr(args, self.dest):
                    if values[1] not in getattr(args, self.dest)[values[0]]:
                        # Append the argument with the file
                        getattr(args, self.dest)[values[0]].append(values[1])
                else:
                    # Set the argument with the file (as list)
                    arg_dict = defaultdict(list)
                    arg_dict[values[0]].append(values[1])
                    setattr(args, self.dest, arg_dict)
        return customAction

    def parser_dict_list_extend ():
        '''Custom action to house data as defaultdict list'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, values, option_string=None):
                if getattr(args, self.dest):
                    for value in values[1:]:
                        if value not in getattr(args, self.dest)[values[0]]:
                            # Append the argument with the file
                            getattr(args, self.dest)[values[0]].append(value)
                else:
                    # Set the argument with the file (as list)
                    arg_dict = defaultdict(list)
                    arg_dict[values[0]].extend(values[1:])
                    setattr(args, self.dest, arg_dict)
        return customAction

    def parser_dict_str ():
        '''Custom action to house data as defaultdict int'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, values, option_string=None):
                # Assign the passed value
                if getattr(args, self.dest):
                    # Append the argument with the file
                    getattr(args, self.dest)[values[0]] = values[1]

                else:
                    # Set the argument with the file (as list)
                    arg_dict = defaultdict(str)
                    arg_dict[values[0]] = values[1]
                    setattr(args, self.dest, arg_dict)
        return customAction

    def parser_dict_int ():
        '''Custom action to house data as defaultdict int'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, values, option_string=None):

                # Check the passed value is an int
                try:
                    int(values[1])
                except:
                    raise Exception('--%s only accepts integers' % self.dest)

                # Assign the passed value
                if getattr(args, self.dest):
                    # Append the argument with the file
                    getattr(args, self.dest)[values[0]] = int(values[1])

                else:
                    # Set the argument with the file (as list)
                    arg_dict = defaultdict(int)
                    arg_dict[values[0]] = int(values[1])
                    setattr(args, self.dest, arg_dict)
        return customAction

    def parser_dict_file ():
        '''Custom action to house data as defaultdict int'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, values, option_string=None):

                if not os.path.isfile(values[1]):
                    raise IOError('%s not found' % values[1])

                # Assign the passed value
                if getattr(args, self.dest):
                    # Append the argument with the file
                    getattr(args, self.dest)[values[0]] = values[1]

                else:
                    # Set the argument with the file (as list)
                    arg_dict = defaultdict(str)
                    arg_dict[values[0]] = values[1]
                    setattr(args, self.dest, arg_dict)
        return customAction

    def metavar_list (var_list):
        '''Create a formmated metavar list for the help output'''
        return '{' + ', '.join(var_list) + '}'

    model_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    # General model arguments
    model_parser.add_argument('--model-file', help = 'Defines the name of a model file', type = str, action = parser_confirm_file())
    model_parser.add_argument('--model', help = 'Defines the name of a model', type = str, action = 'append', default = [])
    model_parser.add_argument('--update-model', help = "Specifies a model to update. Allows for: i) tree update and ii) populations to be added", type = str, default = [], action = 'append')
    model_parser.add_argument('--update-pop', help = "Specifies a population to update. Allows for individuals to be added", type = str, default = [], action = 'append')

    # Tree arguments
    model_parser.add_argument('--model-tree', metavar = ('MODEL', 'NEWICK_TREE'), help = 'Assign population tree string to a model, in Newick format', type = str, nargs = 2, action = parser_dict_str())
    model_parser.add_argument('--model-tree-file',  metavar = ('MODEL', 'NEWICK_FILE'), help = 'Assign population tree file to a model, in Newick format', type = str, nargs = 2, action = parser_dict_file())

    # Population arguments
    model_parser.add_argument('--model-pop', metavar = ('MODEL', 'POP'), dest = 'pops', help = 'Assign a population name to a model', type = str, nargs = 2, action = parser_dict_list_append())
    model_parser.add_argument('--model-pops', metavar = ('MODEL', 'POP'), dest = 'pops', help = 'Assign multiple population names to a model', type = str, nargs = '+', action = parser_dict_list_extend())
    model_parser.add_argument('--model-pop-file', metavar = ('MODEL', 'POP_FILE'), help = 'Assign multiple population names to a model usign a file', type = str, nargs = 2, action = parser_dict_file())

    model_parser.add_argument('--model-rm-pop', metavar = ('MODEL', 'POP'), dest = 'rm_pops', help = 'Assign a population name to a model', type = str, nargs = 2, action = parser_dict_list_append())
    model_parser.add_argument('--model-rm-pops', metavar = ('MODEL', 'POP'), dest = 'rm_pops', help = 'Assign multiple population names to a model', type = str, nargs = '+', action = parser_dict_list_extend())
    model_parser.add_argument('--model-rm-pop-file', metavar = ('MODEL', 'POP_FILE'), help = 'Assign multiple population names to a model usign a file', type = str, nargs = 2, action = parser_dict_file())

    # Individual arguments
    model_parser.add_argument('--pop-ind', metavar = ('POP', 'IND'), dest = 'inds', help = 'Assign an individual name to a population', type = str, nargs = 2, default = {}, action = parser_dict_list_append())
    model_parser.add_argument('--pop-inds', metavar = ('POP', 'IND'), dest = 'inds', help = 'Assign multiple individual names to a population', type = str, nargs = '+', default = {}, action = parser_dict_list_extend())
    model_parser.add_argument('--pop-ind-file', metavar = ('POP', 'IND_FILE'), help = 'Assign multiple individual names to a population using a file', type = str, nargs = 2, default = {}, action = parser_dict_file())

    model_parser.add_argument('--pop-rm-ind', metavar = ('POP', 'IND'), dest = 'rm_inds', help = 'Assign an individual name to a population', type = str, nargs = 2, action = parser_dict_list_append())
    model_parser.add_argument('--pop-rm-inds', metavar = ('POP', 'IND'), dest = 'rm_inds', help = 'Assign multiple individual names to a population', type = str, nargs = '+', action = parser_dict_list_extend())
    model_parser.add_argument('--pop-rm-ind-file', metavar = ('POP', 'IND_FILE'), help = 'Assign multiple individual names to a population using a file', type = str, nargs = 2, action = parser_dict_file())

    # Output Arguments
    model_parser.add_argument('--out', help = 'Specifies the complete output filename.', type = str, default = 'out.model')
    model_parser.add_argument('--overwrite', help = "Specifies if previous output files should be overwritten", action = 'store_true')

    if passed_arguments:
        return vars(model_parser.parse_args(passed_arguments))
    else:
        return vars(model_parser.parse_args())

def incompatible_duplicates_check (term, *arguments_to_test):

    # Duplicate count
    duplicates_count = 0

    # Loop the passed arguments
    for argument_to_test in arguments_to_test:

        # Check if the argument has been populated
        if argument_to_test:

            # Check if the argument has the term
            if term in argument_to_test:

                # Add to the duplicate count if the incompatible term is found
                duplicates_count += 1

    # Check if there were duplicates_count
    if duplicates_count > 1:
        return True
    else:
        return False

def assigment_check (term, *arguments_to_be_assigned):

    # Term assignment boolean
    term_assigned = False

    # Loop the passed arguments
    for argument_to_be_assigned in arguments_to_be_assigned:

        # Check if the argument has been populated
        if argument_to_be_assigned:

            # Check if the argument has the term
            if term in argument_to_be_assigned:

                # Mark the term assigned
                term_assigned = True

    # Return the boolean
    return term_assigned

def dict_argument_to_list (term, list_argument):

    # Check if either the argument or term was not assigned
    if not list_argument or term not in list_argument:
        return []

    # Return data as list
    return list_argument[term]

def dict_argument_to_str (term, str_argument):

    # Check if either the argument or term was not assigned
    if not str_argument or term not in str_argument:
        return ''

    # Return data as str
    return str_argument[term]

def file_dict_argument_to_list (term, file_argument):

    # Check if either the argument or term was not assigned
    if not file_argument or term not in file_argument:
        return []

    # List to hold the file data
    data_list = []

    # Check if running via python 3
    if sys.version_info[0] == 3:

        # Open the file and save the data
        with open(file_argument[term], 'r') as data_file:
            for data_line in data_file:
                data_list.append(data_line.strip())

    # Check if running via python 2
    else:

        # Open the file and save the data
        with open(file_argument[term], 'rb') as data_file:
            for data_line in data_file:
                data_list.append(data_line.strip())

    # Return data as list
    return data_list

def file_dict_argument_to_str (term, file_argument):

    # Check if either the argument or term was not assigned
    if not file_argument or term not in file_argument:
        return ''

    # str to hold the file data
    data_str = ''

    # Check if running via python 3
    if sys.version_info[0] == 3:

        # Open the file and save the data
        with open(file_argument[term], 'r') as data_file:
            for data_line in data_file:
                data_str = data_line.strip()

    # Check if running via python 2
    else:

        # Open the file and save the data
        with open(file_argument[term], 'rb') as data_file:
            for data_line in data_file:
                data_str = data_line.strip()

    # Return data as str
    return data_str

def run (**kwargs):
    '''

    Parameters
    ----------
    --model : str
        Name of the model
    --model-file : str
        Name of a model file
    --update-model : str
        Name of a model to update. Allows for: i) tree update 
        and ii) populations to be added
    --update-pop : str
        Name of a population to update. Allows for individuals 
        to be added
    --model-tree : str  str
        Assign newick-formatted population tree to model
        Example: --model-tree model_name newick_tree
    --model-tree-file : str  str
        Assign newick-formatted population tree file to model
        Example: --model-tree-file model_name newick_file
    --model-pop : str  str
        Assign single population name to a model
        Example: --model-pop model_name pop_name
    --model-pops : str  list
        Assign multiple population names to a model
        Example: --model-pops model_name pop_name1 pop_name2 etc.
    --model-pop-file : str  str
        Assign population names to a model using a file
        Example: --model-pop-file model_name pop_file
    --pop-ind : str  str
        Assign single individual name to a population
        Example: --pop-ind pop_name ind_name
    --pop-inds : str  list
        Assign multiple population names to a model
        Example: --pop-inds pop_name ind_name1 ind_name2 etc.
    --pop-ind-file : str  str
        Assign population names to a model using a file
        Example: --pop-ind-file pop_name ind_file
    --out : str
        Filename of the model output

    Returns
    -------

    Raises
    ------
    Exception
        No tree assigned to model
    Exception
        Multiple trees assigned to model
    Exception
        No population assigned to model
    Exception
        No individuals assigned to population
    '''

    # Update kwargs with defaults
    if __name__ != "__main__":
        kwargs = argprase_kwargs(kwargs, model_creator_parser)

    # Assign arguments
    creator_args = argparse.Namespace(**kwargs)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(creator_args, 'model_creator')

    # Create population dictionary
    pop_dict = {}

    # Check if a model file was given
    if creator_args.model_file:

        # Load the model file
        model_file = read_model_file(creator_args.model_file)

        # Add the pops to the dict, return error is different populations exist with the same name
        for model in model_file.values():

            for pop, ind_list in model.ind_dict.items():

                if pop not in pop_dict:

                    pop_dict[pop] = ind_list

                else:

                    if pop_dict[pop] != ind_list:
                         raise Exception('Population (%s) differ between models. Please check input' % pop)

        # Combine all models that need to be created/edited
        models = list(set(creator_args.model + list(model_file)))

    else:

        # Create a new model file
        model_file = ModelFile()

        # Assign the models
        models = creator_args.model

    # Loop each model specified to confirm parameters are valid
    for model in models:

        # Create an empty list
        combined_pops = []

        # Check if the model is being updated or created
        if model in creator_args.model or model in creator_args.update_model:

            # Check that a tree has been assigned to the model
            if model not in model_file and not assigment_check(model, creator_args.model_tree_file, creator_args.model_tree):
                logging.warning('No tree assigned to model (%s)' % model)

            # Check that a tree has been assigned to the model
            if model in model_file and not model_file[model].tree and not assigment_check(model, creator_args.model_tree_file, creator_args.model_tree) and not assigment_check(model, creator_args.update_model):
                logging.warning('No tree assigned to model (%s)' % model)

            # Check that multiple trees have not been assigned to the model
            if incompatible_duplicates_check(model, creator_args.model_tree_file, creator_args.model_tree):
                raise Exception('Multiple trees assigned to model (%s)' % model)

            # Check if the model exists, but isn't being updated
            if model in model_file and not assigment_check(model, creator_args.update_model):
                raise Exception('%s already exists. Please use --update-model to update' % model)

            # Check if a pop assignment arguments has been used
            if model in creator_args.model and model in model_file and not assigment_check(model, creator_args.pops):
                raise Exception('No populations %s' % model)

            # Convert the pop dict into a list
            combined_pops = dict_argument_to_list(model, creator_args.pops)

            # Convert the pop file into a list
            pop_file_list = file_dict_argument_to_list(model, creator_args.model_pop_file)

            # Combine the pop arguments
            combined_pops.extend([file_pop for file_pop in pop_file_list if file_pop not in combined_pops])

            # Assign the model if it already exists
            if model in model_file:

                # Open the model
                model_data = model_file[model]

                # Convert the pop dict into a list
                rm_combined_pops = dict_argument_to_list(model, creator_args.rm_pops)

                # Convert the pop file into a list
                rm_pop_file_list = file_dict_argument_to_list(model, creator_args.model_rm_pop_file)

                # Combine the pop arguments
                rm_combined_pops.extend([rm_file_pop for rm_file_pop in rm_pop_file_list if rm_file_pop not in rm_combined_pops])

                # Check if a population being removed is also being added
                if set(combined_pops) & set(rm_combined_pops):
                    logging.warning('Population being added and removed in same model (%s)' % model)

                # Check if there are populations to remove
                if rm_combined_pops:

                    # Remove the populations from the model
                    for rm_pop in rm_combined_pops:
                        model_data.remove_pop(rm_pop)

                        logging.info('Population (%s) removed for model (%s)' % (rm_pop, model))

            else:

                # Create the model
                model_data = Model(name = model)

                logging.info('Model (%s) created' % model)

            # Convert the inds dict into a list
            tree_str = dict_argument_to_str(model, creator_args.model_tree)

            # Convert the pop file into a list
            tree_file_str = file_dict_argument_to_str(model, creator_args.model_tree_file)

            # Assign the tree
            model_data.assign_tree(tree = tree_str + tree_file_str)

        else:

            # Open the model for possible pop editing
            model_data = model_file[model]

        # Check if pop_dict is populated and assign the pops
        if pop_dict:

            pops = list(pop_dict)

        else:
            # Assign the pops from the model otherwise
            pops = combined_pops

        # Loop each population in the model
        for pop in pops:

            # Check if the population is being added or updated
            if pop in combined_pops or pop in creator_args.update_pop or pop in list(creator_args.inds) or pop in list(creator_args.pop_ind_file):

                # Check if a individual assignment arguments has been used
                if not pop_dict and not assigment_check(pop, creator_args.inds, creator_args.pop_ind_file):
                    raise Exception('No individuals assigned to pop (%s)' % pop)

                # Check if the pop exists, but isn't being updated
                if pop in model_data.pop_list and not assigment_check(pop, creator_args.update_pop):
                    raise Exception('Pop (%s) already exists. Please use --update-pop to update' % pop)

                # Convert the inds dict into a list
                combined_inds = dict_argument_to_list(pop, creator_args.inds)

                # Convert the inds file into a list
                ind_file_list = file_dict_argument_to_list(pop, creator_args.pop_ind_file)

                # Combine the inds arguments
                combined_inds.extend([file_ind for file_ind in ind_file_list if file_ind not in combined_inds])

                # Update the pop if it already exists
                if pop in creator_args.update_pop:

                    # Convert the inds dict into a list
                    rm_combined_inds = dict_argument_to_list(pop, creator_args.rm_inds)

                    # Convert the inds file into a list
                    rm_ind_file_list = file_dict_argument_to_list(pop, creator_args.pop_rm_ind_file)

                    # Combine the pop arguments
                    rm_combined_inds.extend([rm_file_ind for rm_file_ind in rm_ind_file_list if rm_file_ind not in rm_combined_inds])

                    # Check if a population being removed is also being added
                    if set(combined_inds) & set(rm_combined_inds):
                        logging.warning('Individual being added and removed in same population (%s)' % pop)

                    # Check if there are populations to remove
                    if rm_combined_inds:

                        # Update the population
                        model_data.update_pop(pop = pop, rm_inds = rm_combined_inds)

                        logging.info('Populations removed for population (%s)' % pop)

                    # Update the population
                    model_data.update_pop(pop = pop, inds = combined_inds)

                    logging.info('Individuals updated in population (%s)' % pop)

                else:

                    # Check if individuals were assigned
                    if combined_inds:

                        # Assign the population
                        model_data.assign_pop(pop = pop, inds = combined_inds)

                    else:

                        # Assign the population using the pop_dict
                        model_data.assign_pop(pop = pop, inds = pop_dict[pop])

                    logging.info('Individuals assigned to population (%s)' % pop)

                logging.info('Population (%s) assigned to model (%s)' % (pop, model))


        # If the model has a tree, check that is has all pops
        if model_data.tree:
            if not all([pop in model_data.tree for pop in model_data.pop_list]):
                raise Exception('Tree (%s) does not contain all populations' % model_data.tree)

        # Add the model to the model file
        model_file[model] = model_data

        logging.info('Model (%s) assigned to model file' % model)

    # Create the model file
    write_model_file(model_file, creator_args.out, overwrite = creator_args.overwrite)


if __name__ == "__main__":
    initLogger()
    run(**model_creator_parser())
