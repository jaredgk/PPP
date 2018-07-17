import os
import sys
import json
import subprocess
import argparse
import logging
import itertools

from collections import defaultdict

# Insert Jared's directory path, required for calling Jared's functions. Change when directory structure changes.
sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'jared')))

from logging_module import initLogger

class ModelFile(dict):
    def __init__(self, *arg, **kw):
        super(ModelFile, self).__init__(*arg, **kw)
        self.inds = []
        self.ind_file = ''
        self.exclude_file = ''

    def assign_inds (self, inds = []):
        # Return error if inds is empty
        if not inds:
            raise IOError('No individuals found in the model file.')
        # Store the individuals
        self.inds = [str(ind) for ind in inds]

    def create_ind_file (self, file_ext = '', file_path = '', overwrite = False):
        # Assign the filename for the population file
        ind_filename = 'unique_individuals' + file_ext

        # If a path is assigned, create the file at the specified location
        if file_path:
            ind_filename = os.path.join(file_path, ind_filename)

        # Check if previous files should be overwriten
        if not overwrite:
            # Check if the file already exists
            if os.path.isfile(ind_filename):
                raise IOError('Individuals file exists.')

        # Create the population file
        ind_file = open(ind_filename, 'w')
        ind_file.write('%s\n' %'\n'.join(self.inds))
        ind_file.close()

        # Save the individuals filename
        self.ind_file = ind_filename

    def delete_ind_file (self):
        # Check if an individuals file was created
        if self.ind_file:

            # Delete the individuals file
            os.remove(self.ind_file)

            # Remove the filename
            self.ind_file = ''

    def create_exclude_ind_file (self, inds_to_include = [], file_ext = '', file_path = '', overwrite = False):
        # Assign the filename for the population file
        ind_filename = 'exclude_individuals' + file_ext

        # If a path is assigned, create the file at the specified location
        if file_path:
            ind_filename = os.path.join(file_path, ind_filename)

        # Check if previous files should be overwriten
        if not overwrite:
            # Check if the file already exists
            if os.path.isfile(ind_filename):
                raise IOError('Individuals file exists.')

        # Create exclude list by removing included individuals
        exclude_inds = list(set(self.inds) - set(inds_to_include))

        # Create the population file
        ind_file = open(ind_filename, 'w')
        ind_file.write('%s\n' %'\n'.join(exclude_inds))
        ind_file.close()

        # Save the individuals filename
        self.exclude_file = ind_filename

    def delete_ind_file (self):
        # Check if an individuals file was created
        if self.exclude_file:

            # Delete the individuals file
            os.remove(self.exclude_file)

            # Remove the filename
            self.exclude_file = ''

class Model:
    def __init__ (self, name):
        self.name = name
        self.tree = ''
        self.npop = 0
        self.pop_list = []
        self.nind = defaultdict(int)
        self.ind_dict = defaultdict(list)
        self.pop_files = []
        self.ind_file = ''

    @property
    def inds(self):
        return list(itertools.chain.from_iterable(self.ind_dict.values()))

    def assign_tree (self, tree):
        self.tree = str(tree)

    def assign_pop (self, pop, inds = []):
        self.npop += 1
        self.pop_list.append(str(pop))
        if inds:
            self.nind[pop] = len(inds)
            self.ind_dict[pop] = [str(ind) for ind in inds]

    def create_pop_files (self, file_ext = '', file_path = '', overwrite = False):
        for pop in self.pop_list:
            # Assign the filename for the population file
            pop_filename = pop + file_ext

            # If a path is assigned, create the file at the specified location
            if file_path:
                pop_filename = os.path.join(file_path, pop_filename)

            # Check if previous files should be overwriten
            if not overwrite:
                # Check if the file already exists
                if os.path.isfile(pop_filename):
                    raise IOError('Population file exists.')

            # Create the population file
            pop_file = open(pop_filename, 'w')
            pop_file.write('%s\n' %'\n'.join(self.ind_dict[pop]))
            pop_file.close()

            # Save the population filename
            self.pop_files.append(pop_filename)

    def delete_pop_files (self):
        # Check if pop files were created
        if len(self.pop_files) != 0:

            # Loop the created pop files
            for pop_file in self.pop_files:
                # Delete the pop file
                os.remove(pop_file)

            # Remove the filenames
            self.pop_files = []

    def create_ind_file (self, file_ext = '', file_path = '', overwrite = False):
        # Assign the filename for the population file
        ind_filename = 'individual.keep' + file_ext

        # If a path is assigned, create the file at the specified location
        if file_path:
            ind_filename = os.path.join(file_path, ind_filename)

        # Check if previous files should be overwriten
        if not overwrite:
            # Check if the file already exists
            if os.path.isfile(ind_filename):
                raise IOError('Individuals file exists.')

        # Create the population file
        ind_file = open(ind_filename, 'w')
        ind_file.write('%s\n' %'\n'.join(self.inds))
        ind_file.close()

        # Save the individuals filename
        self.ind_file = ind_filename

    def delete_ind_file (self):
        # Check if an individuals file was created
        if self.ind_file:

            # Delete the individuals file
            os.remove(self.ind_file)

            # Remove the filename
            self.ind_file = ''

def read_model_file (model_filename):

    # Check that the file exists
    if not os.path.isfile(model_filename):
        raise IOError

    # Create ModelFile object
    models_to_return = ModelFile()

    # Check if using python 2 or 3
    if sys.version_info[0] == 2:
        # Open the model file in python 2
        model_file = open(model_filename, 'rU')
    else:
        # Open the model file in python 3
        model_file = open(model_filename, 'r', newline=None)

    # Parse the model file using the json reader
    models_dict = json.load(model_file)

    # List to store all unique individuals (i.e. individuals in all models)
    individual_list = []

    # Loop the parsed models
    for model_dict in models_dict:

        # Create the model
        model = Model(model_dict['name'])

        # Loop the populations in the model
        for pop, pop_dict in model_dict['pops'].items():

            # Assign the population ans it's individuals to the model
            model.assign_pop(pop, pop_dict['inds'])
            # Assign the individuals to the unique individual list
            individual_list.extend(pop_dict['inds'])

        # Remove duplicates from the unique individual list
        individual_list = list(set(individual_list))

        # Save the model
        models_to_return[str(model.name)] = model

    # Store the  unique individuals within the ModelFile object
    models_to_return.assign_inds(individual_list)

    # Return the models
    return models_to_return
