import os
import sys
import json
import subprocess
import argparse
import logging
import itertools
import copy

import numpy as np

from collections import defaultdict, OrderedDict

# Import logging module
from pgpipe.logging_module import initLogger

class ModelFile(dict):
    def __init__(self, *arg, **kw):
        super(ModelFile, self).__init__(*arg, **kw)
        self.inds = []
        self.ind_file = ''
        self.exclude_file = ''

        if arg and self.confirm_model_instance(arg[1]):
            self.update_inds(arg[1])

    def __setitem__(self, *arg, **kw):
        super(ModelFile, self).__setitem__(*arg, **kw)

        if arg and self.confirm_model_instance(arg[1]):
            self.update_inds(model = arg[1])

    def __delitem__(self, key):
        super(ModelFile, self).__delitem__(key)
        self.update_inds()

    def confirm_model_instance (self, unknown):

        if isinstance(unknown, Model):

            return True

        else:

            return False

    def copy_model (self, src_model_name, new_model_name):

        src_model = super(ModelFile, self).__getitem__(src_model_name)

        src_model_copy = copy.deepcopy(src_model)

        src_model_copy.name = new_model_name

        super(ModelFile, self).__setitem__(new_model_name, src_model_copy)

    def rename_model (self, src_model_name, new_model_name):

        src_model = super(ModelFile, self).pop(src_model_name)

        src_model.name = new_model_name

        super(ModelFile, self).__setitem__(new_model_name, src_model)

    def update_inds (self, model = None):

        if self.confirm_model_instance(model):

            # Return error if inds is empty
            if not model.inds:
                raise IOError('No individuals found in %s.' % model.name)

            # Create a list of the unique individuals
            unique_inds = list(set(self.inds + model.inds))

        else:

            # Create an empty list for the unique individuals
            unique_inds = []

            # Loop the models in the file
            for model_in_file in super(ModelFile, self).values():

                # Create a list of the unique individuals
                unique_inds = list(set(unique_inds + model_in_file.inds))


        # Store the individuals
        self.inds = unique_inds

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

    def delete_exclude_ind_file (self):
        # Check if an individuals file was created
        if self.exclude_file:

            # Delete the individuals file
            os.remove(self.exclude_file)

            # Remove the filename
            self.exclude_file = ''

    def to_json (self):

        model_file_json = []

        for model_name, model_data in super(ModelFile, self).items():
            model_file_json.append(model_data.to_json())

        return model_file_json

class Model:
    def __init__ (self, name):
        self.name = name
        self.tree = ''
        self.pop_list = []
        self.ind_dict = defaultdict(list)
        self.nind = defaultdict(int)
        self.pop_files = []
        self.ind_file = ''

    def __str__ (self):
        return self.name

    @property
    def npop (self):
        return len(self.pop_list)

    @property
    def inds (self):
        return list(itertools.chain.from_iterable(self.ind_dict.values()))

    def assign_tree (self, tree):
        self.tree = str(tree)

    def assign_pop (self, pop, inds = []):
        self.pop_list.append(str(pop))
        if inds:
            self.ind_dict[pop] = [str(ind) for ind in inds]
        self.nind[pop] = len(self.ind_dict[pop])

    def remove_pop (self, pop):

        # Confirm the pop is in the model
        if str(pop) not in self.pop_list:

            # Raise error if pop not found
            raise Exception('%s not found' % pop)
            
        self.pop_list.remove(str(pop))
        del self.ind_dict[pop]
        del self.nind[pop]

    def update_pop (self, pop, inds = [], rm_inds = []):
        
        # Confirm the pop is in the model
        if str(pop) not in self.pop_list:

            # Raise error if pop not found
            raise Exception('%s not found' % pop)

        if inds:
            self.ind_dict[pop].extend([str(ind) for ind in inds])

        if rm_inds:
            for rm_ind in rm_inds:
                if str(rm_ind) in self.ind_dict[pop]: 
                    self.ind_dict[pop].remove(str(rm_ind))
            
        self.nind[pop] = len(self.ind_dict[pop])

    def sample_pop (self, pop, sample_size, with_replacements = False):

        # Confirm the pop is in the model
        if str(pop) not in self.pop_list:

            # Raise error if pop not found
            raise Exception('%s not found' % pop)

        # Confirm the sample size is an int
        try:

            sample_size = int(sample_size)

        except:

            # Raise error if sample_size not an int
            raise Exception('%s not int' % sample_size)

        # Check if the sample size is larger than the pop
        if int(sample_size) > self.nind[pop]:

            # Raise error if sample_size is larger
            raise Exception('%s is larger than %s' % (sample_size, pop))

        # Use numpy choice to randomly sample the pop
        sampled_inds = np.random.choice(self.ind_dict[pop], sample_size, replace = with_replacements)

        # Save the sampled inds as a list
        self.ind_dict[pop] = list(sampled_inds)

    def sample_pops (self, sample_size, with_replacements = False):

        # Confirm the sample size is an int
        try:

            sample_size = int(sample_size)

        except:

            # Raise error if sample_size not an int
            raise Exception('%s not int' % sample_size)

        # Loop each pop in the pop list
        for pop in self.pop_list:

            # Check if the sample size is larger than the pop
            if int(sample_size) > self.nind[pop]:

                # Raise error if sample_size is larger
                raise Exception('%s is larger than %s' % (sample_size, pop))

        # Loop each pop in the pop list, if no error raised
        for pop in self.pop_list:

            # Use numpy choice to randomly sample the pop
            sampled_inds = np.random.choice(self.ind_dict[pop], sample_size, replace = with_replacements)

            # Save the sampled inds as a list
            self.ind_dict[pop] = list(sampled_inds)

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

        logging.info('Population files created for %s' % self)

    def delete_pop_files (self):
        # Check if pop files were created
        if len(self.pop_files) != 0:

            # Loop the created pop files
            for pop_file in self.pop_files:
                # Delete the pop file
                os.remove(pop_file)

            # Remove the filenames
            self.pop_files = []

            logging.info('Population files deleted for %s' % self)

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

        logging.info('Individuals file created for %s' % self)

    def delete_ind_file (self):
        # Check if an individuals file was created
        if self.ind_file:

            # Delete the individuals file
            os.remove(self.ind_file)

            # Remove the filename
            self.ind_file = ''

            logging.info('Individuals file deleted for %s' % self)

    def to_json (self):

        model_json = OrderedDict()

        model_json['name'] = self.name

        pop_json = OrderedDict()

        model_json['tree'] = self.tree

        for pop in self.pop_list:

            pop_json[pop] = OrderedDict()

            pop_json[pop]['indv'] = self.ind_dict[pop]

        model_json['pops'] = pop_json

        return model_json

    def return_pop (self, ind):

        # Loop the pops within the model
        for pop, inds in self.ind_dict.items():

            # Check if the ind belongs to the pop
            if ind in inds:

                return pop

        raise Exception('Individual (%s) not found within %s' % (ind, self.name))

def pops_not_in_model (model, pops):

    # Check if a pop in pops is not within model.pop_list
    pops_not_found = set(pops) - set(model.pop_list)

    # Check any pops were not found
    if pops_not_found:
        return pops_not_found

    return None

def read_model_file (filename):
    """
    returns a ModelFile, which is a dictionary    
    """

    # Check that the file exists
    if not os.path.isfile(filename):
        raise IOError

    # Create ModelFile object
    models_to_return = ModelFile()

    # Check if using python 2 or 3
    if sys.version_info[0] == 2:
        # Open the model file in python 2
        model_file = open(filename, 'rU')
    else:
        # Open the model file in python 3
        model_file = open(filename, 'r', newline=None)

    # Parse the model file using the json reader
    models_dict = json.load(model_file)
    model_file.close()

    # List to store all unique individuals (i.e. individuals in all models)
    individual_list = []

    # Loop the parsed models
    for model_dict in models_dict:

        # Create the model
        model = Model(str(model_dict['name']))

        # Check if model has a tree
        if 'tree' in model_dict:

            # Assign the model tree
            model.assign_tree(model_dict['tree'])

        # Loop the populations in the model
        for pop, pop_dict in model_dict['pops'].items():

            # Convert all individuals names to str
            ind_list = [str(pop_ind) for pop_ind in pop_dict['inds']]

            # Assign the population ans it's individuals to the model
            model.assign_pop(str(pop), ind_list)

            # Assign the individuals to the unique individual list
            individual_list.extend(ind_list)

        # Remove duplicates from the unique individual list
        individual_list = list(set(individual_list))

        # Save the model
        models_to_return[str(model.name)] = model

    logging.info('Finished reading model file (%s)' % filename)

    # Return the models
    return models_to_return

def read_single_model(filename, popname = None):
    #Returns single model from file, assumes either only one
    #in file or popname is provided
    popmodels = read_model_file(filename)
    if popname is None:
        if len(popmodels) > 1:
            raise Exception("Model file %s has %d models, must specify which to use" % (filename,len(popmodels)))
        pp = list(popmodels.keys())
        return popmodels[pp[0]]
    else:
        if popname not in popmodels.keys():
            raise Exception("Model %s not found in model file %s" % (popname,filename))
        return popmodels[popname]

def write_model_file (model_file, filename, overwrite = False):

    # Check if the file is to be overwritten
    if not overwrite:

        # Check if the file exists
        if os.path.exists(filename):
            raise Exception('%s already exists' % filename)

    # Open the output file
    output_file = open(filename, 'w')

    # Write the json-formmated data to the output file
    output_file.write(json.dumps(model_file.to_json(), indent = 4))

    # Close the output file
    output_file.close()

    logging.info('Finished writing model file (%s)' % filename)
