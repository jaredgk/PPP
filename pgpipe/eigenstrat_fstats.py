import os
import sys
import argparse
import logging
import subprocess
import random
import string
import shutil

# Import PPP modules and scripts
from admixtools import *
from pgpipe.model import read_model_file, pops_not_in_model
from pgpipe.logging_module import initLogger, logArgs

def admix_parser (passed_arguments):
    '''admix Argument Parser - Assigns arguments from command line'''

    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError('%s not found' % value)
                setattr(args, self.dest, value)
        return customAction

    def parser_add_to_list ():
        '''Custom action to add items to a list'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):

                # Clean up any commas
                value = [item.strip(',') for item in value]

                if not getattr(args, self.dest):
                    setattr(args, self.dest, value)
                else:
                    getattr(args, self.dest).extend(value)
        return customAction

    def metavar_list (var_list):
        '''Create a formmated metavar list for the help output'''
        return '{' + ', '.join(var_list) + '}'

    admix_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Input BED arguments
    admix_parser.add_argument("--geno", dest = 'geno_filename', help = "Defines the filename of the GENO file. Called alongside --ind and --snp", type = str, action = parser_confirm_file())
    admix_parser.add_argument("--ind", dest = 'ind_filename', help = "Defines the filename of the IND file. Called alongside --geno and --snp", type = str, action = parser_confirm_file())
    admix_parser.add_argument("--snp", dest = 'snp_filename', help = "Defines the filename of the SNP file. Called alongside --geno and --ind", type = str, action = parser_confirm_file())

    # Input Eigenstrat prefix arguments
    admix_parser.add_argument("--eigenstrat-prefix", help = "Defines the filename prefix of the GENO, IND, and SNP files", type = str)

    # Model file arguments.
    admix_parser.add_argument('--model-file', help = 'Defines the model file', type = str, action = parser_confirm_file())
    admix_parser.add_argument('--model', help = 'Defines the model and the individual(s) to include', type = str)

    admix_parser.add_argument('--out-prefix', help = 'Defines the output prefix (i.e. filename without file extension)', default = 'out')

    # General arguments.
    admix_parser.add_argument('--overwrite', help = "Overwrite previous output files", action = 'store_true')

    # Statistic based arguments
    stat_list = ['D', 'F4', 'F4-ratio', 'F3']
    admix_parser.add_argument('--calc-admix-statistic', metavar = metavar_list(stat_list), help = 'The admixture statistic to calculate', required = True, type = str, choices = stat_list)
    
    # Admix analysis pops
    admix_parser.add_argument('--admix-w-pop', dest = 'w_pops', help = 'W population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admix_parser.add_argument('--admix-w-pops', dest = 'w_file', help = 'W populations for admixure analysis', type = str, action = parser_confirm_file())
    
    admix_parser.add_argument('--admix-x-pop', dest = 'x_pops', help = 'X population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admix_parser.add_argument('--admix-x-pops', dest = 'x_file', help = 'X populations for admixure analysis', type = str, action = parser_confirm_file())

    admix_parser.add_argument('--admix-y-pop', dest = 'y_pops', help = 'Y population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admix_parser.add_argument('--admix-y-pops', dest = 'y_file', help = 'Y populations for admixure analysis', type = str, action = parser_confirm_file())

    admix_parser.add_argument('--admix-z-pop', dest = 'z_pops', help = 'Z population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admix_parser.add_argument('--admix-z-pops', dest = 'z_file', help = 'Z populations for admixure analysiss', type = str, action = parser_confirm_file())

    admix_parser.add_argument('--admix-a-pop', dest = 'a_pops', help = 'A population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admix_parser.add_argument('--admix-a-pops', dest = 'a_file', help = 'A populations for admixure analysis', type = str, action = parser_confirm_file())

    admix_parser.add_argument('--admix-b-pop', dest = 'b_pops', help = 'B population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admix_parser.add_argument('--admix-b-pops', dest = 'b_file', help = 'B populations for admixure analysis', type = str, action = parser_confirm_file())

    admix_parser.add_argument('--admix-c-pop', dest = 'c_pops', help = 'C population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admix_parser.add_argument('--admix-c-pops', dest = 'c_file', help = 'C populations for admixure analysis', type = str, action = parser_confirm_file())

    admix_parser.add_argument('--admix-o-pop', dest = 'o_pops', help = 'O population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admix_parser.add_argument('--admix-o-pops', dest = 'o_file', help = 'O populations for admixure analysis', type = str, action = parser_confirm_file())

    if passed_arguments:
        return admix_parser.parse_args(passed_arguments)
    else:
        return admix_parser.parse_args()

def read_admix_pops_file (pop_filename):

    # List to hold admix pops 
    admix_pops = []

    # Open the ind file
    with open(pop_filename, 'r') as pop_file:

        # Loop the ind file
        for pop_line in pop_file:

            # Append the population
            admix_pops.append(pop_line.strip())

    # Return the admix pops
    return admix_pops

def run(passed_arguments = []):

    # Grab admixtools arguments from command line
    admix_args = admix_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(admix_args, func_name='admixtools')

    # Read in the models
    models_in_file = read_model_file(admix_args.model_file)

    # Check that the selected model was not found in the file
    if admix_args.model not in models_in_file:
        raise IOError('Selected model "%s" not found in: %s' % (admix_args.model, admix_args.model_file))

    # Select model, might change this in future versions
    selected_model = models_in_file[admix_args.model]
      
    # Check if the W populations file has been assigned
    if admix_args.w_file:

        # Save the pops from the file
        w_file_pops = read_admix_pops_file(admix_args.w_file)

        # Extend the pops list
        admix_args.w_pops.extend(w_file_pops)

    # Check if the X populations file has been assigned
    if admix_args.x_file:

        # Save the pops from the file
        x_file_pops = read_admix_pops_file(admix_args.x_file)

        # Extend the pops list
        admix_args.x_pops.extend(x_file_pops)

    # Check if the Y populations file has been assigned
    if admix_args.y_file:

        # Save the pops from the file
        y_file_pops = read_admix_pops_file(admix_args.y_file)

        # Extend the pops list
        admix_args.y_pops.extend(y_file_pops)

    # Check if the Z populations file has been assigned
    if admix_args.z_file:

        # Save the pops from the file
        z_file_pops = read_admix_pops_file(admix_args.z_file)

        # Extend the pops list
        admix_args.z_pops.extend(z_file_pops)

    # Check if the A populations file has been assigned
    if admix_args.a_file:

        # Save the pops from the file
        a_file_pops = read_admix_pops_file(admix_args.a_file)

        # Extend the pops list
        admix_args.a_pops.extend(a_file_pops)

    # Check if the B populations file has been assigned
    if admix_args.b_file:

        # Save the pops from the file
        b_file_pops = read_admix_pops_file(admix_args.b_file)

        # Extend the pops list
        admix_args.b_pops.extend(b_file_pops)

    # Check if the C populations file has been assigned
    if admix_args.c_file:

        # Save the pops from the file
        c_file_pops = read_admix_pops_file(admix_args.c_file)

        # Extend the pops list
        admix_args.c_pops.extend(c_file_pops)

    # Check if the O populations file has been assigned
    if admix_args.o_file:

        # Save the pops from the file
        o_file_pops = read_admix_pops_file(admix_args.o_file)

        # Extend the pops list
        admix_args.o_pops.extend(o_file_pops)

    # Check if W pop(s) were assigned
    if admix_args.w_pops:

        # Assign pops not found within the model
        pops_not_found = pops_not_in_model(selected_model, admix_args.w_pops)

        # Check if any pops were not found
        if pops_not_found:
            raise Exception('W population(s) not found in model: %s' % ', '.join(pops_not_found))

    # Check if X pop(s) were assigned
    if admix_args.x_pops:

        # Assign pops not found within the model
        pops_not_found = pops_not_in_model(selected_model, admix_args.x_pops)

        # Check if any pops were not found
        if pops_not_found:
            raise Exception('X population(s) not found in model: %s' % ', '.join(pops_not_found))

    # Check if Y pop(s) were assigned
    if admix_args.y_pops:

        # Assign pops not found within the model
        pops_not_found = pops_not_in_model(selected_model, admix_args.y_pops)

        # Check if any pops were not found
        if pops_not_found:
            raise Exception('Y population(s) not found in model: %s' % ', '.join(pops_not_found))

    # Check if Z pop(s) were assigned
    if admix_args.z_pops:

        # Assign pops not found within the model
        pops_not_found = pops_not_in_model(selected_model, admix_args.z_pops)

        # Check if any pops were not found
        if pops_not_found:
            raise Exception('Z population(s) not found in model: %s' % ', '.join(pops_not_found))

    # Check if A pop(s) were assigned
    if admix_args.a_pops:

        # Assign pops not found within the model
        pops_not_found = pops_not_in_model(selected_model, admix_args.a_pops)

        # Check if any pops were not found
        if pops_not_found:
            raise Exception('A population(s) not found in model: %s' % ', '.join(pops_not_found))

    # Check if B pop(s) were assigned
    if admix_args.b_pops:

        # Assign pops not found within the model
        pops_not_found = pops_not_in_model(selected_model, admix_args.b_pops)

        # Check if any pops were not found
        if pops_not_found:
            raise Exception('B population(s) not found in model: %s' % ', '.join(pops_not_found))

    # Check if C pop(s) were assigned
    if admix_args.c_pops:

        # Assign pops not found within the model
        pops_not_found = pops_not_in_model(selected_model, admix_args.c_pops)

        # Check if any pops were not found
        if pops_not_found:
            raise Exception('C population(s) not found in model: %s' % ', '.join(pops_not_found))

    # Check if O pop(s) were assigned
    if admix_args.o_pops:

        # Assign pops not found within the model
        pops_not_found = pops_not_in_model(selected_model, admix_args.o_pops)

        # Check if any pops were not found
        if pops_not_found:
            raise Exception('O population(s) not found in model: %s' % ', '.join(pops_not_found))

    # Assign ind filename, update with better method
    admix_args.ind_filename = admix_args.eigenstrat_prefix + '.ind'

    # Assign a ind filename with the model
    model_ind_filename = assign_model_ind_filename(admix_args.ind_filename, selected_model, overwrite = admix_args.overwrite)
    
    # Create a population filename
    create_ind_w_pops(admix_args.ind_filename, selected_model, model_ind_filename)

    # String to hold the expected statistic output
    expected_statistic_output = ''

    # Check if the D statistic was specified
    if admix_args.calc_admix_statistic == 'D':

        # Check that the all the requried pops were specified
        if not admix_args.w_pops or not admix_args.x_pops or not admix_args.y_pops or not admix_args.z_pops:
            raise Exception('D requires the the following populations to be defined: W, X, Y, and Z')

        # Assign the expected statistic output 
        expected_statistic_output = admix_args.out_prefix + '.d'

        # Run the D Statistic
        r_admixr_d(admix_args.eigenstrat_prefix, admix_args.out_prefix, admix_args.w_pops, admix_args.x_pops, admix_args.y_pops, admix_args.z_pops)

    # Check if the F4 statistic was specified
    elif admix_args.calc_admix_statistic == 'F4':

        # Check that the all the requried pops were specified
        if not admix_args.w_pops or not admix_args.x_pops or not admix_args.y_pops or not admix_args.z_pops:
            raise Exception('F4 requires the the following populations to be defined: W, X, Y, and Z')

        # Assign the expected statistic output
        expected_statistic_output = admix_args.out_prefix + '.f4'
        
        # Run the F4 Statistic
        r_admixr_f4(admix_args.eigenstrat_prefix, admix_args.out_prefix, admix_args.w_pops, admix_args.x_pops, admix_args.y_pops, admix_args.z_pops)

    elif admix_args.calc_admix_statistic == 'F4-ratio':

        # Check that the all the requried pops were specified
        if not admix_args.a_pops or not admix_args.b_pops or not admix_args.c_pops or not admix_args.x_pops or not admix_args.o_pops:
            raise Exception('F4-ratio requires the the following populations to be defined: A, B, C, X, and O')

        # Assign the expected statistic output
        expected_statistic_output = admix_args.out_prefix + '.f4ratio'

        # Run the F4ratio Statistic
        r_admixr_f4ratio(admix_args.eigenstrat_prefix, admix_args.out_prefix, admix_args.a_pops, admix_args.b_pops, admix_args.c_pops, admix_args.x_pops, admix_args.o_pops)

    elif admix_args.calc_admix_statistic == 'F3':

        # Check that the all the requried pops were specified
        if not admix_args.a_pops or not admix_args.b_pops or not admix_args.c_pops:
            raise Exception('F3 requires the the following populations to be defined: A, B, and C')

        # Assign the expected statistic output
        expected_statistic_output = admix_args.out_prefix + '.f3'

        # Run the F3 Statistic
        r_admixr_f3(admix_args.eigenstrat_prefix, admix_args.out_prefix, admix_args.a_pops, admix_args.b_pops, admix_args.c_pops)

    # Remove the model ind file
    os.remove(model_ind_filename)

if __name__ == "__main__":
    initLogger()
    run()
