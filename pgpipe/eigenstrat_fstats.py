#!/usr/bin/env python
'''
    Automates the calculation of multiple admixture statistics, including: Patterson's D, 
    F4 statistic, F4-ratio statistic, and F3 statistic.

    ##################
    Command-line Usage
    ##################
    The admixture statistics automater may be called using the following command:

    .. code-block:: bash
        
        eigenstrat_fstats.py

    *************
    Example usage
    *************
    Command-line to calculate Patterson's D:

    .. code-block:: bash
        
        eigenstrat_fstats.py --eigenstrat-prefix snps --calc-admix-statistic D --admix-w-pop French --admix-x-pop Yoruba --admix-y-pop Vindija --admix-z-pop Chimp 
    
    Command-line to calculate the F4-ratio: 

    .. code-block:: bash
        
        eigenstrat_fstats.py --eigenstrat-prefix snps --calc-admix-statistic F4-ratio --admix-a-pop Altai --admix-b-pop Vindija --admix-c-pop Yoruba --admix-x-pop French --admix-o-pop Chimp  
    
    ############
    Dependencies 
    ############
    * `AdmixTools <https://github.com/DReichLab/AdmixTools>`_
    * `admixr <https://github.com/bodkan/admixr>`_


    ############################
    Input Command-line Arguments
    ############################
    **--eigenstrat-prefix** *<input_prefix>*
        Argument used to define the filename prefix shared by the genotype file (.geno), 
        the individual file (.ind), and the SNP file (.snp). Should not be used alongside
        the specific file arguments (e.g. --geno).
    **--geno** *<geno_filename>*
        Argument used to define the filename of the eigenstrat genotype file (.geno). 
        Must be called alongside --ind and --snp. Cannot be called alongside 
        --eigenstrat-prefix.
    **--ind** *<ind_filename>*
        Argument used to define the filename of the eigenstrat individual file (.ind). 
        Must be called alongside --geno and --snp. Cannot be called alongside 
        --eigenstrat-prefix.
    **--snp** *<snp_filename>*
        Argument used to define the filename of the eigenstrat SNP file (.snp). Must be
        called alongside --geno and --ind. Cannot be called alongside --eigenstrat-prefix.
    **--model-file** *<model_filename>*
        Argument used to define the model file. Please note that this argument cannot be 
        used with the individual-based filters.
    **--model** *<model_str>*
        Argument used to define the model (i.e. the individual(s) to include and/or the 
        populations for relevant statistics). May be used with any statistic. Please note 
        that this argument cannot be used with **--pop-file** argument or the 
        individual-based filters.
    
    #############################
    Output Command-line Arguments
    #############################
    **--out** *<output_filename>*
        Argument used to define the complete output filename, overrides **--out-prefix**.
        Cannot be used if multiple output files are created.
    **--out-prefix** *<output_prefix>*
        Argument used to define the output prefix (i.e. filename without file extension)
    **--overwrite**
        Argument used to define if previous output should be overwritten.
    
    ####################################
    Statistic Command-line Specification
    ####################################
    **--calc-admix-statistic** *<D, F4, F4-ratio, F3>*
        Argument used to define the admix statistic to be calculated. Patterson's D (D), 
        F4 statistic (F4), F4-ratio statistic (F4-ratio), and F3 statistic (F3). See below
        for details on the arguments requried by each statistic .

    ***********************************
    Statistic Command-line Requirements
    ***********************************
    It should be noted that each admix statistic has a specific set of population labels
    arguments. These labels are used to specify a representive population. For instance, 
    the argument '--admix-w-pop CEU' will replace the W label of Patterson's D and the 
    F4 statistic with the CEU population. These arguments may be found in the next section.
 
    **--calc-admix-statistic** *D*
        Requires: **--admix-w-pop**/**--admix-w-pop-file**, 
        **--admix-x-pop**/**--admix-x-pop-file**,
        **--admix-y-pop**/**--admix-y-pop-file**, 
        and **--admix-z-pop**/**--admix-z-pop-file**.

    **--calc-admix-statistic** *F4*
        Requires: **--admix-w-pop**/**--admix-w-pop-file**, 
        **--admix-x-pop**/**--admix-x-pop-file**,
        **--admix-y-pop**/**--admix-y-pop-file**, 
        and **--admix-z-pop**/**--admix-z-pop-file**.

    **--calc-admix-statistic** *F4-ratio*
        Requires: **--admix-a-pop**/**--admix-a-pop-file**,
        **--admix-b-pop**/**--admix-b-pop-file**,
        **--admix-c-pop**/**--admix-c-pop-file**, 
        **--admix-x-pop**/**--admix-x-pop-file**,
        and **--admix-o-pop**/**--admix-o-pop-file**.
    **--calc-admix-statistic** *F3*
        Requires: **--admix-a-pop**/**--admix-a-pop-file**,
        **--admix-b-pop**/**--admix-b-pop-file**,
        and **--admix-c-pop**/**--admix-c-pop-file**.

    *******************************************
    Additional Statistic Command-line Arguments
    *******************************************
    **--admix-w-pop** *<w_pop_str>* *<w_pop1_str, w_pop2_str, etc.>*
        Argument used to define the population(s) to represent W in the supported admixure 
        statistic. This argument may be used multiple times if desired. If multiple 
        populations the statistic will be repeated until each population has represented W.
    **--admix-w-pop-file** *<w_pop_filename>*
        Argument used to define a file of population(s) to represent W in the supported 
        admixure statistic. If multiple populations the statistic will be repeated until each 
        population has represented W.
    **--admix-x-pop** *<x_pop_str>* *<x_pop1_str, x_pop2_str, etc.>*
        Argument used to define the population(s) to represent X in the supported admixure 
        statistic. This argument may be used multiple times if desired. If multiple 
        populations the statistic will be repeated until each population has represented X.
    **--admix-x-pop-file** *<x_pop_filename>*
        Argument used to define a file of population(s) to represent X in the supported 
        admixure statistic. If multiple populations the statistic will be repeated until each 
        population has represented X.
    **--admix-y-pop** *<y_pop_str>* *<y_pop1_str, y_pop2_str, etc.>*
        Argument used to define the population(s) to represent Y in the supported admixure 
        statistic. This argument may be used multiple times if desired. If multiple 
        populations the statistic will be repeated until each population has represented Y.
    **--admix-y-pop-file** *<y_pop_filename>*
        Argument used to define a file of population(s) to represent Y in the supported 
        admixure statistic. If multiple populations the statistic will be repeated until each 
        population has represented Y.
    **--admix-z-pop** *<z_pop_str>* *<z_pop1_str, z_pop2_str, etc.>*
        Argument used to define the population(s) to represent Z in the supported admixure 
        statistic. This argument may be used multiple times if desired. If multiple 
        populations the statistic will be repeated until each population has represented Z.
    **--admix-z-pop-file** *<z_pop_filename>*
        Argument used to define a file of population(s) to represent Z in the supported 
        admixure statistic. If multiple populations the statistic will be repeated until each 
        population has represented Z.
    **--admix-a-pop** *<a_pop_str>* *<a_pop1_str, a_pop2_str, etc.>*
        Argument used to define the population(s) to represent A in the supported admixure 
        statistic. This argument may be used multiple times if desired. If multiple 
        populations the statistic will be repeated until each population has represented A.
    **--admix-a-pop-file** *<a_pop_filename>*
        Argument used to define a file of population(s) to represent A in the supported 
        admixure statistic. If multiple populations the statistic will be repeated until each 
        population has represented A.
    **--admix-b-pop** *<b_pop_str>* *<b_pop1_str, b_pop2_str, etc.>*
        Argument used to define the population(s) to represent B in the supported admixure 
        statistic. This argument may be used multiple times if desired. If multiple 
        populations the statistic will be repeated until each population has represented B.
    **--admix-b-pop-file** *<b_pop_filename>*
        Argument used to define a file of population(s) to represent B in the supported 
        admixure statistic. If multiple populations the statistic will be repeated until each 
        population has represented B.
    **--admix-c-pop** *<c_pop_str>* *<c_pop1_str, c_pop2_str, etc.>*
        Argument used to define the population(s) to represent C in the supported admixure 
        statistic. This argument may be used multiple times if desired. If multiple 
        populations the statistic will be repeated until each population has represented C.
    **--admix-c-pop-file** *<c_pop_filename>*
        Argument used to define a file of population(s) to represent C in the supported 
        admixure statistic. If multiple populations the statistic will be repeated until each 
        population has represented C.
'''

import os
import sys
import argparse
import logging
import subprocess
import random
import string
import shutil

# Import PPP modules and scripts
from pgpipe.eigenstrat_wrapper import *
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

    # Input Eigenstrat prefix arguments
    admix_parser.add_argument("--eigenstrat-prefix", help = "Defines the filename prefix of the GENO, IND, and SNP files", type = str)

    # Input Eigenstrat file arguments
    admix_parser.add_argument("--geno", dest = 'geno_filename', help = "Defines the filename of the GENO file. Called alongside --ind and --snp", type = str, action = parser_confirm_file())
    admix_parser.add_argument("--ind", dest = 'ind_filename', help = "Defines the filename of the IND file. Called alongside --geno and --snp", type = str, action = parser_confirm_file())
    admix_parser.add_argument("--snp", dest = 'snp_filename', help = "Defines the filename of the SNP file. Called alongside --geno and --ind", type = str, action = parser_confirm_file())

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
    admix_parser.add_argument('--admix-w-pop-file', dest = 'w_file', help = 'W populations for admixure analysis', type = str, action = parser_confirm_file())
    
    admix_parser.add_argument('--admix-x-pop', dest = 'x_pops', help = 'X population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admix_parser.add_argument('--admix-x-pop-file', dest = 'x_file', help = 'X populations for admixure analysis', type = str, action = parser_confirm_file())

    admix_parser.add_argument('--admix-y-pop', dest = 'y_pops', help = 'Y population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admix_parser.add_argument('--admix-y-pop-file', dest = 'y_file', help = 'Y populations for admixure analysis', type = str, action = parser_confirm_file())

    admix_parser.add_argument('--admix-z-pop', dest = 'z_pops', help = 'Z population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admix_parser.add_argument('--admix-z-pop-file', dest = 'z_file', help = 'Z populations for admixure analysiss', type = str, action = parser_confirm_file())

    admix_parser.add_argument('--admix-a-pop', dest = 'a_pops', help = 'A population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admix_parser.add_argument('--admix-a-pop-file', dest = 'a_file', help = 'A populations for admixure analysis', type = str, action = parser_confirm_file())

    admix_parser.add_argument('--admix-b-pop', dest = 'b_pops', help = 'B population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admix_parser.add_argument('--admix-b-pop-file', dest = 'b_file', help = 'B populations for admixure analysis', type = str, action = parser_confirm_file())

    admix_parser.add_argument('--admix-c-pop', dest = 'c_pops', help = 'C population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admix_parser.add_argument('--admix-c-pop-file', dest = 'c_file', help = 'C populations for admixure analysis', type = str, action = parser_confirm_file())

    admix_parser.add_argument('--admix-o-pop', dest = 'o_pops', help = 'O population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admix_parser.add_argument('--admix-o-pop-file', dest = 'o_file', help = 'O populations for admixure analysis', type = str, action = parser_confirm_file())

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

    # Bool to indicate if file names were changed
    eigenstrat_filenames_changed = False

    # Check if an eigenstrat prefix was assigned
    if admix_args.eigenstrat_prefix:

        # Confirm all eigenstrat files from prefix, will raise exception upon failure
        confirm_eigenstrat_files_from_prefix(admix_args.eigenstrat_prefix)

    # Check if the eigenstrat files were assigned separately
    elif admix_args.geno_filename or admix_args.ind_filename or admix_args.snp_filename:

        # Confirm that all eigenstrat files exist, will raise exception upon failure
        confirm_eigenstrat_files(admix_args.geno_filename, admix_args.ind_filename, admix_args.snp_filename)

        # Confirm that an eigenstrat prefix may be assigned
        if check_eigenstrat_prefix(admix_args.geno_filename, admix_args.ind_filename, admix_args.snp_filename):

            # Assign the prefix
            admix_args.eigenstrat_prefix = os.path.splitext(admix_args.geno_filename)[0]

        else:

            # Confirm the eigenstrat filenames were changed
            eigenstrat_filenames_changed = True

            # Rename files and assign the prefix
            admix_args.eigenstrat_prefix = assign_eigenstrat_prefix(admix_args.geno_filename, admix_args.ind_filename, admix_args.snp_filename)

    # Assign a ind filename with the model
    model_ind_filename = assign_model_ind_filename(admix_args.eigenstrat_prefix, selected_model, overwrite = admix_args.overwrite)
    
    # Create a population filename
    create_ind_w_pops(admix_args.eigenstrat_prefix, selected_model, model_ind_filename)

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

    # Check if the eigenstrat filenames were changed
    if eigenstrat_filenames_changed:

        # Restore the filename if they were changed
        restore_eigenstrat_files(admix_args.geno_filename, admix_args.ind_filename, admix_args.snp_filename)

if __name__ == "__main__":
    initLogger()
    run()
