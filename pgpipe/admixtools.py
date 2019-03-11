import os
import sys
import argparse
import logging
import subprocess
import random
import string
import shutil

# Import PPP modules and scripts
from pgpipe.plink import *
from pgpipe.model import read_model_file
from pgpipe.logging_module import initLogger, logArgs

def check_convertf_for_errors (convertf_stderr):
    '''
        Checks the convertf stderr for errors

        Parameters
        ----------
        convertf_stderr : str
            convertf stderr

        Raises
        ------
        IOError
            If convertf stderr returns an error
    '''
    
    # Check if IDs are too long
    if 'too long' in convertf_stderr:
        raise Exception('Individual and/or family names too long. Eigenstrat supports a maximum of 39 characters')

    # Expand as needed
    elif convertf_stderr:
        raise Exception(convertf_stderr)

def call_convertf (convertf_call_args):
    '''
        Calls convertf

        The function calls convertf. Returns the stderr of convertf to
        create log file of the call.

        Parameters
        ----------
        convertf_call_args : list
            convertf arguments

        Returns
        -------
        convertf_err : str
            convertf log output
    '''

    # convertf subprocess call without stdout
    convertf_call = subprocess.Popen(['convertf'] + list(map(str, convertf_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Wait for convertf to finish
    convertf_stdout, convertf_stderr = convertf_call.communicate()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        convertf_stderr = convertf_stderr.decode()

    logging.info('convertf call complete')

    # Check that the log file was created correctly
    check_convertf_for_errors(convertf_stderr)

    return convertf_stderr

def assign_unique_filename (out_prefix, out_suffix, str_size = 10):

    # Generate a random string for par filename
    random_str = ''.join(random.choice(string.ascii_uppercase + string.digits) for digit in range(str_size))

    # Save the unique filename
    unique_filename = '%s.%s.%s' % (out_prefix, random_str, out_suffix)

    # Loop until filename is unique
    while os.path.isfile(unique_filename):

        # Generate a random string for par filename
        random_str = ''.join(random.choice(string.ascii_uppercase + string.digits) for digit in range(str_size))

        # Save the unique filename
        unique_filename = '%s.%s.%s' % (out_prefix, random_str, out_suffix)

    return unique_filename

def par_assign_input (gt_filename, snp_filename, ind_filename):
    
    # Assign the genotype filename
    par_input_str =  'genotypename:    %s\n' % gt_filename

    # Assign the snpname filename
    par_input_str += 'snpname:         %s\n' % snp_filename

    # Assign the indivname filename
    par_input_str += 'indivname:       %s\n' % ind_filename

    return par_input_str

def par_assign_eigenstrat (out_prefix):

    # Confirm an output prefix has been assigned
    if not out_prefix:
        raise Exception('No output prefix assigned. Unable to complete conversion')

    # Assign the output filenames
    par_eigenstrat_str =  'outputformat:    EIGENSTRAT\n'
    par_eigenstrat_str += 'genotypeoutname: %s.geno\n' % out_prefix
    par_eigenstrat_str += 'snpoutname:      %s.snp\n' % out_prefix
    par_eigenstrat_str += 'indivoutname:    %s.ind\n' % out_prefix

    return par_eigenstrat_str

def par_assign_family_ids (family_id):

    # Check if the user has requested family IDs
    if family_id:

        # Add the family ID option
        family_id_str =  'familynames:     YES'
        
    else:

        # Add the family ID option
        family_id_str = 'familynames:     NO'

    return family_id_str

def ped_to_eigenstrat (ped_filename = None, map_filename = None, ped_prefix = None, family_id = False, out_prefix = None, overwrite = False, keep_original = True):
    
    # Check if the input is a specified by a ped prefix, and check that the files exist
    if ped_prefix and assign_ped_from_prefix(ped_prefix):

        # Assign the files
        ped_filename = ped_prefix + '.ped'
        map_filename = ped_prefix + '.map'

    # Check that all the input files exist
    elif assign_ped_from_files(ped_filename, map_filename):
        pass

    # Assign the par filename
    par_filename = assign_unique_filename(out_prefix, 'par')

    # Assign the input files
    par_file_str = par_assign_input(ped_filename, map_filename, ped_filename)

    # Add the Eigenstrat assignment lines
    par_file_str += par_assign_eigenstrat(out_prefix)

    # Add the family ID line
    par_file_str += par_assign_family_ids(family_id)

    # Create the par file
    par_file = open(par_filename, 'w')
    par_file.write(par_file_str)
    par_file.close()

    # Call convertf
    call_convertf(['-p', par_filename])

    # Remove the par file
    os.remove(par_filename)

def bed_to_eigenstrat (bed_filename = None, bim_filename = None, fam_filename = None, bed_prefix = None, family_id = False, out_prefix = None, overwrite = False, keep_original = True):

    logging.info('Beginning eigenstrat conversion')

    # Check if the input is a specified by a bed prefix, and check that the files exist
    if bed_prefix and assign_bed_from_prefix(bed_prefix):

        # Assign the files
        bed_filename = bed_prefix + '.bed'
        bim_filename = bed_prefix + '.bim'
        fam_filename = bed_prefix + '.fam'

    # Check that all the input files exist
    elif assign_bed_from_files(bed_filename, bim_filename, fam_filename):
        pass

    # Rename the fam file
    shutil.copy(fam_filename, fam_filename + '.pedind')

    # Update the fam filename
    fam_filename += '.pedind'

    # Assign the par filename
    par_filename = assign_unique_filename(out_prefix, 'par')

    # Assign the input files
    par_file_str = par_assign_input(bed_filename, bim_filename, fam_filename)

    # Add the Eigenstrat assignment lines
    par_file_str += par_assign_eigenstrat(out_prefix)

    # Add the family ID line
    par_file_str += par_assign_family_ids(family_id)

    # Create the par file
    par_file = open(par_filename, 'w')
    par_file.write(par_file_str)
    par_file.close()

    # Call convertf
    call_convertf(['-p', par_filename])

    # Remove the par file
    os.remove(par_filename)

    # Remove the copy
    os.remove(fam_filename)

    # Update the fam filename
    fam_filename = fam_filename[:-7]

def update_eigenstrat_pops (ind_filename, model):

    # Create a filename for the temporary ind file
    tmp_ind_filename = assign_unique_filename(ind_filename, 'tmp')

    # Create a temporary file to replace the ind file
    tmp_ind_file = open(tmp_ind_filename, 'w')

    # Open the ind file
    with open(ind_filename, 'r') as ind_file:

        # Loop the ind file
        for ind_line in ind_file:

            # Split the ind line
            ind_data = ind_line.strip().split()

            # Assign the current ind
            ind_name = ind_data[0]

            # Assign the pop from the ind
            pop_name = model.return_pop(ind_name)

            # Update the ind data
            ind_data[2] = pop_name
            
            # Write the updated line to the temporary file
            tmp_ind_file.write(' '.join(ind_data) + '\n')

    tmp_ind_file.close()

    # Rename the ind file
    shutil.move(tmp_ind_filename, ind_filename)

def r_admixr (list1, list2):
    import rpy2.robjects.numpy2ri
    from rpy2.robjects.packages import importr
    import rpy2.robjects as robjects
    import numpy as np

    rpy2.robjects.numpy2ri.activate()
    rc = robjects.r['c']

    r_grdevices = importr('grDevices')
    r_cor = robjects.r['cor.test']
    r_summary = robjects.r.summary

    cor_output = r_cor(np.array(list1), np.array(list2), method = "pearson", alternative = "two.sided")
    return cor_output[cor_output.names.index('estimate')][0], cor_output[cor_output.names.index('p.value')][0]

def admixtools_parser(passed_arguments):
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
                if not getattr(args, self.dest):
                    setattr(args, self.dest, value)
                else:
                    getattr(args, self.dest).extend(value)
        return customAction

    admixtools_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Input PED arguments
    admixtools_parser.add_argument("--ped", dest = 'ped_filename', help = "Defines the filename of the PED file. Called alongside --map", type = str, action = parser_confirm_file())
    admixtools_parser.add_argument("--map", dest = 'map_filename', help = "Defines the filename of the MAP file. Called alongside --ped", type = str, action = parser_confirm_file())

    # Input BED arguments
    admixtools_parser.add_argument("--binary-ped", dest = 'bed_filename', help = "Defines the filename of the Binary-PED (i.e. BED) file. Called alongside --fam and --bim", type = str, action = parser_confirm_file())
    admixtools_parser.add_argument("--fam", dest = 'fam_filename', help = "Defines the filename of the FAM file. Called alongside --binary-ped and --bim", type = str, action = parser_confirm_file())
    admixtools_parser.add_argument("--bim", dest = 'bim_filename', help = "Defines the filename of the BIM file. Called alongside --binary-ped and --fam", type = str, action = parser_confirm_file())

    # Input PED/BED prefix arguments
    admixtools_prefix = admixtools_parser.add_mutually_exclusive_group()
    admixtools_prefix.add_argument("--ped-prefix", help = "Defines the filename prefix of both PED and MAP files", type = str)
    admixtools_prefix.add_argument("--binary-ped-prefix", dest = 'bed_prefix', help = "Defines the filename prefix of the Binary-PED, FAM, and BIM files", type = str)

    # Model file arguments.
    admixtools_parser.add_argument('--model-file', help = 'Defines the model file', type = str, action = parser_confirm_file())
    admixtools_parser.add_argument('--model', help = 'Defines the model and the individual(s) to include', type = str)

    admixtools_prefix.add_argument('--out-prefix', help = 'Defines the output prefix (i.e. filename without file extension)', default = 'out')

    # Admix Pops
    admixtools_prefix.add_argument('--admix-w-pop', help = 'W population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admixtools_prefix.add_argument('--admix-w-pops', help = 'W populations for admixure analysis', type = str, action = parser_confirm_file())
    
    admixtools_prefix.add_argument('--admix-x-pop', help = 'X population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admixtools_prefix.add_argument('--admix-x-pops', help = 'X populations for admixure analysis', type = str, action = parser_confirm_file())

    admixtools_prefix.add_argument('--admix-y-pop', help = 'Y population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admixtools_prefix.add_argument('--admix-y-pops', help = 'Y populations for admixure analysis', type = str, action = parser_confirm_file())

    admixtools_prefix.add_argument('--admix-z-pop', help = 'Z population(s) for admixure analysiss', nargs = '+', type = str, action = parser_add_to_list())
    admixtools_prefix.add_argument('--admix-z-pops', help = 'Z populations for admixure analysiss', type = str, action = parser_confirm_file())

    admixtools_prefix.add_argument('--admix-a-pop', help = 'A population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admixtools_prefix.add_argument('--admix-a-pops', help = 'A populations for admixure analysis', type = str, action = parser_confirm_file())

    admixtools_prefix.add_argument('--admix-b-pop', help = 'B population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admixtools_prefix.add_argument('--admix-b-pops', help = 'B populations for admixure analysis', type = str, action = parser_confirm_file())

    admixtools_prefix.add_argument('--admix-c-pop', help = 'C population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admixtools_prefix.add_argument('--admix-c-pops', help = 'C populations for admixure analysis', type = str, action = parser_confirm_file())

    admixtools_prefix.add_argument('--admix-o-pop', help = 'Outgroup population(s) for admixure analysis', nargs = '+', type = str, action = parser_add_to_list())
    admixtools_prefix.add_argument('--admix-o-pops', help = 'Outgroup populations for admixure analysis', type = str, action = parser_confirm_file())

    if passed_arguments:
        return admixtools_parser.parse_args(passed_arguments)
    else:
        return admixtools_parser.parse_args()

def run(passed_arguments = []):

    # Grab admixtools arguments from command line
    admixtools_args = admixtools_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(admixtools_args, func_name='admixtools')

    # Read in the models
    models_in_file = read_model_file(admixtools_args.model_file)

    # Check that the selected model was not found in the file
    if admixtools_args.model not in models_in_file:
        raise IOError('Selected model "%s" not found in: %s' % (admixtools_args.model, admixtools_args.model_file))

    # Select model, might change this in future versions
    selected_model = models_in_file[admixtools_args.model]

    # Check if ped-based files were assigned
    if admixtools_args.ped_prefix or (admixtools_args.ped_filename and admixtools_args.map_filename):

        ped_to_eigenstrat(ped_filename = admixtools_args.ped_filename, map_filename = admixtools_args.map_filename, ped_prefix = admixtools_args.ped_prefix, out_prefix = admixtools_args.out_prefix)

    elif admixtools_args.bed_prefix or (admixtools_args.bed_filename and admixtools_args.bim_filename and admixtools_args.fam_filename):

        bed_to_eigenstrat(bed_filename = admixtools_args.bed_filename, bim_filename = admixtools_args.bim_filename, fam_filename = admixtools_args.fam_filename, bed_prefix = admixtools_args.bed_prefix, out_prefix = admixtools_args.out_prefix)
        update_eigenstrat_pops(admixtools_args.out_prefix + '.ind', selected_model)
    

if __name__ == "__main__":
    initLogger()
    run()
