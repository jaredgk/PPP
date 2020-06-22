import os
import sys
import argparse
import logging
import subprocess
import random
import string
import shutil

import rpy2.robjects as robjects

from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

# Import PPP modules and scripts
from pgpipe.plink import *
from pgpipe.model import read_model_file, pops_not_in_model
from pgpipe.logging_module import initLogger, logArgs
from pgpipe.misc import confirm_executable

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

    # Confirm where the specifed executable is located
    convertf_path = confirm_executable('convertf')
    
    # Check if executable is installed
    if not convertf_path:
        raise Exception('convertf not found. Please confirm the executable is installed')

    # convertf subprocess call without stdout
    convertf_call = subprocess.Popen([convertf_path] + list(map(str, convertf_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

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

def ped_to_eigenstrat (ped_filename = None, map_filename = None, ped_prefix = None, out_prefix = None, family_id = False, overwrite = False, **kwargs):

    # Check if ped-based files were assigned
    if ped_prefix and confirm_ped_prefix(ped_prefix):

        # Assign the files
        ped_filename = ped_prefix + '.ped'
        map_filename = ped_prefix + '.map'

    # Check that all the ped files are assigned
    elif ped_filename and confirm_ped_files(ped_filename, map_filename):

        # No action requried, just needs to pass
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

def bed_to_eigenstrat (bed_filename = None, bim_filename = None, fam_filename = None, bed_prefix = None, out_prefix = None, family_id = False, overwrite = False, **kwargs):

    # Check if the ped prefix files exist
    if bed_prefix and confirm_bed_prefix(bed_prefix):

        # Assign the files
        bed_filename = bed_prefix + '.bed'
        bim_filename = bed_prefix + '.bim'
        fam_filename = bed_prefix + '.fam'

    # Check that all the bed files are assigned
    elif bed_filename and confirm_bed_files(bed_filename, bim_filename, fam_filename):
        
        # No action requried, just needs to pass
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

def confirm_eigenstrat_files_from_prefix (eigenstrat_prefix):

    # Assign expected filenames
    geno_filename = eigenstrat_prefix + '.geno'
    ind_filename = eigenstrat_prefix + '.ind'
    snp_filename = eigenstrat_prefix + '.snp'

    # Check if a genotype file is assigned
    if not geno_filename:
        raise IOError('Unable to assign genotype file (i.e. --geno). Please confirm the file is named correctly')

    # Check if the ind and snp files are not assigned
    if not ind_filename and not snp_filename:
        raise IOError('Unable to assign the ind and snp files. Please confirm the files (i.e. --ind, --snp) are assigned')

    # Check if the ind is not assigned
    if not ind_filename:
        raise IOError('Unable to assign ind file. Please confirm the ind file (i.e. --ind) is assigned')

    # Check if the snp file is not assigned
    if not snp_filename:
        raise IOError('Unable to assign snp file. Please confirm the snp file (i.e. --snp) is assigned')

def confirm_eigenstrat_files (geno_filename, ind_filename, snp_filename):

    # Check if a genotype file is assigned
    if not geno_filename:
        raise IOError('Unable to assign genotype file (i.e. --geno). Please confirm the file is named correctly')

    # Check if the ind and snp files are not assigned
    if not ind_filename and not snp_filename:
        raise IOError('Unable to assign the ind and snp files. Please confirm the files (i.e. --ind, --snp) are assigned')

    # Check if the ind is not assigned
    if not ind_filename:
        raise IOError('Unable to assign ind file. Please confirm the ind file (i.e. --ind) is assigned')

    # Check if the snp file is not assigned
    if not snp_filename:
        raise IOError('Unable to assign snp file. Please confirm the snp file (i.e. --snp) is assigned')

def check_eigenstrat_prefix (geno_filename, ind_filename, snp_filename):

    # Assign the base name of each file
    geno_filename_wo_ext = os.path.splitext(geno_filename)[0]
    ind_basename_wo_ext = os.path.splitext(ind_filename)[0]
    snp_basename_wo_ext = os.path.splitext(snp_filename)[0]

    # Check if a prefix may be assigned without any changes
    if geno_filename_wo_ext == ind_basename_wo_ext == snp_basename_wo_ext:
        
        # Return True if a prefix may be assigned
        return True

    return False

def assign_eigenstrat_prefix (geno_filename, ind_filename, snp_filename):

    # Assign the base name of each file
    geno_filename_wo_ext = os.path.splitext(geno_filename)[0]

    # Check that the ind file needs to be changed
    if ind_filename != geno_filename_wo_ext + '.ind':

        # Rename the file
        shutil.move(ind_filename, geno_filename_wo_ext + '.ind')

    # Check that the ind file needs to be changed
    if snp_filename != geno_filename_wo_ext + '.snp':

        # Rename the file
        shutil.move(snp_filename, geno_filename_wo_ext + '.snp')

    # Return the prefix
    return geno_filename_wo_ext

def restore_eigenstrat_files (geno_filename, ind_filename, snp_filename):

    # Assign the base name of each file
    geno_filename_wo_ext = os.path.splitext(geno_filename)[0]

    # Check that the ind file needs to be changed
    if ind_filename != geno_filename_wo_ext + '.ind':

        # Rename the file
        shutil.move(geno_filename_wo_ext + '.ind', ind_filename)

    # Check that the ind file needs to be changed
    if snp_filename != geno_filename_wo_ext + '.snp':

        # Rename the file
        shutil.move(geno_filename_wo_ext + '.snp', snp_filename)

def assign_model_ind_filename (eigenstrat_prefix, model, overwrite = False):

    # Assign a model-based ind filename
    model_ind_filename = '%s.%s.ind' % (eigenstrat_prefix, model.name)

    # Check if the file should be overwritten
    if not overwrite:

        # Check if the output filename exists
        if os.path.isfile(model_ind_filename):
            raise Exception('Model-based individual file exists. Please use --overwrite to ignore.')

    return model_ind_filename

def create_ind_w_pops (eigenstrat_prefix, model, out_filename):

    # Create a temporary file to replace the ind file
    ind_pop_file = open(out_filename, 'w')

    # Open the ind file
    with open(eigenstrat_prefix + '.ind', 'r') as ind_file:

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
            ind_pop_file.write(' '.join(ind_data) + '\n')

    ind_pop_file.close()

def update_ind_w_pops (eigenstrat_prefix, model):

    # Create a filename for the temporary ind file
    tmp_ind_filename = assign_unique_filename(eigenstrat_prefix + '.ind', 'tmp')

    # Create a temporary file that will replace the ind file
    ind_pop_file = open(tmp_ind_filename, 'w')

    # Open the ind file
    with open(eigenstrat_prefix + '.ind', 'r') as ind_file:

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
            ind_pop_file.write(' '.join(ind_data) + '\n')

    ind_pop_file.close()

    # Rename the ind file
    shutil.move(tmp_ind_filename, eigenstrat_prefix + '.ind')

def r_admixr_d (eigenstrat_prefix, out_prefix, w_pops, x_pops, y_pops, z_pops):

    # Activate converions between pandas and R 
    pandas2ri.activate()

    # Activate the admixr library
    try:
        r_admixr_lib = importr('admixr')
    except:
        raise Exception('Please install the admixr R library')

    # Create the c() R object
    rc = robjects.r['c']

    # Create the eigenstrat file object from admixr
    r_eigenstrat = robjects.r.eigenstrat

    # Create the statistic objects from admixr
    r_d = robjects.r.d
    
    # Load the eigenstrat file into R
    r_eigenstrat_input = r_eigenstrat(eigenstrat_prefix)

    # Run the D statistic
    r_d_results = r_d(W = rc(w_pops), X = rc(x_pops), Y = rc(y_pops), Z = rc(z_pops), data = r_eigenstrat_input)

    # Convert the results into a Pandas dataframe
    d_results = pandas2ri.ri2py(r_d_results)

    # Change the datatype of the nsnps, ABBA, and BABA columns
    d_results[['nsnps', 'BABA', 'ABBA']] = d_results[['nsnps', 'BABA', 'ABBA']].astype(int)

    # Write the dataframe as a file
    d_results.to_csv(out_prefix + '.d', sep = '\t', index = False)

def r_admixr_f4 (eigenstrat_prefix, out_prefix, w_pops, x_pops, y_pops, z_pops):

    # Activate converions between pandas and R 
    pandas2ri.activate()

    # Activate the admixr library
    try:
        r_admixr_lib = importr('admixr')
    except:
        raise Exception('Please install the admixr R library')

    # Create the c() R object
    rc = robjects.r['c']

    # Create the eigenstrat file object from admixr
    r_eigenstrat = robjects.r.eigenstrat

    # Create the statistic objects from admixr
    r_f4 = robjects.r.f4
    
    # Load the eigenstrat file into R
    r_eigenstrat_input = r_eigenstrat(eigenstrat_prefix)

    # Run the F4 statistic
    r_f4_results = r_f4(W = rc(w_pops), X = rc(x_pops), Y = rc(y_pops), Z = rc(z_pops), data = r_eigenstrat_input)

    # Convert the results into a Pandas dataframe
    f4_results = pandas2ri.ri2py(r_f4_results)

    # Change the datatype of the nsnps, ABBA, and BABA columns
    f4_results[['nsnps', 'BABA', 'ABBA']] = f4_results[['nsnps', 'BABA', 'ABBA']].astype(int)

    # Write the dataframe as a file
    f4_results.to_csv(out_prefix + '.f4', sep = '\t', index = False)

def r_admixr_f4ratio (eigenstrat_prefix, out_prefix, a_pops, b_pops, c_pops, x_pops, o_pops):

    # Activate converions between pandas and R 
    pandas2ri.activate()

    # Activate the admixr library
    try:
        r_admixr_lib = importr('admixr')
    except:
        raise Exception('Please install the admixr R library')

    # Create the c() R object
    rc = robjects.r['c']

    # Create the eigenstrat file object from admixr
    r_eigenstrat = robjects.r.eigenstrat

    # Create the statistic objects from admixr
    r_f4ratio = robjects.r.f4ratio
    
    # Load the eigenstrat file into R
    r_eigenstrat_input = r_eigenstrat(eigenstrat_prefix)

    # Run the f4ratio statistic
    r_f4ratio_results = r_f4ratio(A = rc(a_pops), B = rc(b_pops), C = rc(c_pops), X = rc(x_pops), O = rc(o_pops), data = r_eigenstrat_input)

    # Convert the results into a Pandas dataframe
    f4ratio_results = pandas2ri.ri2py(r_f4ratio_results)

    # Write the dataframe as a file
    f4ratio_results.to_csv(out_prefix + '.f4ratio', sep = '\t', index = False)

def r_admixr_f3 (eigenstrat_prefix, out_prefix, a_pops, b_pops, c_pops):

    # Activate converions between pandas and R 
    pandas2ri.activate()

    # Activate the admixr library
    try:
        r_admixr_lib = importr('admixr')
    except:
        raise Exception('Please install the admixr R library')

    # Create the c() R object
    rc = robjects.r['c']

    # Create the eigenstrat file object from admixr
    r_eigenstrat = robjects.r.eigenstrat

    # Create the statistic objects from admixr
    r_f3 = robjects.r.f3
    
    # Load the eigenstrat file into R
    r_eigenstrat_input = r_eigenstrat(eigenstrat_prefix)

    # Run the f3 statistic
    r_f3_results = r_f3(A = rc(a_pops), B = rc(b_pops), C = rc(c_pops), data = r_eigenstrat_input)

    # Convert the results into a Pandas dataframe
    f3_results = pandas2ri.ri2py(r_f3_results)

    # Change the datatype of the nsnps column
    f3_results['nsnps'] = f3_results['nsnps'].astype(int)

    # Write the dataframe as a file
    f3_results.to_csv(out_prefix + '.f3', sep = '\t', index = False)