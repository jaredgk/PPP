import sys
import os
import filecmp
import gzip
import random 
import shutil
import string
import subprocess

from pgpipe.misc import confirm_executable

def fileComp (test_output, expected_output):

    # Compare two files, return bool
    return filecmp.cmp(test_output, expected_output)

def gzFileComp (test_output, expected_output, tmp_dir):

    # Create the tmp paths
    tmp_test_path = os.path.join(tmp_dir, 'Test')
    tmp_expected_path = os.path.join(tmp_dir, 'Expected')

    # Create test tmp directories, if needed
    if not os.path.exists(tmp_test_path):
        os.makedirs(tmp_test_path)

    # Create test expected directories, if needed
    if not os.path.exists(tmp_expected_path):
        os.makedirs(tmp_expected_path)

    # Assign the tmp output files
    tmp_test_output = os.path.join(tmp_test_path, os.path.basename(test_output))
    tmp_expected_output = os.path.join(tmp_expected_path, os.path.basename(expected_output))

    # Open the gzip file
    with gzip.open(test_output, 'rb') as test_file:

        # Open the gunzip file
        with open(tmp_test_output, 'wb') as tmp_test_file:
            
            # Copy the file
            shutil.copyfileobj(test_file, tmp_test_file)

    # Open the gzip file
    with gzip.open(expected_output, 'rb') as expected_file:

        # Open the gunzip file
        with open(tmp_expected_output, 'wb') as tmp_expected_file:
            
            # Copy the file
            shutil.copyfileobj(expected_file, tmp_expected_file)

    # Check if the files have the same content
    file_compare_results = fileComp(tmp_test_output, tmp_expected_output)

    # Remove the tmp dirs
    shutil.rmtree(tmp_test_path)
    shutil.rmtree(tmp_expected_path)

    # Return the results
    return file_compare_results

def vcfFileComp (test_output, test_format, expected_output, tmp_dir):

    # Confirm where the specifed executable is located
    bcftools_path = confirm_executable('bcftools')

    # Check if the executable was found
    if not bcftools_path:
        raise IOError('bcftools not found. Please confirm the executable is installed')

    # Index the vcf.gz test file
    if test_format == 'vcf.gz':

        # Confirm where the specifed executable is located
        tabix_path = confirm_executable('tabix')

        # Check if the executable was found
        if not bcftools_path:
            raise IOError('tabix not found. Please confirm the executable is installed')

        index_call = subprocess.Popen([tabix_path, '-p', 'vcf', test_output], stdout = subprocess.PIPE, stderr = subprocess.PIPE)

    # Index the bcf test file
    elif test_format == 'bcf':

        index_call = subprocess.Popen([bcftools_path, 'index', test_output], stdout = subprocess.PIPE, stderr = subprocess.PIPE)

    # Return error if not vcf.gz or bcf
    else:
        raise Exception('Format must be vcf.gz or bcf')

    # Save the stdout and stderr from index call
    index_stdout, index_stderr = index_call.communicate()

    # Decode if running in python 3
    if sys.version_info[0] == 3:
        index_stderr = index_stderr.decode()

    # Check if error occured
    if index_stderr:
        raise Exception('Error indexing')

    # Assign the isec dir
    isec_path = os.path.join(tmp_dir, 'isec_' + randomGenerator())

    # Run isec
    isec_call = subprocess.Popen([bcftools_path, 'isec', test_output, expected_output, '-p', isec_path], stdout = subprocess.PIPE, stderr = subprocess.PIPE)

    # Save the stdout and stderr from isec call
    isec_stdout, isec_stderr = isec_call.communicate()

    # Decode if running in python 3
    if sys.version_info[0] == 3:
        isec_stderr = isec_stderr.decode()

    if isec_stderr:
        raise Exception(isec_stderr)

    # Create bool to store the file comparison
    file_compare_results = True

    # Check the vcf output for differences
    for isec_filename in os.listdir(isec_path):
        if '0000' in isec_filename or '0001' in isec_filename:
            isec_file_path = os.path.join(isec_path, isec_filename)
            with open(isec_file_path, 'r') as isec_file:
                for isec_line in isec_file:
                    if isec_line.strip() and not isec_line.startswith('#'):
                        file_compare_results = False

    # Delete the isec dir
    shutil.rmtree(isec_path)

    return file_compare_results

def randomGenerator (length = 10, char_set = string.ascii_uppercase + string.digits):

    # Return a random string of letters and numbers
    return ''.join(random.choice(char_set) for i in range(length))
