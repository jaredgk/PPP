#!/usr/bin/env python
'''
    Automates the estimation of individual ancestries using Admixture. The
    functions allows for input as: i) Binary-PED files or ii) PED 12-formatted
    files. The function is also capable of configuring the optional arguments
    of Admixture.

    ##################
    Command-line Usage
    ##################
    The admixture automater may be called using the following command:

    .. code-block:: bash
        
        admixture.py

    *************
    Example usage
    *************
    Estimating individual ancestries for each sample within *hapmap3.bed*
    for three ancestral populations.

    .. code-block:: bash
        
        admixture.py --binary-ped-prefix hapmap3 --pop 3

    ############
    Dependencies 
    ############
    * `Admixture <http://software.genetics.ucla.edu/admixture/>`_

    ############################
    Input Command-line Arguments
    ############################
    **--ped-12-prefix** *<input_prefix>*
        Argument used to define the filename prefix shared by the 12-formatted ped 
        file (.ped) and the map file (.map). Should not be used alongside the 
        specific file arguments (e.g. --ped).
    **--ped-12** *<ped_filename>*
        Argument used to define the filename of the plink 12-formatted ped file 
        (.ped). Must be called alongside --map. Cannot be called alongside 
        --ped-prefix.
    **--map** *<map_filename>*
        Argument used to define the filename of the plink map file (.map). Must be 
        called alongside --ped. Cannot be called alongside --ped-prefix.
    **--binary-ped-prefix** *<input_prefix>*
        Argument used to define the filename prefix shared by the binary ped file
        (.bed), the fam file (.fam), and the bim file (.bim). Should not be used 
        alongside the specific file arguments (e.g. --binary-ped).
    **--binary-ped** *<binary_ped_filename>*
        Argument used to define the filename of the plink binary ped file (.bed). 
        Must be called alongside --fam and --bim. Cannot be called alongside 
        --binary-ped-prefix.
    **--fam** *<fam_filename>*
        Argument used to define the filename of the plink fam file (.fam). Must be 
        called alongside --binary-ped and --bim. Cannot be called alongside 
        --binary-ped-prefix.
    **--bim** *<bim_filename>*
        Argument used to define the filename of the plink bim file (.bim). Must be 
        called alongside --binary-ped and --fam. Cannot be called alongside 
        --binary-ped-prefix.

    #############################
    Output Command-line Arguments
    #############################
    **--overwrite**
        Argument used to define if previous output should be overwritten.

    ###############################
    Required Command-line Arguments
    ###############################
    **--pop** *<K_int>*
        Argument used to defines the number of ancestral populations.
    **--admix-method** *<em, block>*
        Argument used to define the algorithm to use. Two algorithm are supported: 
        Block relaxation algorithm (block) or EM algorithm (em). By default, the 
        Block relaxation algorithm is used.

    ###############################
    Optional Command-line Arguments
    ###############################
    **--acceleration** *<acceleration_int>*
        Argument used to defines the value of quasi-Newton acceleration method.
    **--major-converge-likelihood** *<likelihood_float>*
        Argument used to define the major terminaton criterion. Halt when the
        log-likelihood increases by less than the specified value between
        iterations.
    **--major-converge-iter** *<iter_int>*
        Argument used to define the major terminaton criterion. Defines the
        maximum number of iterations.
    **--minor-converge-likelihood** *<likelihood_float>*
        Argument used to define the minor terminaton criterion. Halt when the
        log-likelihood increases by less than the specified value between
        iterations.
    **--minor-converge-iter** *<iter_int>*
        Argument used to define the minor terminaton criterion. Defines the
        maximum number of iterations.
    **--bootstrap** *<bootstrap_int>*
        Argument used to define the number of bootstrap replicates.
    **--random-seed** *<seed_int>*
        Argument used to define the seed value for the random number generator.
    **--threads** *<thread_int>*
        Argument used to define the number of threads to be used for computation.
'''

import os
import argparse
import logging
import subprocess

from pgpipe.logging_module import initLogger, logArgs
from pgpipe.plink import confirm_ped_prefix, confirm_bed_prefix, confirm_ped_files, confirm_bed_files
from pgpipe.misc import confirm_executable

def admix_parser(passed_arguments):
    '''admix Argument Parser - Assigns arguments from command line'''

    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError('%s not found' % value)
                setattr(args, self.dest, value)
        return customAction

    def metavar_list (var_list):
        '''Create a formmated metavar list for the help output'''
        return '{' + ', '.join(var_list) + '}'

    admix_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Input PED/BED prefix arguments
    admix_prefix = admix_parser.add_mutually_exclusive_group()
    admix_prefix.add_argument("--ped-12-prefix", dest = 'ped_prefix', help = "Defines the filename prefix of both PED-12 and MAP files", type = str)
    admix_prefix.add_argument("--binary-ped-prefix", dest = 'bed_prefix', help = "Defines the filename prefix of the Binary-PED, FAM, and BIM files", type = str)

    # Input PED arguments
    admix_parser.add_argument("--ped-12", dest = 'ped_filename', help = "Defines the filename of the PED-12 file. Called alongside --map", type = str, action = parser_confirm_file())
    admix_parser.add_argument("--map", dest = 'map_filename', help = "Defines the filename of the MAP file. Called alongside --ped-12", type = str, action = parser_confirm_file())

    # Input BED arguments
    admix_parser.add_argument("--binary-ped", dest = 'bed_filename', help = "Defines the filename of the Binary-PED (i.e. BED) file. Called alongside --fam and --bim", type = str, action = parser_confirm_file())
    admix_parser.add_argument("--fam", dest = 'fam_filename', help = "Defines the filename of the FAM file. Called alongside --binary-ped and --bim", type = str, action = parser_confirm_file())
    admix_parser.add_argument("--bim", dest = 'bim_filename', help = "Defines the filename of the BIM file. Called alongside --binary-ped and --fam", type = str, action = parser_confirm_file())

    # Population argument
    admix_parser.add_argument('--pop', help = 'Defines the number of populations', required = True, type = int)

    # Method argument
    method_list = ['em', 'block']
    method_default = 'block'
    admix_parser.add_argument('--admix-method', metavar = metavar_list(method_list), help = 'Defines the algorithm to use. Default: block', type = str, choices = method_list, default = method_default)

    # Acceleration arguments
    admix_parser.add_argument('--acceleration', help = 'Defines the quasi-Newton acceleration value', type = int)

    # Convergence arguments
    admix_major = admix_parser.add_mutually_exclusive_group()
    admix_major.add_argument('--major-converge-likelihood', dest = 'major_converge', help = 'Defines the major terminaton criterion for the log-likelihood', type = float)
    admix_major.add_argument('--major-converge-iter', dest = 'major_converge', help = 'Defines the major terminaton criterion for the maximum number of iterations.', type = int)
    admix_minor = admix_parser.add_mutually_exclusive_group()
    admix_minor.add_argument('--minor-converge-likelihood', dest = 'minor_converge', help = 'Defines the minor convergence criterion for the log-likelihood', type = float)
    admix_minor.add_argument('--minor-converge-iter', dest = 'minor_converge', help = 'Defines the minor convergence criterion for the maximum number of iterations', type = int)

    # General arguments
    admix_parser.add_argument('--bootstrap', help = 'Defines the number of bootstrap replicates', type = int)
    admix_parser.add_argument('--threads', help = 'Defines the number of threads to be used for computation', type = int)
    admix_parser.add_argument('--random-seed', help = "Defines the seed value for the random number generator", type = int)
    admix_parser.add_argument('--overwrite', help = "Defines if previous output should be overwritten", action = 'store_true')


    if passed_arguments:
        return admix_parser.parse_args(passed_arguments)
    else:
        return admix_parser.parse_args()


def run(passed_arguments = []):

    # Grab admixture arguments from command line
    admix_args = admix_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(admix_args, func_name = 'admixture')    

    # List to store arguments to be passed to admixture
    admix_call_args = []

    # Check if a method not has been assigned. Error should not be seen
    if not admix_args.admix_method:
        raise Exception('Algorithm assignment error. Please contact development team')

    # Assign the algorithm to the argument list
    admix_call_args.append('--method=' + admix_args.admix_method)

    # Check if the user has assigned an acceleration value
    if admix_args.acceleration:

        # Assign the acceleration value to the argument list
        admix_call_args.extend(['-a', 'qn' + str(admix_args.acceleration)])

    # Check if the user assigned a major terminaton criterion
    if admix_args.major_converge:

        # Assign the criterion to the argument list
        admix_call_args.append('-C=' + str(admix_args.major_converge))

    # Check if the user assigned a minor terminaton criterion
    if admix_args.minor_converge:

        # Assign the criterion to the argument list
        admix_call_args.append('-c=' + str(admix_args.minor_converge))

    # Check if a random seed has been assigned
    if admix_args.random_seed:

        # Assign the seed value to the argument list
        admix_call_args.append('--seed=' + str(admix_args.random_seed))

    # Check if the bootstrap replicates have been assigned 
    if admix_args.bootstrap:

        # Assign the number of bootstrap replicates to the argument list
        admix_call_args.append('-B' + str(admix_args.bootstrap))
    
    # Check if the user has defined the number of threads to use
    if admix_args.threads:

        # Assign the number of threads to the argument list
        admix_call_args.append('-j' + str(admix_args.threads))

    # Check if the a ped prefix was assigned
    if admix_args.ped_prefix and confirm_ped_prefix(admix_args.ped_prefix):

        # Add the ped-12 filename to the argument list
        admix_call_args.append(admix_args.ped_prefix + '.ped')

    # Check if the a bed prefix was assigned
    elif admix_args.bed_prefix and confirm_bed_prefix(admix_args.bed_prefix):

        # Add the binary-ped filename to the argument list
        admix_call_args.append(admix_args.bed_prefix + '.bed')

    # Check if ped files were assigned
    elif confirm_ped_files (admix_args.ped_filename, admix_args.map_filename):

    	# Check if the files have the same prefix. Update later
    	if admix_args.ped_filename[:-4] != admix_args.map_filename[:-4]:
    		raise Exception('PED and MAP files must share the same prefix - i.e. input.ped/input.map')

    	# Add the ped-12 filename to the argument list
    	admix_call_args.append(admix_args.ped_filename)

    # Check if bed files were assigned
    elif confirm_bed_files (admix_args.bed_filename, admix_args.bim_filename, admix_args.fam_filename):

    	# Check if the files have the same prefix. Update later
    	if admix_args.bed_filename[:-4] != admix_args.bim_filename[:-4] or admix_args.bed_filename[:-4] != admix_args.fam_filename[:-4]:
    		raise Exception('Binary-PED, BIM, and FAM files must share the same prefix - i.e. input.bed/input.bim/input.fam')

    	# Add the ped-12 filename to the argument list
    	admix_call_args.append(admix_args.bed_filename)

    else:

    	raise Exception('No input specified. Please check command-line')

    # Add the number of populations to the argument list
    admix_call_args.append(admix_args.pop)

    logging.info('admixture parameters assigned')

    # Confirm where the specifed executable is located
    admixture_path = confirm_executable('admixture')

    # Check if the executable was found
    if not admixture_path:
        raise IOError('admixture not found. Please confirm the executable is installed')

    #admixture_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'bin','admixture')

    # Run 'admixture' executable file with options provided by user
    admixture_call = subprocess.Popen([admixture_path] + list(map(str, admix_call_args)), stdout = subprocess.PIPE, stderr = subprocess.PIPE)

    # Store command output and/or error to variables
    admix_stdout, admix_stderr = admixture_call.communicate()

    # Check if admixture returned an error
    if admix_stderr:
        raise Exception(admix_stderr)

if __name__ == "__main__":
    initLogger()
    run()
