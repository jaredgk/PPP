import os
import argparse
import logging
import subprocess


from pgpipe.logging_module import initLogger, logArgs
from pgpipe.plink import confirm_ped_prefix, confirm_bed_files, confirm_ped_files, confirm_bed_files

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

    admix_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Input PED arguments
    convert_parser.add_argument("--ped-12", dest = 'ped_filename', help = "Defines the filename of the PED-12 file. Called alongside --map", type = str, action = parser_confirm_file())
    convert_parser.add_argument("--map", dest = 'map_filename', help = "Defines the filename of the MAP file. Called alongside --ped-12", type = str, action = parser_confirm_file())

    # Input BED arguments
    convert_parser.add_argument("--binary-ped", dest = 'bed_filename', help = "Defines the filename of the Binary-PED (i.e. BED) file. Called alongside --fam and --bim", type = str, action = parser_confirm_file())
    convert_parser.add_argument("--fam", dest = 'fam_filename', help = "Defines the filename of the FAM file. Called alongside --binary-ped and --bim", type = str, action = parser_confirm_file())
    convert_parser.add_argument("--bim", dest = 'bim_filename', help = "Defines the filename of the BIM file. Called alongside --binary-ped and --fam", type = str, action = parser_confirm_file())

    # Input PED/BED prefix arguments
    convert_prefix = convert_parser.add_mutually_exclusive_group()
    convert_prefix.add_argument("--ped-12-prefix", help = "Defines the filename prefix of both PED-12 and MAP files", type = str)
    convert_prefix.add_argument("--binary-ped-prefix", dest = 'bed_prefix', help = "Defines the filename prefix of the Binary-PED, FAM, and BIM files", type = str)

    admix_parser.add_argument('--pop', help = 'Number of populations', required = True, type = int)

    admix_parser.add_argument('--random-seed', help='random seed', type=int)
    admix_parser.add_argument('--threads', help='No. of threads to be used for computation', type=int)
    admix_parser.add_argument('--method', help='Algorithm to be used (em or block)', type=str, choices={"em", "block"})
    admix_parser.add_argument('--acceleration', help='Set acceleration(sqs<X> or qn<X>)', type=str)
    admix_parser.add_argument('--major-converge', help='set major convergence criterion (for point estimation)', type=int)
    admix_parser.add_argument('--minor-converge', help='set minor convergence criterion (for bootstrap and CV reestimates)', type=int)
    admix_parser.add_argument('--bootstrap', help='Bootstrapping [with X replicates]', type=int)


    if passed_arguments:
        return admix_parser.parse_args(passed_arguments)
    else:
        return admix_parser.parse_args()


def run(passed_arguments = []):

	# Assign location of admixture file
    admixture_exec = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'bin','admixture')

    # Check that admixture exists at specified path
    if not os.path.isfile(admixture_exec):
        raise IOError('admixture executable not found in Path specified: %s' % admixture_exec)


    # Grab admixture arguments from command line
    admix_args = admix_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(admix_args, func_name='admixture')    

    # List to store arguments to be passed to admixture
    admix_call_args = []

    # Check if ped files were assigned
    if confirm_ped_files (admix_args.ped_filename, admix_args.map_filename):

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

    # Check if the a ped prefix was assigned
    elif admix_args.ped_12_prefix and confirm_ped_prefix(admix_args.ped_12_prefix):

    	# Add the ped-12 filename to the argument list
    	admix_call_args.append(admix_args.ped_12_prefix + '.ped')

    # Check if the a bed prefix was assigned
    elif admix_args.binary_ped_prefix and confirm_bed_prefix(admix_args.binary_ped_prefix):

    	# Add the binary-ped filename to the argument list
    	admix_call_args.append(admix_args.binary_ped_prefix + '.bed')

    else:

    	raise Exception('No input specified. Please check command-line')

    # Add the number of populations to the argument list
    admix_call_args.append(admix_args.pop)

    # Add random seed(int) to the list
    if admix_args.random_seed:
        admix_call_args.extend(['--seed=' + str(admix_args.random_seed)])

    # Add threads(int) to the list
    if admix_args.threads:
        admix_call_args.extend(['-j' + str(admix_args.threads)])

    # Add algorithm(em or block) to the list
    if admix_args.method:
        admix_call_args.extend(['--method=' + admix_args.method])

    # Add acceleration(int) to the list
    if admix_args.acceleration:
        admix_call_args.extend(['-a', admix_args.acceleration])

    # Add convergence(int) to the list
    if admix_args.major_converge:
        admix_call_args.extend(['-C=' + str(admix_args.major_converge)])

    if admix_args.minor_converge:
        admix_call_args.extend(['-c=' + str(admix_args.minor_converge)])

    # Add bootstrap(int) to the list
    if admix_args.bootstrap:
        admix_call_args.extend(['-B[' + str(admix_args.bootstrap) + ']'])

    logging.info('admixture parameters assigned')

    # Run 'admixture' executable file with options provided by user
    admixture_call = subprocess.Popen(admixture_exec + " " + ' '.join(map(str, admix_call_args)),
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

    # Store command output and/or error to variables
    admix_stdout, admix_stderr = admixture_call.communicate()

    if admixture_call.returncode == 0:
        logging.info("Successful")
    else:
        logging.error(admix_stderr)


if __name__ == "__main__":
    initLogger()
    run()
