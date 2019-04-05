import os
import argparse
import logging
import subprocess


from pgpipe.logging_module import initLogger, logArgs

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

    # Input arguments.
    admix_parser.add_argument('--file', help="Input file", type=str, required=True, action=parser_confirm_file())
    admix_parser.add_argument('--pop', help='Number of populations', required=True, type=int)
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
    # Grab admixture arguments from command line
    admix_args = admix_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(admix_args, func_name='admixture')

    # Assign location of admixture file
    #admixture_exec = os.path.join(os.pardir, 'bin', 'admixture')
    admixture_exec = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'bin','admixture')
    #admixture_exec = 'admixture'
    # Check that admixture exists at specified path
    if not os.path.isfile(admixture_exec):
        raise IOError('admixture executable not found in Path specified: %s' % admixture_exec)

    # New list to store user input and options
    admix_call_args = []

    # Add filename to the list
    if admix_args.file:
        admix_call_args.extend([admix_args.file])

    # Add population size(int) to the list
    if admix_args.pop:
        admix_call_args.extend([admix_args.pop])

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
