import os
import argparse
import logging
import subprocess


from logging_module import initLogger, logArgs

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

    admixtools_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Input arguments.
    admixtools_parser.add_argument('--file', help="Input file", type=str, required=True, action=parser_confirm_file())

    if passed_arguments:
        return admixtools_parser.parse_args(passed_arguments)
    else:
        return admixtools_parser.parse_args()


def run(passed_arguments = []):
    # Grab admixtools arguments from command line
    admixtools_args = admixtools_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(admixtools_args, func_name='admixtools')

    # Assign location of admixtools file
    admixtools_exec = os.path.join(os.getcwd(), '', 'AdmixTools/bin/convertf')

    # Check that beagle.jar exists
    if not os.path.isfile(admixtools_exec):
        raise IOError('executable not found in Path specified: %s' % admixtools_exec)

    # New list to store user input and options
    admixtools_call_args = []

    # Add formatted input format type and output format type to the list (e.g "par.EIGENSTRAT.PED")

    # Compulsory parameter argument
    admixtools_call_args.extend(['-p'])

    if admixtools_args.file:
        admixtools_call_args.extend([admixtools_args.file])

    logging.info('admixtools parameters assigned')

    # Run 'admixtools' executable file with options provided by user
    admixture_call = subprocess.Popen(admixtools_exec + " " + ' '.join(map(str, admixtools_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

    # Store command output and/or error to variables
    admixtools_stdout, admixtools_stderr = admixture_call.communicate()

    if admixture_call.returncode == 0:
        logging.info("Successful")
    else:
        logging.error(admixtools_stderr)


if __name__ == "__main__":
    initLogger()
    run()
