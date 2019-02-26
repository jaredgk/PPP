import os
import sys
import subprocess
import shutil
import argparse
import glob
import logging
import pybedtools

#sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'pppipe')))

from pgpipe.logging_module import initLogger

def bed_argument_parser(passed_arguments):
    '''Phase Argument Parser - Assigns arguments for vcftools from command line.
    Depending on the argument in question, a default value may be specified'''

    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError # File not found
                setattr(args, self.dest, value)
        return customAction

    def parser_confirm_append_file ():
        '''Custom action to house data as defaultdict list'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if getattr(args, self.dest):
                    # Check the file exists
                    if not os.path.isfile(value):
                        raise IOError
                    # Append the argument with the file
                    getattr(args, self.dest).append(value)
                else:
                    # Check the file exists
                    if not os.path.isfile(value):
                        raise IOError
                    # Create the argument list with the file
                    setattr(args, self.dest, [value])
        return customAction

    def metavar_list (var_list):
        '''Create a formmated metavar list for the help output'''
        return '{' + ', '.join(var_list) + '}'

    bed_parser = argparse.ArgumentParser()

    # Input arguments.
    bed_parser.add_argument("--bed", help = "Input BED filename", type = str, action = parser_confirm_append_file())

    task_list = ['cat']
    task_default = 'cat'
    bed_parser.add_argument('--bed-task', metavar = metavar_list(task_list), help = 'Specifies the phase algorithm to be used', type = str, choices = task_list, default = task_default)

    # Other basic arguments. Expand as needed
    bed_parser.add_argument('--out', help = 'Defines the output filename', default = 'out.bed')

    if passed_arguments:
        return bed_parser.parse_args(passed_arguments)
    else:
        return bed_parser.parse_args()

def logArgs(args, pipeline_function):
    '''Logs arguments from argparse system. May be replaced with a general function in logging_module'''
    logging.info('Arguments for %s:' % pipeline_function)
    for k in vars(args):
        logging.info('Arguments %s: %s' % (k, vars(args)[k]))

def return_BedTools (bed_files):

    # List to hold BedTool objects
    bedtool_list = []

    # Loop BED files
    for bed_file in bed_files:
        # Save BedTool object
        bedtool_list.append(pybedtools.BedTool(bed_file))

    # Return BedTool objects
    return bedtool_list


def run (passed_arguments = []):
    '''
        Summary Line

        Complete Summary

        Parameters
        ----------

        Returns
        -------

        Raises
        ------

    '''

    # Grab BED arguments from command line
    bed_args = bed_argument_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(bed_args, 'bed_tasks')

    # Load BED files
    bedtool_files = return_BedTools(bed_args.bed)

    # Check if the cat function was selected
    if bed_args.bed_task == 'cat':

        # Run the cat function and save the results
        bedtool_files[0].cat(*bedtool_files[1:], output = bed_args.out)

if __name__ == "__main__":
    initLogger()
    run()
