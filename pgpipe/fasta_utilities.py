import os
import sys
import argparse
import itertools
import copy
import shutil
import logging
import pandas as pd

# Insert Jared's directory path, required for calling Jared's functions. Change when directory structure changes.
#sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'pppipe')))

from pgpipe.logging_module import initLogger, logArgs

from pgpipe.picard import call_picard

def fasta_utility_parser(passed_arguments):
    '''fasta Argument Parser - Assigns arguments from command line'''

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

    fasta_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    # Input arguments.
    fasta_parser.add_argument('--fasta', help = "Input FASTA filename", type = str, nargs = '+', required = True, action = parser_confirm_file())

    # Utility based arguments
    utility_list = ['ref-dictionary']
    fasta_parser.add_argument('--utility', metavar = metavar_list(utility_list), help = 'The utility to use', type = str, choices = utility_list)

    # Output file arguments
    out_format_list = ['fasta', 'fasta.gz']
    out_format_default = 'fasta.gz'

    fasta_parser.add_argument('--out-format', metavar = metavar_list(out_format_list), help = 'Specifies the format of FASTA-based output', type = str, choices = out_format_list, default = out_format_default)


    fasta_parser.add_argument('--out', help = 'Output filename. If used, overrides --out-prefix', type = str)
    fasta_parser.add_argument('--out-prefix', help = 'Defines the output prefix', default = 'out')

    fasta_parser.add_argument('--overwrite', help = "Overwrite previous output files", action = 'store_true')

    fasta_parser.add_argument('--picard-path', help = "Defines path to locate picard.jar", type = str)

    if passed_arguments:
        return fasta_parser.parse_args(passed_arguments)
    else:
        return fasta_parser.parse_args()

def run (passed_arguments = []):
    '''
    Liftover for VCF-formatted files

    Automates the liftover function from picard for VCF-formatted files.

    Parameters
    ----------
    --fasta : str
        Input FASTA filename
	--utility : str
        Specifies the utility to be used. Choices: ref-dictionary
    --out-format : str
        Format of the output file
    --out-prefix : str
        Output filename prefix
    --out : str
        Output filename. Overrides --out-prefix
    --overwrite
        If previous output should be overwriten

    Returns
    -------
    output : file
        Statistic file output
    log : file
        Log file output

    Raises
    ------
    IOError
        Input FASTA file does not exist
    IOError
        Output file already exists
    '''

    # Grab FASTA arguments from command line
    fasta_args = fasta_utility_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(fasta_args, func_name = 'fasta_utilities')

    # Check if a reference dictionary was not assigned
    if fasta_args.utility == 'ref-dictionary':

        # Create arg list for sequence dictionary
        picard_seqdict_call_args = ['CreateSequenceDictionary']

        # Add the reference file to the sequence dictionary args
        picard_seqdict_call_args.append('R=' + fasta_args.fasta)

        # Create dictionary for fasta reference
        call_picard(fasta_args.picard_path, picard_seqdict_call_args)


if __name__ == "__main__":
    #initLogger()
    run()