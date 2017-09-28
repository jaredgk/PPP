import os
import sys
import subprocess
import argparse
import glob
import logging

sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'jared')))

import vcf_reader_func
from logging_module import initLogger

def phase_argument_parser(passed_arguments):
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

    def parser_confirm_no_file ():
        '''Custom action to confirm file does not exist'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if os.path.isfile(value):
                    raise IOError # File found
                setattr(args, self.dest, value)
        return customAction

    def metavar_list (var_list):
        '''Create a formmated metavar list for the help output'''
        return '{' + ', '.join(var_list) + '}'

    phase_parser = argparse.ArgumentParser()

    # Input arguments.
    phase_parser.add_argument("vcfname", metavar='VCF_Input', help = "Input VCF filename", type = str, action = parser_confirm_file())

    phasing_list = ['beagle', 'shapeit']
    phasing_default = 'beagle'
    phase_parser.add_argument('--phase-algorithm', metavar = metavar_list(phasing_list), help = 'Specifies the phase algorithm to be used', type = str, choices = phasing_list, default = phasing_default)

    # Other basic arguments. Expand as needed
    phase_parser.add_argument('--out', help = 'Defines the output filename', default = 'out', action = parser_confirm_no_file())

    if passed_arguments:
        return phase_parser.parse_args(passed_arguments)
    else:
        return phase_parser.parse_args()

def logArgs(args, pipeline_function):
    '''Logs arguments from argparse system. May be replaced with a general function in logging_module'''
    logging.info('Arguments for %s:' % pipeline_function)
    for k in vars(args):
        logging.info('Arguments %s: %s' % (k, vars(args)[k]))

def possible_beagle_paths ():
    possible_paths = ['beagle*.jar']
    if os.name == 'posix':
        possible_paths.extend(['\usr\local\bin\beagle*.jar', '\usr\bin\beagle*.jar'])

    for current_path in possible_paths:
        print [n for n in glob.glob(current_path) if os.path.isfile(n)]

def assign_vcf_extension (filename):
    # Checks if the file is unzipped, bgzipped, or gzipped
    vcfname_format = vcf_reader_func.checkFormat(filename)

    if vcfname_format == 'nozip':
        return '.vcf'
    elif vcfname_format == 'gzip' or vcfname_format == 'bgzip':
        return '.vcf.gz'
    else:
        sys.exit('Unknown file format')

def run (passed_arguments = []):
    '''
        Phaser for VCF files.

        Automates the phasing process for a specified VCF file. The function
        allows users to select between multiple phasing algorithms: beagle
        (default) and shapit.

        Parameters
        ----------
        VCF_Input : str
            Specifies the input VCF filename
        --phase-algorithm : str
            Specifies the algorithm to be used. Choices: beagle (default) and
            shapit
        --out : str
            Specifies the output filename
        Returns
        -------
        output : file
            Phased VCF file

        Raises
        ------
        IOError
            Input VCF file does not exist
        IOError
            Output file already exists

    '''

    # Grab VCF arguments from command line
    phase_args = phase_argument_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(phase_args, 'vcf_phase')

    # Assign file extension for VCF input file
    vcfname_ext = assign_vcf_extension(phase_args.vcfname)

    logging.info('Input file assigned')

    # Used to confirm if the VCF input was renamed
    vcfname_renamed = False

    # Confirm input has correct file extension
    if vcfname_ext not in phase_args.vcfname:
        vcfname_renamed = True
        os.rename(phase_args.vcfname, phase_args.vcfname + vcfname_ext)
        phase_args.vcfname += vcfname_ext

    if phase_args.phase_algorithm == 'beagle':
        # Assign the algorithm
        algorithm_call_args = ['java', '-jar', 'bin/beagle.jar']

        # Assigns the arguments for phasing
        phase_call_args = ['gt=' + phase_args.vcfname, 'out=' + phase_args.out]

        logging.info('beagle phasing parameters assigned')

        # Phasing subprocess call
        phase_call = subprocess.Popen(algorithm_call_args + phase_call_args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        phase_out, phase_err = phase_call.communicate()
        print phase_out, phase_err

        logging.info('beagle phasing complete')

    elif phase_args.phase_algorithm == 'shapeit':
        # Assign the algorithm
        algorithm_call_args = ['./bin/shapeit']

        # Assigns the arguments for phasing
        phase_call_args = ['--input-vcf', phase_args.vcfname, '-O', phase_args.out]

        logging.info('shapeit phasing parameters assigned')

        # Phasing subprocess call
        phase_call = subprocess.Popen(algorithm_call_args + phase_call_args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        phase_out, phase_err = phase_call.communicate()

        # Confirms that shapeit finished without error
        if phase_err:
            logging.error('Error occured in phasing. Please check input file.')
            sys.exit('Error occured in phasing. Please check input file.')

        logging.info('shapeit phasing complete (non-VCF format)')

        # Assigns the arguments for converting the phased output into a vcf file
        convert_call = ['-convert', '--input-haps', phase_args.out, '--output-vcf', phase_args.out + vcfname_ext]

        logging.info('shapeit conversion parameter assigned')

        # Convert subprocess call
        convert_call = subprocess.Popen(algorithm_call_args + convert_call, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        convert_out, convert_err = convert_call.communicate()
        print convert_out, convert_err

        logging.info('shapeit conversion complete (non-VCF to VCF)')


    # Reverts the VCF input file
    if vcfname_renamed:
        os.rename(phase_args.vcfname, phase_args.vcfname[:-len(vcfname_ext)])

if __name__ == "__main__":
    initLogger()
    run()
