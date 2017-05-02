import os, sys, subprocess, argparse, glob

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

    phase_parser = argparse.ArgumentParser()

    # Input arguments.
    phase_parser.add_argument("vcfname", metavar='VCF_Input', help = "Input VCF filename", type = str, action = parser_confirm_file())

    phasing_list = ['beagle', 'shapeit']
    phasing_default = 'beagle'
    phase_parser.add_argument('--phase-algorithm', metavar = '{' + ', '.join(phasing_list) + '}', help = 'Specifies the phase algorithm to be used', type = str, choices = phasing_list, default = phasing_default)

    # Other basic arguments. Expand as needed
    phase_parser.add_argument('--out', help = 'Defines the output filename', default = 'out', action = parser_confirm_no_file())

    if passed_arguments:
        return phase_parser.parse_args(passed_arguments)
    else:
        return phase_parser.parse_args()

def possible_beagle_paths ():
    possible_paths = ['beagle*.jar']
    if os.name == 'posix':
        possible_paths.extend(['\usr\local\bin\beagle*.jar', '\usr\bin\beagle*.jar'])

    for current_path in possible_paths:
        print [n for n in glob.glob(current_path) if os.path.isfile(n)]

def run (passed_arguments = []):
    ''' Wrapper code for Phasing. Commands are assigned using argparse.'''

    # Grab VCF arguments from command line
    phase_args = phase_argument_parser(passed_arguments)

    # Container for the phasing subprocess command
    phase_call_args = ['out=' + phase_args.out]

    if phase_args.phase_algorithm:
        # Assign the algorithm
        phase_call_algorithm = ['java', '-jar', 'bin/beagle.jar']
        # Assign the arguments for the algorithm
        phase_call_args.append('gl=' + phase_args.vcfname)

    # Assigns the vcftools arguments to a single command line
    #beagle_command =  [split_arg for arg in vars(beagle_args) for split_arg in getattr(beagle_args, arg)]

    # vcftools subprocess call
    phase_call = subprocess.Popen(phase_call_algorithm + phase_call_args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    phase_out, phase_err = phase_call.communicate()
    print phase_out, phase_err

if __name__ == "__main__":
    run()
