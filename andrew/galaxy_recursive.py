import os
import sys
import subprocess
import argparse
import logging

sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'jared')))

from logging_module import initLogger


def galaxy_parser():
    '''Galaxy Argument Parser - Assigns arguments for from command line.
    Depending on the argument in question, a default value may be specified'''

    galaxy_parser = argparse.ArgumentParser()

    # Arguments

    # Python script to run (e.g. vcf_calc.py)
    galaxy_parser.add_argument('--script-to-run', help = "The script to be automated", type = str, required = True)

    # List of input files
    galaxy_parser.add_argument('--input-list', help = "Input List", type = str, action='append', required = True)
    # List of output files
    galaxy_parser.add_argument('--output-list', help = "Output List", type = str, action='append')
    # List of links, for creating symlinks of input
    galaxy_parser.add_argument('--link-list', help = "Link List", type = str, action='append')
    # Output Directory, if option is not in script
    galaxy_parser.add_argument('--galaxy-out-arg', help = "Output Arg", type = str, default = 'out')
    # Output call, if option if not --out in script
    galaxy_parser.add_argument('--galaxy-out-dir', help = "Output Directory", type = str)

    return galaxy_parser.parse_known_args()


def run ():

    # Grab arguments from command line
    galaxy_args, command_args = galaxy_parser()

    # Check that there is an eual number of input and links
    if galaxy_args.link_list:
        if len(galaxy_args.link_list) != len(galaxy_args.input_list):
            raise Exception('Link and input number do not match')

    # Check that there is an eual number of input and output
    if galaxy_args.output_list:
        if len(galaxy_args.output_list) != len(galaxy_args.input_list):
            raise Exception('Output and input number do not match')

    # Create the output dir, if required
    if galaxy_args.galaxy_out_dir:
        if not os.path.exists(galaxy_args.galaxy_out_dir):
            os.makedirs(galaxy_args.galaxy_out_dir)

    for input_pos, input_file in enumerate(galaxy_args.input_list):

        if galaxy_args.link_list:
            if os.path.islink(galaxy_args.link_list[input_pos]):
                os.unlink(galaxy_args.link_list[input_pos])
            os.symlink(input_file, galaxy_args.link_list[input_pos])
            input_arg = galaxy_args.link_list[input_pos]
        else:
            input_arg = input_file

        output_arg = []
        if galaxy_args.output_list:
            if galaxy_args.galaxy_out_dir:
                output_path = os.path.join(galaxy_args.galaxy_out_dir, galaxy_args.output_list[input_pos])
            else:
                output_path = galaxy_args.output_list[input_pos]

            output_arg.extend(['--' + galaxy_args.galaxy_out_arg, output_path])

        print ['python', galaxy_args.script_to_run, input_arg] + output_arg + command_args

        python_call = subprocess.Popen(['python', galaxy_args.script_to_run, input_arg] + output_arg + command_args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        python_out, python_err = python_call.communicate()

        if python_err:
            logging.error(python_err)
            raise Exception(python_out)

        if galaxy_args.link_list:
            os.unlink(galaxy_args.link_list[input_pos])


if __name__ == "__main__":
    initLogger()
    run()
