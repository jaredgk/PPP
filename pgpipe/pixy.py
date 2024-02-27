#!/usr/bin/env python
"""
    Automates estimating nucleotide diversity within (π) and between (dxy) 
    populations from a VCF using pixy

    ##################
    Command-line Usage
    ##################
    The pixy statistic calculator may be called using the following command:

     .. code-block:: bash
        
        pixy.py
"""
import os
import argparse

from pgpipe.logging_module import initLogger, logArgs
from pgpipe.misc import argprase_kwargs


def pixy_parser(passed_arguments=[]):
    """
    Pixy Argument Parser

    Assign the parameters for Pixy using argparse.

    Parameters
    ----------
    passed_arguments : list, optional
        Parameters passed by another function. sys.argv is used if
        not given. 

    Raises
    ------
    IOError
        If the input, or other specified files do not exist
    """

    def parser_confirm_file():
        """Custom action to confirm file exists"""

        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError("%s not found" % value)
                setattr(args, self.dest, value)

        return customAction

    def parser_confirm_file_list():
        """Custom action to confirm file exists in list"""

        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                # Loop the list
                for value_item in value:
                    # Check if the file exists
                    if not os.path.isfile(value_item):
                        raise IOError("%s not found" % value_item)
                if not getattr(args, self.dest):
                    setattr(args, self.dest, value)
                else:
                    getattr(args, self.dest).extend(value)

        return customAction

    def parser_add_to_list():
        """Custom action to add items to a list"""

        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):

                # Clean up any commas
                value = [item.strip(",") for item in value]

                if not getattr(args, self.dest):
                    setattr(args, self.dest, value)
                else:
                    getattr(args, self.dest).extend(value)

        return customAction

    def metavar_list(var_list):
        """Create a formmated metavar list for the help output"""
        return "{" + ", ".join(var_list) + "}"

    pixy_parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaulsHelpFormatter
    )

    # Required arguments
    pixy_parser.add_argument(
        "-vcf",
        help="Defines the path and/or filename of the VCF",
        type=str,
        required=True,
        action=parser_confirm_file,
    )
    pixy_parser.add_argument(
        "--pop-file",
        help="Defines the path and/or filename of the populations file",
        type=str,
        required=True,
        action=parser_confirm_file,
    )
    pixy_parser.add_argument(
        "--stats",
        help="Defines which statistic(s) pixy should perform [pi|fst|dxy]",
        type=str,
        required=True,
        nargs="+",
        choices=["pi", "dxy", "fst"],
    )

    # One of either of these is required
    pixy_parser.add_argument(
        "--statistic-window-size",
        type=int,
        nargs="?",
        help="Defines window size in base pairs over which to calculate stats",
        required=False,
    )
    pixy_parser.add_argument(
        "--bed_file",
        type=str,
        nargs="?",
        help="Path to a headerless .BED file with columns [chrom chromStart chromEnd]",
        required=False,
    )

    # Optional arguments
    pixy_parser.add_argument(
        "--out-dir",
        type=str,
        nargs="?",
        default="",
        help="Defines the output directory",
        required=False,
    )
    pixy_parser.add_argument(
        "--out",
        type=str,
        nargs="?",
        default="pixy",
        help="Defines prefix for output file(s)",
        required=False,
    )

    if passed_arguments:
        return vars(pixy_parser.parse_args(passed_arguments))
    else:
        return vars(pixy_parser.parse_args())
    
