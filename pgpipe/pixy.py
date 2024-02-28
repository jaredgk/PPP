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


def pixy_calc_parser(passed_arguments=[]):
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
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required arguments
    pixy_parser.add_argument(
        "--vcf",
        help="Defines the path and/or filename of the VCF",
        type=str,
        required=True,
        action=parser_confirm_file(),
    )
    pixy_parser.add_argument(
        "--pop-file",
        help="Defines the path and/or filename of the populations file",
        type=str,
        required=True,
        action=parser_confirm_file(),
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
        "--bed-file",
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
        help="Defines the output directory",
        required=False,
    )
    pixy_parser.add_argument(
        "--out",
        type=str,
        nargs="?",
        help="Defines prefix for output file(s)",
        required=False,
    )

    if passed_arguments:
        return vars(pixy_parser.parse_args(passed_arguments))
    else:
        return vars(pixy_parser.parse_args())


def run(**kwargs):
    """
    Statistic Calculations using Pixy.

    This function uses the argparse-based function :py:func:`pixy_parser`
    to parse either sys.argv or passed_arguments to obtain the parameters below.
    The parameters are then translated to their pixy equivalent. Once all the
    parameters are assigned, pixy is called.

    Parameters
    ----------
    --vcf : str
        Input VCF filename and/or path
    --pop-file : str
        Input populations filename and/or path
    --stats : str
        Specifies statistic(s) to calculate
    --statistic-window-size : int
        Specifices window size
    --bed-file : str
        Path to a .BED file
    --out-dir : str
        Output directory
    --out  : str
        Output prefix

    Raises
    ------
    IOError
        Output file already exists and --overwrite is not specified
    Exception
        Incompatible arguments
    """
    # Update kwargs with defaults
    if __name__ != "__main__":
        kwargs = argprase_kwargs(kwargs, pixy_calc_parser)

    # Assign arguments
    pixy_args = argparse.Namespace(**kwargs)

    # Adds arguments to the log file
    logArgs(pixy_args, func_name="pixy")

    # Argument container for pixy
    pixy_call_args = []

    # Add required arguments
    if pixy_args.vcf:
        pixy_call_args.extend(["--vcf", pixy_args.vcf])

    if pixy_args.pop_file:
        pixy_call_args.extend(["--populations", pixy_args.pop_file])

    if pixy_args.stats:
        pixy_call_args.extend(["--stats"])
        pixy_call_args.extend(pixy_args.stats)

    # Assign one of either
    if pixy_args.bed_file and pixy_args.statistic_window_size:
        raise Exception(
            "Cannot have both bed-file and statistic-window-size parameters. Must be one or the other"
        )

    elif not pixy_args.bed_file and not pixy_args.statistic_window_size:
        raise Exception("Must specificy a bed file or a window size")

    if pixy_args.bed_file and not pixy_args.statistic_window_size:
        pixy_call_args.extend(["--bed_file", pixy_args.bed_file])

    elif pixy_args.statistic_window_size and not pixy_args.bed_file:
        pixy_call_args.extend(["--window_size", str(pixy_args.statistic_window_size)])

    # Add optional arguments
    if pixy_args.out:
        pixy_call_args.extend(["--output_prefix", pixy_args.out])

    if pixy_args.out_dir:
        pixy_call_args.extend(["--output_folder", pixy_args.out_dir])

    # TODO for testing
    print(pixy_call_args)


if __name__ == "__main__":
    initLogger()
    run(**pixy_calc_parser())
