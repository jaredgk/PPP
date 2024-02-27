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
    
    pixy_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaulsHelpFormatter)

