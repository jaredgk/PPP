import os
import sys
import argparse
import itertools
import copy
import shutil
import logging
import subprocess



# Import basic structure functions
#from structure import *
from logging_module import initLogger

def run ():
    
    structure_call = subprocess.Popen(["./structure"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    global structure_err 
    structure_err= structure_call.communicate()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        structure_out = structure_out.decode()
        structure_err = structure_err.decode()

    check_structure_for_errors(structure_err)

    return structure_err


def check_structure_for_errors (structure_stderr):

    # Returns True if the job completed without error
    if 'Final results' in str(structure_err):
        pass

    # Print output for structure if error is detected
    elif 'error' in str(structure_err):
          logging.error(structure_stderr)
          raise Exception(structure_stderr)
  




run()
    

if __name__ == "__main__":
    initLogger()
    run()

