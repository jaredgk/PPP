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

def call_structure ():
    
    structure_call = subprocess.Popen(["./structure"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    structure_err= structure_call.communicate()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        
        structure_err = structure_err.decode()

    check_structure_for_errors(structure_err)

    return structure_err


def check_structure_for_errors (structure_stderr):


    # Print output for structure if error is detected
    if 'error' in str(structure_stderr):
       raise Exception(structure_stderr)
  
def produce_structure_log (output, append_mode = False):

    if append_mode:
        structure_log_file = open('structurefile' + '.log','a')
    else:
        structure_log_file = open('strcuturefile' + '.log','w')

    structure_log_file.write(str(output))
    structure_log_file.close()

    


