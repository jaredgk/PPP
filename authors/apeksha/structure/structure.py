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
    structure_out, structure_err= structure_call.communicate()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        
        structure_out = structure_out.decode()
        structure_err = structure_err.decode()
       

    check_structure_for_errors(structure_err)

    return structure_out, structure_err


def check_structure_for_errors (structure_stderr):


    # Print output for structure if error is detected
    if 'error' in str(structure_stderr):
       #print 'error'
       raise Exception(structure_stderr)
      # structure_log_file = open('structureError' + '.log','w')
      # structure_log_file.write('error')
      # structure_log_file.close()
       
  
def produce_structure_log (output, append_mode = False):

    if append_mode:
        structure_log_file = open('structureOutput' + '.log','a')
    else:
        structure_log_file = open('structureOutput' + '.log','w')

    structure_log_file.write(str(output))
    structure_log_file.close()

    


