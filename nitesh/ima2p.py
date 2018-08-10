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

def call_ima2p (ima2p_call_args):
     

    #ima2p_call = subprocess.Popen(['mpirun']+['-np']+['5']+['IMa2p'] + list(map(str,ima2p_call_args)) , stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    a = ['mpirun','-np','5','IMa2p','-i','Sim1_5loci.u',
                        '-o','Sim1_5loci.out',
                        '-q', '2',
                        '-m','1',
                        '-t','3',  
                        '-hf','g',                      
                        '-ha', '0.98',
                        '-hb', '0.75',
                        '-r', '245',
                        '-b', '1000',
                        '-l', '100']

    ima2p_call = subprocess.Popen(list(map(str,ima2p_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    b = subprocess.list2cmdline(a)
    print b



    #print ( ['IMa2p']  + list(map(str,ima2p_call_args)))

    ima2p_out, ima2p_err= ima2p_call.communicate()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        
        ima2p_out = ima2p_out.decode()
        ima2p_err = ima2p_err.decode()
       

    check_ima2p_for_errors(ima2p_err)

    return ima2p_out, ima2p_err


def check_ima2p_for_errors (ima2p_stderr):


    if 'Error' in str(ima2p_stderr):
        # Splits log into list of lines
        ima2p_stderr_lines = vcftools_stderr.splitlines()
        # Prints the error(s)
        logging.error('\n'.join((output_line for output_line in ima2p_stderr_lines if output_line.startswith('Error'))))
        raise Exception('\n'.join((output_line for output_line in ima2p_stderr_lines if output_line.startswith('Error'))))



def check_for_ima2p_output (ima2p_output):
 
    # Check if output file already exists
    if os.path.isfile(ima2p_output):
        logging.error('IMa2p output file already exists')
        raise IOError('IMa2p output file already exists')

    logging.info('Output file assigned')

    # Check if log file already exists
    if os.path.isfile(ima2p_output + '.log'):
        logging.error('Log file already exists')
        raise IOError('Log file already exists')

    logging.info('Log file assigned')
       
 
def produce_ima2p_output (output, filename, append_mode = False):


    # Check if single log file is required from multiple calls
    if append_mode:
        ima2p_log_file = open(filename + '.output.log','a')
    else:
        ima2p_log_file = open(filename + '.output.log','w')

    ima2p_log_file.write(str(output))
    ima2p_log_file.close()


def produce_ima2p_log (output, filename, append_mode = False):

    if append_mode:
        ima2p_log_file = open(filename + '.log','a')
    else:
        ima2p_log_file = open(filename + '.log','w')

    ima2p_log_file.write(str(output))
    ima2p_log_file.close()

    


