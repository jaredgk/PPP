import os
import sys
import subprocess
import argparse
import logging
import itertools

# Import basic vcftools functions
from ima2p import *

# Insert Jared's directory path, required for calling Jared's functions. Change when directory structure changes.
#sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'jared')))

from logging_module import initLogger

def ima2p_filter_parser(passed_arguments):
    '''Ima2p Argument Parser - Assigns arguments from command line'''

    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError # File not found
                setattr(args, self.dest, value)
        return customAction

    def metavar_list (var_list):
        '''Create a formmated metavar list for the help output'''
        return '{' + ', '.join(var_list) + '}'

    ima2p_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)



    ima2p_parser.add_argument("-i", metavar = 'Ima2p_Input', help = "Specifies the name of input data file", type = str)
    ima2p_parser.add_argument("-o", help = "Specifies the name for program output files", type = str)

    ima2p_parser.add_argument("-b", help = "Specifies the duration of burn : number of burnin steps")

    ima2p_parser.add_argument("-c", help = "Specifies the calculation options", type = int)

    ima2p_parser.add_argument("-d", help = "Specifies the number of steps between tree saving", type = int)

    ima2p_parser.add_argument("-f", help = "Specifies the name of file with saved Markov chain state generated in previous run -use with -r3", type = str)

    ima2p_parser.add_argument("-g", help = "Specifies the name of file with parameter priors  (requires -c3)", type = str)

    #ima2p_parser.add_argument("-hfl", help = "Specifies the heating model : linear")

    #ima2p_parser.add_argument("-hfg", help = "Specifies the heating model : geometric") 
    ima2p_parser.add_argument("-hf", help = "Specifies the heating model : geometric", type = str) 

    ima2p_parser.add_argument("-hn", help = "Specifies the number of chains", type = int)

    ima2p_parser.add_argument("-hk", help = "Specifies the number of chain swap attempts per step", type = int)

    ima2p_parser.add_argument("-ha", help = "Specifies the first heating parameter, effect depends on heating model", type = float)

    ima2p_parser.add_argument("-hb", help = "Specifies the second heating parameter, effect depends on heating model", type = float)

    ima2p_parser.add_argument("-j", help = "Specifies the model options", type = int)

    ima2p_parser.add_argument("-l", help = "Specifies the run duration", type = int)

    ima2p_parser.add_argument("-m", help = "Specifies the migration prior value", type = float)

    ima2p_parser.add_argument("-p", help = "Specifies the output options", type = int)

    ima2p_parser.add_argument("-q", help = "Specifies the maximum for population size parameters", type = int)

    ima2p_parser.add_argument("-r", help = "Specifies the run options", type = int)

    ima2p_parser.add_argument("-s", help = "Specifies the random number seed", type = int)

    ima2p_parser.add_argument("-t", help = "Specifies the maximum time of population splitting", type = int)

    ima2p_parser.add_argument("-u", help = "Specifies the Generation time in years -for use with -p3", type = int)

    ima2p_parser.add_argument("-v", help = "Specifies the base name", type = str)

    ima2p_parser.add_argument("-w", help = "Specifies the name of file with nested models to be tested (LOAD-GENEALOGY mode only), invokes -c2)", type = str)

    ima2p_parser.add_argument("-y", help = "Mutation rate scalar for relevant loci - for use with -p3")

    ima2p_parser.add_argument("-z", help = "Specifies the number of stepsbetween screen output", type = int)
 
    if passed_arguments:
      return ima2p_parser.parse_args(passed_arguments)     
    else:
      return ima2p_parser.parse_args()


def logArgs(args, pipeline_function):
    '''Logs arguments from argparse system. May be replaced with a general function in logging_module'''
    logging.info('Arguments for %s:' % pipeline_function)
    for k in vars(args):
        logging.info('Arguments %s: %s' % (k, vars(args)[k]))

def run (passed_arguments = []):
    
    ima2p_args = ima2p_filter_parser(passed_arguments)
    
    c = []
 
    for arg in vars(ima2p_args):
        if getattr(ima2p_args, arg) != None:
           a = '-'+arg, getattr(ima2p_args, arg)
           c.append(a)

    
    d = list(itertools.chain.from_iterable(c))
    e = list(map(str,d))
    #print e
    

    logArgs(ima2p_args, 'ima2p_func')
    

    ima2p_out, ima2p_err = call_ima2p(['mpirun','-np','5','IMa2p','-i','Sim1_5loci.u',
                        '-o','Sim1_5loci.out',
                        '-q', '2',
                        '-m','1',
                        '-t','3',  
                        '-hf','g',                      
                        '-ha', '0.98',
                        '-hb', '0.75',
                        '-r', '245',
                        '-b', '1000',
                        '-l', '100'])
   # check_for_ima2p_output(ima2p_args.o)
    produce_ima2p_output(ima2p_out, ima2p_args.o)
    produce_ima2p_log(ima2p_err, ima2p_args.o)

run(['-i','Sim1_5loci.u',
                        '-o','Sim1_5loci.out',
                        '-q', '2',
                        '-m','1',
                        '-t','3',
                        '-hf','g',
                        '-hn','1',                        
                        '-ha', '0.98',
                        '-hb', '0.75',
                        '-r', '245',
                        '-b', '1000',
                        '-l', '100'])

if __name__ == "__main__":
    initLogger()
    


