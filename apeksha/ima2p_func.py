import os
import sys
import subprocess
import argparse
import logging
import itertools
from distutils.spawn import find_executable

# Import basic vcftools functions
from ima2p import *

# Insert Jared's directory path, required for calling Jared's functions. Change when directory structure changes.
#sys.path.insert(0, os.path.abspath(os.path.join(os.pardir, 'jared')))

#from logging_module import initLogger

def imaParser():
    '''Ima2p Argument Parser - Assigns arguments from command line'''

    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError # File not found
                setattr(args, self.dest, value)
        return customAction

    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter, add_help=False)

    parser.add_argument("-i", metavar = 'Ima2p_Input', help = "Specifies the name of input data file", type = str)
    parser.add_argument("-o", help = "Specifies the name for program output files", type = str)
    parser.add_argument("-b", help = "Specifies the duration of burn : number of burnin steps")
    parser.add_argument("-c", help = "Specifies the calculation options", type = int)
    parser.add_argument("-d", help = "Specifies the number of steps between tree saving", type = int)
    parser.add_argument("-f", help = "Specifies the name of file with saved Markov chain state generated in previous run -use with -r3", type = str)
    parser.add_argument("-g", help = "Specifies the name of file with parameter priors  (requires -c3)", type = str)
    parser.add_argument("-hf", help = "Specifies the heating model : geometric", type = str) 
    parser.add_argument("-hn", help = "Specifies the number of chains", type = int)
    parser.add_argument("-hk", help = "Specifies the number of chain swap attempts per step", type = int)
    parser.add_argument("-ha", help = "Specifies the first heating parameter, effect depends on heating model", type = float)
    parser.add_argument("-hb", help = "Specifies the second heating parameter, effect depends on heating model", type = float)
    parser.add_argument("-j", help = "Specifies the model options", type = int)
    parser.add_argument("-l", help = "Specifies the run duration", type = int)
    parser.add_argument("-m", help = "Specifies the migration prior value", type = float)
    parser.add_argument("-p", help = "Specifies the output options", type = int)
    parser.add_argument("-q", help = "Specifies the maximum for population size parameters", type = int)
    parser.add_argument("-r", help = "Specifies the run options", type = int)
    parser.add_argument("-s", help = "Specifies the random number seed", type = int)
    parser.add_argument("-t", help = "Specifies the maximum time of population splitting", type = int)
    parser.add_argument("-u", help = "Specifies the Generation time in years -for use with -p3", type = int)
    parser.add_argument("-v", help = "Specifies the base name", type = str)
    parser.add_argument("-w", help = "Specifies the name of file with nested models to be tested (LOAD-GENEALOGY mode only), invokes -c2)", type = str)
    parser.add_argument("-y", help = "Mutation rate scalar for relevant loci - for use with -p3")
    parser.add_argument("-z", help = "Specifies the number of stepsbetween screen output", type = int)
    return parser

def wrapperParser():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--threads",type=int,help="Use multi-threaded mode with this thread count")
    parser.add_argument("--ima-path",type=str,help="Path to IMa3 exe if not on system path")
    return parser

def bothParser(p1,p2):
    parser = argparse.ArgumentParser(parents=[p1,p2])
    return parser


def getIMaMT():
    for exe in ['IMa3','IMa3_xml']:
        if find_executable(exe) is not None:
            return exe
    return None

def getIMaST():
    for exe in ['IMa3_singlecpu','IMa3_xml_singlecpu']:
        if find_executable(exe) is not None:
            return exe
    return None

def getIMaExe(args):
    if args.ima_path is not None:
        if not os.path.isfile(args.ima_path):
            raise Exception("IMa exe at %s does not exist"%args.ima_path)
        return args.ima_path
    imamt = getIMaMT()
    if args.threads is not None:
        return imamt
    if imamt is not None:
        return imamt
    return getIMaST()
    


def call_ima3(args,ima_arglist):
    run_args = []
    #if find_executable('mpirun') is not None:
    if args.threads is not None:
        if find_executable('mpirun') is None:
            raise Exception('mpirun required on path for multi-threading')
        run_args.extend(['mpirun','-np',str(ima3_args.threads)])
    ima_exe = getIMaExe(args)
    if ima_exe is None:
        raise Exception("No IMa executable found on system path")
    run_args.append(ima_exe)
    run_args.extend(map(str,ima_arglist))
    ima_call = subprocess.Popen(run_args,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    ima_stdout, ima_stderr = ima_call.communicate()
    if sys.version_info[0] == 3:
        ima_stdout = ima_stdout.decode()
        ima_stderr = ima_stderr.decode()
    if ima_call.returncode != 0:
        raise Exception("Error with IMa call:%d\n%s"%(ima_call.returncode,ima_stderr))
    sys.stdout.write(ima_stdout+'\n')
    sys.stdout.write(ima_stderr+'\n')
    


def logArgs(args, pipeline_function):
    '''Logs arguments from argparse system. May be replaced with a general function in logging_module'''
    logging.info('Arguments for %s:' % pipeline_function)
    for k in vars(args):
        logging.info('Arguments %s: %s' % (k, vars(args)[k]))

def run (sysargv):
    
    ima_parser = imaParser()
    misc_parser = wrapperParser()
    main_parser = bothParser(ima_parser,misc_parser)
    if len(sysargv) == 0:
        main_parser.print_help()
        exit()
    misc_args, ima_arglist = misc_parser.parse_known_args(sysargv)

    #Check that arguments don't cause error, since list is passed to 
    #popen the returned namespace isn't used
    ima_parser.parse_known_args(ima_arglist)
    call_ima3(misc_args,ima_arglist)
    logging.info("Run successful")


if __name__ == "__main__":
    #initLogger()
    run(sys.argv[1:])


