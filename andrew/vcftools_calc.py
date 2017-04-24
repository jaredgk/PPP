import os
import sys
import subprocess
import argparse
import logging

# Insert Jared's directory path, required for calling Jared's functions. Change when directory structure changes.
#sys.path.insert(0, os.path.abspath(os.path.join(os.pardir,'jared')))

#from logging_module import initLogger

def vcf_argument_parser(passed_arguments):
    '''VCF Argument Parser - Assigns arguments for vcftools from command line.
    Depending on the argument in question, a default value may be specified
    ([--fst-window-step, 2000]) or not (i.e. [])'''
    
    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError # File not found
                setattr(args, self.dest, value)
        return customAction
    
    def parser_confirm_no_file ():
        '''Custom action to confirm file does not exist'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if os.path.isfile(value):
                    raise IOError # File found
                setattr(args, self.dest, value)
        return customAction
    
    def parser_confirm_files ():
        '''Custom action to confirm multiple file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError # File not found
                getattr(args, self.dest).append(value)                
        return customAction
        
    vcf_parser = argparse.ArgumentParser()
    
    # Input arguments. 
    vcf_parser.add_argument("vcfname", metavar='VCF_Input', help = "Input VCF filename", type = str, action = parser_confirm_file())
    vcf_parser.add_argument("--ext", help = "Format for variant file if filename doesn't contain extension")
    
    # Other file arguments. Expand as needed
    vcf_parser.add_argument('--out', help = 'Specifies the output filename', type = str, action = parser_confirm_no_file())
    vcf_parser.add_argument('--pop-file', help = 'Defines the population files for calculating specific statistics', type = str, action='append')

    
    # Statistic based arguments.
    statistic_list = ['weir-fst', 'TajimaD', 'pi', 'freq', 'het']
    vcf_parser.add_argument('--calc-statistic', metavar = '{{{0}}}'.format(', '.join(statistic_list)) , help = "Specifies the statistic to calculate", type=str, choices = statistic_list, default = 'weir-fst')
    
    # Statistic window options
    vcf_parser.add_argument('--statistic-window-size', help = 'Specifies the size of window calculations', type = int, default = 10000)
    vcf_parser.add_argument('--statistic-window-step', help = 'Specifies step size between windows', type = int, default = 20000)

          
    if passed_arguments:
        return vcf_parser.parse_args(passed_arguments)
    else:
        return vcf_parser.parse_args()

def check_vcftools_for_errors (vcftools_output):
    '''Checks the vcftools stderr for reported errors'''
    if 'Run Time' in str(vcftools_output):
        return True
    else:
        sys.exit('Error with vcftools')
    
def produce_vcftools_log (output, filename, function):
    '''Creates a log file for the vcftools run. Also reports if the log file already exits'''
    if not os.path.isfile(filename + function + '.log'):
        vcftools_log_file = open(filename + function + '.log','w')
        vcftools_log_file.write(str(output))
        vcftools_log_file.close()
    else:
        sys.exit('Error: Log file already exits')
                
def assign_log_suffix (passed_command):
    '''Assigns a suffix for the log file to reduce overwriting'''
    if '--weir-fst-pop' in passed_command and '--fst-window-size' in passed_command and '--fst-window-step' in passed_command:
        return '.windowed.weir.fst'
    elif '--weir-fst-pop' in passed_command:
        return '.weir.fst'
    elif '--TajimaD' in passed_command:
        return '.Tajima.D'
    elif '--window-pi' in passed_command:
        return '.windowed.pi'
    elif '--het' in passed_command:
        return '.het'
    elif '--freq' in passed_command or '--freq2' in passed_command:
        return '.frq'
    

def vcftools_subprocess_call (vcfname, vcftools_args):
    
    # vcftools subprocess call
    vcftools_call = subprocess.Popen(['vcftools'] + vcfname + vcftools_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    vcftools_out, vcftools_err = vcftools_call.communicate()
    
    # Check that the log file was created correctly, get the suffix for the log file, and create the file
    if check_vcftools_for_errors(vcftools_err):
        log_suffix = assign_log_suffix(vcftools_command)
        produce_vcftools_log(vcftools_err, vcf_args.out[1], log_suffix)
    
    
def run (passed_arguments = []):
    '''
    PPP script for automating VCFTools.
    
    The VCFTools wrapper was designed to simplify the use of VCFTools within the PPP.
    The function is capable of operating within the pipeline (i.e. imported & called)
    or individually. Currently, the functions uses argparse to correctly assign the
    the passed parameters to VCFTools.
        
    '''
    
    # Grab VCF arguments from command line
    vcf_args = vcf_argument_parser(passed_arguments)
    
    # Argument container for vcftools
    vcftools_call_args = []
    
    if vcf_args.calc_statistic == 'weir-fst':
        if not vcf_args.pop_file or len(vcf_args.pop_file) < 2:
            sys.exit('Two or more population files requried. Please assign using --pop-file')
        
        # Assigns specific vcftools arguments for calculating fst
        vcftools_pop_args = [population_args for population_file in vcf_args.pop_file for population_args in ['--weir-fst-pop', population_file]]
        vcftools_window_args = ['--fst-window-size', vcf_args.statistic_window_size, '--fst-window-step', vcf_args.statistic_window_step]
        # Assigns all the vcftools arguments for calculating fst
        vcftools_call_args = vcftools_pop_args + vcftools_window_args
        
    elif vcf_args.calc_statistic == 'TajimaD':
        # Assigns all the vcftools arguments for calculating TajimaD
        vcftools_call_args = ['--TajimaD', vcf_args.statistic_window_size]
        
    elif vcf_args.calc_statistic == 'pi':
        # Assigns all the vcftools arguments for calculating pi
        vcftools_call_args = ['--window-pi', vcf_args.statistic_window_size, '--window-pi-step', vcf_args.statistic_window_step]
        
    elif vcf_args.calc_statistic == 'freq':
        # Assigns all the vcftools arguments for the allele frequency
        vcftools_call_args = ['--freq']
        
    elif vcf_args.calc_statistic == 'het':
        # Assigns all the vcftools arguments for calculating heterozygosity
        vcftools_call_args = ['--het']        
    
    # Temp file assignment
    vcfname_arg = ['--gzvcf', vcf_args.vcfname]
    
    print ['vcftools'] + vcfname_arg + vcftools_call_args
    
    # vcftools subprocess call
    '''
    vcftools_call = subprocess.Popen(['vcftools'] + vcfname_arg + vcftools_call_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    vcftools_out, vcftools_err = vcftools_call.communicate()
    
    # Check that the log file was created correctly, get the suffix for the log file, and create the file
    if check_vcftools_for_errors(vcftools_err):
        log_suffix = assign_log_suffix(vcftools_command)
        produce_vcftools_log(vcftools_err, vcf_args.out[1], log_suffix)
    '''
if __name__ == "__main__":
    #initLogger()
    run()