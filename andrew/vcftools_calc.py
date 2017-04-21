import os
import sys
import subprocess
import argparse
import logging

# Insert Jared's directory path, required for calling Jared's functions. Change when directory structure changes.
sys.path.insert(0, os.path.abspath(os.path.join(os.pardir,'jared')))

from logging_module import initLogger

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
        
    vcf_parser = argparse.ArgumentParser()
    
    # Input arguments. 
    vcf_parser.add_argument("vcfname", metavar='VCF_Input', help = "Input VCF filename", type = str, action = parser_confirm_file())
    vcf_parser.add_argument("--ext", help = "Format for variant file if filename doesn't contain extension")
    
    # Other basic arguments. Expand as needed
    vcf_parser.add_argument('--out', help = 'Specifies the output filename', type = str, action = parser_confirm_no_file())
    
    # Statistic based arguments.
    statistic_list = ['weir-fst', 'TajimaD', 'pi', 'freq', 'het']
    vcf_parser.add_argument('--calc-statistic', help = "Specifies the statistic to calculate", type=str, choices = statistic_list)
    
    vcf_parser.add_argument('--statistic-window-size', help = 'Specifies the size of window calculations', type = int, default = 10000)
    vcf_parser.add_argument('--statistic-window-step', help = 'Specifies step size between windows', type = int, default = 20000)
       
    '''
    # Fst arguments
    vcf_parser.add_argument('--weir-fst-pop', help = 'Defines the population (list of individuals) files for calculating Fst', default = [], type = str, action = vcf_argument_append_file('--weir-fst-pop'))
    vcf_parser.add_argument('--fst-window-size', help = 'Defines the size of the Fst calculation windows (rather than Fst calculations per site)', nargs='?', default = [], type=int, const = 10000, action = vcf_argument_attribute('--fst-window-size'))
    vcf_parser.add_argument('--fst-window-step', help = 'Defines the step size between Fst windows', nargs='?', default = [], type=int, const = 20000, action = vcf_argument_attribute('--fst-window-step'))
    
    # Tajima's D arguments   
    vcf_parser.add_argument('--TajimaD', help = "Defines the Tajima's D bin size", nargs='?', default = [], const = 10000, type=int, action = vcf_argument_attribute('--TajimaD'))
    
    # Nucleotide Diversity (Pi) arguments
    vcf_parser.add_argument('--window-pi', help = 'Calcualtes nucleotide diversity for the defined window size', nargs='?', default = [], type=int, const = 10000, action = vcf_argument_attribute('--window-pi'))
    vcf_parser.add_argument('--window-pi-step', help = 'Calcualtes nucleotide diversity by site', nargs='?', default = [], type=int, const = 20000, action = vcf_argument_attribute('--window-pi-step'))
    
    # Allele frequency arguments. Currently mutually exclusive to only allow a single allele frequency reporting method
    vcf_allele_freq = vcf_parser.add_mutually_exclusive_group()
    vcf_allele_freq.add_argument('--freq', help = 'Outputs the allele frequency', default = [], action = 'store_const', const = ['--freq'])
    vcf_allele_freq.add_argument('--freq2', help = 'Outputs the allele frequency without information on the alleles', default = [], action = 'store_const', const = ['--freq2']) #vcf_argument_append_attribute
    
    # Heterozygosity arguments
    vcf_parser.add_argument('--het', help = 'Outputs the heterozygosity', default = [], action = 'store_const', const = ['--het'])
    '''       
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
    
    '''
         
    # Only allows a single statistic to be run at a time, may update to loop statistics
    stats_count = sum(1 for data in (vcf_args.weir_fst_pop, vcf_args.TajimaD, vcf_args.window_pi, vcf_args.het, vcf_args.freq, vcf_args.freq2) if data)
    if not stats_count:
        sys.exit('No statistic selected')
    if stats_count > 1:
        sys.exit('Statistic limit')
    
    # Checks that two or more population files are assigned if calculating Fst 
    if vcf_args.weir_fst_pop and (len(vcf_args.weir_fst_pop) / 2) < 2:
        sys.exit('Two or more population files requried. Please assign using --weir-fst-pops')
    
    # Assigns the vcftools arguments to a single command line
    vcftools_command = [split_arg for arg in vars(vcf_args) for split_arg in getattr(vcf_args, arg)]
    
    # vcftools subprocess call
    vcftools_call = subprocess.Popen(['vcftools'] + vcftools_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    vcftools_out, vcftools_err = vcftools_call.communicate()
    
    # Check that the log file was created correctly, get the suffix for the log file, and create the file
    if check_vcftools_for_errors(vcftools_err):
        log_suffix = assign_log_suffix(vcftools_command)
        produce_vcftools_log(vcftools_err, vcf_args.out[1], log_suffix)
    '''

if __name__ == "__main__":
    #initLogger()
    run()