import os, sys, subprocess, argparse, types

def vcf_argument_parser(passed_arguments):
    '''VCF Argument Parser - Assigns arguments for vcftools from command line.
    Depending on the argument in question, a default value may be specified
    ([--fst-window-step, 2000]) or not (i.e. [])'''
    
    
    def vcf_argument_file (vcf_argument):
        '''Custom action to transform input files into a list with it's
        respective argument prefix.'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not isinstance(value, list):
                    if os.path.isfile(value):
                        setattr(args, self.dest, [vcf_argument, value])
                    else:
                        setattr(args, self.dest, [])
                else:
                    setattr(args, self.dest, [current_argument_element for value_in_list in value for current_argument_element in [vcf_argument, value_in_list] if os.path.isfile(value_in_list)])
        return customAction
    
    def vcf_argument_attribute (vcf_argument):
        '''Custom action to transform input attributes into a list with it's
        respective argument prefix.'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                setattr(args, self.dest, [vcf_argument, str(value)])
        return customAction
    
    vcf_parser = argparse.ArgumentParser()
       
    # Input arguments. Currently mutually exclusive to only allow a single input type
    vcf_input = vcf_parser.add_mutually_exclusive_group(required=True)
    vcf_input.add_argument('--vcf', dest = 'input', help='Defines the VCF file to be processed', default = [], action=vcf_argument_file('--vcf'))
    vcf_input.add_argument('--gzvcf', dest = 'input', help='Defines the compressed (gzipped) VCF file to be processed', default = [], action=vcf_argument_file('--gzvcf'))
    vcf_input.add_argument('--bcf', dest = 'input', help='Defines the BCF file to be processed', default = [], action=vcf_argument_file('--bcf'))
    
    # Other basic arguments. Expand as needed
    vcf_parser.add_argument('--out', help='Defines the output filename', default = ['--out', 'out'], action=vcf_argument_attribute('--out'))
    
    # Fst arguments
    vcf_parser.add_argument('--weir-fst-pops', help='Defines the population (list of individuals) files for calculating Fst', default = [], nargs='+', action=vcf_argument_file('--weir-fst-pop'))
    vcf_parser.add_argument('--fst-window-size', help='Defines the size of the Fst calculation windows (rather than Fst calculations per site)', default = ['--fst-window-size', '10000'], action=vcf_argument_attribute('--fst-window-size'))
    vcf_parser.add_argument('--fst-window-step', help='Defines the step size between Fst windows', default = ['--fst-window-step', '20000'], action=vcf_argument_attribute('--fst-window-step'))
    
    # Tajima's D arguments
    vcf_parser.add_argument('--TajimaD', help="Defines the Tajima's D bin size", type = int, default = ['--TajimaD', '10000'], action=vcf_argument_attribute('--TajimaD'))
    
    # Nucleotide Diversity (Pi) arguments
    vcf_parser.add_argument('--window-pi', help="Calcualtes nucleotide diversity for the defined window size", type = int, default = ['--window-pi', '10000'], action=vcf_argument_attribute('--window-pi'))
    vcf_parser.add_argument('--window-pi-step', help="Calcualtes nucleotide diversity by site", type = int, default = ['--window-pi-step', '20000'], action=vcf_argument_attribute('--window-pi-step'))
    
    # Allele frequency arguments. Currently mutually exclusive to only allow a single allele frequency reporting method
    vcf_allele_freq = vcf_parser.add_mutually_exclusive_group()
    vcf_allele_freq.add_argument('--freq', dest = 'allele_freq', help="Outputs the allele frequency", action = 'store_const', const = ['--freq'])
    vcf_allele_freq.add_argument('--freq2', dest = 'allele_freq', help="Outputs the allele frequency without information on the alleles", action = 'store_const', const = ['--freq2'])
    
    # Heterozygosity arguments
    vcf_parser.add_argument('--het', help="Outputs the heterozygosity", action = 'store_const', const = ['--het'], default = ['--het'])
    
    # Sets the --freq argument as default for allele frequency
    vcf_parser.set_defaults(allele_freq = ['--freq'])
    
    if passed_arguments:
        return vcf_parser.parse_args(passed_arguments)
    else:
        return vcf_parser.parse_args()

def check_vcftools_for_errors (vcftools_output):
    '''Checks the vcftools stderr for reported errors'''
    if 'Run Time' in vcftools_output:
        return True
    else:
        print vcftools_output
        sys.exit('Error with vcftools')
    
def produce_vcftools_log (output, filename, function):
    '''Creates a log file for the vcftools run. Also reports if the log file already exits'''
    if not os.path.isfile(filename + function + '.log'):
        vcftools_log_file = open(filename + function + '.log','w')
        vcftools_log_file.write(output)
        vcftools_log_file.close()
    else:
        sys.exit('Error: Log file already exits')
                
def return_basic_args (vcf_arguments):
    '''Returns a list with the general arguments for VCFtools. Will be expanded as needed.'''
    return ['vcftools'] + vcf_arguments.input + vcf_arguments.out

def calculate_Fst (passed_arguments = []):
    '''Calculates Fst using VCFTools'''
    
    # Grab VCF arguments from command line
    vcf_args = vcf_argument_parser(passed_arguments)
    
    if (len(vcf_args.weir_fst_pops) / 2) < 2:
        sys.exit('Two or more population files requried. Please assign using --weir-fst-pops')
    
    vcftools_command = return_basic_args(vcf_args) + vcf_args.weir_fst_pops + vcf_args.fst_window_size + vcf_args.fst_window_step
    vcftools_Fst_call = subprocess.Popen(vcftools_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    vcftools_Fst_out, vcftools_Fst_err = vcftools_Fst_call.communicate()

    if check_vcftools_for_errors(vcftools_Fst_err):
        produce_vcftools_log(vcftools_Fst_err, vcf_args.out[1], '.windowed.weir.fst')
  
def calculate_tajimasD (passed_arguments = []):
    '''Calculates Tajima's D using VCFTools'''
    
    # Grab VCF arguments from command line
    vcf_args = vcf_argument_parser(passed_arguments)
    
    vcftools_command = return_basic_args(vcf_args) + vcf_args.TajimaD
    vcftools_tajimasD_call = subprocess.Popen(vcftools_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    vcftools_tajimasD_out, vcftools_tajimasD_err = vcftools_tajimasD_call.communicate()
    if check_vcftools_for_errors(vcftools_tajimasD_err):
        produce_vcftools_log(vcftools_tajimasD_err, vcf_args.out[1], '.Tajima.D')  
        
def calculate_pi (passed_arguments = []):
    '''Calculates nucleotide diversity (pi) using VCFTools'''
    
    # Grab VCF arguments from command line
    vcf_args = vcf_argument_parser(passed_arguments)
    
    vcftools_command = return_basic_args(vcf_args) + vcf_args.window_pi + vcf_args.window_pi_step
    vcftools_pi_call = subprocess.Popen(vcftools_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    vcftools_pi_out, vcftools_pi_err = vcftools_pi_call.communicate()
    if check_vcftools_for_errors(vcftools_pi_err):
        produce_vcftools_log(vcftools_pi_err, vcf_args.out[1], '.windowed.pi')
        
def calculate_af (passed_arguments = []):
    '''Calculates the allele frequency using VCFTools. Currently works with --freq as default
    but the user can specify --freq2 if that is prefered.'''
    
    # Grab VCF arguments from command line
    vcf_args = vcf_argument_parser(passed_arguments)
      
    vcftools_command = return_basic_args(vcf_args) + vcf_args.allele_freq
    vcftools_af_call = subprocess.Popen(vcftools_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    vcftools_af_out, vcftools_af_err = vcftools_af_call.communicate()
    
    if check_vcftools_for_errors(vcftools_af_err):
        produce_vcftools_log(vcftools_af_err, vcf_args.out[1], '.frq')
        
def calculate_heterozygosity (passed_arguments = []):
    '''Calculates the heterozygosity using VCFTools.'''
    
    # Grab VCF arguments from command line
    vcf_args = vcf_argument_parser(passed_arguments)
       
    vcftools_command = return_basic_args(vcf_args) + vcf_args.het
    vcftools_het_call = subprocess.Popen(vcftools_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    vcftools_het_out, vcftools_het_err = vcftools_het_call.communicate()
    if check_vcftools_for_errors(vcftools_het_err):
        produce_vcftools_log(vcftools_het_err, vcf_args.out[1], '.het')
        
if __name__ == "__main__":
    cmd_parser = argparse.ArgumentParser()
    cmd_choices = ['calculate_Fst', 'calculate_tajimasD', 'calculate_pi', 'calculate_af', 'calculate_heterozygosity']
    cmd_parser.add_argument('command', help='Defines the command to be invoked', choices = cmd_choices)
    cmd_args = cmd_parser.parse_args()
    
    locals()[cmd_args.command]()