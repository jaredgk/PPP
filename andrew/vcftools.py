import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.pardir,'jared')))

import vcf_reader_func

def check_vcftools_for_errors (vcftools_output):
    '''Checks the vcftools stderr for reported errors'''
    
    # Returns True if the job completed without error
    if 'Run Time' in str(vcftools_output):
        return True
    
    # Print output for vcftools if error is detected
    elif 'Error' in str(vcftools_output):
        # Splits log into list of lines
        vcftools_output_lines = vcftools_output.splitlines()
        # Prints the error(s)
        sys.exit('\n'.join((output_line for output_line in vcftools_output_lines if output_line.startswith('Error'))))
    
    # Print output if not completed and no error found. Unlikely to be used, but included. 
    else:
        sys.exit(vcftools_output)
        
def produce_vcftools_log (output, filename, function):
    '''Creates a log file for the vcftools run. Also reports if the log file already exits'''
    if not os.path.isfile(filename + function + '.log'):
        vcftools_log_file = open(filename + function + '.log','w')
        vcftools_log_file.write(str(output))
        vcftools_log_file.close()
    else:
        sys.exit('Error: Log file already exits')
        

def assign_vcftools_input_arg (filename):
    # True if file extensions is recognized by vcftools
    if filename.split('.', 1)[-1] in ['vcf', 'vcf.gz', 'bcf']:
        
        if filename.endswith('.vcf'):
            return ['--vcf', filename]
        elif filename.endswith('.vcf.gz'):
            return ['--gzvcf', filename]
        elif filename.endswith('.bcf'):
            return ['--bcf', filename]
    
    # True if file extension is unknown or not recognized
    else:
        
        # Checks if the file is unzipped, bgzipped, or gzipped
        vcfname_format = vcf_reader_func.checkFormat(filename)
    
        if vcfname_format == 'nozip':
            return ['--vcf', filename]
        elif vcfname_format == 'gzip':
            return ['--gzvcf', filename]
        elif vcfname_format == 'bgzip':
            return ['--bcf', filename]
        else:
            sys.exit('Unknown file format')