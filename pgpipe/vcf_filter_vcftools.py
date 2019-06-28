'''
    Automates VCF file filering using VCFtools. Filters may be used independently or
    combined to create complex filtering operations.        

    Please note, `VCFtools <https://vcftools.github.io/index.html>`_ should be cited
    alonside the PPP when using this function. Please see the log file for the relevant 
    citation.

    ############################
    Input Command-line Arguments
    ############################
    **--vcf** *<input_filename>*
        Argument used to define the filename of the VCF file to be filtered.
    **--model-file** *<model_filename>*
    	Argument used to define the model file. Please note that this argument cannot be 
    	used with the individual-based filters.
    **--model** *<model_str>*
    	Argument used to define the model (i.e. the individual(s) to include). Please note
    	that this argument cannot be used with the individual-based filters.

    #############################
    Output Command-line Arguments
    #############################
    **--out** *<output_filename>*
        Argument used to define the complete output filename, overrides **--out-prefix**
    **--out-prefix** *<output_prefix>*
        Argument used to define the output prefix (i.e. filename without file extension)
    **--out-format** *<vcf, vcf.gz, bcf, removed_sites, kept_sites, removed_bed, kept_bed>*
        Argument used to define the desired output format. Formats include: uncompressed 
        VCF (vcf); compressed VCF (vcf.gz) [default]; BCF (bcf); removed or kept variants
        in sites format (removed_sites, kept_sites); or removed or kept variants in bed
        format (removed_bed, kept_bed).
    **--overwrite**
        Argument used to define if previous output should be overwritten.

    #############################
    Filter Command-line Arguments
    #############################
    The filtering arguments below are roughly seperated into catagoires. Please not that
    mulitple filters are seperated into two opposing function types **include** and 
    **exclude**.

    **************************
    Individual-Based Arguments
    **************************
    Please note that all individual-based arguments are not compatible with either the 
    **--model** or **--model-file** command-line arguments.

    **--filter-include-indv** *<indv_str>* *<indv1_str, indv2_str, etc.>*
        Argument used to define the individual(s) to include. This argument may be used 
        multiple times if desired.
    **--filter-exclude-indv** *<indv_str>* *<indv1_str, indv2_str, etc.>*
        Argument used to define the individual(s) to exclude. This argument may be used 
        multiple times if desired.
    **--filter-include-indv-file** *<indv_filename>*
        Argument used to define a file of individuals to include.
    **--filter-exclude-indv-file** *<indv_filename>*
        Argument used to define a file of individuals to exclude.

    *******************************
    Allele/Genotype-Based Arguments
    *******************************
    **--filter-min-alleles** *min_int*
    	Argument used to include positions with a number of allele >= to the 
    	given number.
    **--filter-max-alleles** *max_int*
        Argument used to include positions with a number of allele <= to the 
        given number.
    **--filter-maf-min** *maf_proportion*
    	Argument used to include sites with equal or greater MAF values.
    **--filter-maf-max** *maf_proportion*
    	Argument used to include sites with equal or lesser MAF values.
    **--filter-mac-min** *mac_int*
    	Argument used to include sites with equal or greater MAC values.
    **--filter-mac-max** *mac_int*
    	Argument used to include sites with equal or lesser MAC values.
    	**--filter-include-indels**
    	Argument used to include positions if they contain an insertion or a deletion. 
    **--filter-exclude-indels**
    	Argument used to exclude positions if they contain an insertion or a deletion. 
    **--filter-include-snp** *<rs#>* *<rs#1, rs#2, etc.>*
    	Argument used to include SNP(s) with the matching ID. This argument may be used 
        multiple times if desired.
    **--filter-include-snp-file** *<snp_filename>*
    	Argument used to define a file of SNP IDs to include.
    **--filter-exclude-snp-file** *<snp_filename>*
    	Argument used to define a file of SNP IDs to exclude.
    **--filter-max-missing** *proportion_float*
    	Argument used to filter positions by their proportion of missing data, a value of 
    	0.0 will ignore missing data whereas a value of 1.0 will allow no missing data.

    ************************
    Position-Based Arguments
    ************************
    **--filter-include-chr** *<chr>* *<chr1, chr2, etc.>*
        Argument used to define the chromosome(s) to include. This argument may be used 
        multiple times if desired.
    **--filter-exclude-chr** *<chr>* *<chr1, chr2, etc.>*
        Argument used to define the chromosome(s) to exclude. This argument may be used 
        multiple times if desired.
    **--filter-from-bp** *<bp_int>*
        Argument used to define the lower bound of positions to include. May only be used 
        with a single chromosome.
    **--filter-to-bp** *<bp_int>*
        Argument used to define the upper bound of positions to include. May only be used 
        with a single chromosome.
    **--filter-distance** *<bp_int>*
    	Argument used to define the distance that no two sites may be within
    **--filter-include-positions** *<position_filename>*
        Argument used to define a file of positions to include within a tsv file 
        (chromosome and position).
    **--filter-exclude-positions** *<position_filename>*
        Argument used to define a file of positions to exclude within a tsv file 
        (chromosome and position).
    **--filter-include-bed** *<position_bed_filename>*
        Argument used to define a BED file of positions to include.
    **--filter-exclude-bed** *<position_bed_filename>*
        Argument used to define a BED file of positions to exclude.

    ********************
    Flag-Based Arguments
    ********************
    **--filter-include-passed**
        Argument used to include positions with the 'PASS' filter flag.
    **--filter-include-filtered** *<filter_flag>*
        Argument used to include positions with the specified filter flag.
    **--filter-exclude-filtered** *<filter_flag>*
        Argument used to exclude positions with the specified filter flag.
    **--filter-include-info** *<info_flag>*
        Argument used to include positions with the specified info flag.
    **--filter-exclude-info** *<info_flag>*
        Argument used to exclude positions with the specified info flag.

    ##########################
    Example Command-line Usage
    ##########################
    Command-line to only include biallelic sites:

    .. code-block:: bash
        
        python vcf_filter.py --vcf input.bcf --filter-min-alleles 2 --filter-max-alleles 2

    Command-line to only include chromosome chr1:

    .. code-block:: bash
        
        python vcf_filter.py --vcf input.vcf.gz --filter-include-chr chr1

    Command-line to remove indels and ouput a BCF file:

    .. code-block:: bash
        
        python vcf_filter.py --vcf input.vcf --filter-exclude-indels --out-format bcf
'''

import os
import sys
import subprocess
import argparse
import logging

# Import basic vcftools functions
from pgpipe.vcftools import *

# Model file related functions
from pgpipe.model import read_model_file

# Import logging functions
from pgpipe.logging_module import initLogger, logArgs

def vcf_filter_parser(passed_arguments):
    '''
    VCF Filter Argument Parser

    Assign the parameters for the VCF Filter using argparse.

    Parameters
    ----------
	passed_arguments : list, optional
		Parameters passed by another function. sys.argv is used if
		not given. 

	Raises
    ------
    IOError
        If the input, or other specified files do not exist
    '''

    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError('%s not found' % value)
                setattr(args, self.dest, value)
        return customAction

    def parser_add_to_list ():
        '''Custom action to add items to a list'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):

                # Clean up any commas
                value = [item.strip(',') for item in value]
                
                if not getattr(args, self.dest):
                    setattr(args, self.dest, value)
                else:
                    getattr(args, self.dest).extend(value)
        return customAction

    def metavar_list (var_list):
        '''Create a formmated metavar list for the help output'''
        return '{' + ', '.join(var_list) + '}'

    vcf_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    # Input arguments.
    vcf_parser.add_argument('--vcf', help = "Defines the filename of the VCF", type = str, required = True, action = parser_confirm_file())

    # Model file arguments.
    vcf_parser.add_argument('--model-file', help = 'Defines the model file', type = str, action = parser_confirm_file())
    vcf_parser.add_argument('--model', help = 'Defines the model and the individual(s) to include', type = str)

    # Other file arguments. Expand as needed
    vcf_parser.add_argument('--out', help = 'Defines the complete output filename, overrides --out-prefix', type = str)
    vcf_parser.add_argument('--out-prefix', help = 'Defines the output prefix (i.e. filename without file extension)', type = str,  default = 'out')
    #vcf_parser.add_argument('--log', help = "Specifies if the vcftools log should be saved", action = 'store_false')

    out_format_list = ['vcf', 'vcf.gz', 'bcf', 'removed_sites', 'kept_sites', 'removed_bed', 'kept_bed']
    out_format_default = 'vcf.gz'

    vcf_parser.add_argument('--out-format', metavar = metavar_list(out_format_list), help = 'Defines the desired output format', type = str, choices = out_format_list, default = out_format_default)

    # General arguments.
    vcf_parser.add_argument('--overwrite', help = "Overwrite previous output", action = 'store_true')

    # Galaxy Option to pipe log to stdout
    vcf_parser.add_argument('--log-stdout', help = argparse.SUPPRESS, action = 'store_true')

    ### Filters

    # Non-model file arguments
    vcf_parser.add_argument('--filter-include-indv', help = 'Defines the individual(s) to include. This argument may be used multiple times if desired', nargs = '+', type = str, action = parser_add_to_list())
    vcf_parser.add_argument('--filter-exclude-indv', help = 'Defines the individual(s) to exclude. This argument may be used multiple times if desired', nargs = '+', type = str, action = parser_add_to_list())
    vcf_parser.add_argument('--filter-include-indv-file', help = 'Defines a file of individuals to include', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-indv-file', help = 'Defines a file of individuals to exclude', action = parser_confirm_file())

    # Allele count filters
    vcf_parser.add_argument('--filter-min-alleles', help = 'Include positions with a number of allele >= to the given number', type = int)
    vcf_parser.add_argument('--filter-max-alleles', help = 'Include positions with a number of allele <= to the given number', type = int)

    # Missing data filter
    vcf_parser.add_argument('--filter-max-missing', help = 'Filter positions by their proportion of missing data (0.0: include all, 1.0: no missing data)', type = float)

    # Indel filters
    indel_filters = vcf_parser.add_mutually_exclusive_group()
    indel_filters.add_argument('--filter-include-indels', help = 'Include positions if they contain an insertion or a deletion', action = 'store_true')
    indel_filters.add_argument('--filter-exclude-indels', help = 'Exclude positions if they contain an insertion or a deletion', action = 'store_true')

    # Chromosome filters
    vcf_parser.add_argument('--filter-include-chr', help = 'The chromosome(s) to include. This argument may be used multiple times', nargs = '+', type = str, action = parser_add_to_list())
    vcf_parser.add_argument('--filter-exclude-chr', help = 'The chromosome(s) to exclude. This argument may be used multiple times', nargs = '+', type = str, action = parser_add_to_list())

    # Basic position filters
    vcf_parser.add_argument('--filter-from-bp', help = 'Defines the lower bound of positions to include. May only be used with a single chromosome', type = int)
    vcf_parser.add_argument('--filter-to-bp', help = 'Defines the upper bound of positions to include. May only be used with a single chromosome', type = int)

    # Position-based position filters
    vcf_parser.add_argument('--filter-include-positions', help = 'Defines a file of positions to include within a tsv file (chromosome and position)', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-positions', help = 'Defines a file of positions to exclude within a tsv file (chromosome and position)', action = parser_confirm_file())

    # BED-based position filters
    vcf_parser.add_argument('--filter-include-bed', help = 'Defines a BED file of positions to include', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-bed', help = 'Defines a BED file of positions to exclude', action = parser_confirm_file())

    # Filter-flag filters
    vcf_parser.add_argument('--filter-include-passed', help = "Include positions with the 'PASS' filter flag", action = 'store_true')
    vcf_parser.add_argument('--filter-include-flag', help = 'Include positions with the specified filter flag', nargs = '+', type = str, action = parser_add_to_list())
    vcf_parser.add_argument('--filter-exclude-flag', help = 'Exclude positions with the specified filter flag', nargs = '+', type = str, action = parser_add_to_list())

    # Info-flag filters
    vcf_parser.add_argument('--filter-include-info', help = 'Include positions with the specified info flag', nargs = '+', type = str, action = parser_add_to_list())
    vcf_parser.add_argument('--filter-exclude-info', help = 'Exclude positions with the specified info flag', nargs = '+', type = str, action = parser_add_to_list())

    # SNP filters
    vcf_parser.add_argument('--filter-include-snp', help = 'Include SNP(s) with the matching ID. This argument may be used multiple times.', nargs = '+', type = str, action = parser_add_to_list())
    vcf_parser.add_argument('--filter-include-snp-file', help = 'Defines a file of SNP IDs to include', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-snp-file', help = 'Defines a file of SNP IDs to exclude', action = parser_confirm_file())

    # MAF Filter
    vcf_parser.add_argument('--filter-maf-min', help = 'Include sites with equal or greater MAF values', type = float)
    vcf_parser.add_argument('--filter-maf-max', help = 'Include sites with equal or lesser MAF values', type = float)

    # MAC Filter
    vcf_parser.add_argument('--filter-mac-min', help = 'Include sites with equal or greater MAC values', type = int)
    vcf_parser.add_argument('--filter-mac-max', help = 'Include sites with equal or lesser MAC values', type = int)

    # Additional Filters
    vcf_parser.add_argument('--filter-distance', help = 'Distance that no two sites may be within', type = int)

    if passed_arguments:
        return vcf_parser.parse_args(passed_arguments)
    else:
        return vcf_parser.parse_args()

def run (passed_arguments = []):
    '''
    Filter function for the PPP

    This function uses the argparse-based function :py:func:`vcf_filter_parser`
    to parse either sys.argv or passed_arguments to obtain the parameters below. 
    The parameters are then translated to their VCFtools equivalent. Once all the
    VCFtools-based parameters are assigned, VCFtools is called.

    Parameters
    ----------
    --vcf : str
        Input VCF filename
    --out : str
        Complete output filename, overrides --out-prefix
    --out-prefix : str
        Output filename prefix
    --out-format : str
        Output format
	--model-file : str
		Model filename
	--model : str
		Model to use
    --filter-include-indv : list or str
        Individual(s) to include. May be used multiple times. Not usable w/ --model
    --filter-exclude-indv : list or str
        Individual(s) to exclude. May be used multiple times. Not usable w/ --model
    --filter-include-indv-file : str
        File of individuals to include. Not usable w/ --model
    --filter-exclude-indv-file : str
        File of individuals to exclude.  Not usable w/ --model
    --filter-min-alleles : int
        Include positions with a number of allele >= to the given number
    --filter-max-alleles : int
        Include positions with a number of allele <= to the given number
    --filter-maf-min
		Include sites with equal or greater MAF values
    --filter-maf-max
    	Include sites with equal or lesser MAF values
    --filter-mac-min
    	Include sites with equal or greater MAC values
    --filter-mac-max
    	Include sites with equal or lesser MAC values
   	--filter-include-indels
    	Include positions if they contain an insertion or a deletion
    --filter-exclude-indels
    	Exclude positions if they contain an insertion or a deletion
 	--filter-include-snp
 		Include SNP(s) with the matching ID. May be used multiple times
    --filter-include-snp-file
    	File of SNP IDs to include
    --filter-exclude-snp-file
    	File of SNP IDs to exclude
    --filter-max-missing
		Filter positions by their proportion of missing data. (0.0: include all, 1.0: no missing data)
    --filter-include-chr
        Chromosome(s) to include. May be used multiple times.
    --filter-exclude-chr
        Chromosome(s) to exclude. May be used multiple times.
    --filter-from-bp : int
        Lower bound of positions to include. Only usable with a single chromosome
    --filter-to-bp : int
        Upper bound of positions to include. Only usable with a single chromosome
    --filter-include-positions : str
        File of positions to include (tsv: chromosome and position)
    --filter-exclude-positions : str
        File of positions to exclude (tsv: chromosome and position)
    --filter-include-bed : str
        BED file of positions to include 
    --filter-exclude-bed : str
        BED file of positions to exclude
    --filter-include-passed : bool
        Include positions with the 'PASS' filter flag
    --filter-include-filtered : list or str
        Include positions with the specified filter flag
    --filter-exclude-filtered : list or str
        Exclude positions with the specified filter flag
    --filter-include-info : list or str
        Include positions with the specified info flag
    --filter-exclude-info : list or str
        Exclude positions with the specified info flag

    Raises
    ------
    IOError
        Output file already exists and --overwrite is not specified
    Exception
        Incompatible arguments

    '''

    # Grab VCF arguments from command line
    vcf_args = vcf_filter_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(vcf_args, func_name = 'vcf_filter')

    # Argument container for vcftools
    vcftools_call_args = []

    # Check if the user has specified a model file
    if vcf_args.model_file and vcf_args.model:
        # Read in the models
        models_in_file = read_model_file(vcf_args.model_file)

        # Check that the selected model was not found in the file
        if vcf_args.model not in models_in_file:
            raise IOError('Selected model "%s" not found in: %s' % (vcf_args.model, vcf_args.model_file))

        # Check that individual-based filters are not also being used
        if vcf_args.filter_include_indv or vcf_args.filter_exclude_indv:
            if vcf_args.filter_include_indv:
                raise Exception('--model and --filter-include-indv arguments are incompatible')
            if vcf_args.filter_exclude_indv:
                raise Exception('--model and --filter-exclude-indv arguments are incompatible')

        # Check that individuals-based filters are not also being used
        if vcf_args.filter_include_indv_file or vcf_args.filter_exclude_indv_file:
            if vcf_args.filter_include_indv_file:
                raise Exception('--model and --filter_include_indv-file arguments are incompatible')
            if vcf_args.filter_exclude_indv_file:
                raise Exception('--model and --filter-exclude-indv-file arguments are incompatible')

        # Select model, might change this in future versions
        selected_model = models_in_file[vcf_args.model]

        # Create individuals file
        selected_model.create_ind_file(overwrite = vcf_args.overwrite)

        # Assign the individuals file to vcftools
        vcftools_call_args.extend(['--keep', selected_model.ind_file])

    # Holds the filename suffix vcftools assigns to the filtered output
    vcftools_out_suffix = None

    # Used to assign the output format to the vcftools call, and assign vcftools output suffix
    if vcf_args.out_format == 'removed_sites':
        vcftools_call_args.append('--removed-sites')
        vcftools_out_suffix = '.removed.sites'
    elif vcf_args.out_format == 'kept_sites':
        vcftools_call_args.append('--kept-sites')
        vcftools_out_suffix = '.kept.sites'
    elif vcf_args.out_format == 'removed_bed':
        vcftools_call_args.append('--removed-sites')
        vcftools_out_suffix = '.removed.bed'
    elif vcf_args.out_format == 'kept_bed':
        vcftools_call_args.append('--kept-sites')
        vcftools_out_suffix = '.kept.bed'
    elif vcf_args.out_format == 'bcf':
        vcftools_call_args.append('--recode')
        vcftools_out_suffix = '.recode.bcf'
    elif vcf_args.out_format == 'vcf':
        vcftools_call_args.append('--recode')
        vcftools_out_suffix = '.recode.vcf'
    elif vcf_args.out_format == 'vcf.gz':
        vcftools_call_args.append('--recode')
        vcftools_out_suffix = '.recode.vcf.gz'

    # Assign expected vcftools output filename
    vcftools_output_filename = vcf_args.out_prefix + vcftools_out_suffix

    # Check if the user has specified a output filename
    if vcf_args.out:

        # Rename the vcftools output filename, using the output filename
        vcftools_output_filename = vcf_args.out

    # Check if previous output should be overwritten
    if not vcf_args.overwrite:

        # Confirm the vcftools output and log file do not exist
        check_for_vcftools_output(vcftools_output_filename)

    # Individuals-based filters
    if vcf_args.filter_include_indv_file or vcf_args.filter_exclude_indv_file:
        # Used to include a file of individuals to keep
        if vcf_args.filter_include_indv_file:
            vcftools_call_args.extend(['--keep', vcf_args.filter_include_indv_file])

        # Used to include a file of individuals to remove
        if vcf_args.filter_exclude_indv_file:
            vcftools_call_args.extend(['--remove', vcf_args.filter_exclude_indv_file])

    # Individual-based filters
    if vcf_args.filter_include_indv or vcf_args.filter_exclude_indv:
        if vcf_args.filter_include_indv:
            for indv_to_include in vcf_args.filter_include_indv:
                vcftools_call_args.extend(['--indv', indv_to_include])
        if vcf_args.filter_exclude_indv:
            for indv_to_exclude in vcf_args.filter_exclude_indv:
                vcftools_call_args.extend(['--remove-indv', indv_to_exclude])

    # Indel-based filters
    if vcf_args.filter_include_indels or vcf_args.filter_exclude_indels:
        if vcf_args.filter_include_indels:
            vcftools_call_args.append('--keep-only-indels')
        if vcf_args.filter_exclude_indels:
            vcftools_call_args.append('--remove-indels')

    # Chromosome-based filters
    if vcf_args.filter_include_chr or vcf_args.filter_exclude_chr:
        if vcf_args.filter_include_chr:
            for chr_to_include in vcf_args.filter_include_chr:
                vcftools_call_args.extend(['--chr', chr_to_include])
        if vcf_args.filter_exclude_chr:
            for chr_to_exclude in vcf_args.filter_exclude_chr:
                vcftools_call_args.extend(['--not-chr', chr_to_exclude])

    # Site (i.e. basepair) filters
    if vcf_args.filter_from_bp or vcf_args.filter_to_bp:
        if vcf_args.filter_include_chr:
            vcftools_call_args.extend(['--from-bp', vcf_args.filter_from_bp])
        if vcf_args.filter_exclude_chr:
            vcftools_call_args.extend(['--to-bp', vcf_args.filter_to_bp])

    # Position (vcftools output file) filters
    if vcf_args.filter_include_positions or vcf_args.filter_exclude_positions:
        if vcf_args.filter_include_positions:
            vcftools_call_args.extend(['--positions', vcf_args.filter_include_positions])
        if vcf_args.filter_exclude_positions:
            vcftools_call_args.extend(['--exclude-positions', vcf_args.filter_exclude_positions])

    # Position (BED format file) filters
    if vcf_args.filter_include_bed or vcf_args.filter_exclude_bed:
        if vcf_args.filter_include_bed:
            vcftools_call_args.extend(['--bed', vcf_args.filter_include_bed])
        if vcf_args.filter_exclude_bed:
            vcftools_call_args.extend(['--exclude-bed', vcf_args.filter_exclude_bed])

    # Flag-based Filters
    if vcf_args.filter_include_passed or vcf_args.filter_include_flag or vcf_args.filter_exclude_flag:
        if vcf_args.filter_include_passed:
            vcftools_call_args.append('--remove-filtered-all')
        if vcf_args.filter_include_flag:
            for filtered_to_include in vcf_args.filter_include_flag:
                vcftools_call_args.extend(['--keep-filtered', filtered_to_include])
        if vcf_args.filter_exclude_flag:
            for filtered_to_exclude in vcf_args.filter_exclude_flag:
                vcftools_call_args.extend(['--remove-filtered', filtered_to_exclude])

    # Infor-based filters
    if vcf_args.filter_include_info or vcf_args.filter_exclude_info:
        if vcf_args.filter_include_info:
            for info_to_include in vcf_args.filter_include_info:
                vcftools_call_args.extend(['--keep-INFO', info_to_include])
        if vcf_args.filter_exclude_info:
            for info_to_exclude in vcf_args.filter_exclude_info:
                vcftools_call_args.extend(['--remove-INFO', info_to_exclude])

    # Allele-based filters
    if vcf_args.filter_min_alleles or vcf_args.filter_max_alleles:
        if vcf_args.filter_min_alleles:
            vcftools_call_args.extend(['--min-alleles', vcf_args.filter_min_alleles])
        if vcf_args.filter_max_alleles:
            vcftools_call_args.extend(['--max-alleles', vcf_args.filter_max_alleles])

    # Missing data filters
    if vcf_args.filter_max_missing:
        vcftools_call_args.extend(['--max-missing', vcf_args.filter_max_missing])

    # SNP-based filters
    if vcf_args.filter_include_snp or vcf_args.filter_include_snp_file or vcf_args.filter_exclude_snp_file:
        if vcf_args.filter_include_snp:
            for snp_to_include in vcf_args.filter_include_snp:
                vcftools_call_args.extend(['--snp', snp_to_include])
        if vcf_args.filter_include_snp_file:
            vcftools_call_args.extend(['--snps', vcf_args.filter_include_snp_file])
        if vcf_args.filter_exclude_snp_file:
            vcftools_call_args.extend(['--exclude', vcf_args.filter_exclude_snp_file])

    # MAF-based filters
    if vcf_args.filter_maf_min or vcf_args.filter_maf_max:
        if vcf_args.filter_maf_min:
            vcftools_call_args.extend(['--maf', vcf_args.filter_maf_min])
        if vcf_args.filter_maf_max:
            vcftools_call_args.extend(['--max-maf', vcf_args.filter_maf_max])

    # MAC-based filters
    if vcf_args.filter_mac_min or vcf_args.filter_mac_max:
        if vcf_args.filter_mac_min:
            vcftools_call_args.extend(['--mac', vcf_args.filter_mac_min])
        if vcf_args.filter_mac_max:
            vcftools_call_args.extend(['--max-mac', vcf_args.filter_mac_max])

    # Distance (between sites) filters
    if vcf_args.filter_distance:
        vcftools_call_args.extend(['--thin', vcf_args.filter_distance])

    logging.info('vcftools parameters assigned')

    # Assigns the file argument for vcftools
    vcfname_arg = assign_vcftools_input_arg(vcf_args.vcf)
    logging.info('Input file assigned')

    # Check if the output format is vcf
    if vcf_args.out_format == 'vcf':

        # Call vcftools with the specifed arguments
        vcftools_err = call_vcftools(vcfname_arg + vcftools_call_args, output_format = vcf_args.out_format, output_filename = vcftools_output_filename)

    else:

        # Call vcftools with the specifed arguments
        vcftools_err = call_vcftools(vcfname_arg + vcftools_call_args, output_format = vcf_args.out_format, output_filename = vcftools_output_filename)

    # Check if the log should be piped to the stdout
    if vcf_args.log_stdout:

        # Write the log to stdout
        sys.stdout.write(vcftools_err)

    # Check if log should be saved as a file
    else:

        # Create the log file
        produce_vcftools_log(vcftools_err, vcftools_output_filename)

    # Delete any files that were created for vcftools
    if vcf_args.model_file:
        selected_model.delete_ind_file()

if __name__ == "__main__":
    initLogger()
    run()
