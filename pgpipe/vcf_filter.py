#!/usr/bin/env python
'''
    Depending on the analysis being conducted, a number of variant sites and/or 
    samples may be unsuitable and must be removed. Given an unfiltered VCF and 
    the desired filters, vcf_filter will apply the filters and produce a filtered 
    VCF. Filters may be used independently or combined as needed. In addition, a
    number of the filters are seperated into two types: include (include/keep 
    all relevant variant sites or samples) and exclude (exclude/remove all 
    relevant variant sites or samples).

    .. image:: ../../PPP_assets/PPP_Filter.png
        :width: 100 %
        :align: center

    In this illustration of the filtering process (within a locus of interest), variant 
    sites were kept only if they: i) were biallelic and ii) passed all filters. These 
    requirements resulted in the removal of two variant sites (i.e. 197557 and 198510) 
    within the given locus.

    ##################
    Command-line Usage
    ##################
    The VCF file filter may be called using the following command:

    .. code-block:: bash
        
        vcf_filter.py

    *************
    Example usage
    *************
    Command-line to create a BCF with only biallelic sites:

    .. code-block:: bash
        
        vcf_filter.py --vcf examples/files/merged_chr1_10000.vcf.gz --filter-only-biallelic --out-format bcf

    Command-line to only include variants on chr1 from 1 to 1509546:

    .. code-block:: bash
        
        vcf_filter.py --vcf examples/files/merged_chr1_10000.bcf --filter-include-pos chr1:1-1509546

    Command-line to remove indels and ouput a BCF file:

    .. code-block:: bash
        
        vcf_filter.py --vcf examples/files/merged_chr1_10000.indels.vcf.gz --filter-exclude-indels --out-format bcf

    ############
    Dependencies 
    ############
    * `BCFtools <https://samtools.github.io/bcftools/bcftools.html>`_    

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
    **--out-format** *<vcf, vcf.gz, bcf, bed, sites>*
        Argument used to define the desired output format. Formats include: uncompressed 
        VCF (vcf); compressed VCF (vcf.gz) [default]; BCF (bcf); variants in bed format; 
        or variants in sites format.
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
    **--filter-only-biallelic**
        Argument used to only include variants that are biallelic.
    **--filter-min-alleles** *min_int*
        Argument used to include variants with a number of allele >= to the 
        given number.
    **--filter-max-alleles** *max_int*
        Argument used to include variants with a number of allele <= to the 
        given number.
    **--filter-maf-min** *maf_proportion*
        Argument used to include variants with equal or greater MAF values.
    **--filter-maf-max** *maf_proportion*
        Argument used to include variants with equal or lesser MAF values.
    **--filter-mac-min** *mac_int*
        Argument used to include variants with equal or greater MAC values.
    **--filter-mac-max** *mac_int*
        Argument used to include variants with equal or lesser MAC values.
    **--filter-include-indels**
        Argument used to include variants if they contain an insertion or a deletion. 
    **--filter-exclude-indels**
        Argument used to exclude variants if they contain an insertion or a deletion.
    **--filter-include-snps**
        Argument used to include variants if they contain a SNP. 
    **--filter-exclude-snps**
        Argument used to exclude variants if they contain a SNP. 
    **--filter-include-snp** *<rs#>* *<rs#1, rs#2, etc.>*
        Argument used to include SNP(s) with the matching ID. This argument may be used 
        multiple times if desired.
    **--filter-exclude-snp** *<rs#>* *<rs#1, rs#2, etc.>*
        Argument used to exclude SNP(s) with the matching ID. This argument may be used 
        multiple times if desired.
    **--filter-include-snp-file** *<snp_filename>*
        Argument used to define a file of SNP IDs to include.
    **--filter-exclude-snp-file** *<snp_filename>*
        Argument used to define a file of SNP IDs to exclude.
    **--filter-max-missing** *proportion_float*
        Argument used to filter positions by their proportion of missing data, a value of
        0.0 allows for no missing whereas a value of 1.0 ignores missing data. 
    **--filter-max-missing-count** *count_int*
        Argument used to filter positions by the number of samples with missing data, a 
        value of 0 allows for no samples to have missing data.

    ************************
    Position-Based Arguments
    ************************
    **--filter-include-pos** *<chr, chr:pos, chr:start-end, chr:start->*
        Argument used to include matching positions. May be used to include: an
        entire chromosome (i.e. chr); a single position (i.e. chr:pos); a 
        chromosomal locus (i.e. chr:start-end); or a chromosomal span (i.e.
        chr:start-/chr:0-end). This argument may be used multiple times if desired.
    **--filter-exclude-pos** *<chr, chr:pos, chr:start-end, chr:start->*
        Argument used to exclude matching positions. May be used to exclude: an
        entire chromosome (i.e. chr); a single position (i.e. chr:pos); a 
        chromosomal locus (i.e. chr:start-end); or a chromosomal span (i.e.
        chr:start-/chr:0-end). This argument may be used multiple times if desired.
    **--filter-include-pos-file** *<position_filename>*
        Argument used to define a file of positions to include within a tsv file 
        (chromosome and position).
    **--filter-exclude-pos-file** *<position_filename>*
        Argument used to define a file of positions to exclude within a tsv file 
        (chromosome and position).
    **--filter-include-bed** *<position_bed_filename>*
        Argument used to define a BED file of positions to include. Please note that
        filename must end in .bed.
    **--filter-exclude-bed** *<position_bed_filename>*
        Argument used to define a BED file of positions to exclude. Please note that
        filename must end in .bed.

    ********************
    Flag-Based Arguments
    ********************
    **--filter-include-passed**
        Argument used to include positions with the 'PASS' filter flag.
    **--filter-exclude-passed**
        Argument used to exclude positions with the 'PASS' filter flag.
    **--filter-include-filtered** *<filter_flag>*
        Argument used to include positions with the specified filter flag.
    **--filter-exclude-filtered** *<filter_flag>*
        Argument used to exclude positions with the specified filter flag.
'''

import os
import sys
import subprocess
import argparse
import logging

# Import basic bcftools functions
from pgpipe.bcftools import *

# Model file related functions
from pgpipe.model import read_model_file

# Import logging functions
from pgpipe.logging_module import initLogger, logArgs
from pgpipe.misc import argprase_kwargs

def vcf_filter_parser(passed_arguments = []):
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

    out_format_list = ['vcf', 'vcf.gz', 'bcf', 'bed', 'sites']
    out_format_default = 'vcf.gz'

    vcf_parser.add_argument('--out-format', metavar = metavar_list(out_format_list), help = 'Defines the desired output format', type = str, choices = out_format_list, default = out_format_default)

    # General arguments.
    vcf_parser.add_argument('--overwrite', help = "Overwrite previous output", action = 'store_true')

    # Galaxy Option to pipe log to stdout
    vcf_parser.add_argument('--log-stdout', help = argparse.SUPPRESS, action = 'store_true')

    ### Filters

    # Individual list arguments
    indv_filters = vcf_parser.add_mutually_exclusive_group()
    indv_filters.add_argument('--filter-include-indv', help = 'Defines the individual(s) to include. May be used multiple times', nargs = '+', type = str, action = parser_add_to_list())
    indv_filters.add_argument('--filter-exclude-indv', help = 'Defines the individual(s) to exclude. May be used multiple times', nargs = '+', type = str, action = parser_add_to_list())
    
    # Individual file arguments
    indv_file_filters = vcf_parser.add_mutually_exclusive_group()
    indv_file_filters.add_argument('--filter-include-indv-file', help = 'Defines a file of individuals to include', action = parser_confirm_file())
    indv_file_filters.add_argument('--filter-exclude-indv-file', help = 'Defines a file of individuals to exclude', action = parser_confirm_file())

    # Allele count filters
    vcf_parser.add_argument('--filter-only-biallelic', help = 'Only include variants that are biallelic', action = 'store_true')
    vcf_parser.add_argument('--filter-min-alleles', help = 'Include variants with a number of allele >= to the given number', type = int)
    vcf_parser.add_argument('--filter-max-alleles', help = 'Include variants with a number of allele <= to the given number', type = int)

    # Missing data filters
    missing_filters = vcf_parser.add_mutually_exclusive_group()
    missing_filters.add_argument('--filter-max-missing', help = 'Max proportion of missing data allowed (0.0: no missing data, 1.0: include all data)', type = float)
    missing_filters.add_argument('--filter-max-missing-count', help = 'Max number of sample with missing data allowed', type = int)

    # Indel variant-type filters
    indel_filters = vcf_parser.add_mutually_exclusive_group()
    indel_filters.add_argument('--filter-include-indels', help = 'Include variants if they contain an insertion or a deletion', action = 'store_true')
    indel_filters.add_argument('--filter-exclude-indels', help = 'Exclude variants if they contain an insertion or a deletion', action = 'store_true')
    
    # SNP variant-type filter
    snps_filters = vcf_parser.add_mutually_exclusive_group()
    snps_filters.add_argument('--filter-include-snps', help = 'Include variants if they contain a SNP', action = 'store_true')
    snps_filters.add_argument('--filter-exclude-snps', help = 'Exclude variants if they contain a SNP', action = 'store_true')

    # Position filters
    vcf_parser.add_argument('--filter-include-pos', help = 'Defines comma seperated positions (i.e. CHROM:START-END) to include. START and END are optional. May be used multiple times', nargs = '+', type = str, action = parser_add_to_list())
    vcf_parser.add_argument('--filter-exclude-pos', help = 'Defines comma seperated positions (i.e. CHROM:START-END) to exclude. START and END are optional. May be used multiple times', nargs = '+', type = str, action = parser_add_to_list())

    # Position file filters
    vcf_parser.add_argument('--filter-include-pos-file', help = 'Defines a file of positions to include within a file', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-pos-file', help = 'Defines a file of positions to exclude within a file', action = parser_confirm_file())

    # BED position filters
    vcf_parser.add_argument('--filter-include-bed', help = 'Defines a BED file of positions to include', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-bed', help = 'Defines a BED file of positions to exclude', action = parser_confirm_file())

    # Filter-flag filters
    vcf_parser.add_argument('--filter-include-passed', help = "Include variants with the 'PASS' filter flag", action = 'store_true')
    vcf_parser.add_argument('--filter-exclude-passed', help = "Exclude variants with the 'PASS' filter flag", action = 'store_true')
    vcf_parser.add_argument('--filter-include-flag', help = 'Include variants with the specified filter flag', nargs = '+', type = str, action = parser_add_to_list())
    vcf_parser.add_argument('--filter-exclude-flag', help = 'Exclude variants with the specified filter flag', nargs = '+', type = str, action = parser_add_to_list())

    # Info-flag filters
    #vcf_parser.add_argument('--filter-include-info', help = 'Include positions with the specified info flag', nargs = '+', type = str, action = parser_add_to_list())
    #vcf_parser.add_argument('--filter-exclude-info', help = 'Exclude positions with the specified info flag', nargs = '+', type = str, action = parser_add_to_list())

    # SNP filters
    vcf_parser.add_argument('--filter-include-snp', help = 'Include SNP(s) with the matching ID. This argument may be used multiple times', nargs = '+', type = str, action = parser_add_to_list())
    vcf_parser.add_argument('--filter-exclude-snp', help = 'Exclude SNP(s) with the matching ID. This argument may be used multiple times', nargs = '+', type = str, action = parser_add_to_list())

    # SNP file fitlers
    vcf_parser.add_argument('--filter-include-snp-file', help = 'Defines a file of SNP IDs to include', action = parser_confirm_file())
    vcf_parser.add_argument('--filter-exclude-snp-file', help = 'Defines a file of SNP IDs to exclude', action = parser_confirm_file())

    # MAF Filter
    vcf_parser.add_argument('--filter-maf-min', help = 'Include variants with equal or greater MAF values', type = float)
    vcf_parser.add_argument('--filter-maf-max', help = 'Include variants with equal or lesser MAF values', type = float)

    # MAC Filter
    vcf_parser.add_argument('--filter-mac-min', help = 'Include variants with equal or greater MAC values', type = int)
    vcf_parser.add_argument('--filter-mac-max', help = 'Include variants with equal or lesser MAC values', type = int)

    if passed_arguments:
        return vars(vcf_parser.parse_args(passed_arguments))
    else:
        return vars(vcf_parser.parse_args())

def run (**kwargs):
    '''
    Filter function for the PPP

    This function uses the argparse-based function :py:func:`vcf_filter_parser`
    to parse either sys.argv or passed_arguments to obtain the parameters below. 
    The parameters are then translated to their bcftools equivalent. Once all the
    bcftools-based parameters are assigned, bcftools is called.

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
    --filter-only-biallelic
        Only include variants that are biallelic
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
        Filter positions by their missing data proportion. Values should be given
        as max proportion of missing data to accept. For example, a value of 0.0
        allows no missing data whereas a value of 0.1 allows for 10% of samples to 
        have missig data.
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
    --filter-include-passed
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

    # Update kwargs with defaults
    if __name__ != "__main__":
        kwargs = argprase_kwargs(kwargs, vcf_filter_parser)

    # Assign arguments
    vcf_args = argparse.Namespace(**kwargs)

    # Check if the input file is indexed
    vcf_is_indexed = check_for_index(vcf_args.vcf)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(vcf_args, func_name = 'vcf_filter')

    # Argument container for bcftools
    bcftools_call_args = ['view']

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

        # Assign the individuals file to bcftools
        bcftools_call_args.extend(['--samples-file', selected_model.ind_file])

    # Check for incompatible arguments if the user specified only biallelic alleles
    if vcf_args.filter_only_biallelic:

        # Check type arguments
        if vcf_args.filter_include_indels:
            raise Exception('--filter-only-biallelic and --filter-include-indel arguments are incompatible')
        if vcf_args.filter_exclude_indels:
            raise Exception('--filter-only-biallelic and --filter-exclude-indel arguments are incompatible')
        if vcf_args.filter_include_snps:
            raise Exception('--filter-only-biallelic and --filter-include-snps arguments are incompatible')
        if vcf_args.filter_exclude_snps:
            raise Exception('--filter-only-biallelic and --filter-exclude-snps arguments are incompatible')

        # Check allele arguments
        if vcf_args.filter_min_alleles:
            raise Exception('--filter-only-biallelic and --filter-min-alleles arguments are incompatible')
        if vcf_args.filter_max_alleles:
            raise Exception('--filter-only-biallelic and --filter-max-alleles arguments are incompatible')

    # Holds the filename suffix bcftools assigns to the filtered output
    bcftools_out_suffix = None

    # Assign expected bcftools output filename
    bcftools_output_filename = '%s.%s' % (vcf_args.out_prefix, vcf_args.out_format)

    # Check if the user has specified a output filename
    if vcf_args.out:

        # Rename the bcftools output filename, using the output filename
        bcftools_output_filename = vcf_args.out

    # Check if previous output should be overwritten
    if not vcf_args.overwrite:

        # Confirm the bcftools output and log file do not exist
        check_for_bcftools_output(bcftools_output_filename)

    # Individuals-based filters
    if vcf_args.filter_include_indv_file or vcf_args.filter_exclude_indv_file:
        # Used to include a file of individuals to keep
        if vcf_args.filter_include_indv_file:
            bcftools_call_args.extend(['--samples-file', vcf_args.filter_include_indv_file])

        # Used to include a file of individuals to remove
        if vcf_args.filter_exclude_indv_file:
            bcftools_call_args.extend(['--samples-file', '^' + vcf_args.filter_exclude_indv_file])

    # Individual-based filters
    if vcf_args.filter_include_indv or vcf_args.filter_exclude_indv:
        if vcf_args.filter_include_indv:
            bcftools_call_args.extend(['--samples', ','.join(vcf_args.filter_include_indv)])
        if vcf_args.filter_exclude_indv:
            bcftools_call_args.extend(['--samples', '^' + ','.join(vcf_args.filter_exclude_indv)])

    # Position-list filter
    if vcf_args.filter_include_pos or vcf_args.filter_exclude_pos:

        # Check if positions are to be included
        if vcf_args.filter_include_pos:
            # Check if the input is indexed
            if vcf_is_indexed == True:
                bcftools_call_args.extend(['--regions', ','.join(vcf_args.filter_include_pos)])

            # Check if the input is not indexed
            elif vcf_is_indexed != True:
                bcftools_call_args.extend(['--targets', ','.join(vcf_args.filter_include_pos)])
        
        # Check if positions are to be excluded
        if vcf_args.filter_exclude_pos:
            bcftools_call_args.extend(['--targets', '^' + ','.join(vcf_args.filter_exclude_pos)])

    # Position file filters that do not require indexed input
    if (vcf_args.filter_include_pos_file or vcf_args.filter_exclude_pos_file) and vcf_is_indexed != True:
        if vcf_args.filter_include_pos_file:
            bcftools_call_args.extend(['--targets-file', vcf_args.filter_include_pos_file])
        if vcf_args.filter_exclude_pos_file:
            bcftools_call_args.extend(['--targets-file', '^' + vcf_args.filter_exclude_pos_file])

    # Position file filters that requires indexed input
    if (vcf_args.filter_include_pos_file or vcf_args.filter_exclude_pos_file) and vcf_is_indexed == True:
        if vcf_args.filter_include_pos_file:
            bcftools_call_args.extend(['--regions-file', vcf_args.filter_include_pos_file])
        if vcf_args.filter_exclude_pos_file:
            bcftools_call_args.extend(['--targets-file', '^' + vcf_args.filter_exclude_pos_file])

    # Position BED filters that do not require indexed input
    if (vcf_args.filter_include_bed or vcf_args.filter_exclude_bed) and vcf_is_indexed != True:
        if vcf_args.filter_include_bed:
            bcftools_call_args.extend(['--targets-file', vcf_args.filter_include_bed])
        if vcf_args.filter_exclude_bed:
            bcftools_call_args.extend(['--targets-file', '^' + vcf_args.filter_exclude_bed])

    # Position BED filters that requires indexed input
    if (vcf_args.filter_include_bed or vcf_args.filter_exclude_bed) and vcf_is_indexed == True:
        if vcf_args.filter_include_bed:
            bcftools_call_args.extend(['--regions-file', vcf_args.filter_include_bed])
        if vcf_args.filter_exclude_bed:
            bcftools_call_args.extend(['--targets-file', '^' + vcf_args.filter_exclude_bed])

    # Include variants (indel, snp, or both)
    if vcf_args.filter_include_indels or vcf_args.filter_include_snps:

        # List to store variants to include
        include_variant_types = []

        # Check if indels should be included
        if vcf_args.filter_include_indels:
            
            # Add the indels type
            include_variant_types.append('indels')

        # Check if SNPs should be included
        if vcf_args.filter_include_snps:
            
            # Add the snps type
            include_variant_types.append('snps')

        # Store the types to include argument
        bcftools_call_args.extend(['--types', '%s' % ','.join(include_variant_types)])

    # Exclude variants (indel, snp, or both)
    if vcf_args.filter_exclude_indels or vcf_args.filter_exclude_snps:

        # List to store variants to exclude
        exclude_variant_types = []

        # Check if indels should be excluded
        if vcf_args.filter_exclude_indels:
            
            # Add the indels type
            exclude_variant_types.append('indels')

        # Check if SNPs should be excluded
        if vcf_args.filter_exclude_snps:
            
            # Add the snps type
            exclude_variant_types.append('snps')

        # Store the types to include argument
        bcftools_call_args.extend(['--exclude-types', '%s' % ','.join(exclude_variant_types)])

    # Allele-based filters, min and max 
    if vcf_args.filter_only_biallelic or vcf_args.filter_min_alleles or vcf_args.filter_max_alleles:
        if vcf_args.filter_only_biallelic:
            bcftools_call_args.extend(['--min-alleles', '2', '--max-alleles', '2', '--types', 'snps'])
        if vcf_args.filter_min_alleles:
            bcftools_call_args.extend(['--min-alleles', vcf_args.filter_min_alleles])
        if vcf_args.filter_max_alleles:
            bcftools_call_args.extend(['--max-alleles', vcf_args.filter_max_alleles])

    # MAF-based filters
    if vcf_args.filter_maf_min != None or vcf_args.filter_maf_max != None:
        if vcf_args.filter_maf_min != None:
            bcftools_call_args.extend(['--min-af', vcf_args.filter_maf_min])
        if vcf_args.filter_maf_max != None:
            bcftools_call_args.extend(['--max-af', vcf_args.filter_maf_max])

    # MAC-based filters
    if vcf_args.filter_mac_min != None or vcf_args.filter_mac_max != None:
        if vcf_args.filter_mac_min != None:
            bcftools_call_args.extend(['--min-ac', vcf_args.filter_mac_min])
        if vcf_args.filter_mac_max != None:
            bcftools_call_args.extend(['--max-ac', vcf_args.filter_mac_max])

    # List to store expression for BCFtools include expression functionality
    bcftools_expressions = []

    # SNP filters
    if vcf_args.filter_include_snp or vcf_args.filter_exclude_snp:

        # Check for SNPs to include
        if vcf_args.filter_include_snp:

            # Create list to store include SNP expression
            include_snp_expression_list = []

            # Loop the SNPs
            for snp_to_include in vcf_args.filter_include_snp:
                
                # Save the include SNP expression       
                include_snp_expression_list.append('%%ID=="%s"' % snp_to_include)

            # Merge the include SNP expression
            include_snp_expression_str = '%s' % ' || '.join(include_snp_expression_list)

            # Check if parentheses are needed
            if len(include_snp_expression_list) > 1:

                # Add parentheses
                include_snp_expression_str = '(%s)' % include_snp_expression_str

            # Add the include expression to the list
            bcftools_expressions.append(include_snp_expression_str)

        # Check for SNPs to exclude
        if vcf_args.filter_exclude_snp:

            # Create list to store exclude SNP expression
            exclude_snp_expression_list = []

            # Loop the SNPs
            for snp_to_exlude in vcf_args.filter_exclude_snp:
                
                # Save the include SNP expression       
                exclude_snp_expression_list.append('%%ID!="%s"' % snp_to_exlude)

            # Merge the exclude SNP expression
            exclude_snp_expression_str = '%s' % ' && '.join(exclude_snp_expression_list)

            # Check if parentheses are needed
            if len(exclude_snp_expression_list) > 1:

                # Add parentheses
                exclude_snp_expression_str = '(%s)' % exclude_snp_expression_str

            # Add the exclude expression to the list
            bcftools_expressions.append(exclude_snp_expression_str)

    # SNP-file filters
    if vcf_args.filter_include_snp_file or vcf_args.filter_exclude_snp_file:
        if vcf_args.filter_include_snp_file:
            bcftools_expressions.append('%%ID==@%s' % vcf_args.filter_include_snp_file)
        if vcf_args.filter_exclude_snp_file:
            bcftools_expressions.append('%%ID!=@%s' % vcf_args.filter_exclude_snp_file)

    # Missing data filters
    if vcf_args.filter_max_missing != None or vcf_args.filter_max_missing_count != None:
        if vcf_args.filter_max_missing != None:
            bcftools_expressions.append('F_MISSING <= %s' % vcf_args.filter_max_missing)
        if vcf_args.filter_max_missing_count  != None:
            bcftools_expressions.append('N_MISSING <= %s' % vcf_args.filter_max_missing_count)

    # PASS flag Filters
    if vcf_args.filter_include_passed or vcf_args.filter_exclude_passed:

        if vcf_args.filter_include_passed:
            bcftools_expressions.append('%FILTER=="PASS"')

        if vcf_args.filter_exclude_passed:
            bcftools_expressions.append('%FILTER!="PASS"')

    # General Flag filters
    if vcf_args.filter_include_flag or vcf_args.filter_exclude_flag:

        # Check for flags to include
        if vcf_args.filter_include_flag:

            # List to store variants to include
            include_flag_list = []

            # Loop flags to include
            for filtered_to_include in vcf_args.filter_include_flag:
                
                # Save the flag to include
                include_flag_list.append('%%FILTER=="%s"' % filtered_to_include)

            # Merge the flags to include
            include_flag_str = '%s' % ' || '.join(include_flag_list)

            # Check if parentheses are needed
            if len(include_flag_list) > 1:

                # Add parentheses
                include_flag_str = '(%s)' % include_flag_str

            # Append the include flag expression
            bcftools_expressions.append(include_flag_str)

        # Check for flags to exclude
        if vcf_args.filter_exclude_flag:

            # List to store variants to exclude
            exclude_flag_list = []

            # Loop flags to exclude
            for filtered_to_exclude in vcf_args.filter_exclude_flag:
                
                # Save the flag to exclude
                exclude_flag_list.append('%%FILTER!="%s"' % filtered_to_exclude)

            # Merge the flags to exclude
            exclude_flag_str = '%s' % ' && '.join(exclude_flag_list)

            # Check if parentheses are needed
            if len(exclude_flag_list) > 1:

                # Add parentheses
                exclude_flag_str = '(%s)' % exclude_flag_str

            # Append the exclude flag expression
            bcftools_expressions.append(exclude_flag_str)

    '''
    # Info-based filters
    if vcf_args.filter_include_info or vcf_args.filter_exclude_info:
        if vcf_args.filter_include_info:
            for info_to_include in vcf_args.filter_include_info:
                bcftools_call_args.extend(['--keep-INFO', info_to_include])
        if vcf_args.filter_exclude_info:
            for info_to_exclude in vcf_args.filter_exclude_info:
                bcftools_call_args.extend(['--remove-INFO', info_to_exclude])
    '''

    # Save the BCFtools expression argument
    bcftools_expression_arg = '(%s)' % ' && '.join(bcftools_expressions)

    # Add the expression argument to the argument list
    bcftools_call_args.extend(['--include', bcftools_expression_arg])

    logging.info('bcftools parameters assigned')

    # Assign the input file
    bcftools_call_args.append(vcf_args.vcf)

    logging.info('Input file assigned')

    if vcf_args.out_format in ['vcf', 'vcf.gz', 'bcf']:

        # Assign the output format
        bcftools_call_args.extend(return_output_format_args(vcf_args.out_format))

        # Assign the output filename
        bcftools_call_args.extend(['-o', bcftools_output_filename])

        # Call bcftools with the specifed arguments
        bcftools_err = call_bcftools(bcftools_call_args)

    elif vcf_args.out_format in ['bed', 'sites']:

        # Create list for the conversion call
        bcftools_cvt_args = ['query', '-f']

        # Check if the output format is bed
        if vcf_args.out_format == 'bed':

            # Add the arguments to create a bed file
            bcftools_cvt_args.append('%CHROM\t%POS0\t%END\t%ID\n')

        # Check if the output format is sites
        elif vcf_args.out_format == 'sites':

            # Add the arguments to create a bed file
            bcftools_cvt_args.append('%CHROM\t%POS\n')

        # Assign the output filename
        bcftools_cvt_args.extend(['-o', bcftools_output_filename])

        print(bcftools_call_args)

        # Call BCFtools twice
        pipe_bcftools_bcftools(bcftools_call_args, bcftools_cvt_args)

    else:
        raise Exception('Unknown file format: %s' % vcf_args.out_format)

    # Log the bcftools reference
    log_bcftools_reference()

    # Delete any files that were created for bcftools
    if vcf_args.model_file:
        selected_model.delete_ind_file()

if __name__ == "__main__":
    initLogger()
    run(**vcf_filter_parser())

