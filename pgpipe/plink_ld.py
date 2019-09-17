#!/usr/bin/env python
'''
    Automates the calculation of multiple LD statistics using Plink.

    ##################
    Command-line Usage
    ##################
    The LD statistics automater may be called using the following command:

    .. code-block:: bash
        
        plink_ld.py

    *************
    Example usage
    *************
    Command-line to calculate Lewontin's D-prime statistic

    .. code-block:: bash
        
        plink_ld.py --ped-prefix hapmap1 --ld-format table --ld-statistic r2 --table-d-statistic dprime

    ############
    Dependencies 
    ############
    * `plink 1.9 <https://www.cog-genomics.org/plink2/>`_
    * `plink 2.0 <https://www.cog-genomics.org/plink/2.0/>`

    ############################
    Input Command-line Arguments
    ############################
    **--ped-prefix** *<input_prefix>*
        Argument used to define the filename prefix shared by the ped file (.ped) and
        the map file (.map). Should not be used alongside the specific file arguments 
        (e.g. --ped).
    **--ped** *<ped_filename>*
        Argument used to define the filename of the plink ped file (.ped). Must be 
        called alongside --map. Cannot be called alongside --ped-prefix.
    **--map** *<map_filename>*
        Argument used to define the filename of the plink map file (.map). Must be 
        called alongside --ped. Cannot be called alongside --ped-prefix.
    **--binary-ped-prefix** *<input_prefix>*
        Argument used to define the filename prefix shared by the binary ped file
        (.bed), the fam file (.fam), and the bim file (.bim). Should not be used 
        alongside the specific file arguments (e.g. --binary-ped).
    **--binary-ped** *<binary_ped_filename>*
        Argument used to define the filename of the plink binary ped file (.bed). 
        Must be called alongside --fam and --bim. Cannot be called alongside 
        --binary-ped-prefix.
    **--fam** *<fam_filename>*
        Argument used to define the filename of the plink fam file (.fam). Must be 
        called alongside --binary-ped and --bim. Cannot be called alongside 
        --binary-ped-prefix.
    **--bim** *<bim_filename>*
        Argument used to define the filename of the plink bim file (.bim). Must be 
        called alongside --binary-ped and --fam. Cannot be called alongside 
        --binary-ped-prefix.
    
    #############################
    Output Command-line Arguments
    #############################
    **--out-format** *<output_format>*
        Argument used to define the output format. Supported formats include: gzip
        compressed (gzipped); standard uncompresed (standard); single-precision binary
        (bin32); and double-precision binary (bin64). Please note that both binary formats
        are only supported when called with the square **--lf-format**. By default gzip
        compressed files are produced.
    **--out** *<output_filename>*
        Argument used to define the complete output filename, overrides **--out-prefix**.
        Cannot be used if multiple output files are created.
    **--out-prefix** *<output_prefix>*
        Argument used to define the output prefix (i.e. filename without file extension)
    **--overwrite**
        Argument used to define if previous output should be overwritten.
    
    ###############################
    Basic LD Command-line Arguments
    ###############################
    **--ld-statistic** *<r, r2>*
        Argument used to define the correlation statistic to report. Two options are 
        supported: the raw inter-variant allele count correlations (r) and squared 
        correlations (r2).
    **--ld-format** *<square, square-zero, triangle, table>*
        Argument used to define the matrix result format. Four formats are supported:
        A symmetric matrix (square); a square matrix in which the cells of the upper 
        right triangle are zeroed out (square-zero); only the lower-triangular of the 
        matrix (triangle); and the matrix as a table (table). 
	**--ld-window-snps** *<snp_int>*
		Argument used to define the maximum number of SNPs between LD comparisons.
	**--ld-window-kb** *<snp_int>*
		Argument used to define the maximum distance in bp between LD comparisons.
	**--ld-window-cm** *<snp_int>*
		Argument used to define the maximum distance in cM between LD comparisons.
   
    ****************************
    Table Command-line Arguments
    ****************************
	Please note that the following arguments may only be used with **--ld-format** table.

	**--table-d-statistic** *<dprime, dprime-signed, d>*
        Argument used to add the specified D statistic to table-formatted results. Three options
        are supported: the absolute value of Lewontin's D-prime statistic (dprime); Lewontin's 
        D-prime statistic (dprime-signed); and the value of D prior to division by D\ :sub:`max` 
        (d).  
	**--table-in-phase**
		Argument used to add in-phase allele pairs to table-formatted results.
	**--table-maf**
		Argument used to add MAF values to table-formatted results.
	**--table-r2-threshold** *<r2_float>*
		Argument used to define the threshold for filtering pairs of r2 values.
	**--table-snp** *<snp_str>* *<snp1_str, snp2_str, etc.>*
		Argument used to define one or more SNP(s) for LD analysis. This argument may be used 
        multiple times if desired.
    **--table-snps** *<snp_filename>*
		Argument used to define a file with one or more SNP(s) for LD analysis.
'''

import os
import sys
import subprocess
import shutil
import argparse
import glob
import copy
import logging

# Call PPP-based scripts
#sys.path.insert(0, os.path.abspath(os.path.join(os.pardir,'pppipe')))

from pgpipe.logging_module import initLogger, logArgs
from pgpipe.plink import *

def plink_argument_parser(passed_arguments):
    '''Phase Argument Parser - Assigns arguments for vcftools from command line.
    Depending on the argument in question, a default value may be specified'''

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
                if not getattr(args, self.dest):
                    setattr(args, self.dest, value)
                else:
                    getattr(args, self.dest).extend(value)
        return customAction

    def metavar_list (var_list):
        '''Create a formmated metavar list for the help output'''
        return '{' + ', '.join(var_list) + '}'

    plink_parser = argparse.ArgumentParser()

    plink_prefix = plink_parser.add_mutually_exclusive_group()

    # Input PED arguments
    plink_parser.add_argument("--ped", dest = 'ped_filename', help = "Input PED filename", type = str, action = parser_confirm_file())
    plink_parser.add_argument("--map", dest = 'map_filename', help = "Input MAP filename. Called alongside --ped", type = str, action = parser_confirm_file())
    plink_prefix.add_argument("--ped-prefix", help = "Input PED filename prefix", type = str)

    # Input BED arguments
    plink_parser.add_argument("--binary-ped", dest = 'bed_filename', help = "Input Binary-PED (i.e. BED) filename", type = str, action = parser_confirm_file())
    plink_parser.add_argument("--fam", dest = 'fam_filename', help = "Input FAM filename. Called alongside --bed", type = str, action = parser_confirm_file())
    plink_parser.add_argument("--bim", dest = 'bim_filename', help = "Input BIM filename. Called alongside --bed", type = str, action = parser_confirm_file())
    plink_prefix.add_argument("--binary-ped-prefix", dest = 'bed_prefix', help = "Input Binary-PED (i.e. BED) filename prefix", type = str)

    # Method arguments
    ld_statistics = ['r', 'r2']
    ld_default = 'r2'
    plink_parser.add_argument('--ld-statistic', metavar = metavar_list(ld_statistics), help = 'Specifies the LD statistic', type = str, choices = ld_statistics, default = ld_default)

    shape_formats = ['triangle', 'square', 'square-zero', 'table']
    shape_default = 'table'
    plink_parser.add_argument('--ld-format', metavar = metavar_list(shape_formats), help = 'Specifies the format of the LD results', type = str, choices = shape_formats, default = shape_default)

    # LD-based arguments
    plink_parser.add_argument('--ld-window-snps', help = 'Sets the maximum number of SNPs between LD comparisons', type = int)
    plink_parser.add_argument('--ld-window-kb', help = 'Sets the maximum distance in bp between LD comparisons', type = int)
    plink_parser.add_argument('--ld-window-cm', help = 'Sets the maximum distance in cM between LD comparisons', type = float)

    # Table specific options
    d_statistics = ['d', 'dprime', 'dprime-signed']
    plink_parser.add_argument('--table-d-statistic', metavar = metavar_list(d_statistics), help = 'Adds the specified D statistic to table-formatted results', type = str, choices = d_statistics)

    plink_parser.add_argument('--table-in-phase', help = 'Adds in-phase allele pairs to table-formatted results', action = 'store_true')
    plink_parser.add_argument('--table-maf', help = 'Adds MAF to table-formatted results', action = 'store_true')
    plink_parser.add_argument('--table-r2-threshold', help = "Sets the threshold for filtering pairs of r2 values", type = float)
    plink_parser.add_argument('--table-snp', help = "Specifies one or more SNP(s) for LD analysis", type = str, action = parser_add_to_list())
    plink_parser.add_argument('--table-snps', help = "Specifies a file with one or more SNP(s) for LD analysis", type = str, action = parser_confirm_file())

    # Output arguments
    output_formats = ['standard', 'gzipped', 'bin64', 'bin32']
    output_default = 'gzipped'
    plink_parser.add_argument('--out-format', metavar = metavar_list(output_formats), help = 'Specifies the output format.', type = str, choices = output_formats, default = output_default)

    plink_parser.add_argument('--out', help = 'Defines the output filename. Only usable with vcf-based output')
    plink_parser.add_argument('--out-prefix', help = 'Defines the output filename prefix', default = 'out')

    # General arguments.
    plink_parser.add_argument('--overwrite', help = "Overwrite previous output files", action = 'store_true')
    plink_parser.add_argument('--threads', help = "Set the number of threads", type = int, default = 1)

    if passed_arguments:
        return plink_parser.parse_args(passed_arguments)
    else:
        if len(sys.argv) == 1:
            plink_parser.print_help(sys.stderr)
            sys.exit(1)
        return plink_parser.parse_args()

def run (passed_arguments = []):
    '''
        PLINK-based LD Analysis

        Automates plink-based LD analysis. The function has arguments for
        selecting the LD statistic to be used, the output format of the
        analysis, in addition to many of the optional arguments from plink.

        Parameters
        ----------
        --vcf : str
            Specifies the input VCF filename
        --vcf-fid : str
            Specifies the family ID for all samples
        --ped : str
            Specifies the input PED filename
        --map : str
            Specifies the input MAP filename. Called alongside --ped
        --ped-prefix : str
            Specifies the input PED filename prefix
        --bed : str
            Specifies the input BED filename
        --map : str
            Specifies the input FAM filename. Called alongside --bed
        --bim : str
            Specifies the input BIM filename. Called alongside --bed
        --bed-prefix : str
            Specifies the input BED filename prefix
        --ld-statistic : str
            Specifies the LD statistic {r, r2} (Default: r2)
        --lf-format : str
            Specifies the format of the LD results {triangle, square,
            square-zero, table} (Default: table)
        --ld-window-snps : int
            Sets the maximum number of SNPs between LD comparisons
        --ld-window-kb : int
            Sets the maximum distance in bp between LD comparisons
        --ld-window-cm : float
            Sets the maximum distance in cM between LD comparisons
        --table-d-statistic : str
            Adds the specified D statistic to table-formatted results {d,
            dprime, dprime-signed}
        --table-in-phase
            Adds in-phase allele pairs to table-formatted results
        --table-maf
            Adds MAF to table-formatted results
        --table-r2-threshold : float
            Sets the threshold for filtering pairs of r2 values"
        --table-snp : str, list
            Specifies one or more SNP(s) for LD analysis
        --table-snps : str
            Specifies a file with one or more SNP(s) for LD analysis
        --out-format : str
            Specifies the output format {standard, gzipped, bin64, bin32}
            (Default: gzipped)
        --out : str
            Specifies the output filename. Only usable with vcf-based output
        --out-prefix : str
            Specifies the output filename prefix
        --overwrite
            Overwrite previous output files
        --threads : int
            Set the number of threads

        Returns
        -------
        output : file
            Statistic file output
        log : file
            Log file output

        Raises
        ------
        IOError
            Input file does not exist

    '''

    # Grab plink arguments from command line
    plink_args = plink_argument_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(plink_args, func_name = 'plink_ld')

    # List to hold general plink arguments
    plink_general_args = []

    # Assign the number of threads
    plink_general_args.extend(['--threads', plink_args.threads])

    # List to hold the input argument(s)
    plink_input_args = []

    # Check if the user has selected the table ld format
    if plink_args.ld_format == 'table':

        # Check if the output format isn't supported
        if plink_args.out_format not in ['standard', 'gzipped']:
            raise Exception('The %s output only supports the square ld format (i.e. --lf-format square)' % plink_args.out_format)

        # Check if a r2 threshold was specified alonside the r2 statistic
        if plink_args.table_r2_threshold and plink_args.ld_statistic != 'r2':
            raise Exception('The --table-r2-threshold argument is only compatible with the r2 ld statistic')

    # Check if the user has selected an ld format that isn't table
    else:

        # Check if the output format isn't supported
        if plink_args.out_format != 'standard':
            raise Exception('The %s ld format is incompatible with the standard output format' % plink_args.ld_format)

        # Check if a d statistic was specified, which is incompatible
        if plink_args.table_d_statistic:
            raise Exception('The --table-d-statistic argument is incompatible with the %s ld format' % plink_args.ld_format)

        # Check if in-phase allele pairs were specified, which are incompatible
        if plink_args.table_in_phase:
            raise Exception('The --table-in-phase argument is incompatible with the %s ld format' % plink_args.ld_format)

        # Check if in-phase allele pairs were specified, which are incompatible
        if plink_args.table_maf:
            raise Exception('The --table-maf argument is incompatible with the %s ld format' % plink_args.ld_format)

        # Check if a r2 threshold was specified, which is incompatible
        if plink_args.table_r2_threshold:
            raise Exception('The --table-r2-threshold argument is incompatible with the %s ld format' % plink_args.ld_format)

        # Check if a r2 threshold was specified, which is incompatible
        if plink_args.table_snp:
            raise Exception('The --table-snp argument is incompatible with the %s ld format' % plink_args.ld_format)

        # Check if a r2 threshold was specified, which is incompatible
        if plink_args.table_snps:
            raise Exception('The --table-snps argument is incompatible with the %s ld format' % plink_args.ld_format)

    # Confirm if a prefix has been specified and the files exist
    if plink_args.ped_prefix and confirm_ped_prefix(plink_args.ped_prefix):

        # Assign bed input arguments from a prefix
        plink_input_args = ['--file', plink_args.ped_prefix]

    # Confirm if a prefix has been specified and the files exist
    elif plink_args.bed_prefix and confirm_ped_prefix(plink_args.bed_prefix):

        # Assign the bed input from a prefix
        plink_input_args = ['--bfile', plink_args.bed_prefix]

    # Confirm if required ped files exist
    elif plink_args.ped_filename and confirm_ped_files(plink_args.ped_filename, plink_args.map_filename):

        # Assign the ped input
        plink_input_args = ['--ped', plink_args.ped_filename, '--map', plink_args.map_filename]

    # Confirm if required binary-ped files exist
    elif plink_args.bed_filename and confirm_bed_files(plink_args.bed_filename, plink_args.bim_filename, plink_args.fam_filename):

        # Assign bed input
        plink_input_args = ['--bed', plink_args.bed_filename, '--bim', plink_args.bim_filename, '--fam', plink_args.fam_filename]

    else:
        raise Exception('No PED or Binary-PED input assigned')


    logging.info('Input parameters assigned')

    # List to hold the ld argument(s)
    plink_ld_args = []

    # Add the ld statistic method to the ld argument list
    plink_ld_args.append('--%s' % plink_args.ld_statistic)

    # Check if the format is either triangle or square
    if plink_args.ld_format in ['triangle', 'square']:

        # Add the ld format to the ld argument list
        plink_ld_args.append(plink_args.ld_format)

    # Check if the format is square-zero
    elif plink_args.ld_format == 'square-zero':

        # Add the ld format to the ld argument list
        plink_ld_args.append('square0')

    # Check if the format is a table
    elif plink_args.ld_format == 'table':

        # Add the ld format to the ld argument list
        plink_ld_args.append('inter-chr')

    # Check if an unsupported format was found and raise an exception
    else:
        raise Exception('The format specified (%s) is not supported' % plink_args.ld_format)

    # String to hold output suffix
    plink_output_suffix = 'ld'

    # Check if the output format is gzipped
    if plink_args.out_format == 'gzipped':

        # Add the output format to the ld argument list
        plink_ld_args.append('gz')

        # Add gz to the output suffix
        plink_output_suffix += '.gz'

    # Check if the output format is double-percision binary
    elif plink_args.out_format == 'bed64':

        # Add the output format to the ld argument list
        plink_ld_args.append('bed')

        # Add bed to the output suffix
        plink_output_suffix += '.bed'

    # Check if the output format is single-percision binary
    elif plink_args.out_format == 'bed32':

        # Add the output format to the ld argument list
        plink_ld_args.append('bed4')

        # Add bed to the output suffix
        plink_output_suffix += '.bed'

    elif plink_args.out_format == 'standard':

        # No output format to add to the ld argument list, but not error
        pass

    # Check if an unsupported format was found and raise an exception
    else:
        raise Exception('The format specified (%s) is not supported' % plink_args.out_format)

    # Check if the user has specified the table-specific in-phase allele pairs argument
    if plink_args.table_in_phase:

        # Add the table-specific in-phase allele pairs to the ld argument list
        plink_ld_args.append('in-phase')

    # Check if the user has specified the table-specific d statistic argument
    if plink_args.table_d_statistic:

        # Add the table-specific d statistic to the ld argument list
        plink_ld_args.append(plink_args.table_d_statistic)

    # Check if the user has specified the table-specific maf argument
    if plink_args.table_maf:

        # Add the table-specific maf to the ld argument list
        plink_ld_args.append('with-freqs')

    # Check if the user set a maximum number of SNPs between LD comparisons
    if plink_args.ld_window_snps:

        # Assign the maximum SNP arguments
        plink_ld_args.extend(['--ld-window', plink_args.ld_window_snps])

    # Check if the user set a maximum distance in bp between LD comparisons
    if plink_args.ld_window_kb:

        # Assign the maximum bp distance arguments
        plink_ld_args.extend(['--ld-window-kb', plink_args.ld_window_kb])

    # Check if the user set a maximum distance in cM between LD comparisons
    if plink_args.ld_window_cm:

        # Assign the maximum cm distance arguments
        plink_ld_args.extend(['--ld-window-cm', plink_args.ld_window_cm])

    # Check if the user set a r2 threshold
    if plink_args.table_r2_threshold:

        # Assign the r2 threshold
        plink_ld_args.extend(['--ld-window-r2', plink_args.table_r2_threshold])


    # Check if the user has specified SNP(s) to analyze
    if plink_args.table_snp:

        # Check if there is only a single snp
        if len(plink_args.table_snp) == 1 and '-' not in plink_args.table_snp[0]:

            # Assign the snp
            plink_ld_args.extend(['--ld-snp', plink_args.table_snp[0]])

        # Check if there is more than a single snp
        else:

            # Assign the snps
            plink_ld_args.extend(['--ld-snps', ', '.join(plink_args.table_snp)])

    # Check if the user has specified a SNP file
    if plink_args.table_snps:

        # Assign the snp file
        plink_ld_args.extend(['--ld-snp-list', plink_args.table_snps])

    logging.info('LD parameters assigned')

    # Check if the has specified the output filename, without a prefix
    if plink_args.out and '--out-prefix' not in passed_arguments and '--out-prefix' not in sys.argv:

        # Assign a prefix based on the output filename
        phase_args.out_prefix = copy.deepcopy(phase_args.out)

    # Assign the prefix argument
    plink_output_args = ['--out', plink_args.out_prefix]

    # Assign the expected output filename
    plink_expected_output = '%s.%s' % (plink_args.out_prefix, plink_output_suffix)

    # Check if previous output should not be overwritten
    if not plink_args.overwrite:

        # Check if previous output exists
        if os.path.exists(plink_expected_output):
            raise Exception('%s already exists. Add --overwrite to ignore' % plink_expected_output)

    logging.info('Output parameters assigned')

    # Call plink with the selected arguments
    call_plink(plink_input_args + plink_ld_args + plink_output_args + plink_general_args)

    # Rename output to plink_args.out, if specified
    if plink_args.out:
        shutil.move(plink_expected_output, plink_args.out)
        shutil.move(plink_args.out_prefix + '.log', plink_args.out + '.log')

    # Rename log using plink_args.out_prefix
    else:
        shutil.move(plink_args.out_prefix + '.log', '%s.log' % plink_expected_output)


if __name__ == "__main__":
    initLogger()
    run()
