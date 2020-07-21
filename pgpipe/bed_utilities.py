#!/usr/bin/env python
'''
    Automates various utilites for BED-formatted files. This currently includes: i)
    sample a BED file; ii) subtract from a BED that overlap with a second BED file; 
    iii) extend a BED upstream, downstream, or both upstream and downstream; iv) sort 
    a single BED; v) merge features within one or more BED files; vi) create a BED of 
    complementary features.

    ******************
    Command-line Usage
    ******************
    The BED utilites function may be called using the following command:

    .. code-block:: bash
        
        bed_utilities.py


    ########
    Utilites
    ########

    **************
    Sample Utility
    **************

    .. image:: ../../PPP_assets/PPP_BED_Sample.png
        :width: 100 %
        :align: center

    Given a BED file and a sample size, the sample utility will generate a pseudorandomly 
    sampled BED. Please note that the random seed may be used to reproduced the sample.

    =============
    Example usage
    =============
    Sample 20 features from a BED file:

    .. code-block:: bash
        
        bed_utilities.py --utility sample --bed examples/files/chr1_sites.bed --sample-size 20

    ************
    Sort Utility
    ************

    Given an unsorted BED file, the sort utility will generate a sorted BED file. 

    =============
    Example usage
    =============
    Sort an unsorted BED file:

    .. code-block:: bash
        
        bed_utilities.py --utility sort --bed examples/files/chr1_sites.unsorted.bed

    **************
    Extend Utility
    **************

    Given a BED file and an extend length, the extend utility will increase the length of each 
    feature upstream, downstream, or both upstream and downstream.

    =============
    Example usage
    =============
    Extend upstream by 1kb:

    .. code-block:: bash
        
        bed_utilities.py --utility extend --bed examples/files/chr1_sites.bed --chrom-file examples/files/chr_sizes.txt --extend-upstream 1000

    Extend downstream by 1kb:

    .. code-block:: bash
        
        bed_utilities.py --utility extend --bed examples/files/chr1_sites.bed --chrom-file examples/files/chr_sizes.txt --extend-downstream 1000
    
    Extend flanks (i.e. both upstream and downstream) by 1kb:

    .. code-block:: bash
        
        bed_utilities.py --utility extend --bed examples/files/chr1_sites.bed --chrom-file examples/files/chr_sizes.txt--extend-flanks 1000

    ****************
    Subtract Utility
    ****************

    Given two BED files, the subtract utility will remove BED features from a BED file
    if they overlap with the features from a second BED file.

    =============
    Example usage
    =============
    Remove BED features if they overlap features within the subtract BED file:

    .. code-block:: bash
        
        bed_utilities.py --utility subtract --bed examples/files/chr1_sites.bed --subtract-bed examples/files/chr1_sites.1.bed --subtract-entire-feature


    ******************
    Complement Utility
    ******************

    Given a BED file, the complementary utility will generate a BED file of complementary 
    features.

    =============
    Example usage
    =============
    Return a BED with features that do not overlap within the given file:

    .. code-block:: bash
        
        bed_utilities.py --utility complement --bed examples/files/chr1_sites.bed --chrom-file examples/files/chr_sizes.txt

    *************
    Merge Utility
    *************

    Given one or more BED files, the merge utility will generate a single sorted BED file of
    merged BED features.

    =============
    Example usage
    =============
    Merge BED features from a single BED file:

    .. code-block:: bash
        
        bed_utilities.py --utility merge --bed examples/files/chr1_sites.bed 

    Merge BED features from multiple BED files:

    .. code-block:: bash
        
        bed_utilities.py --utility merge --beds examples/files/chr1_sites.1.bed examples/files/chr1_sites.2.bed examples/files/chr1_sites.3.bed examples/files/chr1_sites.4.bed
    
    ############
    Dependencies 
    ############
    * `BEDtools <https://bedtools.readthedocs.io/en/latest/>`_
    
    ############################
    Input Command-line Arguments
    ############################
    **--bed** *<input_filename>*
        Argument used to define the filename of the BED file.
    **--beds** *<input_filename>* *<input1_filename, input2_filename, etc.>*
        Argument used to define the filename of the BED file(s). May be used multiple 
        times.
    **--chrom-file** *<chrom_filename>*
        Argument used to define the filename of a file with the sizes of each
        chromosome file. Appropriate files may be downloaded from the UCSC Genome
        Browser.
    
    #############################
    Output Command-line Arguments
    #############################
    **--out** *<output_filename>*
        Argument used to define the complete output filename.
    **--overwrite**
        Argument used to define if previous output should be overwritten.
    
    ##################################
    Utility Command-line Specification
    ##################################
    **--utility** *<sample, subtract, extend, sort, merge, complement>*
        Argument used to define the desired utility. Current utilities include: sample
        features from a BED file (sample); subtract features from a BED file that overlap 
        with features within a second BED file (subtract); extend the flanks of features 
        upstream, downstream, or both within a single BED file (extend); sort the 
        features within a single BED file (sort); merge features within one or more BED 
        files (merge); create a BED file of complementary features - i.e. features that 
        do not overlap - from a BED file (complement).

    *************************************
    Sample Utility Command-line Arguments
    *************************************
    **--sample-size** *<sample_size_int>*
        Argument used to define the total sample size.
    **--random-seed** *<seed_int>*
        Argument used to define the seed value for the random number generator.

    ***************************************
    Subtract Utility Command-line Arguments
    ***************************************
    **--subtract-bed** *<subtract_file_filename>*
        Argument used to define the BED file used for removing features/positions.
    **--subtract-entire-feature**
        Argument used to define if entire features within the input BED should be 
        removed if they overlap with features in subtract-bed.
    **--min-reciprocal-overlap** *<overlap_float>*
        Argument used to define the minimum reciprocal overlap of features required 
        for removal (e.g. 0.1 indicates 10% overlap).
    **--min-input-overlap** *<overlap_float>*
        Argument used to define the minimum overlap of input features required for 
        removal.
    **--min-subtract-overlap** *<overlap_float>*
        Argument used to define the minimum overlap of subtract-bed features required 
        for removal.
    **--subtract-entire-feature**
        Argument used to define that features should be removed from the input BED if
        the minimum overlap of **--min-input-overlap** or **--min-subtract-overlap** 
        is reached.

    *************************************
    Extend Utility Command-line Arguments
    *************************************
    **--extend-flanks** *<bp_int>*
        Argument used to define the length of base pairs (bp) to extend both upstream
        and downstream of features.
    **--extend-upstream** *<bp_int>*
        Argument used to define the length of base pairs (bp) to extend upstream of 
        features.
    **--extend-downstream** *<bp_int>*
        Argument used to define the length of base pairs (bp) to extend downstream of 
        features.

    ************************************
    Merge Utility Command-line Arguments
    ************************************
    **--max-merge-distance** *<bp_int>*
        Argument used to define the maximum distance allowed between features to be 
        merged.
'''

import os
import sys
import subprocess
import shutil
import argparse
import glob
import logging
import csv
import numpy as np

from pgpipe.bedtools_wrapper import merge_bed_files, standard_bedtools_call, log_bedtools_reference
from pgpipe.logging_module import initLogger, logArgs
from pgpipe.misc import argprase_kwargs

def bed_argument_parser(passed_arguments = []):
    '''Phase Argument Parser - Assigns arguments for vcftools from command line.
    Depending on the argument in question, a default value may be specified'''

    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError # File not found
                setattr(args, self.dest, value)
        return customAction

    def parser_confirm_file_list ():
        '''Custom action to confirm file exists in list'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                # Loop the list
                for value_item in value:
                    # Check if the file exists
                    if not os.path.isfile(value_item):
                        raise IOError('%s not found' % value_item)
                if not getattr(args, self.dest):
                    setattr(args, self.dest, value)
                else:
                    getattr(args, self.dest).extend(value)
        return customAction

    def metavar_list (var_list):
        '''Create a formmated metavar list for the help output'''
        return '{' + ', '.join(var_list) + '}'

    bed_parser = argparse.ArgumentParser()

    # Input arguments
    bed_input = bed_parser.add_mutually_exclusive_group(required = True)
    bed_input.add_argument('--bed', help = 'Defines the filename of the BED file', type = str, action = parser_confirm_file())
    bed_input.add_argument('--beds', help = 'Defines the filenames of the BED files (may be used multiple times)', type = str, nargs = '+', action = parser_confirm_file_list())
    bed_parser.add_argument('--chrom-file', help = 'File of chromosome sizes', type = str, action = parser_confirm_file())

    # Utility based arguments
    utility_list = ['sample', 'subtract', 'extend', 'sort', 'merge', 'complement']
    bed_parser.add_argument('--utility', metavar = metavar_list(utility_list), help = 'Specifies the utility to be used', type = str, choices = utility_list, required = True)

    # Sample-specific arguments
    sample_group = bed_parser.add_argument_group('Sample Utility Arguments')
    sample_group.add_argument('--sample-size', help="Defines the total sample size", type = int)
    sample_group.add_argument('--random-seed', help="Defines the seed value for the random number generator", type = int, default = np.random.randint(1, high = 1000000000))

    # Subtract-specific arguments
    subtract_group = bed_parser.add_argument_group('Subtract Utility Arguments')
    subtract_group.add_argument('--subtract-bed', help = 'BED file used for removing features/positions', type = str, action = parser_confirm_file())
    subtract_group.add_argument('--subtract-entire-feature', help = 'Remove entire features within the input if they overlap with features in subtract-bed', action = 'store_true')
    subtract_group.add_argument('--min-reciprocal-overlap', help = 'Minimum reciprocal overlap of features required for removal', type = float)
    subtract_group.add_argument('--min-input-overlap', help = 'Minimum overlap of input features required for removal', type = float)
    subtract_group.add_argument('--min-subtract-overlap', help = 'Minimum overlap of subtract-bed features required for removal', type = float)
    subtract_group.add_argument('--either-min-overlap', help = 'Remove a feature if either minimum overlap is reached ', action = 'store_true')

    # Extend-specific arguments
    extend_group = bed_parser.add_argument_group('Extend Utility Arguments')
    extend_group.add_argument('--extend-flanks', help = 'Base pairs to extend both flanks', type = int)
    extend_group.add_argument('--extend-upstream', help = 'Base pairs to extend upstream', type = int)
    extend_group.add_argument('--extend-downstream', help = 'Base pairs to extend downstream', type = int)

    # Extend-specific arguments
    merge_group = bed_parser.add_argument_group('Merge Utility Arguments')
    merge_group.add_argument('--max-merge-distance', help = 'Max distance allowed between features to be merged', type = int)

    # Other basic arguments. Expand as needed
    bed_parser.add_argument('--out', help = 'Defines the output filename', default = 'out.bed')
    bed_parser.add_argument('--overwrite', help = "Defines that previous output files should be overwritten", action = 'store_true')

    if passed_arguments:
        return vars(bed_parser.parse_args(passed_arguments))
    else:
        return vars(bed_parser.parse_args())

def random_bed_sampler (bed_filename, out_filename, sample_size, with_replacements = False):
    '''
        Random Sampler.

        Returns a list of randomly selected elements from stat_file_data. The
        length of this list is defined by sample_size.

        Parameters
        ----------
        bed_filename : str
            Filename of the bed file
        out_filename : str
            Filename of the output file
        sample_size : int
            Size of the total sample
        with_replacements : bool, optional
            If the function should sample with replacements
    '''

    # Open the BED file
    with open(bed_filename, 'r') as bed_file:

        # Start a CSV reader, split by tab
        bed_reader = csv.reader(bed_file, delimiter = '\t')

        # Store the split contents within file
        bed_list = list(bed_reader)

    # Use numpy choice to randomly sample the bed list by index
    random_samples = np.random.choice(len(bed_list), int(sample_size), replace = with_replacements)

    # Convert the sample indices to a list
    random_samples = [bed_list[random_sample] for random_sample in random_samples]

    # Open the output file
    with open(out_filename,'w') as out_file:

        # Create the output writer
        out_writer = csv.writer(out_file, delimiter = '\t')

        # Write the output
        out_writer.writerows(random_samples)

def run (**kwargs):
    '''
        Utilites for BED-formatted files

        This function uses the argparse-based function :py:func:`bed_argument_parser`
        to parse either sys.argv or passed_arguments to obtain the parameters below. 
        The parameters are then translated to their BEDtools equivalent. Once all the 
        parameters are assigned, BEDtools is called.

        Parameters
        ----------
        --bed : str
            Filename of the BED
        --beds : list
            List of BED filenames
        --chrom-file : str
            File of chromosome sizes
        --utility : str
            Utility to be used
        --out : str
            Complete output filename
        --overwrite
            Species if previous output should be overwriten
        --sample-size : int
            Specifies the total sample size.
        --random-seed : int
            Specifies the random seed value for the random number generator
        --subtract-bed : str
            Filename of the BED to subtract
        --subtract-entire-feature : bool
            Remove entire feature if there is an overlap in --subtract-bed
        --min-reciprocal-overlap : float
            Minimum reciprocal overlap for feature removal
        --min-input-overlap : float
            Minimum input (i.e. --bed) overlap for feature removal
        --min-subtract-overlap : float
            Minimum --subtract-bed overlap for feature removal
        --either-min-overlap : bool
            Remove a feature if either minimum overlap is reached
        --extend-flanks : int
            Base pairs to extend both flanks
        --extend-upstream : int
            Base pairs to extend upstream
        --extend-downstream : int
            Base pairs to extend downstream
        --max-merge-distance : int
            Max distance allowed between features to be merged

        Raises
        ------
        IOError
            Output file already exists
        Exception
            More than one BED given - unless using the merge utility
        Exception
            No extension arguments given with extend utility
        Exception
            No --chrom-file given with extend/complement utilities

    '''

    # Update kwargs with defaults
    if __name__ != "__main__":
        kwargs = argprase_kwargs(kwargs, bed_argument_parser)

    # Assign arguments
    bed_args = argparse.Namespace(**kwargs)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(bed_args, 'bed_utilities')

    # Check if previous output should not be overwriten
    if not bed_args.overwrite:

        # Check if previous output exists
        if os.path.exists(bed_args.out):
            raise IOError('Utility output already exists')

    # Check if the utility only functions with a single BED file
    if bed_args.beds and bed_args.utility not in ['merge']:
        raise Exception('%s only supports a single BED file' % bed_args.utility)

    # Check if the user has selected the sample utility
    if bed_args.utility == 'sample':

        # Check that a sample size was specified
        if not bed_args.sample_size:
            raise Exception('%s requires a sample size. Please use --sample-size' % bed_args.utility)

        # Set the random seed
        np.random.seed(bed_args.random_seed)

        # Randomly sample the BED file
        random_bed_sampler(bed_args.bed, bed_args.out, bed_args.sample_size, with_replacements = False)

    # Check if the user has selected the merge utility
    elif bed_args.utility == 'merge':

        # Create list of BED files to merge
        bed_files_to_merge = []

        # Check if multiple BED files have been specified
        if bed_args.beds:

            # Assign the specified BED files 
            bed_files_to_merge = bed_args.beds

        # Check if a single BED file has been specified
        if bed_args.bed:

            # Assign the specified BED files 
            bed_files_to_merge = [bed_args.bed]

        # Create list of optional merge args
        optional_merge_args = []

        # Check if --max-merge-distance has been specified 
        if bed_args.max_merge_distance:
            optional_merge_args.extend(['-d', bed_args.max_merge_distance])

        # Merge the BED files
        merge_bed_files(bed_files_to_merge, bed_args.out, optional_merge_args)

        # Log the bedtools reference
        log_bedtools_reference()

    # Catch the rest of the utilites 
    else:

        # Check if the utility requires a chromosome size file
        if bed_args.utility in ['complement', 'extend']:

            # Check if a chromosome sizes file has been defined
            if not bed_args.chrom_file:
                raise Exception('%s requires a file of chromosome sizes. Please use --chrom-file' % bed_args.utility)

        # Create list of bedtools
        bedtools_arg_list = []
        
        # Check if the subtract utility has been selected
        if bed_args.utility == 'subtract':

            # Check that a subtract BED file has been specified
            if not bed_args.subtract_bed:
                raise Exception('%s requires a subtract BED file. Please use --subtract-bed' % bed_args.utility)

            # Assign the utility
            bedtools_arg_list.append(bed_args.utility)

            # Assign the input arguments
            bedtools_arg_list.extend(['-a', bed_args.bed])

            # Assign the --subtract-bed argument
            bedtools_arg_list.extend(['-b', bed_args.subtract_bed])

            # Check if --subtract-entire-feature has been specified 
            if bed_args.subtract_entire_feature:
                bedtools_arg_list.append('-A')

            # Check if --min-reciprocal-overlap has been specified 
            if bed_args.min_reciprocal_overlap:
                bedtools_arg_list.extend(['-f', bed_args.min_reciprocal_overlap, '-r'])

            # Check if --min-input-overlap has been specified 
            if bed_args.min_input_overlap:
                bedtools_arg_list.extend(['-f', bed_args.min_input_overlap])

            # Check if --min-subtract-overlap has been specified 
            if bed_args.min_subtract_overlap:
                bedtools_arg_list.extend(['-F', bed_args.min_subtract_overlap])

            # Check if --either-min-overlap has been specified 
            if bed_args.either_min_overlap:
                bedtools_arg_list.append('-e')

        # Check if the extend utility has been selected
        elif bed_args.utility == 'extend':

            # Confirm at least one extension argument has been assigned
            if not bed_args.extend_flanks and not bed_args.extend_upstream and not bed_args.extend_downstream:
                raise Exception('%s requires at least one extension argument (e.g. --extend-flanks)' % bed_args.utility)

            # Confirm extend flanks is not being called alone
            if bed_args.extend_flanks and (bed_args.extend_upstream or bed_args.extend_downstream):
                raise Exception('--extend-flanks cannot be used with other extend commands')

            # Assign the utility
            bedtools_arg_list.append('slop')

            # Assign the input arguments
            bedtools_arg_list.extend(['-i', bed_args.bed])

            # Assign the --chrom-file argument
            bedtools_arg_list.extend(['-g', bed_args.chrom_file])

            # Check if --extend-flanks has been specified 
            if bed_args.extend_flanks:
                bedtools_arg_list.extend(['-b', bed_args.extend_flanks])

            # Check if --extend-upstream/--extend-downstream have been specified
            else:
                bedtools_arg_list.extend(['-l', bed_args.extend_upstream if bed_args.extend_upstream is not None else 0])
                bedtools_arg_list.extend(['-r', bed_args.extend_downstream if bed_args.extend_downstream is not None else 0])

        # Check if the sort utility has been selected
        elif bed_args.utility == 'sort':

            # Assign the utility
            bedtools_arg_list.append(bed_args.utility)

            # Assign the input arguments
            bedtools_arg_list.extend(['-i', bed_args.bed])

        # Check if the complement utility has been selected
        elif bed_args.utility == 'complement':

            # Assign the utility
            bedtools_arg_list.append(bed_args.utility)

            # Assign the input arguments
            bedtools_arg_list.extend(['-i', bed_args.bed])

            # Assign the --chrom-file argument
            bedtools_arg_list.extend(['-g', bed_args.chrom_file])

        # Catch unknown utilities    
        else:
            raise Exception('%s is an unknown utility. Please check input' % bed_args.utility)

        # Call bedtools
        standard_bedtools_call(bedtools_arg_list, bed_args.out)

        # Log the bedtools reference
        log_bedtools_reference()

if __name__ == "__main__":
    initLogger()
    run(**bed_argument_parser())
