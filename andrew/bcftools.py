import os
import sys
import logging
import subprocess

sys.path.insert(0, os.path.abspath(os.path.join(os.pardir,'jared')))

from vcf_reader_func import checkFormat

def return_output_format_args (output_format):
    '''
        Return bcftools arguments for output format

        Parameters
        ----------
        output_format : str
            The specified output format

        Raises
        ------
        Exception
            If output format is unsupported by bcftools
    '''

    # Return the output format arguments
    if output_format == 'vcf':
        return ['-O', 'v']
    elif output_format == 'bcf':
        return ['-O', 'b']
    elif output_format == 'vcf.gz':
        return ['-O', 'z']
    else:
        raise Exception('Unsupported file format')

def check_bcftools_for_errors (bcftools_stderr):
    '''
        Checks the bgzip stderr for errors

        Parameters
        ----------
        bcftools_stderr : str
            bcftools stderr

        Raises
        ------
        Exception
            If bcftools stderr returns an error
    '''

    # Expand as errors are discovered

    # Log warning messages
    if 'W::' in bcftools_stderr:
        logging.warning(bcftools_stderr.strip())

    # Report errors that are not warnings
    elif bcftools_stderr:
        raise Exception(bcftools_stderr)

def pipe_bcftools (bcftools_call_args):
    '''
        Calls bcftools with pipe output

        The output of this function is the stdout and stderr of bcftools. This
        function should only be used if bcftools is being used as the stdin of
        another function. Please note that this function does not check the for
        errors in the bcftools call. Please check for errors after the call is
        closed using check_bcftools_for_errors.

        Parameters
        ----------
        bcftools_stderr : str
            bcftools stderr

        Returns
        -------
        bcftools_call : PIPE
            Pipe of subprocess call, including both stdout and stderr

    '''

    # bcftools subprocess call
    bcftools_call = subprocess.Popen(['bcftools'] + list(map(str, bcftools_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    return bcftools_call

def pipe_bcftools_to_list (bcftools_call_args):
    '''
        Pipes output of bcftools to a set

        The purpose of this function is to return a set of the output of
        various bcftools commands, this will result in the removal of duplicate
        entries.

        Parameters
        ----------
        bcftools_call_args : list
            bcftools arguments

        Returns
        -------
        set_to_return : set
            Set of unique elements from the bcftools output
    '''

    # Open bcftools pipe
    bcftools_call = pipe_bcftools(bcftools_call_args)

    # Create a list to hold unique elements
    list_to_return = []

    try:

        # Iterate the bcftools stdout unless error occurs
        for bcftools_stdout_line in iter(bcftools_call.stdout.readline, b''):

            # Check if code is running in python 3
            if sys.version_info[0] == 3:
                # Convert bytes to string
                bcftools_stdout_line = bcftools_stdout_line.decode()

            # Remove the newline character
            bcftools_line = bcftools_stdout_line.strip()
            # Save the line to the list
            list_to_return.append(bcftools_line)

    except:
        raise Exception('bcftools call error')

    # Close the bcftools stdout
    bcftools_call.stdout.close()

    # Wait for bctools to finish
    bcftools_call.wait()

    # Read the bcftools stderr
    bcftools_stderr = bcftools_call.stderr.read()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        bcftools_stderr = bcftools_stderr.decode()

    # Check that the log file was created correctly
    check_bcftools_for_errors(bcftools_stderr)

    logging.info('bcftools call complete')

    return list_to_return

def pipe_bcftools_to_set (bcftools_call_args):
    '''
        Pipes output of bcftools to a set

        The purpose of this function is to return a set of the output of
        various bcftools commands, this will result in the removal of duplicate
        entries.

        Parameters
        ----------
        bcftools_call_args : list
            bcftools arguments

        Returns
        -------
        set_to_return : set
            Set of unique elements from the bcftools output
    '''

    # Open bcftools pipe
    bcftools_call = pipe_bcftools(bcftools_call_args)

    # Create a set to hold unique elements
    set_to_return = set()

    try:

        # Current element, reduces duplicates if VCF is sorted
        previous_element = None

        # Iterate the bcftools stdout unless error occurs
        for bcftools_stdout_line in iter(bcftools_call.stdout.readline, b''):

            # Check if code is running in python 3
            if sys.version_info[0] == 3:
                # Convert bytes to string
                bcftools_stdout_line = bcftools_stdout_line.decode()

            # Remove the newline character
            bcftools_line = bcftools_stdout_line.strip()
            # Check if the element is different from the previous element
            if bcftools_line != previous_element:
                # Store the new element for comparisons to reduce duplicates
                previous_element = bcftools_line
                # Save the element in the set
                set_to_return.add(bcftools_line)

    except:
        raise Exception('bcftools call error')

    # Close the bcftools stdout
    bcftools_call.stdout.close()

    # Wait for bctools to finish
    bcftools_call.wait()

    # Read the bcftools stderr
    bcftools_stderr = bcftools_call.stderr.read()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        bcftools_stderr = bcftools_stderr.decode()

    # Check that the log file was created correctly
    check_bcftools_for_errors(bcftools_stderr)

    logging.info('bcftools call complete')

    return set_to_return

def call_bcftools (bcftools_call_args):
    '''
        Calls bcftools

        The function calls bcftools.

        Parameters
        ----------
        bcftools_call_args : list
            bcftools arguments

        Returns
        -------
        vcftools_err : str
            vcftools log output

        Raises
        ------
        Exception
            If bcftools stderr returns an error
    '''

    # bcftools subprocess call
    bcftools_call = subprocess.Popen(['bcftools'] + list(map(str, bcftools_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Wait for bcftools to finish
    bcftools_stdout, bcftools_stderr = bcftools_call.communicate()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        bcftools_stderr = bcftools_stderr.decode()

    check_bcftools_for_errors(bcftools_stderr)

    logging.info('bcftools call complete')

    return bcftools_stderr

def get_unique_chrs (filename):

    # Get set of the chromosomes
    chromosome_set = pipe_bcftools_to_set(['query', '-f', '%CHROM\n', filename])

    # Return list of the unique chromosomes
    return list(chromosome_set)

def get_samples (filename):

    # Get list of the samples
    sample_list = pipe_bcftools_to_list(['query', '-l', filename])

    # Return list of the samples
    return sample_list

def check_for_index (filename):
    '''
        Checks for index file

        If the file is capable of having an index (i.e. bgzipped-VCF or BCF) the
        function will return either True (i.e. index found) or False. However,
        if the file is a VCF the function will return None (as VCF files cannot
        have an index). An error is returned if the file is either a
        gzipped-VCF file or not a VCF-based format.

        Parameters
        ----------
        filename : str
            Filename of VCF-formatted file

        Returns
        -------
        bool, None
            Returns bool for VCF.GZ and BCF files. Returns None for VCF files

        Raises
        ------
        Exception
            If the file is a gzipped-VCF or of an unknown format
    '''

    # Assign the file format
    file_format = checkFormat(filename)

    # Check if the file to be indexed is a vcf.gz
    if file_format == 'bgzip':
        # Check if the index (.tbi) exists
        if os.path.isfile(filename + '.tbi'):
            return True

    # Check if the file to be indexed is a bcf
    elif file_format == 'bcf':
        # Check if the index (.csi) exists
        if os.path.isfile(filename + '.csi'):
            return True

    # Check if the file is vcf (does not need an index)
    elif file_format == 'vcf':
        return None

    # Check if the file is gzip-compressed vcf (cannot have an index)
    elif file_format == 'gzip':
        raise Exception('GZIP-compressed VCF files do not support index files.')

    # Check if the file is an unknown format
    else:
        raise Exception('Unknown file format')

    # Return false if no index is found
    return False

def delete_index (filename):
    '''
        Deletes an index file

        If the file is capable of having an index (i.e. bgzipped-VCF or BCF)
        this function will delete the index. However, if the file is either a
        VCF or a gzip-compressed VCF the function will return an error. The
        function also results in an error if the index cannot be found. This
        function should be used following check_for_index.

        Parameters
        ----------
        filename : str
            Filename of VCF-formatted file

        Raises
        ------
        Exception
            No index file could be found
        Exception
            If the file is a gzipped-VCF or a VCF
    '''

    # Assign the file format
    file_format = checkFormat(filename)

    # Check if the file to be indexed is a vcf.gz
    if file_format == 'bgzip':
        # Check if the index (.tbi) exists
        if os.path.isfile(filename + '.tbi'):
            # Delete the index
            os.remove(filename + '.tbi')
            return

    # Check if the file to be indexed is a bcf
    elif file_format == 'bcf':
        # Check if the index (.csi) exists
        if os.path.isfile(filename + '.csi'):
            # Delete the index
            os.remove(filename + '.csi')
            return

    # Check if the file is vcf (cannot have an index)
    elif file_format == 'vcf':
        raise Exception('VCF format does not support index files.')

    # Check if the file is gzip-compressed vcf (cannot have an index)
    elif file_format == 'gzip':
        raise Exception('GZIP-compressed VCF files do not support index files.')

    # Return error if no index is found
    raise Exception('No index file found.')

def create_index (filename):
    '''
        Creates an index file

        If the file is capable of having an index (i.e. bgzipped-VCF or BCF)
        this function will create an index file. However, if the file is a
        different format the function will return an error.

        Parameters
        ----------
        filename : str
            Filename of VCF-formatted file

        Raises
        ------
        Exception
            If the file is not a bgzipped-VCF or BCF
    '''

    # Assign the file format
    file_format = checkFormat(filename)

    # Check if the file to be indexed is a vcf.gz
    if file_format == 'bgzip':
        # Create a index (.tbi)
        call_bcftools(['index', '-t', filename])

    # Check if the file to be indexed is a bcf
    elif file_format == 'bcf':
        # Create a index (.csi)
        call_bcftools(['index', '-c', filename])

    # Report if file cannot be indexed
    else:
        raise Exception('Error creating index for: %s. Only .bcf and .vcf.gz (bgzip) files are supported.' % filename)

def chr_subset_file (filename, chromosome, output_prefix, output_format, from_bp = None, to_bp = None):
    '''
        Creates chromosome subset

        This function is used to create a VCF-formatted subset with only
        the data from a single chromosome.

        Parameters
        ----------
        filename : str
            Filename of VCF-formatted input
        chromosome : str
            Chromosome to subset
        output_prefix : str
            Prefix of the VCF-formatted output (i.e. without file extension)
        output_format : str
            The format of the output (e.g. vcf, bcf, vcf.gz)
        from_bp : int, optional
            Lower bound of sites to include
        to_bp : int, optional
            Upper bound of sites to include
    '''

    # Creates a list to the arguments and store the bcftools call
    subset_args = ['view']

    # Assign the output format arguments
    output_format_args = return_output_format_args(output_format)

    # Store the output format arguments
    subset_args.extend(output_format_args)

    # Stores the specified output filename
    vcf_output = '%s.%s' % (output_prefix, output_format)

    # Assigns the output file to the arguments
    subset_args.extend(['-o', vcf_output])

    # Holds the subset argument
    chr_subet_arg = chromosome

    # Check if either bp position arguments were specified
    if from_bp or to_bp:

        # List of the position arguments, in their required order
        position_args = [':', from_bp, '-', to_bp]

        # Filter the position arguments to remove empty values
        filttered_position_args = filter(None, position_args)

        # Map the arguments to str and add them to the chromosome argument
        chr_subet_arg += ''.join(map(str, filttered_position_args))

    # Checks if the input file has an index, then subset to the arguments
    if check_for_index(filename):
        # Subsets using the index
        subset_args.extend(['-r', chr_subet_arg])
    else:
        # Subsets using the stdout
        subset_args.extend(['-t', chr_subet_arg])

    # Assigns the input file to the arguments
    subset_args.append(filename)

    # Call bcftools
    call_bcftools(subset_args)

def concatenate (filenames, output_prefix, output_format, keep_original = False):
    '''
        Concatenate multiple VCF-formatted files

        This function will concatenate multiple VCF-formatted files into a
        single VCF-formatted file of the specifed format.

        Parameters
        ----------
        filenames : list
            List of VCF-formatted input filenames
        output_prefix : str
            Prefix of the VCF-formatted output (i.e. without file extension)
        output_format : str
            The format of the output (e.g. vcf, bcf, vcf.gz)
    '''

    # Holds the arguments to convert to VCF format
    concat_args = ['concat']

    # Assign the output format arguments
    output_format_args = return_output_format_args(output_format)

    # Store the output format arguments
    concat_args.extend(output_format_args)

    # Stores the specified output filename
    vcf_output = '%s.%s' % (output_prefix, output_format)

    # Assigns the output file to the arguments
    concat_args.extend(['-o', vcf_output])

    # Assigns the input files to merge
    concat_args.extend(filenames)

    # Call bcftools
    call_bcftools(concat_args)

    # Delete the original files once the merged file is created
    if not keep_original:
        for filename in filenames:
            if check_for_index(filename) == True:
                delete_index(filename)
            os.remove(filename)

def convert_to_bcf (filename, output_prefix, keep_original = False):
    '''
        Converts a VCF-formatted file to BCF

        This function will convert a VCF-formatted file to BCF with the
        specified filename prefix. The function also has the option to keep or
        delete the input file once the BCF file has been created.

        Parameters
        ----------
        filename : str
            Filename of VCF-formatted input
        output_prefix : str
            Prefix of the BCF output (i.e. without file extension)
        keep_original : bool, optional
            If the input file should be kept once converted
    '''

    # Holds the arguments to convert to BCF format
    convert_args = ['convert', '-O', 'b']

    # Stores the specified output_prefix to the BCF file
    bcf_output = '%s.bcf' % output_prefix

    # Assigns the output file to the arguments
    convert_args.extend(['-o', bcf_output])

    # Assigns the specified input to the arguments
    convert_args.append(filename)

    # Call bcftools
    call_bcftools(convert_args)

    # Delete the original file once the bcf file is created
    if not keep_original:
        if check_for_index(filename) == True:
            delete_index(filename)
        os.remove(filename)

def convert_to_vcf (filename, output_prefix, keep_original = False):
    '''
        Converts a VCF-formatted file to VCF

        This function will convert a VCF-formatted file to VCF with the
        specified filename prefix. The function also has the option to keep or
        delete the input file once the VCF file has been created.

        Parameters
        ----------
        filename : str
            Filename of VCF-formatted input
        output_prefix : str
            Prefix of the VCF output (i.e. without file extension)
        keep_original : bool, optional
            If the input file should be kept once converted
    '''

    # Holds the arguments to convert to VCF format
    convert_args = ['view', '-O', 'v']

    # Stores the specified output_prefix to the VCF file
    vcf_output = '%s.vcf' % output_prefix

    # Assigns the output file to the arguments
    convert_args.extend(['-o', vcf_output])

    # Assigns the specified input to the arguments
    convert_args.append(filename)

    # Call bcftools
    call_bcftools(convert_args)

    # Delete the original file once the vcf file is created
    if not keep_original:
        if check_for_index(filename) == True:
            delete_index(filename)
        os.remove(filename)

def convert_to_vcfgz (filename, output_prefix, keep_original = False):
    '''
        Converts a VCF-formatted file to bgzipped-VCF

        This function will convert a VCF-formatted file to bgzipped-VCF with the
        specified filename prefix. The function also has the option to keep or
        delete the input file once the bgzipped-VCF file has been created.

        Parameters
        ----------
        filename : str
            Filename of VCF-formatted input
        output_prefix : str
            Prefix of the bgzipped-VCF output (i.e. without file extension)
        keep_original : bool, optional
            If the input file should be kept once converted
    '''

    # Holds the arguments to convert to VCFGZ format
    convert_args = ['view', '-O', 'z']

    # Stores the specified output_prefix to the VCFGZ file
    vcfgz_output = '%s.vcf.gz' % output_prefix

    # Assigns the output file to the arguments
    convert_args.extend(['-o', vcfgz_output])

    # Assigns the specified input to the arguments
    convert_args.append(filename)

    # Call bcftools
    call_bcftools(convert_args)

    # Delete the original file once the vcfgz file is created
    if not keep_original:
        if check_for_index(filename) == True:
            delete_index(filename)
        os.remove(filename)


import sys
print (get_unique_chrs(sys.argv[1]))