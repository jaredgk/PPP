def check_if_bed (filename):
    # Open file
    check_file = open(filename, 'rb')
    # Go to beginning of file
    check_file.seek(0)
    # Read first two bytes of file
    check_bytes = check_file.read(2)
    # Check if the bytes match the expected output of a bed
    if check_bytes == b'\x6c\x1b':
        return True
    else:
        return False

def check_if_fam (filename):
    # Open file
    check_file = open(filename, 'r')
    # Read the first line of the file
    check_line = check_file.readline()

    print check_line

    sys.exit()


def assign_plink_input_arg (plink_filenames):
    '''
        Confirms file format for vcftools

        Parameters
        ----------
        filename : str
            Specifies the input filename of unknown format

        Returns
        -------
        list
            Returns vcftools input command for `filename`

        Raises
        ------
        IOError
            If filename is an unknown file format
    '''

    # Plink PED and Map files
    plink_ped = None
    plink_map = None

    # Plink BED, BIM, and FAM files
    plink_bed = None
    plink_bim = None
    plink_fam = None

    # Unknown files
    plink_unknown = []

    # Check if expected files are found
    for plink_filename in plink_filenames:
        # Check for PED and MAP files
        if '.ped' in plink_filename:
            plink_ped = plink_filename
        elif '.map' in plink_filename:
            plink_map = plink_filename

        # Check for BED and associated files
        elif '.bed' in plink_filename:
            plink_bed = plink_filename
        elif '.bim' in plink_filename:
            plink_bim = plink_filename
        elif '.fam' in plink_filename:
            plink_fam = plink_filename

        # Assign any file with unknown file formats
        else:
            # Exclude log file
            if not plink_filename.endswith('.log'):
                # Record unknown files
                plink_unknown.append(plink_filename)

    # If file are missing, determine if they may be assigned
    while plink_unknown:

        # Create copy of the files
        unknown_files = copy.deepcopy(plink_unknown)

        # Loop the files
        for unknown_file in unknown_files:

            # Check if any ped-related files are missing
            if plink_ped or plink_map:
                # Check if the ped file is missing
                if not plink_ped:
                    print unknown_file
                # Check if the map file is missing
                if not plink_map:
                    print unknown_file


            # Check if any ped-related files are missing
            if plink_bed or plink_bim or plink_fam:

                # Check if the bed file is missing
                if not plink_bed:
                    if check_if_bed(unknown_file):
                        plink_bed = unknown_file
                        plink_unknown.remove(unknown_file)

                # Check if the bim file is missing
                if not plink_bim:
                    print unknown_file

                # Check if the fam file is missing
                if not plink_fam:
                    if check_if_fam(unknown_file):
                        plink_fam = unknown_file
                        plink_unknown.remove(unknown_file)


        # Break loop if not files could be assigned
        if unknown_files == plink_unknown:
            break


    # Check that input files included a ped or bed file
    if plink_ped or plink_bed:
        # Check that ped and map files are both found
        if plink_ped and plink_map:
            # Check if the basic plink ped call will work
            if plink_ped[:-4] and plink_map[:-4]:
                return ['--file', plink_ped[:-4]]
            else:
                return ['--ped', plink_ped, '--map', plink_map]

        # Check that bed, bim, and fam files were found
        elif plink_bed and plink_bim and plink_fam:
            # Check if the basic plink ped call will work
            if plink_bed[:-4] and plink_bim[:-4]  and plink_fam[:-4]:
                return ['--file', plink_bed[:-4]]
            else:
                return ['--bed', plink_bed, '--bim', plink_bim, '--fam', plink_fam]

        else:
            print 'Temp Error'



    else:
        # Determine format of the input
        pass


    '''
    if filename_prefix + '.ped' in plink_filenames:
        # Confirm the presence of the map file
        if filename_prefix + '.map' in plink_filenames:
            pass
        else:
            print 'Temp Error'

    elif filename_prefix + '.bed' in plink_filenames:
        # Confirm the presence of the other expected files
        # Assign the associated input command
        if filename.endswith('.vcf'):
            return ['--vcf', filename]
        elif filename.endswith('.vcf.gz'):
            return ['--gzvcf', filename]

    # True if file extension is unknown or not recognized
    else:

        # Checks if the file is unzipped, bgzipped, or gzipped
        vcfname_format = vcf_reader_func.checkFormat(filename)

        # Assign the associated input command, or return an error.
        if vcfname_format == 'vcf':
            return ['--vcf', filename]
        elif vcfname_format == 'bgzip':
            return ['--gzvcf', filename]
        elif vcfname_format == 'bcf':
            return ['--bcf', filename]
        else:
            logging.error('Unknown VCF file format')
            raise Exception('Unknown VCF file format')
    '''


def run (passed_arguments = []):
    '''
        Summary Line

        Complete Summary

        Parameters
        ----------

        Returns
        -------

        Raises
        ------

    '''

    # Grab plink arguments from command line
    plink_args = plink_argument_parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(plink_args, 'plink_tasks')

    if plink_args.ped_prefix:
        assign_plink_input_arg(glob.glob(plink_args.ped_prefix + '*'))
    elif plink_args.ped_files:
        assign_plink_input_arg(plink_args.ped_files)
