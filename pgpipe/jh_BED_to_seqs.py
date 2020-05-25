import sys
import os
import subprocess
import logging
import pgpipe.vcf_reader_func as vf
from pgpipe.model import Model, read_model_file
from pgpipe.misc import confirm_executable
from pgpipe.vcf_to_ima import generateSequence
from pgpipe.genome_region import Region, RegionList
import pysam

##given
##    fasta reference file  (e.g. 1 or more chromosomes,  such as a genome reference)
##    bed file 0-based intervals  bedfile chr column must match names in fasta reference
##    vcf filename
##    model filename
##    model name (i.e. one of the models in the model file)
##
##    other args:
##        addchrstr=False,  if True, then the bedfile and fasta ref just use chromosome numbers,  but vcf also uses 'chr'
##        dropchrstr=False, if True then bedfile and fasta ref use 'chr' in chromosome name,  but vcf file does not
##
##    use pysam.FastaFile to get the sequence, uses BED style intervals
##        this requires a faidx index



def check_samtools_for_errors (samtools_stderr):
    '''
        Checks the samtools stderr for errors

        Parameters
        ----------
        samtools_stderr : str
            samtools stderr

        Raises
        ------
        IOError
            If samtools stderr returns an error
    '''

    # Returns True if the job completed without error
    print(samtools_stderr)
    if len(str(samtools_stderr))==0:
        pass

    # Print output for samtools if error is detected
    elif 'Format error' in str(samtools_stderr):
        # Splits log into list of lines
        samtools_stderr_lines = samtools_stderr.splitlines()
        # Prints the error(s)
        raise Exception('\n'.join((output_line for output_line in samtools_stderr_lines if output_line.startswith('Format error'))))

    # Print output if not completed and no error found. Unlikely to be used, but included.
    else:
        raise Exception(samtools_stderr)

def call_samtools_faidx (fastaname,fastaout,rangearg):
    '''
        Calls samtools faidx to generate a fasta file index

        Parameters
        ----------
        fastaname: input fasta file name
    
        Raises
        ------
        Exception
            If fasta file not found
            If samtools faidx returns an error
    '''
    # Confirm where the specifed executable is located
    samtools_path = confirm_executable('samtools')

    # Check if the executable was found
    if not samtools_path:
        raise IOError('samtools not found. Please confirm the executable is installed')
    # check that the fasta file exists
    if os.path.isfile(fastaname) is False:
        raise IOError('fasta file %s does not exist'%fastaname)

    # samtools subprocess call
    samtools_call = subprocess.Popen(cmdargs, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Wait for samtools to finish
    samtools_stdout, samtools_stderr = samtools_call.communicate()

    # Check if code is running in python 3
    if sys.version_info[0] == 3:
        # Convert bytes to string
        samtools_stderr = samtools_stderr.decode()

    # Check that the log file was created correctly
    check_samtools_for_errors(samtools_stderr)

    return samtools_stderr


def get_model_sequences_from_region(vcf_reader=None,popmodel=None,
                                    seq_reference=None, region=None,
                                    useNifmissingdata=None):    
    '''
        gets a list of sequences using a list of vcf records and sequence

 ? can't handle haploidy or mixed ploidy among individuals
        
        Parameters
        ----------
        recordlist
            a list of vcf records as made using VcfReader.getRecordList()
        popmodel
            an instance of class model
            the individuals in the model must also be in the vcf file
        seq_reference
            a string containing the DNA sequence for the region
        region
            an instance of class Region
        useNifmissingdata
            if True  put N's in wherever data is missing 
        
        Returns
        -------
        a list of sequences 

        Raises
        ------
        ?       
    '''
    recordlist = vcf_reader.getRecordList(region=region)
    # argstemp is an undefined object for passings arguments to generateSequence()
    argstemp = type('',(),{})() 
    argstemp.indel_flag = None
    argstemp.trim_seq = True
    argstemp.N_if_missing = useNifmissingdata
    # list that will hold sequences
    seqlist = []
    for pop in popmodel.pop_list:
        for indiv in popmodel.ind_dict[pop]:
            for ai in range(2): # assumes diploidy 
                s = generateSequence(recordlist,seq_reference, region, indiv, ai, argstemp)
                seqlist.append(s)
    return(seqlist)

def get_model_sequences_from_multiple_regions(vcf_filename=None,popmodel=None,
                                     fasta_reference=None,BED_filename = None,
                                    useNifmissingdata = None):

    '''
 ?  not clear how this will work when bed file and vcf file span multiple chromosomes
 
        returns a generator for getting sets of sequences from
        regions given in a BED file for individuals in a model
        e.g. usage:
            a = get_model_sequences_from_multiple_regions(vcf_filename=vcf_filename,
                        popmodel=popmodel,fasta_reference=fasta_reference,
                        BED_filename = BED_file)
            while True:
                try:
                    s = next(a)# Will raise StopIteration if lines is exhausted
                    print(len(s),len(s[0]))
                except StopIteration:
                    break   # end loop 

        Parameters
        ----------
        vcf_filename
            the name of the vcf file
        popmodel
            an instance of class model
            the individuals in the model must also be in the vcf file
        fasta_reference
            a fasta file with one more sequences
            typically these are reference chromosome sequences
        BED_filename
            a sorted BED file giving regions from which to pull sequences
            the first column with the chromosome name must match a name in
            the fasta_reference file

        Returns
        -------
        returns a generator
        subsequent calls using next()  return a list of sequences
    
        
        Raises
        ------
        ?
    ''' 
    if not os.path.isfile(fasta_reference + ".fai"):
        # make an index 
        samtools_stderr = call_samtools_faidx (fasta_reference)
    # pysam.FastaFile object is later used with fetch()
    fasta_access = pysam.FastaFile(fasta_reference)
    vcf_reader = vf.VcfReader(vcf_filename,popmodel=popmodel)
    with open(BED_filename,'r') as bf:
        for line in bf:
            ls = line.split()
            if len(ls) >= 3:  # anything less than 3 is not a valid line
                # make an instand of VcfReader 
##                vcf_reader = vf.VcfReader(vcf_filename,popmodel=popmodel)
                # deal with possibility of chromosome name differences
                BEDchrname = ls[0]
                chr_in_chrome_BED = BEDchrname[0:3] == 'chr'
                if chr_in_chrome_BED:
                    if not vcf_reader.chr_in_chrom:
                        vcfchrname = BEDchrname[3:] # remove 'chr' from ls[0]
                else:
                    if vcf_reader.chr_in_chrom:  ## could use  vcf_reader.info_rec.chrom
                        vcfchrname = 'chr' + BEDchrname
                
                # pysam.fetch can use a samtools faidx style region string 
                regionstr = '%s:%d-%d'%(BEDchrname,int(ls[1])+1,int(ls[2]))
                seqstr = fasta_access.fetch(region=regionstr)
                # make an instance of Region 
                region = Region(int(ls[1]),int(ls[2]),vcfchrname)

                slist = get_model_sequences_from_region(vcf_reader=vcf_reader,
                                                        popmodel=popmodel,
                                                        seq_reference=seqstr,
                                                        region=region,
                                                        useNifmissingdata=useNifmissingdata)                
                yield(slist)

# in  Pan_all_hicov_chr22_decrun_missingasref.vcf chromosome names are 22
vcf_filename = "../jhworkfiles/Pan_all_hicov_chr22_decrun_missingasref.vcf"
##vcf_filename = "../jhworkfiles/Pan_all_hicov_chr22_decrun_missingasref.vcf.gz"
# in chr22_hg18.fa chromosome name is 'chr22'
fasta_reference = "../jhworkfiles/chr22_hg18.fa"
model_file = "../jhworkfiles/panmodels.model"
BED_file = "../jhworkfiles/pan_test.bed"

popmodels = read_model_file(model_file)
popmodel = popmodels['3Pop']

a = get_model_sequences_from_multiple_regions(vcf_filename=vcf_filename,
                        popmodel=popmodel,fasta_reference=fasta_reference,
                        BED_filename = BED_file,useNifmissingdata=True)
while True:
    try:
        s = next(a)# should raise StopIteration if generator is exhausted
        print(len(s),len(s[0]))
    except StopIteration:
        break   # end loop 
