import sys
import os
import subprocess
import logging
import pgpipe.vcf_reader_func as vf
from pgpipe.model import Model
from pgpipe.misc import confirm_executable
from pgpipe.vcf_to_ima import generateSequence
from pgpipe.genome_region import Region, RegionList
import pysam






def get_model_sequences_from_region(vcf=None,popmodel=None,
    seq_reference=None, region=None,return_one_sequence=True,useNifmissingdata=None):
    
    '''
        basically a wrapper for  vcf_to_ima.generateSequence()
        
        gets a list of sequences for individuals in a population model
        for a particular chromosomal region from a vcf file 

        This can operated on its own, or when called by get_model_sequences_from_multiple_regions()

 ? as of 6/2/2020 can't handle haploidy or mixed ploidy among individuals
        
        Parameters
        ----------
        vcf
            is either:
                a vcf_reader  (i.e. see vf.VcfReader())
                or the name of a vcf file
                
        popmodel
            an instance of class model
            the individuals in the model must also be in the vcf file
            
        seq_reference
            a string that either:
                 contains the DNA sequence for the region
                 
                 or a string containting the name of a fasta file
                 if a fasta file, this works better if the chromosome name(s) in the fasta file
                 match those in the vcf file 
            
        region
            either:
                 a samtools/pysam style region string ("chromosome name:start-end")
                     where the chromosome name matches that used in the vcf file
                     where start and end are the 1-based endpoints (closed interval)
                 or an instance of class Region
                     Region uses 0-based open interval on the left
        return_one_sequence
            if True,  return only the first sequence for an individual
            else, return two sequences  (assumes diploid) 
                 
        useNifmissingdata
            if True  put N's in wherever data is missing 
        
        Returns
        -------
        a samtools/pysam-style region string ("chromosome name:start-end")
            if a region string was passed to this function,  this should be the same
            but if a Region was passed,  then it will be new
        and a list of sequences 

        Raises
        ------
        ?       
    '''
    if not isinstance(vcf,vf.VcfReader): # if vcf is a vcf file,  create a vcf_reader from it
        assert isinstance(vcf,str)
        vcf_reader = vf.VcfReader(vcf,popmodel=popmodel)
    else:
        vcf_reader = vcf  # vcf was already a VcfReader

    if not isinstance(region,Region): # if region is not a Region, then make one from it
        assert isinstance(region,str)
        regionstr = region
        chrname = regionstr[:regionstr.find(':')]
        [start,end] = map(int,regionstr[regionstr.find(':')+1:].split(sep='-'))
        region = Region(start-1,end,chrname)
    else: # make regionstr
        regionstr = '%s:%d-%d'%(region.chrom,region.start+1,region.end)# must add 1 to start because Region's use half open 0-based
        
##    print('x',region.start,region.end,region.chrom)
    if os.path.isfile(seq_reference): #seq_reference is the name of a fasta file
        fasta_access = pysam.FastaFile(seq_reference)
        seqstr = fasta_access.fetch(region=regionstr)
    else:  # seq_reference was just a sequence string, but check to see if string looks like DNA
        seqstr = seq_reference
        dnaset = set("agctyrwskmdvhbxnAGCTYRWSKMDVHBXN")
        nondnachar =  set(seqstr)-dnaset
        if len(nondnachar)> 0:
            raise Exception("sequence has non DNA characters:%s"%set(nondnachar))

    
##    print(region.chrom,region.start,region.end)
    recordlist = vcf_reader.getRecordList(region=region)
    
    # argstemp is an undefined object for passings arguments to generateSequence()
    argstemp = type('',(),{})() 
    argstemp.indel_flag = None
    argstemp.trim_seq = True
    argstemp.N_if_missing = useNifmissingdata

    # seqlist is a list that will hold sequences
    seqlist = []
    
    for pop in popmodel.pop_list:
        for indiv in popmodel.ind_dict[pop]:
            numseqs = 1 if return_one_sequence else 2
            for ai in range(numseqs): # assumes diploidy 
                s = generateSequence(recordlist,seqstr, region, indiv, ai, argstemp)
                seqlist.append(s)
    return(regionstr,seqlist)

def get_model_sequences_from_multiple_regions(vcf_filename=None,popmodel=None,
        fasta_reference=None,BED_filename = None,
        return_one_sequence=True,useNifmissingdata = None):

    '''
        returns a generator for getting sets of sequences from
        regions given in a BED file for individuals in a model
        seems to work when vcf and BED span multiple chromosomes
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
            the chromosome name(s) in the fasta file must
            match those in the vcf file
            
        BED_filename
            a sorted BED file giving regions from which to pull sequences
            the first column with the chromosome name must match a name in
            the fasta_reference file
            The chromosome names in the BED file must match those in the fasta file
            and the vcf file
            
        return_one_sequence
            if True,  return only the first sequence for an individual
            else, return two sequences  (assumes diploid) 

        Returns
        -------
        returns a generator
        subsequent calls using next()  return a samtools/pysam-style region string and a list of sequences
    
        
        Raises
        ------
        ?
    ''' 
    if not os.path.isfile(fasta_reference + ".fai"):
        # make an index 
        pysam.faidx(fasta_reference) #easier to call this rather than samtools
    # pysam.FastaFile object is later used with fetch()
    fasta_access = pysam.FastaFile(fasta_reference)
    # make an instance of VcfReader 
    vcf_reader = vf.VcfReader(vcf_filename,popmodel=popmodel)
    with open(BED_filename,'r') as bf:
        for line in bf:
            ls = line.split()
            if len(ls) >= 3:  # anything less than 3 is not a valid line
                # deal with possibility of chromosome name differences
                chrname = BEDchrname = ls[0]
                chr_in_chrome_BED = BEDchrname[0:3] == 'chr'
                if chr_in_chrome_BED:
                    if not vcf_reader.chr_in_chrom:
                        raise Exception("\"chr\" mismatch, check fasta file, BED file and vcf file")
                else:
                    if vcf_reader.chr_in_chrom:  ## could use  vcf_reader.info_rec.chrom
                        raise Exception("\"chr\" mismatch, check fasta file, BED file and vcf file")                        
                
                # pysam.fetch can use a samtools faidx style region string 
                regionstr = '%s:%d-%d'%(BEDchrname,int(ls[1])+1,int(ls[2]))
                seqstr = fasta_access.fetch(region=regionstr)
                # make an instance of Region 
                region = Region(int(ls[1]),int(ls[2]),chrname)

                regionstr,slist = get_model_sequences_from_region(vcf=vcf_reader,
                        popmodel=popmodel,seq_reference=seqstr,
                        region=region,return_one_sequence=return_one_sequence,
                        useNifmissingdata=useNifmissingdata)                
                yield(regionstr,slist)


