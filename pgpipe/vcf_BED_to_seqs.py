#!/usr/bin/env python
'''
    For generating a file with reconstituted sequences for regions of a vcf file.

    Given a vcf file, a corresponding fasta reference file, a population model,
    and a specific interval, this script will return a file containing reconstituted
    sequences - two for each diploid individual.

    This script can be run in stand-alone mode,  or for more flexibility it can
    be imported to give access to more general functions for generating
    reconstituted sequences. 


    ###############
    Required Arguments
    ###############

    **--vcf** *<input_vcf_filename>*
        The name of the vcf file.  This can be a bgzipped vcf file. . 

    **--model-file** *<model_file_name>*
        The name of a PPP model file. 

    **--modelname** *<model_name>*
        The name of a model in the model file.

    **--fasta-reference** *<reference fasta file>*
        The reference genome fasta file is required in order to generate full
        sequences from the SNP data in the vcf file.
    **--region** *<region string>*
        The region of the reference to return sequences from.
        Format: chromosome name,colon (:),first base number,dash(-),last base number.
        e.g. chr1:100392-101391"

    **--out** *<out file name>*
        - the name of an output file. If --out is omitted the default is ppp_vcf_to_sequences.out


    ###############
    Optional Aguments
    ###############
    **--return-single** *<True/False>*
        If true,  only a single sequence is returned for each individual in the model.
        The sequence for a given individual uses the first allele given for that
        individual for each SNP in the vcf file. 


    ###############
    Importing Functions
    ###############   

    This file has two scripts that can be useful for working with site frequency spectra
   
    1. get_model_sequences_from_region()
        Returns a list of sequences for a region in a vcf file.
     
        Parameters
        ----------
        -vcf  either:
            -a vcf_reader  (i.e. see vf.VcfReader())
            -or the name of a vcf file (can be bgzipped)
        -popmodel :an instance of class model. The individuals in the model must also be in the vcf file
        -seq_reference: A string that either contains the DNA sequence for the region
            or is a string containing the name of a fasta file. If a fasta file name,
            the chromosome name(s) in the fasta file must match those in the vcf file 
        - region : either a
            samtools/pysam style region string ("chromosome name:start-end")
            where the chromosome name matches that used in the vcf file
            where start and end are the 1-based endpoints (closed interval)
            Or an instance of class Region (Regions uses 0-based open interval on the left)
        -return_single: if True,  return only the first sequence for an individual
            else, return two sequences  (assumes diploid)
        -out :the name of a file to contain the sequences (optional)

        Returns
        -------
        -a samtools/pysam-style region string ("chromosome name:start-end")
            if a region string was passed to this function,  this should be the same
            but if a Region was passed,  then it will be new
        - a list of sequences as strings

    2. get_model_sequences()
        Returns a generator for getting sets of sequences from
        regions given in a BED file for individuals in a model.
        
        *************
         Example usage
        *************
        Extend flanks (i.e. both upstream and downstream) by 10kb:

        .. code-block:: python
        :linenos:
        
            a = get_model_sequences_from_multiple_regions(vcf_filename=vcf_filename,
                        pppmodel=popmodel,fasta_reference=fasta_reference,
                        BED_filename = BED_file)
            while True:
                try:
                    s = next(a)# Will raise StopIteration if lines is exhausted
                    print(len(s),len(s[0]))
                except StopIteration:
                    break   # end loop 

        Parameters
        ----------
        -vcf_filename: the name of the vcf file (can be bgzipped)
        -model_file : name of a model file
        -modelname: the name of a model in model_file
            Either pppmmodel or both model_file and modelname must be used
        -pppmodel :   an instance of class model
            the individuals in the model must also be in the vcf file
            Either pppmmodel or both model_file and modelname must be used
        -fasta_reference : a fasta file with one more sequences corresponding to the vcf file
            typically these are reference chromosome sequences
            the chromosome name(s) in the fasta file must
            match those in the vcf file
        -BED_filename : a sorted BED file giving regions from which to pull sequences
            the first column with the chromosome name must match a name in
            the fasta_reference file
            The chromosome names in the BED file must match those in the fasta file
            and the vcf file
        -region_string : a pysam style region string,  if only one region is to be returned
        -return_single : if True,  return only the first sequence for an individual
            else, return two sequences  (assumes diploid)

        Returns
        -------
        returns a generator
        subsequent calls using next()  return a samtools/pysam-style region string
        and a list of sequences
            
'''

import sys
import os
import subprocess
import logging
import argparse
import pgpipe.vcf_reader_func as vf
from pgpipe.logging_module import initLogger, logArgs
from pgpipe.model import Model, read_model_file
from pgpipe.misc import confirm_executable
from pgpipe.vcf_to_ima import generateSequence
from pgpipe.genome_region import Region, RegionList
import pysam




def xor(val1,val2):
    """
       makes booleans of val1 and val2 and then return bitwise xor  (i.e. True/False)
       note - bool(None)==False 
    """
    return bool(val1) ^ bool (val2)

def get_model_sequences_from_region(vcf=None,popmodel=None,
    seq_reference=None, region=None,return_single=True,
    out = None, called_from_run = False):
    """
        basically a wrapper for  vcf_to_ima.generateSequence()
        (note - there is also a generateSequence() in vcf_ref_to_seq() ) 
        
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
                 
                 or a string containing the name of a fasta file
                 if a fasta file, this works better if the chromosome name(s) in the fasta file
                 match those in the vcf file 
            
        region
            a valid interval for the vcf file
            either:
                 a samtools/pysam style region string ("chromosome name:start-end")
                     where the chromosome name matches that used in the vcf file
                     where start and end are the 1-based endpoints (closed interval)
                 or an instance of class Region
                     Region uses 0-based open interval on the left
        return_single
            if True,  return only the first sequence for an individual
            else, return two sequences  (assumes diploid)


        out is the name of a file to contain the sequences
            

        called_from_run is true only if build_sfs() was called from run()
            in this case the out file name is set to a default value

        
        Returns
        -------
        a samtools/pysam-style region string ("chromosome name:start-end")
            if a region string was passed to this function,  this should be the same
            but if a Region was passed,  then it will be new
        a list of sequences
        
    """
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
    argstemp.N_if_missing = True


    # seqlist is a list that will hold sequences
    seqlist = []
    
    for pop in popmodel.pop_list:
        for indiv in popmodel.ind_dict[pop]:
            numseqs = 1 if return_single else 2
            for ai in range(numseqs): # assumes diploidy 
                s = generateSequence(recordlist,seqstr, region, indiv, ai, argstemp)
                seqlist.append(s)


    if called_from_run or out != None: #write/append sequences to a file 
        if out == None:
            out = os.path.dirname(vcffile) + "//ppp_vcf_to_sequences.out"
        if os.path.exists(out):
            f = open(out,'a')
        else:
            f = open(out,'w')
        f.write("{}\n".format(regionstr))
        i = 0
        for pop in popmodel.pop_list:
            for indiv in popmodel.ind_dict[pop]:
                numseqs = 1 if return_single else 2
                for ni in range(numseqs):
                    f.write("{}:{}\t{}\n".format(indiv,ni,seqlist[i]))
                    i += 1
        f.close()
    return(regionstr,seqlist)

def get_model_sequences(vcf_filename=None,model_file = None,modelname=None,
        pppmodel=None,fasta_reference=None,BED_filename = None,
        return_single=True):

    

    """
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
        model_file : name of a model file
        modelname: the name of a model in model_file
        pppmodel
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
            
           
        return_single
            if True,  return only the first sequence for an individual
            else, return two sequences  (assumes diploid)

        Returns
        -------
        returns a generator
        subsequent calls using next()  return a samtools/pysam-style region string
        and a list of sequences
    
        
    """ 
    if not os.path.isfile(fasta_reference + ".fai"):
        # make an index 
        pysam.faidx(fasta_reference) #easier to call this rather than samtools
    # pysam.FastaFile object is later used with fetch()
    fasta_access = pysam.FastaFile(fasta_reference)
    #get model
    if pppmodel and (model_file or modelname):
        raise Exception("must specify _either_ pppmodel or both model_file and modelname")
    if xor(model_file,modelname):
        raise Eception("both model_file and modelname must be given, unless pppmodel is being used")
    if not pppmodel:
        popmodels = read_model_file(model_file)
        popmodel = popmodels[model]
    else:
        popmodel = pppmodel
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
                        region=region,return_single=return_single)                
                yield(regionstr,slist)



def parser(passed_arguments):
    """snfs Argument Parser - Assigns arguments from command line"""

    def parser_confirm_file ():
        """Custom action to confirm file exists"""
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError('%s not found' % value)
                setattr(args, self.dest, value)
        return customAction

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--vcf', help = "Defines the filename of the VCF", type = str, required = True, action = parser_confirm_file())
    parser.add_argument('--model-file', help = 'Defines the model filename',required = True, type = str, action = parser_confirm_file())
    parser.add_argument('--modelname', help = 'The name of the model for which sequences are to be returned', required = True,type = str)
    parser.add_argument('--fasta-reference', help = 'Defines a fasta referemce filename', type = str, action = parser_confirm_file())
    parser.add_argument('--region',help='the region of the reference to return sequences from'
                        "format: chromosome name,colon (:),first base number,dash(-),last base number"
                        "e.g. chr1:100392-101391",type=str)
    parser.add_argument('--out', help = 'Optional, the complete output filename of a tab-delimted file', type = str)
    parser.add_argument('--return_single',default=False,help="Optional, return just one sequence per individual")


    if passed_arguments:
        return parser.parse_args(passed_arguments)
    else:
        return parser.parse_args()

def run (passed_arguments = []):
    # get arguments from command line
    args = parser(passed_arguments)

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(args, func_name = 'get_model_sequences_from_region')

    print(args)

    popmodel = read_model_file(args.model_file)[args.modelname]

    if not os.path.isfile(args.fasta_reference + ".fai"):
        # make an index 
        pysam.faidx(args.fasta_reference) #easier to call this rather than samtools
    fasta_access = pysam.FastaFile(args.fasta_reference)
    seqstr = fasta_access.fetch(region=args.region)
       
    get_model_sequences_from_region(vcf=args.vcf,popmodel=popmodel,
        seq_reference=seqstr, region=args.region,return_single=args.return_single,
        called_from_run=True,out=args.out)
    
    
    
if __name__ == "__main__":
    initLogger()
##    run()
    debugargs = ['--vcf','..//jhworkfiles//Pan_chr_21_22_test.vcf.gz','--fasta-reference',
                 "..//jhworkfiles//twochr_test_ref.fa",
            '--model-file',"..//jhworkfiles//panmodels.model",'--modelname',"4Pop",
            '--region',"21:4431001-4499000",
            '--out','..//jhworkfiles//test_vcf_BED_to_seqs.out']    
    run(debugargs)
