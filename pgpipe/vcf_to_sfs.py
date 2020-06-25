"""
   For generating the site frequency spectrum (sfs) from a vcf file.
   
   The sfs is an array with as many dimensions as populations.
   For example, if population samples are in order A,B, C
   then position (i,j,k) of the array refers to the count of SNPs with derived alleles that
   were observed to have a count of i in A, j in B, and k in C
   
   If the sfs is folded then the count in a cell of the sfs is he number
   of SNPs with that combination of minor allele counts.
   
   -can handle an arbitrary number of dimensions (populations)
   -handles downsampling, and reduction of dimensions
   -handles unfolded and folded sfs's
   -handles a BED file
   -can handle an alternative reference genome for rooting, rather than that used in the vcf file
   -primary callable functions:
       build_sfs() - for building an sfs 
       reducesfsdims() - for reducing dimensionality of an sfs by summing across axes
   
"""
import sys
import os
import subprocess
import logging
import argparse
import random
import numpy as np
from functools import reduce
import operator
import pysam
from pgpipe.logging_module import initLogger, logArgs
import pgpipe.vcf_BED_to_seqs as vBs
from pgpipe.model import Model, read_model_file
from pgpipe.genome_region import Region
import pgpipe.vcf_reader_func as vr


def binomial(n,k):
    accum = 1
    for m in range(1,k+1):
        accum = accum*(n-k+m)/m
    return float(accum)

def setup_binomial_coefficients(rsampsize,sampsize):
    b = np.zeros((sampsize+1,sampsize+1),dtype=float)
    for i in range(sampsize+1):
        for j in range(sampsize+1):
            if i>= j:
                b[i][j] = binomial(i,j)
    return b

def odd(num):
    if num % 2:
        return True
    return False

## not used
##def downsampleSFS(rsampsize,sampsize,SFS):
##    """ expect an SFS as a list of length sampsize + 1,  with position i for i gene copies"""
##    b = [[0.0 for i in range(sampsize+1)] for j in range(sampsize+1)]
##    for i in range(sampsize+1):
##        for j in range(sampsize+1):
##            if i>= j:
##                b[i][j] = binomial(i,j)
##    rSFS = [0.0 for i in range(rsampsize+1)]
##    for ai in range(len(SFS)):
##        c = SFS[ai]
##        for ri in range(len(rSFS)):
##            p = 1.0
##            p *= b[ai][ri]*b[sampsize-ai][rsampsize-ri]/ b[sampsize][rsampsize]
##            rSFS[ri]+= c * p
##    # now rescale non-zero values so probability is 1 and 0 cell is 0.0
##    rSFS[0] = 0.0
##    rSFS[rsampsize] = 0.0
##    rsum = sum(rSFS)
##    for ri in range(1,len(rSFS)):
##        rSFS[ri] /= rsum
####    print sum(rSFS)
##    return rsampsize,rSFS

def builddownsamplearrays(rsampsize,sampsize, folded = None):
    
    """ rsampsize is the downsampled size
        sampsize is the full sample size
        considers all possible SNP allele counts (from 0 to sampsize)
        and generates the probability distribution under rsampsize
        returns a list of sampsize+1 np arrays

        downarrays[i] gives the probabilities of different observations,
        in a sample of size rsampsize, given that a SNP was
        observed i times in a sample of size sampsize

        for each population sample, this needs to be called only once 
    """
    b = [[0.0 for i in range(sampsize+1)] for j in range(sampsize+1)]
    for i in range(sampsize+1):
        for j in range(sampsize+1):
            if i>= j:
                b[i][j] = binomial(i,j)
    downarrays = []
    for a in range(sampsize+1):
        temp = np.zeros(rsampsize+1)
        for r in range(rsampsize+1):
            temp[r] = b[a][r]*b[sampsize-a][rsampsize-r]/ b[sampsize][rsampsize]
        if folded:
            dnsim = rsampsize//2 + 1
            foldtemp = []
            for i in range(dnsim):
                foldtemp.append(temp[i] + temp[-(i+1)])
            if not odd(rsampsize): # last element of foldtemp will have been doubled
                foldtemp[rsampsize//2] = temp[rsampsize//2]
            downarrays.append(np.array(foldtemp))
        else:
            downarrays.append(temp)
    return downarrays

def getallelecount(r,popmodel,sampsizes,sfs,altref_access=None,folded=None,downsampsizes = None):

    c = []
    pi = 0
    for pop in popmodel.pop_list:
        adic = vr.getAlleleCountDict(r,popmodel.ind_dict[pop])
        assert r.ref == r.alleles[0]
        temp = [adic[r.alleles[0]],adic[r.alleles[1]]]
        if downsampsizes != None:
            if temp[0] + temp[1] < downsampsizes[pi] : 
                return 1  # too few alleles
        else:
            if temp[0] + temp[1] < sampsizes[pi] : 
                return 1  # too few alleles
        if folded:
            temp = temp if (temp[1] < temp[0]) else [temp[1],temp[0]] # make sure minor allele is in position 1
        pi += 1
        c.append(temp)
    if altref_access != None and not folded:
        altref = altref_access.fetch(r.chrom,r.pos-1,r.pos).upper()
        if not (altref == r.ref.upper() or altref.upper() in r.alts or altref.lower() in r.alts):
            return 2  # altref not found, skip this snp        
        if altref == r.ref.upper() or altref == r.alts[0].upper():
            if altref == r.alts[0]: #switch ref and alt
                altc = [[temp[1],temp[0]] for temp in c]
                c = list(altc)
    return c

def dnprocessSNP(r,popmodel,sampsizes,sfs,altref_access=None,folded=None,downsamplearrays=None,downsamplesizes=None,dndimsrev=None):
    """
      gets an array of counts,c, and then builds a downsampled sfs
        this is then added to the sfs

        the downsampling is tricky 
    """
    def loopfunc(arraylist):
        """
            see  array_product_problem.py
            also https://stackoverflow.com/questions/62313394/incrementing-a-multidimensional-numpy-array-python-with-products-generated-fro
        """
##        m = arraylist[0]
##        for a in arraylist[1:]:
##            ml = []
##            for val in a:
##                ml.append(val*m)
##            m = np.stack(ml)
##        return m
        return reduce(operator.mul, reversed(np.ix_(*reversed(arraylist))))# faster than loop method

    c = getallelecount(r,popmodel,sampsizes,sfs,altref_access=altref_access,folded=folded,downsampsizes = downsamplesizes)
    if c == 1:
        return sfs, 2 # bad alternative reference
    if c == 2:
        return sfs,1 # too few alleles    
    dntemp = []
    i = len(downsamplearrays)-1
    for val in reversed(c):
        dntemp.append(downsamplearrays[i][val[1]]) # val[1] is the alt or derived allele 
        i -= 1
    dntemp = loopfunc(dntemp)
    if sfs.shape != dntemp.shape:
        raise Exception("sfs.shape: {0} dntemp.shape: {1}".format(sfs.shape,dntemp.shape))
    sfs = np.add(sfs,dntemp)
        
    return sfs,0

def processSNP(r,popmodel,sampsizes,sfs,folded = None,altref_access=None):
    """
      gets an array of counts,c, and adds to sfs
         for pop k   c[k] is a list of count of ref alleles and count of alt alleles
         if folded,  ref and alt are ignored and smaller value is listed first

         if altref_access  a new reference allele is obtained from altref_access
         and the counts ordered accordingly

         if the alternate reference allele matches neither the original ref or alt, return 1

         if the total count of alleles is less than the sample size for that population
           (or less than the downsample size if downsamplesizes is not None), return 2 
    """
    c = getallelecount(r,popmodel,sampsizes,sfs,altref_access=altref_access,folded=folded)
    if c == 1:
        return sfs, 2 # bad alternative reference
    if c == 2:
        return sfs,1 # too few alleles
    pos = [val[1] for val in c] # val[1] is the alt or derived allele 
    np.add.at(sfs,tuple(pos),1)
    return sfs,0


def build_sfs(vcffile,popmodel,BEDfilename=None,altreference = None,folded = False,downsamplesizes = None,randomsnpprop = None, seed = None,makeint = True):
    """
        returns a multidimensional numpy array
            for k populations the array has k dimensions
            for populations in order,  the indices of a position in the array are in that order
            e.g. if populations are in order A,B, C
            then position (i,j,k) refers to alleles that were observed to have a count of i in A, j in B, and k in C

        Can handle aribtrary numbers of populations and sample sizes, so long as each
            population has at least a sample size of two
            All vcf handling assumes that inidividuals are diploid at all SNPs.
            Otherwise, if vcf processing were more general then sample sizes could be odd.
            
        vcffile is a vcf, or bgzipped vcf file, with SNP data
        
        popmodel is an instance of Model
            It should specify one or more populations
            sample sizes are determined by the numbers of individuals in the model for each population

        
        BEDfilename is the name of a ucsc-style bedfile with intervals to include

        altreference is the name of a fasta sequence file that contains the reference genome to be used
            this causes the 'ref' allele, as given in the vcf to not be used as the reference
            this can be useful, for example, if an ancestral reference is available to that the reference allele is
            the ancestral allele

            if the base from the alternative references does not match either the vcf reference base or the vcf first alternate base
            the SNP will be ignored 

            using an alternate reference has no effect if folded is True 

        folded indicates that the folded sfs should be returned
            folded causes the count returned for a SNP to be that for the less common base
            ignores alt and ref

        downsamplesizes is an array listing the sample sizes to be used if they are less than given in the model
            2 <= downsamplesizes[i] <= samplesizes[i]
            if None,  then the sample sizes are those given by the popmodel

        randomsnpprop is the proportion of snps to include
            uses random sampling

        seed is a random number seed that can be used with randomsnpprop

        makeint causes the array to be rounded to the nearest integer (dtype remains float)
            == True by default
                        
        
        
    """
    npops = len(list(popmodel.pop_list))
    sampsizes = []
    for pop in popmodel.pop_list:
        sampsizes.append(2*len(popmodel.ind_dict[pop]))

    dodown = downsamplesizes != None
    if dodown:
        if  len(downsamplesizes) != npops:
            raise Exception("downsamplesizes length not the same as # of populations in model")
        sfsdims = [temp + 1 for temp in downsamplesizes]
        # popdnarrays[i] is a list of arrays for population i
        # for each possible count of alleles given the full sample size
        # there is an array that gives the probability distribution of values under the down sampled size
        popdnarrays = []
        for i in range(npops):
            popdnarrays.append(builddownsamplearrays(downsamplesizes[i],sampsizes[i],folded=folded))
    else:
        popdnarrays = None

    if dodown:
        if folded:
            dndims = []
            for s in downsamplesizes:
                dndims.append(s//2 + 1)
        else:
            dndims = [temp + 1 for temp in downsamplesizes]
        dndimsrev = list(dndims)
        dndimsrev.reverse()
        sfs = np.zeros(tuple(dndims),dtype = float)
    else:
        if folded:
            sfsdims = []
            for s in sampsizes:
                sfsdims.append(s//2 + 1)
        else:
            sfsdims = [temp + 1 for temp in sampsizes]        
        sfs = np.zeros(tuple(sfsdims),dtype = float)
        dndimsrev = None
    if altreference != None:
        if not os.path.isfile(altreference + ".fai"):
        # make an index 
            pysam.faidx(altreference)
        altref_access = pysam.FastaFile(altreference)
    else:
        altref_access = None
            
    vcf_reader = vr.VcfReader(vcffile,popmodel=popmodel)

    if randomsnpprop != None and seed != None:
        random.seed(seed)
    numSNPs = 0
    numtoofew = 0
    numbadaltref = 0
    if BEDfilename == None:
        while True:
            r = vcf_reader.getNext()
            if type(r) == type(None):
                break
            if vr.checkRecordPass(r,remove_indels=True,remove_multiallele=True,remove_missing=10):
                if randomsnpprop  == None or random.random() < randomsnpprop :
                    if dodown:
                        sfs,snpcode = dnprocessSNP(r,popmodel,sampsizes,sfs,altref_access=altref_access,folded=folded,downsamplearrays=popdnarrays,downsamplesizes=downsamplesizes,dndimsrev=dndimsrev)
                    else:
                        sfs,snpcode = processSNP(r,popmodel,sampsizes,sfs,altref_access=altref_access,folded=folded)
                    numbadaltref += snpcode==1
                    numtoofew += snpcode==2
                   
    else:
        with open(BEDfilename,'r') as bf:
            for line in bf:
                ls = line.split()
                if len(ls) >= 3:  # anything less than 3 is not a valid line
                    region = Region(int(ls[1]),int(ls[2]),ls[0])
                else:
                    break # get out of loop
                recordlist = vcf_reader.getRecordList(region=region)
                ri = 0
                for r in recordlist:
                    if vr.checkRecordPass(r,remove_indels=True,remove_multiallele=True,remove_missing=10):
                        if randomsnpprop == None or random.random() < randomsnpprop :
                            if dodown:
                                sfs,snpcode = dnprocessSNP(r,popmodel,sampsizes,sfs,altref_access=altref_access,folded=folded,downsamplearrays=popdnarrays,downsamplesizes=downsamplesizes,dndimsrev=dndimsrev)
                            else:
                                sfs,snpcode = processSNP(r,popmodel,sampsizes,sfs,altref_access=altref_access,folded=folded)
                            numbadaltref += snpcode==2
                            numtoofew += snpcode==1
                        ri += 1
    if makeint:
        sfs = np.rint(sfs)
    return sfs
        

def reducesfsdims(sfs,popmodel,keeppops):
    """
        sfs is a full sfs for popmodel that was made using build_sfs() with axes
        ordered the same as populations in popmodel.
        keeppops is a list containing a subset of populations in popmodel 
        returns the sfs after summing along axes not in keeppops
    """
    temp = list(popmodel.pop_list)
    for kp in keeppops:
        temp.remove(kp)
    sumaxes = tuple([popmodel.pop_list.index(kp) for kp in temp])
    sfs = np.sum(sfs,sumaxes)
    return sfs

        
