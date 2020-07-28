#!/usr/bin/env python
'''
For generating the site frequency spectrum (sfs)for a
population model from a vcf file.

The sfs is an array with as many dimensions as populations in the model.
For example, if population samples are in order A,B, C
then position (i,j,k) of the array refers to the count of SNPs with derived alleles that
were observed to have a count of i in A, j in B, and k in C

If the sfs is folded then the count in a cell of the sfs is the number
of SNPs with that combination of minor allele counts.

This script can be run in stand-alone mode,  or for more flexibility it can
be imported to give access to more general functions for building and
manipulating the sfs. 


All vcf handling assumes that all individuals are diploid at all SNPs.

##################
Required Arguments
##################

**--vcf** *<input_vcf_filename>*
    The name of the vcf file.  This can be a bgzipped vcf file. . 

**--model-file** *<model_file_name>*
    The name of a PPP model file. 

**--modelname** *<model_name>*
    The name of a model in the model file.  The treemix file to be
    generated will contain the allele counts for each SNP in each of
    the populations.  The treemix run will estimate the phylogeny
    for the populations in the model.

**--out** *<out file name>*
    The name of an output file. If --out is omitted the default is ppp_sfs.out in
    the same folder as the vcffile This will be a tab-delimited file
    If the number of dimensions is 2, the sfs is contained in the rows and columns,
    otherwise the values are given on the first line of the file

#################
Optional Aguments
#################

  
**--bed-file** *<BED_file_name>*
    The BED file is a sorted UCSC-style bedfile containing chromosome locations of
    the SNPs to be included in the output files. The BED file has no header.
    The first column is the chromosome name (this must match the chromosome
    name in the vcf file).
    The second column is start position (0-based, open interval)
    The third column is end position (closed interval).
    Any other columns are ignored.

**--outgroup_fasta** *<name of alternative reference sequence>*
    This option is used to specify the name of a fasta file to use as an
    alternative reference to that originally used for the vcf file.
    

    This fasta file must have been properly aligned to the reference
    used in the vcf file.  
        
    This option can be useful, for example, if an ancestral or outgroup reference is
    available that more accurately identifies the ancestral (and thus derived)
    allele at each SNP than does the reference used to make the vcf file.

**--downsamplesizes** *<down sample sizes>*
    A sequence of integers,  one for each of the populations in the model in
    the same order as populations listed in the model. The values 
    specify the down sampling to be used for each respective population.
    For a population with k>=1 diploid individuals (2k>=2 genomes) in the model,
    the downsample count d  must be 2<=d<=2k.

**--folded** *<True/False>*
    The folded option indicates that the folded sfs should be returned.
    If folded is False (default) the sfs reports the count of the derived allele.
    If True,  the sfs reports of the count of the minor (less frequent) allele.

**--randomsnpprop** *<floating point value between 0 and 1>*
    This option can be used to randomly sample a subset of SNPs. The default
    is to sample all biallelic SNPs.

**--seed** *<integer>*
    This is used with --randomsnpprop as the seed for the random number generator.

**--makeint** *<True/False>*
    If True, round the counts in the sfs to the nearest integer (default False)


#############
Example usage
#############
Example command-lines:

.. code-block:: bash

    vcf_to_sfs.py -h

   
.. code-block:: bash

    vcf_to_sfs.py --vcf pan_example2.vcf.gz --model-file panmodels.model --modelname 4Pop --downsamplesizes 3 3 3 4 --folded --outgroup-fasta  chr22_pan_example2_ref.fa --out vcf_to_sfs_test1.txt

###################
Importing Functions
###################

    This file has two functions that can be useful for working with site frequency spectra

===========
build_sfs()
===========

The default script that runs when this file is run. This function can also be
accessed directly by importing this file.

    * generates an sfs from a vcf file
    * can handle an arbitrary number of dimensions (populations) and
      sample sizes, so long as each population has at least a
      sample size of two genomes (i.e. one diploid individual)
    * handles downsampling, and reduction of dimensions
    * handles unfolded and folded sfs's
    * can take a BED file name to sample portions of a vcf file
    * can handle an alternative reference genome for rooting, rather
      than that used in the vcf file
    * the arguments closely resemble those used when the function is called
      by running this file

------------------
Required Arguments
------------------

    The first three arguments are, in order, the vcf filename,  the model file name, and the name of the model.
    These are each strings.
    

------------------
Optional arguments
------------------

    Each optional argument requires the use of the argument name.

**BEDfilename**
    The name of a ucsc-style bedfile with intervals to include
**altreference**
    The name of a fasta sequence file that contains the reference genome
**folded**
    True/False. To indicate that the folded sfs should be returned. (False is default)
**downsamplesizes**
    A list of sample sizes to be used if they are
    less than given in the model 2 <= downsamplesizes[i] <= samplesizes[i]
**randomsnpprop**
    The proportion of snps to include using random sampling
**seed**
    A random number seed that can be used with randomsnpprop.
**makeint**
    Causes the array to be rounded to the nearest integer (dtype remains float)
**out**
    The name of a file to contain the sfs
    if out is not None,  this will write a tab-delimited file of the array

-------------
Example usage
-------------
Example python code:

.. code-block:: python
    :linenos:
    
    import vcf_to_sfs as vs
    mysfs = vs.build_sfs(pan_example2.vcf.gz,panmodels.model,'4Pop',
            folded=True,downsamplesizes=[3,3,3,4],
            altreference='chr22_pan_example2_ref.fa',
            out = 'mysfsfile.txt')

=================    
reduce_sfs_dims()
=================
    This function is for reducing the dimensionality of an sfs by summing across axes. It is
    accessed by importing this file. 
    
    There are three required arguments in order:
    
    * the sfs (i.e. a numpy array with as many dimensions as populations)
    * an instance of Class model that specifies the populations and samples to which the
      sfs corresponds.
    * a list of names of the populations to keep in the reduced sfs

    There is one optional argument, 'out', if the reduced sfs is to be written to a file (e.g. out=mysfs.txt)

-------------
Example usage
-------------
Example python code:

.. code-block:: python
    :linenos:

    import vcf_to_sfs as vs
    myreducedsfs = vs.reduce_sfs_dim(mysfs,mypopmodel,
            ['A','B','C'],out="myreducedsfs_file.txt")
            
        
'''

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
from pgpipe.model import Model, read_model_file
from pgpipe.genome_region import Region
import pgpipe.vcf_reader_func as vr
from pgpipe.misc import argprase_kwargs


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

def downsampleSFS(rsampsize,sampsize,SFS):
    """
        standalone function, not called by build_sfs()
        
        expect a 1D SFS as a list of length sampsize + 1,  with position i for i gene copies
        returns a list with the reduced sfs
        
    """
    b = [[0.0 for i in range(sampsize+1)] for j in range(sampsize+1)]
    for i in range(sampsize+1):
        for j in range(sampsize+1):
            if i>= j:
                b[i][j] = binomial(i,j)
    rSFS = [0.0 for i in range(rsampsize+1)]
    for ai in range(len(SFS)):
        c = SFS[ai]
        for ri in range(len(rSFS)):
            p = 1.0
            p *= b[ai][ri]*b[sampsize-ai][rsampsize-ri]/ b[sampsize][rsampsize]
            rSFS[ri]+= c * p
    # now rescale non-zero values so probability is 1 and 0 cell is 0.0
    rSFS[0] = 0.0
    rSFS[rsampsize] = 0.0
    rsum = sum(rSFS)
    for ri in range(1,len(rSFS)):
        rSFS[ri] /= rsum
##    print sum(rSFS)
    return rSFS

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


def reduce_sfs_dims(sfs,popmodel,keeppops, out=None):
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

    if out:
        # detect if array as been rounded to zero decimals, or very close
        temp = sfs.astype(int)
        if np.all(np.isclose(sfs, temp, 1e-8)):
            outfmt = "%d"
        else:
            outfmt = "%.1f"
        if len(sfs.shape) <= 2:
            np.savetxt(out,sfs,fmt = outfmt,delimiter ='\t')
        else:
            sfs.tofile(out,sep='\t',format=outfmt)            
    
    return sfs

def build_sfs(vcffile,model_file,model,BEDfilename=None,altreference = None,folded = False,
              downsamplesizes = None,randomsnpprop = None, seed = None,
              makeint = True, out = None, called_from_run = False):
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
            True by default

        out is the name of a file to contain the sfs
            if out is not None,  this will write a tab-delimited file of the array

        called_from_run is true only if build_sfs() was called from run()
            in this case the out file name is set to a default value 
        
    """


    popmodels = read_model_file(model_file)
    popmodel = popmodels[model]
    
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
        outfmt = "%d"
    else:
        outfmt = "%.1f"
    if called_from_run or out != None:
        if out == None:
            out = os.path.dirname(vcffile) + "//ppp_sfs.out"
        if len(sfs.shape) <= 2:
            np.savetxt(out,sfs,fmt = outfmt,delimiter ='\t')
        else:
            sfs.tofile(out,sep='\t',format=outfmt)
    return sfs
        

def sfs_parser(passed_arguments=[]):
    '''snfs Argument Parser - Assigns arguments from command line'''

    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError('%s not found' % value)
                setattr(args, self.dest, value)
        return customAction

    sfs_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sfs_parser.add_argument('--vcf', help = "Defines the filename of the VCF", type = str, required = True, action = parser_confirm_file())
    sfs_parser.add_argument('--model-file', help = 'Defines the model filename',required = True, type = str, action = parser_confirm_file())
    sfs_parser.add_argument('--modelname', help = 'Defines the model and the individual(s)/population(s) to include', required = True,type = str)
    sfs_parser.add_argument("--bed-file",help="Optional, the name of a BED file", type = str, action = parser_confirm_file())
    sfs_parser.add_argument('--folded',default=False,action="store_true",help="Optional, generate the folded sfs")
    sfs_parser.add_argument('--randomsnpprop',type=float,help="Optinal,the proportion of randomly selected SNPs to include")
    sfs_parser.add_argument('--seed',type=int,help="Optinal,integer random number seed to use with randomsnpprop")
    sfs_parser.add_argument('--downsamplesizes',type=int,nargs='+',help="Optional,sample sizes to use for output "
                                "an array of integers,  one for each population, in the same order "
                                "as populations in popmodel. Values must >=2 and be <= actual number of "
                                "chromosomes in the vcf file for the corresponding population.")    
    sfs_parser.add_argument('--out', help = 'Optional, the complete output filename of a tab-delimted file', type = str)
    sfs_parser.add_argument('--outgroup-fasta',help="Optional, the name of a fasta format file containing"
                                 " the ancestral or outgroup sequence, by default the 'ref' allele "
                                 "of the vcf file is treated as the outgroup", action = parser_confirm_file())
    sfs_parser.add_argument('--makeint',help='Optional, round all values to zero decimal places', default = True)

    if passed_arguments:
        return vars(sfs_parser.parse_args(passed_arguments))
    else:
        return vars(sfs_parser.parse_args())

def run (**kwargs):
    # Update kwargs with defaults
    if __name__ != "__main__":
        kwargs = argprase_kwargs(kwargs, sfs_parser)
    # Assign arguments
    sfs_args = argparse.Namespace(**kwargs)    


    # Adds the arguments (i.e. parameters) to the log file
    logArgs(sfs_args, func_name = 'build_sfs')

##    print(sfs_args)

    build_sfs(sfs_args.vcf,sfs_args.model_file,sfs_args.modelname,BEDfilename=sfs_args.bed_file,
              altreference = sfs_args.outgroup_fasta,folded = sfs_args.folded,
              downsamplesizes = sfs_args.downsamplesizes,
              randomsnpprop = sfs_args.randomsnpprop, seed = sfs_args.seed,
              makeint = sfs_args.makeint, out = sfs_args.out, called_from_run = True)
    
    
if __name__ == "__main__":
    initLogger()
    run(**sfs_parcer())
    exit()

    debugargs=['--vcf',"..//jhtests//pan_example2.vcf.gz",
               '--model-file',"..//jhtests//panmodels.model",'--modelname','4Pop',
               '--downsamplesizes','3','3','3','4',
               '--folded','--outgroup-fasta',"..//jhtests//chr22_pan_example2_ref.fa",
               '--out',"..//jhtests//results//vcf_to_sfs_test1.txt"]
    run(debugargs)
    debugargs=['--vcf',"..//jhtests//pan_example.vcf.gz",
               '--model-file',"..//jhtests//panmodels.model",'--modelname','5Pop',
               '--downsamplesizes','3','3','3','4','2',
               '--folded','--outgroup-fasta',"..//jhtests//pan_example_ref.fa",
               '--out',"..//jhtests//results//vcf_to_sfs_test2.txt"]
    
##    debugargs=['--vcf',"..//jhtests//pan_example2.vcf.gz",
##               '--model-file',"..//jhtests//panmodels.model",'--model','2Pop',
##               '--outgroup-fasta',"..//jhtests//chr22_pan_example2.fa",
##               '--out',"..//jhtests//test_vcf_to_sfs_2pop.txt"]
    
##    debugargs = ['--vcf','pan_example.vcf.gz','--reference',"pan_example2_ref.fa",
##            '--model-file',"panmodels.model",'--model',"4Pop",
##            '--bed-file',"twochr_test.bed",'--out','testgphocsparser.out']#,'--diploid','False','--nloci','4']
    run(debugargs)
