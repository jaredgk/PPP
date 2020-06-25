#!/usr/bin/env python
'''
   Generates a Site Frequency Spectrum (SFS) fils for fastsimcoal
   based on instructions in fastsimcoal ver 2.6 manual.

   Generates one-dimensional (1D), two-dimensional (2D) and multidimensional SFS files

   
   Excoffier, L. and M. Foll. 2011. fastsimcoal: a continuous-time coalescent
   simulator of genomic diversity under arbitrarily complex evolutionary scenarios.
   Bioinformatics 27: 1332-1334.
   


    ###############
    Required Arguments
    ###############
    **--vcf** *<input_vcf_filename>*
        The name of the vcf file.  This can be a bgzipped vcf file. . 

    **--model-file** *<model_file_name>*
        The name of a PPP model file. 

    **--model** *<model_name>*
        The name of a model in the model file.  The treemix file to be
        generated will contain the allele counts for each SNP in each of
        the populations.  The treemix run will estimate the phylogeny
        for the populations in the model.

    **--dim** *<dimension file type signifiers>*
        One or more of '1', '2', or 'm', for 1D, 2D or multidimensional output files.
        
        For 1D files:
           -the filename suffix is _DAFpop#.obs for an array of derived allele counts.
              where '#' is replaced by the population number
           -the filename suffix is _MAFpop#.obs for an array of minor allele counts.
        For 2D files:
           -the filename suffix is _jointDAFpop#_&.obs for an array of derived allele counts.
               where # and & are population numbers, and # is larger than &
           -the filename suffix is _jointDAFpop#_&.obs for an array of minor allele counts.
        For a multidimensional file:
           - the filename suffix is _DSFS.obs for an array of derived allele counts.
           - the filename suffix is _MSFS.obs for an array of minor allele counts.
           
    ###############
    Optional Aguments
    ###############

    **--basename** *<name of outpuf file prefix>*
        This is used to specify the prefix of the output files. The default is
        "ppp_fsc"
        
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
        alternative reference to that used for the vcf file.

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
 
'''
##"""
##   generates Site Frequency Spectrum (SFS) fils for fastsimcoal
##   Excoffier, L. and M. Foll. 2011. fastsimcoal: a continuous-time coalescent
##   simulator of genomic diversity under arbitrarily complex evolutionary scenarios.
##   Bioinformatics 27: 1332-1334.
##
##   Excoffier, L., Dupanloup, I., Huerta-Sánchez, E., and M. Foll (2013) Robust
##   demographic inference from genomic and SNP data. PLOS Genetics 9(10):e1003905.
##
##   Based on instructions in fastsimcoal ver 2.6 manual.
##
##   Generates 1D, 2D and multidimensional SFS files
##   For 1D files:
##       -the filename suffix is _DAFpop#.obs for an array of derived allele counts.
##          where '#' is replaced by the population number
##       -the filename suffix is _MAFpop#.obs for an array of minor allele counts.
##   For 2D files:
##       -the filename suffix is _jointDAFpop#_&.obs for an array of derived allele counts.
##           where # and & are population numbers, and # is larger than &
##       -the filename suffix is _jointDAFpop#_&.obs for an array of minor allele counts.
##   For a multidimensional file:
##       - the filename suffix is _DSFS.obs for an array of derived allele counts.
##       - the filename suffix is _MSFS.obs for an array of minor allele counts.
##       
##"""

import sys
import os
import logging
import argparse
from pgpipe.logging_module import initLogger, logArgs
from pgpipe.model import Model, read_model_file
import pgpipe.vcf_to_sfs as vcfsfs


def writemultidimfile(sfs,popmodel,basefilename,folded,sampsizes):
    dd = sfs.shape
    if folded:
        fn = basefilename + "_MSFS.obs"
    else:
        fn = basefilename + "_DSFS.obs"
    f = open(fn,'w')
    f.write('1 observations. No. of demes and sample sizes are on next line\n')
    f.write("{}\t".format(len(dd)))
    for s in sampsizes:
        f.write("{}\t".format(s))
    f.write("\n")
    temp = ["{:.0f}".format(a) for a in sfs.flat]
    f.write("{}\n".format("\t".join(temp)))
    f.close()
    return fn
                        
    

def write2dimfiles(sfs,popmodel,basename,folded):
    
    dd = sfs.shape
    fns = []
    print(popmodel.pop_list)
    pi = 0
    for popi in popmodel.pop_list[:-1]:
        pj = pi + 1
        for popj in popmodel.pop_list[pi+1:]:
            temp=vcfsfs.reducesfsdims(sfs,popmodel,[popi,popj])
            if folded:
                fn = basename + "_jointMAFpop{}_{}.obs".format(pj,pi)
            else:
                fn = basename + "_jointDAFpop{}_{}.obs".format(pj,pi)
            fns.append(fn)
            f = open(fn,'w')
            f.write('1 observations\n')
            for i in range(dd[pi]):
                f.write('\td{}_{}'.format(pi,i))
            f.write('\n')
            for j in range(dd[pj]):
                f.write('d{}_{}\t'.format(pj,j))
                for i in range(dd[pi]):
                    loc = (i,j)
                    f.write('{:.0f}\t'.format(temp[loc]))
                f.write('\n')
            f.close()
            pj += 1
        pi += 1
    return fns


def write1dimfiles(sfs,popmodel,basename,folded):

    dd = sfs.shape
    fns = []
    for pi,pop in enumerate(popmodel.pop_list):
        temp=vcfsfs.reducesfsdims(sfs,popmodel,[pop])
        if folded:
            fn = basename + "_MAFpop{}.obs".format(pi)
        else:
            fn = basename + "_DAFpop{}.obs".format(pi)
        fns.append(fn)
        f = open(fn,'w')
        for i in range(dd[0]):
            f.write('d0_{}\t'.format(i))
        f.write('\n')
        for c in temp:
            f.write('{:.0f}\t'.format(c))
        f.write('\n')
        f.close()
    return fns

        
def make_fscmsfs_file(vcffile,popmodel, basefilename, dimfiletypes, downsampsizes,
                      folded,outgroup_fasta,BEDfilename,randomsnpprop,seed):
                                                                                                              
    
    """
        calls vcfsfs.build_sfs(), to which most arguments are simply passed
        
        vcffile is a vcf, or bgzipped vcf file, with SNP data
        
        popmodel is an instance of Model
            It should specify one or more populations
        
        basefilename is the basename to be used for the sfs files

        dimfiletypes is a string that takes the value '1', '2', or 'm'

        outgroup_fasta is the name of a fasta sequence file that contains the reference
            genome to be used as the ancestor or root
            this causes the 'ref' allele, as given in the vcf to not be used as the reference
            this can be useful, for example, if an ancestral reference is available to that the reference allele is
            the ancestral allele

            if the base from the alternative references does not match either the vcf reference base or the vcf first alternate base
            the SNP will be ignored
    
        BEDfilename is the name of a BEDfile
            BEDfileis a sorted UCSC-style bedfile containing chromosome locations
            There is no header
            first column is chromosome name (must match chromosome name in vcf file)
            second column is start position (0-based, open interval)
            third column is end position (closed interval)
            other columns are ignored
            
        folded indicates that the folded sfs should be returned
            folded causes the count returned for a SNP to be that for the less common base
            ignores alt and ref

        downsampsizes is an array listing the sample sizes to be used if they are less than given in the model
            2 <= downsamplesizes[i] <= samplesizes[i]
            if None,  then the sample sizes are those given by the popmodel

        randomsnpprop is the proportion of snps to include
            uses random sampling

        seed is a random number seed that can be used with randomsnpprop
            
        
    """
    sfs = vcfsfs.build_sfs(vcffile,popmodel,BEDfilename=BEDfilename,
                    altreference = outgroup_fasta,folded = folded,
                    downsamplesizes = downsampsizes,randomsnpprop =randomsnpprop, seed = seed)
##    print(sfs.shape,sfs.sum())
    numfiles = 0
    if '1' in dimfiletypes:
        num1files = write1dimfiles(sfs,popmodel,basefilename,folded)
        numfiles += len(num1files)
    if '2' in dimfiletypes:
        num2files = write2dimfiles(sfs,popmodel,basefilename,folded)
        numfiles += len(num2files)
    if 'm' in dimfiletypes:
        ss = downsampsizes
        if ss == None:
            ss = []
            for pop in popmodel.pop_list:
                ss.append(2*len(popmodel.ind_dict[pop]))
        multdimfilename = writemultidimfile(sfs,popmodel,basefilename,folded,ss)
        numfiles += 1
        
    infostring = "generated %d. files\n"%(numfiles)
                                                                                                              
    return infostring


def fscmsfs_parser(passed_arguments):
    '''dadisnp Argument Parser - Assigns arguments from command line'''

    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError('%s not found' % value)
                setattr(args, self.dest, value)
        return customAction

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--vcf', help = "Defines the filename of the VCF", type = str, required = True, action = parser_confirm_file())
    parser.add_argument('--model-file', help = 'Defines the model filename',required = True, type = str, action = parser_confirm_file())
    parser.add_argument('--model', help = 'Defines the model and the individual(s)/population(s) to include', required = True,type = str)
    parser.add_argument("--bed-file",help="Defines the BED filename", type = str, action = parser_confirm_file())
    parser.add_argument('--outgroup-fasta',help="The name of a fasta format file containing"
                                 " the ancestral or outgroup sequence, by default the 'ref' allele "
                                 "of the vcf file is treated as the outgroup")
    parser.add_argument('--dim',type=str,nargs='+',required = True,
                        help="1 (single), 2 (two) or m (multi-) dimensional output files. "
                                "Multiple values can be used, e.g.  --dim 1 m ")
    parser.add_argument('--basename',default = 'ppp_fsc',help="base of the file name for sfs obs files")
    parser.add_argument('--folded',default=False,action="store_true",help="Generate the folded sfs")
    parser.add_argument('--randomsnpprop',type=float,help="the proportion of randomly selected SNPs to include")
    parser.add_argument('--seed',type=int,help="integer random number seed to use with randomsnpprop")
    parser.add_argument('--downsamplesizes',type=int,nargs='+',help="sample sizes to use for output "
                                "an array of integers,  one for each population, in the same order "
                                "as populations in popmodel. Values must >=2 and be <= actual number of "
                                "chromosomes in the vcf file for the corresponding population.")
       
    if passed_arguments:
        return parser.parse_args(passed_arguments)
    else:
        return parser.parse_args()


def run (passed_arguments = []):
    # Grab dadisnp arguments from command line
    args = fscmsfs_parser(passed_arguments)
    print(args)
    # Adds the arguments (i.e. parameters) to the log file
    logArgs(args, func_name = 'make_fscmsfs_file')
    popmodels = read_model_file(args.model_file)
    popmodel = popmodels[args.model]


                                                                                                              
    fscmsfs_run_infostring =make_fscmsfs_file(args.vcf,popmodel,args.basename,args.dim,
            args.downsamplesizes,args.folded,args.outgroup_fasta,
            args.bed_file,args.randomsnpprop,args.seed)
    
    logging.info(fscmsfs_run_infostring)

    print(fscmsfs_run_infostring)

if __name__ == "__main__":
    initLogger()
    run()
##    debugargs=['--vcf',"Pan_chr_21_22_test.vcf.gz",'--model-file',"panmodels.model",'--model','4Pop',
##           '--dim','1','2','m']#,'--folded']

##    debugargs=['--vcf',"Pan_all_hicov_chr22_decrun_missingasref.vcf.gz",'--model-file',
##               "panmodels.model",'--model','4Pop','--downsamplesizes','3','3','3','4',
##               '--folded','--dim','1','2','m','--outgroup-fasta',"chr22_hg18.fa"]
##    debugargs=['--vcf',"sfs_test.vcf",'--model-file',
##               "pantest.model",'--model','5Pop','--downsamplesizes','3','3','3','4','3',
##               '--folded','--dim','1','2','m']
##    debugargs=['--vcf',"sfs_test.vcf",'--model-file',
##               "pantest.model",'--model','5Pop','--dim','1','2','m']    
    
##    run(debugargs)
