#!/usr/bin/env python
'''
Generates Site Frequency Spectrum (SFS) files for fastsimcoal
based on instructions in fastsimcoal ver 2.6 manual.

Generates one-dimensional (1D), two-dimensional (2D) and multidimensional SFS files

All generated SFS files are contained in a zip file archive.

Excoffier, L. and M. Foll. 2011. fastsimcoal: a continuous-time coalescent
simulator of genomic diversity under arbitrarily complex evolutionary scenarios.
Bioinformatics 27: 1332-1334.



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

**--dim** *<dimension file type signifiers>*
    One or more of '1', '2', or 'm', for 1D, 2D or multidimensional output files.
    
    For 1D files:
       * the filename suffix is _DAFpop#.obs for an array of derived allele counts.
         where '#' is replaced by the population number
       * the filename suffix is _MAFpop#.obs for an array of minor allele counts.
    For 2D files:
       * the filename suffix is _jointDAFpop#_&.obs for an array of derived allele counts.
         where # and & are population numbers, and # is larger than &
       * the filename suffix is _jointDAFpop#_&.obs for an array of minor allele counts.
    For a multidimensional file:
       * the filename suffix is _DSFS.obs for an array of derived allele counts.
       * the filename suffix is _MSFS.obs for an array of minor allele counts.
       
#################
Optional Aguments
#################

**--basename** *<name of outpuf file prefix>*
    This is used to specify the prefix of the output files and the prefix of the
    zip file archive. The default is
    "ppp_fsc" in the same folder as the vcf file
    
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


#############
Example usage
#############
Example command-lines:

.. code-block:: bash

    vcf_to_fastsimcoal.py -h

   
.. code-block:: bash

    vcf_to_fastsimcoal.py --vcf pan_example2.vcf.gz --model-file panmodels.model --modelname 5Pop --downsamplesizes 3 3 3 4 2  --basename vcf_fsc2 --folded --dim 1 2 m  --outgroup-fasta chr22_pan_example2_ref.fa

'''


import sys
import os
import logging
import argparse
import zipfile
##try:
##    import zlib
##    compression = zipfile.ZIP_DEFLATED
##except:
##    compression = zipfile.ZIP_STORED
from pgpipe.logging_module import initLogger, logArgs
from pgpipe.model import Model, read_model_file
import pgpipe.vcf_to_sfs as vcfsfs
from pgpipe.misc import argprase_kwargs


def writemultidimfile(sfs,popmodel,basename,folded,sampsizes):
    dd = sfs.shape
    if folded:
        fn = basename + "_MSFS.obs"
    else:
        fn = basename + "_DSFS.obs"
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
##    print(popmodel.pop_list)
    pi = 0
    for popi in popmodel.pop_list[:-1]:
        pj = pi + 1
        for popj in popmodel.pop_list[pi+1:]:
            temp=vcfsfs.reduce_sfs_dims(sfs,popmodel,[popi,popj])
            if folded:
                fn = basename + "_jointMAFpop{}_{}.obs".format(pj,pi)
            else:
                fn = basename + "_jointDAFpop{}_{}.obs".format(pj,pi)
##            fns.append(os.path.abspath(fn))
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
        temp=vcfsfs.reduce_sfs_dims(sfs,popmodel,[pop])
        if folded:
            fn = basename + "_MAFpop{}.obs".format(pi)
        else:
            fn = basename + "_DAFpop{}.obs".format(pi)
##        fns.append(os.path.abspath(fn))
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

        
def make_fscmsfs_file(vcffile,model_file,model, basename, dimfiletypes, downsampsizes,
                      folded,outgroup_fasta,BEDfilename,randomsnpprop,seed):
                                                                                                              
    
    """
        calls vcfsfs.build_sfs(), to which most arguments are simply passed
        
        vcffile is a vcf, or bgzipped vcf file, with SNP data

        model_file is a model file
        model is the name of a population in the model
        
        basename is the basename of the zip archive and the basename of the sfs files

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
    
    sfs = vcfsfs.build_sfs(vcffile,model_file,model,BEDfilename=BEDfilename,
                    altreference = outgroup_fasta,folded = folded,
                    downsamplesizes = downsampsizes,randomsnpprop =randomsnpprop, seed = seed)
    popmodels = read_model_file(model_file)
    popmodel = popmodels[model]
##    print(sfs.shape,sfs.sum())
    numfiles = 0
    filenames = []
    if '1' in dimfiletypes:
        filenames += write1dimfiles(sfs,popmodel,basename,folded)
    if '2' in dimfiletypes:
        filenames += write2dimfiles(sfs,popmodel,basename,folded)
    if 'm' in dimfiletypes:
        ss = downsampsizes
        if ss == None:
            ss = []
            for pop in popmodel.pop_list:
                ss.append(2*len(popmodel.ind_dict[pop]))
        multdimfilename = writemultidimfile(sfs,popmodel,basename,folded,ss)
        filenames.append(multdimfilename)
    numfiles = len(filenames)
    if len(os.path.basename(basename)) > 4 and basename[-4:].lower() == ".zip":
        basename = basename[:-4]
    zarch = zipfile.ZipFile(basename + ".zip", mode='w')
    for fn in filenames:
        zarch.write(fn,arcname = os.path.basename(fn),compress_type=zipfile.ZIP_STORED)
    zarch.close()
    for fn in filenames:
        os.remove(fn)
    infostring = "generated %d. files\n"%(numfiles)
                                                                                                              
    return infostring


def fscmsfs_parser(passed_arguments=[]):
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
    parser.add_argument('--modelname', help = 'Defines the model and the individual(s)/population(s) to include', required = True,type = str)
    parser.add_argument("--bed-file",help="Defines the BED filename", type = str, action = parser_confirm_file())
    parser.add_argument('--outgroup-fasta',help="The name of a fasta format file containing"
                                 " the ancestral or outgroup sequence, by default the 'ref' allele "
                                 "of the vcf file is treated as the outgroup", action = parser_confirm_file())
    parser.add_argument('--dim',type=str,nargs='+',required = True,
                        help="1 (single), 2 (two) or m (multi-) dimensional output files. "
                                "Multiple values can be used, e.g.  --dim 1 m ")
    parser.add_argument('--basename',help="base of the file name for sfs obs files")
    parser.add_argument('--folded',default=False,action="store_true",help="Generate the folded sfs")
    parser.add_argument('--randomsnpprop',type=float,help="the proportion of randomly selected SNPs to include")
    parser.add_argument('--seed',type=int,help="integer random number seed to use with randomsnpprop")
    parser.add_argument('--downsamplesizes',type=int,nargs='+',help="sample sizes to use for output "
                                "an array of integers,  one for each population, in the same order "
                                "as populations in popmodel. Values must >=2 and be <= actual number of "
                                "chromosomes in the vcf file for the corresponding population.")
       
    if passed_arguments:
        return vars(parser.parse_args(passed_arguments))
    else:
        return vars(parser.parse_args())


def run (**kwargs):
    
    # Update kwargs with defaults
    if __name__ != "__main__":
        kwargs = argprase_kwargs(kwargs, fscmsfs_parser)
    # Assign arguments
    args = argparse.Namespace(**kwargs)
    
##    print(args)
    # Adds the arguments (i.e. parameters) to the log file
    logArgs(args, func_name = 'make_fscmsfs_file')

    if args.basename == None:
        args.basename = os.path.dirname(args.vcf) + "//ppp_fsc"
                                                                                                              
    fscmsfs_run_infostring =make_fscmsfs_file(args.vcf,args.model_file,args.modelname,args.basename,args.dim,
            args.downsamplesizes,args.folded,args.outgroup_fasta,
            args.bed_file,args.randomsnpprop,args.seed)
    
    logging.info(fscmsfs_run_infostring)

##    print(fscmsfs_run_infostring)

if __name__ == "__main__":
    initLogger()
    run(**fscmsfs_parser)
    exit()
    debugargs=['--vcf',"..//jhtests//pan_example.vcf.gz",
               '--model-file',"..//jhtests//panmodels.model",'--modelname','3Pop',
           '--dim','1','2','m','--basename','../jhtests/results/vcf_fsc1']#,'--folded']
    run(debugargs)
    debugargs=['--vcf',"..//jhtests//pan_example2.vcf.gz",
               '--model-file',"..//jhtests//panmodels.model",'--modelname','5Pop',
               '--downsamplesizes','3','3','3','4','2','--basename','../jhtests/results/vcf_fsc2',
               '--folded','--dim','1','2','m','--outgroup-fasta',"..//jhtests//chr22_pan_example2_ref.fa"]
    run(debugargs)    
##    debugargs=['--vcf',"sfs_test.vcf",'--model-file',
##               "pantest.model",'--model','5Pop','--downsamplesizes','3','3','3','4','3',
##               '--folded','--dim','1','2','m']
##    run(debugargs)    
##    debugargs=['--vcf',"sfs_test.vcf",'--model-file',
##               "pantest.model",'--model','5Pop','--dim','1','2','m']    
##    run(debugargs)    

