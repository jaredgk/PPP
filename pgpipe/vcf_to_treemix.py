#!/usr/bin/env python
'''
The treemix program was developed by Pickrell and Prichard (2012)
to estimate phylogeny and admixture for closely related populations.

Pickrell JK, Pritchard JK (2012) Inference of Population Splits and Mixtures
from Genome-Wide Allele Frequency Data. PLOS Genetics 8(11): e1002967.

The program can make use of very large numbers of SNPs.

vcf_to_treemix.py will generate a treemix input file from a vcf file.


If run using the --bed-file and --kblock options, the resulting treemix
file can be run using the 'linkage disequilibrium' (-k) option.
Under this option each block of kblock SNPs are treated as a linked group
and different groups are treated as unlinked.

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

**--out** *<outpuf_file_name>*
    The name of the treemix file to be generated. The file is bgzipped
    and '.gz' is added to the end of the name
   
#################
Optional Aguments
#################

**--bed-file** *<BED_file_name>*
    The BED file is a sorted UCSC-style bedfile containing chromosome locations of
    the SNPs to be included in the output file. The BED file has no header.
    The first column is the chromosome name (this must match the chromosome
    name in the vcf file).
    The second column is start position (0-based, open interval)
    The third column is end position (closed interval).
    Any other columns are ignored.

    If used with --kblock,  each of the BED file regions is used to generate
    one block of SNPs
    
**--kblock** *<k_block_size>*
    Used with --bed-file, for using treemix runtime option -k.  
    k is the number of SNPs in a block in the treemix file. 
    If the actual number of SNPs in a BED file interval is less than kblock,
    then additional invariant rows are added to the treemix file so
    the total numbers of rows for that and every block is equal to kblock.
    k is set to 1000 by default.  It needs to be increased only when
    one or more BED file regions have more than k snps.
  
#############
Example usage
#############
Example command-lines:

.. code-block:: bash

    vcf_to_treemix.py -h

   
.. code-block:: bash

    vcf_to_treemix.py --vcf pan_example.vcf.gz --model-file panmodels.model --modelname 4Pop --out vcf_treemixtest1 --bed-file pan_example_regions.bed --kblock 1000
 
 .. code-block:: bash

    vcf_to_treemix.py --vcf pan_example.vcf.gz --model-file panmodels.model --modelname 4Pop --out vcf_treemixtest2    


'''


import os
import logging
import argparse
import subprocess
import pgpipe.vcf_bed_to_seq as vBs
from pgpipe.model import Model, read_model_file
from pgpipe.genome_region import Region
import pgpipe.vcf_reader_func as vr
from pgpipe.logging_module import initLogger, logArgs
from pgpipe.misc import argprase_kwargs

def make_treemix_file(vcffile,popmodel, outfilename, BEDfilename=None, kblock = 1000):
    """
        Treemix requires a gzipped input file with SNP allele counts
            Pickrell JK, Pritchard JK: Inference of Population Splits and
            Mixtures from Genome-Wide Allele Frequency Data. PLoS Genet 2012, 8(11):e1002967.

        treemix assumes biallelic states,  this code skips over SNPs with more than 2 alleles

        treemix is fairly foregiving of sample counts:
            counts at different loci do not need to sum to the same values
                this means that program runs regardless of ploidy variation or missing data
            Invariant sites are allowed
            total counts of zero in some populations at some loci
                generate a warning (not clear if the SNP is included) 
        
        vcffile is a vcf, or bgzipped vcf file, with SNP data
        
        popmodel is an instance of Model
            It should specify 3 or more populations
        
        outfilename is the name of the treemix file to be made.
            The file is bgzipped and '.gz' is added to the end of the name 

    
        BEDfile and kblock are used if treemix is to be run
            using the 'linkage disequilibrium' option.
            Under this option each block of kblock SNPs are treated as a linked group
            and different groups are treated as unlinked.

            Each region in the BEDfile designates an interval within which
            the SNPs are to be treated as a block

            kblock (default 1000) is just the length of a set of SNPs.
            If the actual number of SNPs is less than kblock, then additional
            invariant rows are added to the data file so the total numbers of
            rows for that and every block is equal to kblock

            If no variation is found in a region, no rows are added to the output file 

            If any BEDfile interval has a SNP count > kblock, the program is
            halted with an error explaining.  It should be rerun with a higher kblock value
            

            BEDfileis a sorted UCSC-style bedfile containing chromosome locations
            There is no header
            first column is chromosome name (must match chromosome name in vcf file)
            second column is start position (0-based, open interval)
            third column is end position (closed interval)
            other columns are ignored

        returns:
            numblocks - the number of blocks of SNP (1 if BEDfilename==None)
            numSNPs - the number of variable SNPs written to the treemix file
    """



    def bgzipFile(filename,force=False,renameout=None):
        """
            modified from from tabix_wrapper.py
        """
        BGZIP_PATH='bgzip'
        args = [BGZIP_PATH]
        if force:
            args.append('-f')
        args.append(filename)
        if renameout is not None:
            args.append(['-c','>',renameout])
        else: # delete gz file if it exists 
            tempgzfilename = filename + ".gz"
            if os.path.exists(tempgzfilename):
                os.remove(tempgzfilename)            

        bgzip_call = subprocess.Popen(args, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
        bgzip_out, bgzip_err = bgzip_call.communicate()

    def write_treemix_row(r,f,popmodel):
        """
            r is a vcf record
            f is a file handle
            
        """
        buildstr = ""
        checkcount = [0,0]
        spacer = [',',' ']
        for pop in popmodel.pop_list:
            adic = vr.getAlleleCountDict(r,popmodel.ind_dict[pop])
            for ai in range(2):
                allele = r.alleles[ai]
                try:
                    count = adic[allele]
                    checkcount[ai] += count
                    buildstr += "%d%s"%(count,spacer[ai])
                except KeyError:
                    buildstr += "0%s"%(spacer[ai])
        if checkcount[0] > 0 and checkcount[1] > 0:
            f.write("%s\n"%buildstr)
            return True
        else:
            return False
    

    file_handle = open(outfilename, 'w')
    for pop in popmodel.pop_list:
        file_handle.write("%s\t"%pop)
    file_handle.write("\n")

    vcf_reader = vr.VcfReader(vcffile,popmodel=popmodel)
    
    numSNPs = 0
    if BEDfilename == None:
        numblocks = 1
        while True:
            r = vcf_reader.getNext()
            if type(r) == type(None):
                break
            if vr.checkRecordPass(r,remove_indels=True,remove_multiallele=True):
                numSNPs += write_treemix_row(r,file_handle,popmodel)
        infostring = "generating %s.  %d SNPs \n"%(outfilename + ".gz", numSNPs)
    else:
        novarstr = ""
        numblocks = 0
        numnovar = 0
        for pop in popmodel.pop_list:
            numinds = len(popmodel.ind_dict[pop])
            novarstr += "%d,0 "%(2*numinds)
        novarstr = novarstr[:-1] #drop last ' '
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
                    if vr.checkRecordPass(r,remove_indels=True,remove_multiallele=True):
                        numSNPs += write_treemix_row(r,file_handle,popmodel)
                        ri += 1
##                print(len(recordlist),ri)
                if ri > 0:
                    numblocks += 1
                    while ri < kblock: # add rows with no variation
                        file_handle.write("%s\n"%novarstr)
                        ri += 1
                else:
                    numnovar += 1
        infostring ="generating %s.  %d SNPs in %d blocks.  %d BED intervals had 0 SNPs"%(outfilename + ".gz", numSNPs,numblocks,numnovar)
    file_handle.close()
    bgzipFile(outfilename)
    
    return infostring


def treemix_parser(passed_arguments=[]):
    '''treemix Argument Parser - Assigns arguments from command line'''

    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError('%s not found' % value)
                setattr(args, self.dest, value)
        return customAction

    treemix_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    treemix_parser.add_argument('--vcf', help = "Defines the filename of the VCF", type = str, required = True, action = parser_confirm_file())
    treemix_parser.add_argument('--model-file', help = 'Defines the model filename',required = True, type = str, action = parser_confirm_file())
    treemix_parser.add_argument('--modelname', help = 'Defines the model and the individual(s)/population(s) to include', required = True,type = str)
    treemix_parser.add_argument("--bed-file",help="Defines the BED filename", type = str, action = parser_confirm_file())
    treemix_parser.add_argument('--out', help = 'Defines the complete output filename', type = str)
    treemix_parser.add_argument("--kblock",type=int,default = 1000,help="For treemix option of having SNPs in 'blocks'."
                                "The number of SNPs in a block when --bed-file is used (default 1000)"
                                "If the actual number of SNPs is less than kblock, then additional"
                                "invariant rows are added to the data file so the total numbers of"
                                "rows for that and every block is equal to kblock")


    if passed_arguments:
        return vars(treemix_parser.parse_args(passed_arguments))
    else:
        return vars(treemix_parser.parse_args())


def run (**kwargs):
    # Grab treemix arguments from command line

    # Update kwargs with defaults
    if __name__ != "__main__":
        kwargs = argprase_kwargs(kwargs, treemix_parser)

    # Assign arguments
    treemix_args = argparse.Namespace(**kwargs)    
    

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(treemix_args, func_name = 'make_treemix_file')

    if (treemix_args.bed_file != None) and (treemix_args.kblock == None):
        raise Exception("--bed-file and --kblock must be used together")

##    print(treemix_args)
    popmodels = read_model_file(treemix_args.model_file)
    popmodel = popmodels[treemix_args.modelname]
    treemix_run_infostring =make_treemix_file(treemix_args.vcf,popmodel,treemix_args.out,BEDfilename=treemix_args.bed_file, kblock = treemix_args.kblock)
    logging.info(treemix_run_infostring)



if __name__ == "__main__":
    initLogger()
    run(**treemix_parser())
    exit()
    debugargs = ['--vcf','../jhtests/pan_example.vcf.gz','--model-file',
            "../jhtests/panmodels.model",'--modelname',"4Pop",
            '--out','../jhtests/results/vcf_treemixtest1','--bed-file',"../jhtests/pan_example_regions.bed",'--kblock','1000']
    run(**debugargs)
    exit()
    debugargs = ['--vcf','../jhtests/pan_example.vcf.gz','--model-file',
            "../jhtests/panmodels.model",'--modelname',"4Pop",
            '--out','../jhtests/results/vcf_treemixtest2']
    run(debugargs)
    exit()



