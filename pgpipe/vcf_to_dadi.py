#!/usr/bin/env python
'''
Generates a dadi snp file from a vcf file.

The dadi snp file format is described in the dadi manual 

https://dadi.readthedocs.io/en/latest/user-guide/importing-data/#snp-data-format

Gutenkunst RN, Hernandez RD, Williams SH, Bustamante CD (2009) Inferring the
joint demographic history of multiple populations from multidimensional SNP
frequency data. PLoS Genet 5: e1000695. DOI: 10.1371/journal.pgen.1000695

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
    
**--out** *<output file name>*
    Specifies the complete output filename.
       
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
    alternative reference to that used for the vcf file.

    This fasta file must have been properly aligned to the reference
    used in the vcf file.  
        
    This option can be useful, for example, if an ancestral or outgroup reference is
    available that more accurately identifies the ancestral (and thus derived)
    allele at each SNP than does the reference used to make the vcf file.

**--comment** *<comment text>*
    Comment text to be included in the header of the output file.



#############
Example usage
#############
Example command-lines:

.. code-block:: bash

    vcf_to_dadi.py -h

   
.. code-block:: bash

    vcf_to_dadi.py --vcf  pan_example.vcf.gz --model-file  panmodels.model  --modelname 4Pop  --out vcf_dadisnp_bedfile_test.out  --comment testing bedfile --bed-file pan_example_regions.bed  

.. code-block:: bash

    vcf_to_dadi.py --vcf  pan_example2.vcf.gz --model-file  panmodels.model  --modelname 4Pop --out  vcf_dadisnp_test.out  --comment testing comment  

.. code-block:: bash

    vcf_to_dadi.py --vcf pan_example2.vcf.gz --model-file panmodels.model  --modelname 4Pop --out vcf_dadisnp_fasta_test.out  --comment testing outgroup-fasta --outgroup-fasta  chr22_pan_example2_ref.fa 
 

'''

import sys
import os
import subprocess
import logging
import argparse
import pysam
from pgpipe.logging_module import initLogger, logArgs
import pgpipe.vcf_bed_to_seq as vBs
from pgpipe.model import Model, read_model_file
from pgpipe.genome_region import Region
import pgpipe.vcf_reader_func as vr
from pgpipe.misc import argprase_kwargs

def getallelecount(r,popmodel,altref_access=None):

    ref = r.ref.upper()
    alt = r.alts[0].upper()
    c = []
    pi = 0
    for pop in popmodel.pop_list:
        adic = vr.getAlleleCountDict(r,popmodel.ind_dict[pop])
        temp = [adic[r.alleles[0]],adic[r.alleles[1]]]
        pi += 1
        c.append(temp)
    if altref_access != None:
        altref = altref_access.fetch(r.chrom,r.pos-1,r.pos).upper()
##        print(r.chrom,r.pos,r.ref,r.alts,altref)
        if not (altref == ref.upper() or altref.upper() in r.alts or altref.lower() in r.alts):
            print(r.ref,altref,r.alts)
            return '','',-1  # altref not found, skip this snp        
        if altref == r.alts[0]: #switch ref and alt
            print(altref)
            altc = [[temp[1],temp[0]] for temp in c]
            c = list(altc)
            ref = altref
            alt = r.ref.upper()
    counts = [[temp[0] for temp in c],[temp[1] for temp in c]]
    return ref,alt,counts


def write_dadisnp_row(r,f,ref,alt,counts,label):
    """
        r is a vcf record
        f is a file handle
        
    """
    row = ['-{}-'.format(ref),'-{}-'.format(alt)]
    row.append(ref)
    row += list(map(str,counts[0]))
    row.append(alt)
    row += list(map(str,counts[1]))
    row.append(label)
    row.append(str(r.pos))
    f.write("{}\n".format('\t'.join(row)))
        
def make_dadisnp_file(vcffile,popmodel, outfilename, outgroup_fasta = None,
                      BEDfilename=None, comment = None):
    """
        
        vcffile is a vcf, or bgzipped vcf file, with SNP data
        
        popmodel is an instance of Model
            It should specify 2 or more populations
        
        outfilename is the name of the dadi snp file to be made.
            

    
        BEDfile 
            BEDfileis a sorted UCSC-style bedfile containing chromosome locations
            There is no header
            first column is chromosome name (must match chromosome name in vcf file)
            second column is start position (0-based, open interval)
            third column is end position (closed interval)
            other columns are ignored
            
        outgroup_fasta is the name of a fasta sequence file that contains the reference
            genome to be used as the ancestor or root
            this causes the 'ref' allele, as given in the vcf to not be used as the reference
            this can be useful, for example, if an ancestral reference is available to that the reference allele is
            the ancestral allele

            if the base from the alternative references does not match either the vcf reference base or the vcf first alternate base
            the SNP will be ignored

        comment
            text that gets added to the output file 

        
    """
    #open out file
    file_handle = open(outfilename, 'w')

    #make header row
    file_handle.write("# dadi snp file generated using PPP (Weber et al) from {}\n".format(os.path.basename(vcffile)))
    if outgroup_fasta != None:
        file_handle.write("# sequence file {} used to provide ancestral (outgroup) allele\n".format(os.path.basename(outgroup_fasta)))
    else:
        file_handle.write("# reference allele in vcf file treated as ancestral (outgroup) allele\n")
    if BEDfilename != None:
        file_handle.write("# SNPs drawn from regions given in BED file {}\n".format(os.path.basename(BEDfilename)))
    if comment != None:
        file_handle.write("# comment added at runtime: {}\n".format(comment))
        
    header = ['ingroup','outgroup','Allele1']
    for pop in popmodel.pop_list:
        header.append(pop)
    header.append('Allele2')
    for pop in popmodel.pop_list:
        header.append(pop)
    if BEDfilename != None:
        header.append('region')
    else:
        header.append('Chr')
    header.append('Position')
    
    file_handle.write("{}\n".format('\t'.join(header)))

    #make altref_access
    if outgroup_fasta != None:
        if not os.path.isfile(outgroup_fasta + ".fai"):
        # make an index 
            pysam.faidx(outgroup_fasta)
        altref_access = pysam.FastaFile(outgroup_fasta)
    else:
        altref_access = None
                      
    #initiatlize reader 
    vcf_reader = vr.VcfReader(vcffile,popmodel=popmodel)
    
    numSNPs = 0
    if BEDfilename == None:
        while True:
            r = vcf_reader.getNext()
            if type(r) == type(None):
                break
            if vr.checkRecordPass(r,remove_indels=True,remove_multiallele=True):
                ref,alt,counts = getallelecount(r,popmodel,altref_access=altref_access)
                if counts != -1:
                    numSNPs += 1
                    if sum(counts[0]) > 0 and sum(counts[1]) > 0:
                        write_dadisnp_row(r,file_handle,ref,alt,counts,r.chrom)
                    
        
    else:
        with open(BEDfilename,'r') as bf:
            for line in bf:
                ls = line.split()
                if len(ls) >= 3:  # anything less than 3 is not a valid line
                    region = Region(int(ls[1]),int(ls[2]),ls[0])
                    regionstr = '%s:%d-%d'%(region.chrom,region.start+1,region.end)
                else:
                    break # get out of loop
                recordlist = vcf_reader.getRecordList(region=region)
                ri = 0
                for r in recordlist:
                    ref,alt,counts = getallelecount(r,popmodel,altref_access=altref_access)
                    if counts != -1:
                        numSNPs += 1
                        if sum(counts[0]) > 0 and sum(counts[1]) > 0:
                            write_dadisnp_row(r,file_handle,ref,alt,counts,regionstr)
                

        
    file_handle.close()
    infostring = "generating %s.  %d SNPs \n"%(outfilename, numSNPs)
    return infostring


def dadisnp_parser(passed_arguments=[]):
    '''dadisnp Argument Parser - Assigns arguments from command line'''

    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError('%s not found' % value)
                setattr(args, self.dest, value)
        return customAction

    dadisnp_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    dadisnp_parser.add_argument('--vcf', help = "Defines the filename of the VCF", type = str, required = True, action = parser_confirm_file())
    dadisnp_parser.add_argument('--model-file', help = 'Defines the model filename',required = True, type = str, action = parser_confirm_file())
    dadisnp_parser.add_argument('--modelname', help = 'Defines the model and the individual(s)/population(s) to include', required = True,type = str)
    dadisnp_parser.add_argument("--bed-file",help="Defines the BED filename", type = str, action = parser_confirm_file())
    dadisnp_parser.add_argument('--out', help = 'Defines the complete output filename', type = str)
    dadisnp_parser.add_argument('--outgroup-fasta',help="The name of a fasta format file containing"
                                 " the ancestral or outgroup sequence, by default the 'ref' allele "
                                 "of the vcf file is treated as the outgroup")
    dadisnp_parser.add_argument('--comment',help="comment (in quotes) to be added to the snp file ")

    if passed_arguments:
        return vars(dadisnp_parser.parse_args(passed_arguments))
    else:
        return vars(dadisnp_parser.parse_args())


def run (**kwargs):
    # Update kwargs with defaults
    if __name__ != "__main__":
        kwargs = argprase_kwargs(kwargs, dadisnp_parser)
    # Assign arguments
    dadisnp_args = argparse.Namespace(**kwargs)    

    # Adds the arguments (i.e. parameters) to the log file
    logArgs(dadisnp_args, func_name = 'make_dadisnp_file')

##    print(dadisnp_args)
    popmodels = read_model_file(dadisnp_args.model_file)
    popmodel = popmodels[dadisnp_args.modelname]
    dadisnp_run_infostring =make_dadisnp_file(dadisnp_args.vcf,popmodel,dadisnp_args.out,
            outgroup_fasta=dadisnp_args.outgroup_fasta,BEDfilename=dadisnp_args.bed_file,
            comment=dadisnp_args.comment)
    logging.info(dadisnp_run_infostring)

##    print(dadisnp_run_infostring)

if __name__ == "__main__":
    initLogger()
    run(**dadisnp_parser())
    exit()
                      
    debugargs = ['--vcf','..//jhtests//pan_example.vcf.gz','--model-file',
                 "..//jhtests//panmodels.model",
                 '--modelname',"4Pop",'--out','../jhtests/results/vcf_dadisnp_bedfile_test.out',
                 '--comment','testing bedfile','--bed-file',"..//jhtests//pan_example_regions.bed"]
    run(debugargs)
    debugargs = ['--vcf',"..//jhtests//pan_example2.vcf.gz",'--model-file',"..//jhtests//panmodels.model",
                 '--modelname',"4Pop",'--out','..//jhtests//results//vcf_dadisnp_test.out',
                 '--comment','testing comment']
    #'--outgroup-fasta',"..//jhtests//chr22_pan_example2_ref.fa",'--bed-file',"..//jhtests//pan_example_regions.bed",]
    run(debugargs)    
    debugargs = ['--vcf',"..//jhtests//pan_example2.vcf.gz",'--model-file',
                 "..//jhtests//panmodels.model",
                 '--modelname',"4Pop",'--out','..//jhtests//results//vcf_dadisnp_fasta_test.out',
                 '--comment','testing outgroup-fasta','--outgroup-fasta',
                 "..//jhtests//chr22_pan_example2_ref.fa"]

    run(debugargs)



