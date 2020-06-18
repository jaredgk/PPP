##print("hello world")
import sys
import argparse
import pysam
import gzip

##from pgpipe.vcf_reader_func import VcfReader,checkRecordPass,checkPopWithoutMissing,getAlleleCountDict
import pgpipe.vcf_reader_func as vr
from pgpipe.model import Model, read_model_file



def make_treemix_file(vcffile,popmodel, filename=None):
    """
        Makes a SNP data input file from a vcf file for the treemix program
            Pickrell JK, Pritchard JK: Inference of Population Splits and
            Mixtures from Genome-Wide Allele Frequency Data. PLoS Genet 2012, 8(11):e1002967.

        treemix assumes biallelic states,  sites with more than two alleles, or missing data, are skipped
        
        vcffile is a vcf, or bgzipped vcf file
        
        filename is the name of the file to contain the sequences
            if None, results are sent to stdout
        
        popmodel is an instance of Model 
    """

    if filename is None:
        file_handle = sys.stdout
    else:
        file_handle = open(filename, 'w')


    a = vr.VcfReader(vcffile,popmodel=popmodel['popmod1'])

    for pop in popmodel.pop_list:
        file_handle.write("%s\t"%pop)
    file_handle.write("\n")
    numinds = []
    for pop in popmodel.pop_list:
        numinds.append(len(popmodel.ind_dict[pop]))
    ploidy = 2
    while True:
        x = a.getNext()
        if type(x) == type(None):
            break
        if vr.checkRecordPass(x,remove_indels=True,remove_multiallele=True):
            buildstr = ""
            checkcount = [0,0]
            spacer = [',',' ']
            for pop in popmodel.pop_list:
                adic = vr.getAlleleCountDict(x,popmodel.ind_dict[pop])
                
                for ai in range(2):
                    allele = x.alleles[ai]
                    try:
                        count = adic[allele]
                        checkcount[ai] += count
                        buildstr += "%d%s"%(count,spacer[ai])
                    except KeyError:
                        buildstr += "0%s"%(spacer[ai])
            if checkcount[0] > 0 and checkcount[1] > 0 and sum(checkcount) == ploidy*sum(numinds):
                file_handle.write("%s\n"%buildstr)
    file_handle.close()

                
            
        
model_file = r"/home/jodyhey/temp/PPP/tests/input/testmodel.model"

popmodelfile = read_model_file(model_file)
popmodel=popmodelfile['popmod1']

vcffile = r"/home/jodyhey/temp/PPP/tests/input/chr11.subsamples.vcf.gz"
