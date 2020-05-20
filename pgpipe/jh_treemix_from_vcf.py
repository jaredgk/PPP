##print("hello world")
import sys
import argparse
import pysam
import gzip

##from pgpipe.vcf_reader_func import VcfReader,checkRecordPass,checkPopWithoutMissing,getAlleleCountDict
import pgpipe.vcf_reader_func as vr
from pgpipe.model import Model, read_model_file

model_file = r"/home/jodyhey/temp/PPP/tests/input/testmodel.model"

popmodel = read_model_file(model_file)

vcffile = r"/home/jodyhey/temp/PPP/tests/input/chr11.subsamples.vcf.gz"

a = vr.VcfReader(vcffile,popmodel=popmodel['popmod1'])
##x = a.getNext()
##
##
##for pop in popmodel['popmod1'].pop_list:
##    for ind in popmodel['popmod1'].ind_dict[pop]:
##        print(pop,ind,x.samples[ind].alleles)

tfn = "treemixfile.out"
tf = open(tfn,'w')
for pop in popmodel['popmod1'].pop_list:
    tf.write("%s\t"%pop)
tf.write("\n")
numinds = []
for pop in popmodel['popmod1'].pop_list:
    numinds.append(len(popmodel['popmod1'].ind_dict[pop]))
ploidy = 2
while True:
    x = a.getNext()
    if type(x) == type(None):
        break
    if vr.checkRecordPass(x,remove_indels=True,remove_multiallele=True):
        buildstr = ""
        checkcount = [0,0]
        spacer = [',',' ']
        for pop in popmodel['popmod1'].pop_list:
            adic = vr.getAlleleCountDict(x,popmodel['popmod1'].ind_dict[pop])
            
            for ai in range(2):
                allele = x.alleles[ai]
                try:
                    count = adic[allele]
                    checkcount[ai] += count
                    buildstr += "%d%s"%(count,spacer[ai])
                except KeyError:
                    buildstr += "0%s"%(spacer[ai])
        if checkcount[0] > 0 and checkcount[1] > 0 and sum(checkcount) == ploidy*sum(numinds):
            tf.write("%s\n"%buildstr)
tf.close()

                
            
        
