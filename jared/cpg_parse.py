import sys
import pysam
from vcf_reader_func import checkIfCpG, VcfReader
from gene_region import RegionList



vcf_name = str(sys.argv[1])
ref_name = str(sys.argv[2])
reg_name = str(sys.argv[3])


vcfr = VcfReader(vcf_name)
refr = pysam.FastaFile(ref_name)
regr = RegionList(filename=reg_name)

for region in regr.regions:
    rl = vcfr.getRecordList(region)
    print (region.toStr(sep="\t")+'\t'+str(len(rl)))

#for rec in vcfr.reader:
#    print(rec.pos,rec.ref,checkIfCpG(rec,refr))
