import sys
import pysam
from vcf_reader_func import VcfReader






vcf_name = str(sys.argv[1])
subsamp_name = str(sys.argv[2])
outname = str(sys.argv[3])


vcfr = VcfReader(vcf_name, subsamp_fn=subsamp_name)

outvcf = pysam.VariantFile(outname,'w',header=vcfr.reader.header)
#outvcf.header = vcfr.reader.header

for rec in vcfr.reader:
    outvcf.write(rec)
