import sys
from pgpipe.vcf_reader_func import VcfReader, checkPopWithoutMissing, readSinglePopModel
import pysam



vcfname = str(sys.argv[1])
popmodel = readSinglePopModel(sys.argv[2])
vcfr = VcfReader(vcfname,popmodel = popmodel)
rec_list = vcfr.getRecordList()

print (checkPopWithoutMissing(rec_list,popmodel,vcfr.popkeys))

