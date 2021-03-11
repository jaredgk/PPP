import sys
import pysam
from math import ceil,floor
from vcf_reader_func import VcfReader
from gene_region import Region, RegionList

def findPrevLastSnps(vcfr,region,window=10000):
    left_out = 0
    i = region.start
    r = Region(i-window,i,region.chrom)
    tl = vcfr.getRecordList(r)
    if len(tl) != 0:
        left_out = tl[-1].pos
    else:
        left_out = i - window
    right_out = 0
    i = region.end
    r = Region(i,i+window,region.chrom)
    tl = vcfr.getRecordList(r)
    if len(tl) != 0:
        right_out = tl[0].pos
    else:
        right_out = i + window
    return left_out,right_out

#Given a VCF file and a region file, extend regions to
#length halfway between border SNPs inside/outside of the range

vcf_name = str(sys.argv[1])
reg_name = str(sys.argv[2])

vcfr = VcfReader(vcf_name)
regr = RegionList(filename=reg_name)

for region in regr.regions:
    rl = vcfr.getRecordList(region)
    if len(rl) == 0:
        loggin.warning("Region at %s has no variants" % region.toStr())
        print (region.toStr(sep='\t'))
        continue
    left_in = rl[0].pos
    right_in = rl[-1].pos
    left_out,right_out = findPrevLastSnps(vcfr,region)
    left_half = int(ceil((left_in+left_out)/2))
    right_half = int(floor((right_in+right_out)/2))
    print (region.chrom,left_half,right_half,sep='\t')
