from __future__ import print_function
import sys
import pysam
from math import ceil,floor
from pgpipe.vcf_reader_func import VcfReader
from pgpipe.genome_region import Region, RegionList
import logging

def findPrevLastSnps(vcfr,region,window=1000):
    left_out = 0
    i = region.start-1
    r = Region(i-window,i,region.chrom)
    tl = vcfr.getRecordList(r)
    if len(tl) != 0:
        left_out = tl[-1].pos
    else:
        left_out = i - window
    right_out = 0
    i = region.end+1
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
lil_name = str(sys.argv[2])

vcfr = VcfReader(vcf_name)
lilr = VcfReader(lil_name)

lilrecs = lilr.getRecordList()


region = Region(lilrecs[0].pos-1,lilrecs[-1].pos+1,lilrecs[0].chrom)

rl = vcfr.getRecordList(region)
if len(rl) == 0:
    logging.warning("Region at %s has no variants" % region.toStr())
    print (region.toStr(sep='\t'))
    exit()
left_in = rl[0].pos
right_in = rl[-1].pos
left_out,right_out = findPrevLastSnps(vcfr,region)
left_half = int(ceil((left_in+left_out)/2))
right_half = int(floor((right_in+right_out)/2))
print (region.chrom,left_half,right_half,left_out,left_in,right_in,right_out,sep='\t')
