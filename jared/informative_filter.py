import sys
from vcf_reader_func import getPassSites, VcfReader
from gene_region import Region, RegionList
import pysam
import argparse

def createParser():
    parser = argparse.ArgumentParser(description="Given a ")



vcf_name = str(sys.argv[1])
reg_name = str(sys.argv[2])
ref_name = None
fasta_ref = None
if len(sys.argv) > 3:
    ref_name = str(sys.argv[3])
    fasta_ref = pysam.FastaFile(ref_name)

vcf_reader = VcfReader(vcf_name)

regions = RegionList(filename=reg_name)

for region in regions.regions:
    rec_list = vcf_reader.getRecordList(region)
    pass_list = getPassSites(rec_list, remove_cpg=True, fasta_ref=fasta_ref)
    if pass_list.count(True) >= 3:
        print (region.toStr(sep='\t'))
