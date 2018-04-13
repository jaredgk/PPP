import sys
from vcf_reader_func import getPassSites, VcfReader
from gene_region import Region, RegionList
import pysam
import argparse

def createParser():
    parser = argparse.ArgumentParser(description="Given a VCF file and a list of intervals (BED style), will output regions that have over a certain number of qualifying variants in the region")
    parser.add_argument("--vcf", dest="vcfname")
    parser.add_argument("--bed", dest="bedname")
    parser.add_argument("--zero-ho", dest="zeroho", action="store_true")
    parser.add_argument("--zero-closed", dest="zeroclosed", action="store_true")
    parser.add_argument("--parsecpg", dest="refname")
    parser.add_argument("--remove-indels", dest="remove_indels", action="store_true", help=("Removes indels from output VCF files"))
    parser.add_argument("--remove-multi", dest="remove_multiallele", action="store_true")
    parser.add_argument("--remove-missing", dest="remove_missing", default=-1, help=("Will filter out site if more than the given number of individuals (not genotypes) are missing data. 0 removes sites with any missing data, -1 (default) removes nothing"))
    parser.add_argument("--informative-count", dest="informative_count", default=0)
    parser.add_argument("--minsites", dest="minsites", default=3, help=("Regions with at least this many variants passing filters will be output"))
    return parser


#vcf_name = str(sys.argv[1])
#reg_name = str(sys.argv[2])
#ref_name = None
#fasta_ref = None
#if len(sys.argv) > 3:
    #ref_name = str(sys.argv[3])
    #fasta_ref = pysam.FastaFile(ref_name)
def filter_bed_regions(sys_args):
    #parser = argparse.parse_args(sys_args)
    parser = createParser()
    args = parser.parse_args(sys_args)
    vcf_reader = VcfReader(args.vcfname)
    fasta_seq = pysam.FastaFile(args.refname)
    regions = RegionList(filename=args.bedname,zeroho=args.zeroho,zeroclosed=args.zeroclosed)

    for region in regions.regions:
        rec_list = vcf_reader.getRecordList(region)
        pass_list = getPassSites(rec_list, remove_cpg=True, fasta_ref=fasta_seq)
        if pass_list.count(True) >= int(args.minsites):
            print (region.toStr(sep='\t'))


if __name__ == '__main__':
    filter_bed_regions(sys.argv[1:])
