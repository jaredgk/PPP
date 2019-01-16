import sys
from vcf_reader_func import getPassSites, VcfReader
from gene_region import Region, RegionList
import pysam
import argparse
from logging_module import initLogger
import os

sys.path.insert(0,os.path.abspath(os.path.join(os.pardir, 'andrew')))

from model import Model, read_model_file

def createParser():
    parser = argparse.ArgumentParser(description="Given a VCF file and a list of intervals (BED style), will output regions that have over a certain number of qualifying variants in the region")
    parser.add_argument("--vcf", dest="vcfname")
    parser.add_argument("--bed", dest="bedname")
    parser.add_argument("--zero-ho", dest="zeroho", action="store_true")
    parser.add_argument("--zero-closed", dest="zeroclosed", action="store_true")
    parser.add_argument("--parsecpg", dest="refname")
    parser.add_argument("--remove-indels", dest="remove_indels", action="store_true", help=("Removes indels from output VCF files"))
    parser.add_argument("--remove-multi", dest="remove_multiallele", action="store_true")
    parser.add_argument("--remove-missing", dest="remove_missing", default=-1, type=int, help=("Will filter out site if more than the given number of individuals (not genotypes) are missing data. 0 removes sites with any missing data, -1 (default) removes nothing"))
    parser.add_argument("--informative-count", dest="informative_count", type=int, default=2)
    parser.add_argument("--minsites", dest="minsites", default=3, help=("Regions with at least this many variants passing filters will be output"))
    parser.add_argument("--tbi", dest="tabix_index", help="Path to bgzipped file's index if name doesn't match VCF file")
    parser.add_argument("--randcount",dest="randcount",type=int,default=-1,help="If set, will randomly draw from input region file until randcount # of passing BED regions are found")
    parser.add_argument("--no-xy",dest="filter_xy",action="store_true",help="Remove X/Y chromosomes from valid regions")
    parser.add_argument("--min-length",dest="min_length",type=int,default=1000,help="Minimum length of valid regions")
    parser.add_argument("--model",dest="modelname",help="Model file for selecting samples from VCF")
    parser.add_argument("--poptag",dest="poptag",help="Name of pop if model has more than one")
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

    popmodel = None
    if args.modelname is not None:
        popmodels = read_model_file(args.modelname)
        if len(popmodels) != 1:
            popmodel = popmodels[args.poptag]
        else:
            pp = list(popmodels.keys())
            popmodel = popmodels[pp[0]]
    
    vcf_reader = VcfReader(args.vcfname,index=args.tabix_index,popmodel=popmodel)
    fasta_seq = None
    if args.refname is not None:
        fasta_seq = pysam.FastaFile(args.refname)

    #regions = RegionList(filename=args.bedname,zeroho=args.zeroho,zeroclosed=args.zeroclosed,sortlist=(not args.randcoun))
    randomize = False
    if args.randcount != -1:
        randomize = True

    regions = RegionList(filename=args.bedname,zeroho=args.zeroho,
                         zeroclosed=args.zeroclosed,sortlist=(not randomize),
                         randomize=randomize)
    if args.filter_xy:
        regions.filterOutXY()
    regions_output = 0
    for region in regions.regions:
        if len(region) < args.min_length:
            continue
        rec_list = vcf_reader.getRecordList(region)
        pass_list = getPassSites(rec_list, remove_cpg=True,
                    remove_indels=args.remove_indels,
                    remove_multiallele=args.remove_multiallele,
                    remove_missing=args.remove_missing,
                    inform_level=args.informative_count,
                    fasta_ref=fasta_seq)
        if pass_list.count(True) >= int(args.minsites):
            print (region.toStr(sep='\t'))
            regions_output += 1
        if args.randcount != -1 and regions_output == args.randcount:
            break

    if args.randcount != -1 and regions_output != args.randcount:
        sys.stderr.write("Only %d of %d regions found\n"%(regions_output,args.randcount))
        exit(1)


if __name__ == '__main__':
    initLogger()
    filter_bed_regions(sys.argv[1:])
