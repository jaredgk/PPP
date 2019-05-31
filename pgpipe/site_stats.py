import sys
import argparse

from pgpipe.vcf_reader_func import getAlleleStats, VcfReader, getAlleleCountDict, getAlleleStatsTwo
from pgpipe.gene_region import RegionList

class recordCache():
    def __init__(self):
        self.chrom = ''
        self.cache = {}

def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf",dest="vcf",type=str)
    parser.add_argument("--window",dest="window",action="store_true")
    parser.add_argument("--bed",dest="bedfile",type=str)
    return parser

def piSite(record):
    #alleles,total_sites,missing_inds = getAlleleStats(record)
    alleles,total_sites,missing_inds = getAlleleStatsTwo(record)
    #return 2
    site_sum = 0
    for a in record.alleles:
        #print (a)
        a_count = alleles[a]
        site_sum += (a_count)*(total_sites-a_count)
    return float(site_sum)/(float(total_sites)*float(total_sites-1))

def piWindow(recgen,region):
    pi_sum = 0
    for rec in recgen:
        pi_sum += piSite(rec)
    if pi_sum == 0:
        return -1
    return pi_sum/float(len(region))

parser = createParser()
args = parser.parse_args()
#vcfname = str(sys.argv[1])
reglist = None
#if len(sys.argv) == 3:
if args.bedfile is not None:
    regname = str(sys.argv[2])
    reglist = RegionList(filename=args.bedfile,zeroho=True)

vcf_reader = VcfReader(args.vcf)

if reglist is None:
    #region = None
    regcount = 1
else:
    regcount = len(reglist.regions)

for i in range(regcount):
    if reglist is None:
        region = None
    else:
        region = reglist.regions[i]
    if args.window:
        pw = piWindow(vcf_reader.getRegionIter(region=region),region)
        if pw > -.5:
            sys.stdout.write(region.toStr(zeroho=True)+'\t'+str(pw)+'\n')
    else:
        for rec in vcf_reader.getRegionIter(region=region):
            sys.stdout.write(str(rec.pos)+'\t'+str(piSite(rec))+'\n')

