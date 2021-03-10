import sys
import argparse

from pgpipe.vcf_reader_func import getAlleleStats, VcfReader, getAlleleCountDict
from pgpipe.genome_region import RegionList, Region

class recordCache():
    def __init__(self):
        self.chrom = ''
        self.cache = {}

    def getSiteStat(self,record):
        return piSite(record)

    def getWindowStat(self,recgen,region):
        return piWindow(recgen,region)

def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf",dest="vcf",type=str)
    parser.add_argument("--window-stat",dest="window",action="store_true")
    parser.add_argument("--window-size",dest="window_size",type=int)
    parser.add_argument("--step-size",dest="step_size",type=int)
    parser.add_argument("--bed",dest="bedfile",type=str)
    parser.add_argument("--subsamp",dest="subsamp",type=str)
    return parser

def piSite(record):
    alleles,total_sites,missing_inds = getAlleleStats(record)
    #alleles,total_sites,missing_inds = getAlleleStatsTwo(record)
    #return 2
    site_sum = 0
    for a in record.alleles:
        #print (a)
        a_count = alleles[a]
        site_sum += (a_count)*(total_sites-a_count)
    return float(site_sum)/(float(total_sites)*float(total_sites-1))

def getCacheVals(cache,region):
    res_list = []
    for pos in list(cache.keys()):
        if int(pos) < region.start:
            del cache[pos]
        elif int(pos) < region.end:
            res_list.append(cache[pos])
    return res_list


def piWindow(recgen,region,cache=None):
    pi_sum = 0
    if cache is not None:
        cache_list = getCacheVals(cache,region)
        for cv in cache_list:
            pi_sum += cv
    for rec in recgen:
        ps = piSite(rec)
        pi_sum += ps
        if cache is not None:
            cache[str(rec.pos)] = ps
    if pi_sum == 0:
        return -1
    return pi_sum/float(len(region))



parser = createParser()
args = parser.parse_args()
if args.step_size is not None and args.window_size is None:
    raise Exception("--step-size requires --window-size")
if args.step_size is None and args.window_size is not None:
    args.step_size = args.window_size

idx_list = None
if args.subsamp is not None:
    idx_list = [int(i.strip()) for i in open(args.subsamp).readlines()]
cache = None
if args.window_size > args.step_size:
    cache = {}
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
        if args.window_size is not None:
            eof = False
            if region is None:
                current_chrom = vcf_reader.info_rec.chrom
                subregion = Region(0,args.window_size,current_chrom)
            else:
                region_size = len(region)
                subregion = Region(region.start,min(region.start+args.window_size,region.end),region.chrom)
            while True:
                pw = piWindow(vcf_reader.getRegionIter(region=subregion),subregion,cache)
                if pw > -.5:
                    sys.stdout.write(subregion.toStr(zeroho=True)+'\t'+str(pw)+'\n')
                if region is None:
                    if vcf_reader.prev_last_rec is None:
                        eof = True
                        break
                    if vcf_reader.prev_last_rec.chrom != current_chrom:
                        current_chrom = vcf_reader.prev_last_rec.chrom
                        subregion = Region(0,args.window_size,current_chrom)
                    else:
                        subregion.start += args.step_size
                        subregion.end += args.step_size
                else:
                    if vcf_reader.prev_last_rec is None:
                        eof = True
                        break
                    if vcf_reader.prev_last_rec.chrom != region.chrom:
                        break
                    else:
                        subregion.start += args.step_size
                        subregion.end = min(subregion.end+args.step_size,region.end)
                        if subregion.start > region.end:
                            break
            

        else:
            pw = piWindow(vcf_reader.getRegionIter(region=region),region,cache)
            if pw > -.5:
                sys.stdout.write(region.toStr(zeroho=True)+'\t'+str(pw)+'\n')
    else:
        for rec in vcf_reader.getRegionIter(region=region):
            sys.stdout.write(str(rec.pos)+'\t'+str(piSite(rec))+'\n')

