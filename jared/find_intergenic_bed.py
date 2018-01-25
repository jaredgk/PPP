import sys
import pysam
from gene_region import Region, RegionList
import argparse_sets
import argparse

def getIntervalsBetween(region_list, padding=0):
    region_hold = []
    for i in range(len(region_list.regions)-1):
        r1 = region_list.regions[i]
        rn = region_list.regions[i+1]
        if r1.chrom != rn.chrom:
            continue
        new_start = r1.end+padding
        new_end = rn.start-padding
        if new_start < new_end:
            region_hold.append([str(r1.chrom),new_start,new_end])
    out_list = RegionList(reglist=region_hold, oneidx=region_list.oneidx,
                          halfopen=region_list.halfopen)
    return out_list

#Given a list of CDS intervals and optional buffer length,
#Generate full set of intervals between regions.
def createParser():
    parser = argparse.ArgumentParser(description=("Generates list of "
                "intervals within a provided BED file with option "
                "to pad intervals by a fixed amount"))
    argparse_sets.addRegionArgs(parser)
    parser.add_argument('--pad', dest="pad_count", default=0,
                        help="Extend input regions by provided value")
    parser.add_argument('--out', dest="output_name",
                        help="Output filename, default is stdout")
    return parser



def main(sysargs):
    parser = createParser()
    args = parser.parse_args(sysargs[1:])
    reg_list = RegionList(filename=args.region_name, colstr=args.colstr,
                            oneidx=args.region_idx, halfopen=args.halfopen)
    out_list = getIntervalsBetween(reg_list, int(args.pad_count))
    out_list.printList()

if __name__ == "__main__":
    #initLogger()
    main(sys.argv)
