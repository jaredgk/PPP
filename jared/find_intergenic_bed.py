import sys
import pysam
from gene_region import Region, RegionList, getIntervalsBetween
import argparse_sets
import argparse



#Given a list of CDS intervals and optional buffer length,
#Generate full set of intervals between regions.
def createParser():
    parser = argparse.ArgumentParser(description=("Generates list of "
                "intervals within a provided BED file with option "
                "to pad intervals by a fixed amount"))
    argparse_sets.addRegionArgs(parser)
    parser.add_argument('--pad', dest="pad_count", default=0,
                        help="Extend input regions by provided value")
    parser.add_argument('--addchr', dest="add_chr", action="store_true")
    parser.add_argument('--out', dest="output_name",
                        help="Output filename, default is stdout")
    return parser



def main(sysargs):
    parser = createParser()
    args = parser.parse_args(sysargs[1:])
    reg_list = RegionList(filename=args.region_name, colstr=args.colstr,
                          zeroho=args.zeroho, zeroclosed=args.zeroclosed)
    out_list = getIntervalsBetween(reg_list, int(args.pad_count))
    out_list.printList(file_name=args.output_name, add_chr=args.add_chr)

if __name__ == "__main__":
    #initLogger()
    main(sys.argv)
