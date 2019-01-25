import sys
import pysam
from pppipe.gene_region import Region, RegionList, getIntervalsBetween
import pppipe.argparse_sets
import argparse



#Given a list of CDS intervals and optional buffer length,
#Generate full set of intervals between regions.
def createParser():
    parser = argparse.ArgumentParser(description=("Generates list of "
                "intervals within a provided BED file with option "
                "to pad intervals by a fixed amount"))
    parser.add_argument('--bed', dest="region_name",
                        help="Name of gene region file")
    parser.add_argument('--zero-ho', dest="zeroho", action="store_true",
                        help="Region list is 1 indexed, not 0")
    parser.add_argument('--zero-closed', dest="zeroclosed", action="store_true",
                        help="If set, use closed coordinates instead of half-open")
    parser.add_argument('--regcol', dest='colstr', help= (
                        "Comma-separated list of columns for gene region "
                        " data, format is start/end if no chromosome "
                        " data, start/end/chrom if so"))
    parser.add_argument('--pad', dest="pad_count", default=0,
                        help="Extend input regions by provided value")
    cg = parser.add_mutually_exclusive_group()
    cg.add_argument('--addchr', dest="add_chr", action="store_true")
    cg.add_argument('--removechr',dest="remove_chr",action="store_true")
    parser.add_argument('--out', dest="output_name",
                        help="Output filename, default is stdout")
    return parser



def get_intergenic(sysargs):
    parser = createParser()
    args = parser.parse_args(sysargs)
    reg_list = RegionList(filename=args.region_name, colstr=args.colstr,
                          zeroho=args.zeroho, zeroclosed=args.zeroclosed)
    out_list = getIntervalsBetween(reg_list, int(args.pad_count))
    out_list.printList(file_name=args.output_name, add_chr=args.add_chr)

if __name__ == "__main__":
    #initLogger()
    get_intergenic(sys.argv[1:])
