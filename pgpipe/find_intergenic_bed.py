import sys
import pysam
from pgpipe.gene_region import Region, RegionList, getIntervalsBetween
import pgpipe.argparse_sets
import argparse



#Given a list of CDS intervals and optional buffer length,
#Generate full set of intervals between regions.
def createParser():
    parser = argparse.ArgumentParser(description=("Generates list of "
                "intervals within a provided BED file with option "
                "to pad intervals by a fixed amount"))
    parser.add_argument('--bed', dest="region_name",
                        help="Name of gene region file")
    parser.add_argument("--bed-column-index", dest="colstr",
                        help=("Comma-separated list of length 3 with 0-based"
                        " indexes of start, end, and chromosome data in input"
                        " BED file. Default for normal BED is 1,2,0"))
    parser.add_argument('--zero-ho', dest="zeroho", action="store_true",
                        help="Region list is 1 indexed, not 0")
    parser.add_argument('--zero-closed', dest="zeroclosed", action="store_true",
                        help="Use zero-based, closed coordinates")
    parser.add_argument('--pad', dest="pad_count", default=0,
                        help="Extend input regions by provided value")
    cg = parser.add_mutually_exclusive_group()
    cg.add_argument('--addchr', dest="add_chr", action="store_true")
    cg.add_argument('--removechr',dest="remove_chr",action="store_true")
    parser.add_argument('--out', dest="output_name",
                        help="Output filename, default is stdout")
    return parser



def get_intergenic(sysargs):
    """
    Creates a BED file with regions that are not covered in the input BED.

    This function is similar to 'bedtools compliment', except it doesn't 
    require a chromosome end position file (which results in the last region 
    of a chromosome being thrown out). A flanking distance from the ends of 
    input regions can be set with --pad to truncate regions (for example, if 
    trying to find regions that are outside a 5Kb window of a gene).
    Additional options are included for dealing with zero-based, half-open
    based BED files, and for adding/removing 'chr' from chromosome names. 
    By default list will be sorted and overlapping regions will be merged.

    Parameters
    ----------
    --bed : str (required)
        Name of input BED file
    --bed-column-index : ints (3)
        Comma-separated list of length 3 with 0-based indexes of start, end,
        and chromosome data in input BED file. Default for normal BED is
        "1,2,0".
    --zero-ho : bool
        If set, treats regions in BED file as 0-based, (h)alf (o)pen 
        coordinates rather than 1-based, closed
    --pad : int
        Distance to pad input windows. If padding between two 
        consecutive regions cross, will not output region for area.
    --out : str
        Output filename

    Other Parameters
    ----------------
    --addchr : bool
        If set, will add 'chr' to all chromosome names
    --removechr : bool
        If set, will remove 'chr' from all chromosome names
        starting with it
    


    """
    parser = createParser()
    args = parser.parse_args(sysargs)
    if args.region_name is None:
        raise Exception("BED input filename required")
    reg_list = RegionList(filename=args.region_name, colstr=args.colstr,
                          zeroho=args.zeroho, zeroclosed=args.zeroclosed)
    out_list = getIntervalsBetween(reg_list, int(args.pad_count))
    out_list.printList(file_name=args.output_name, add_chr=args.add_chr)

if __name__ == "__main__":
    #initLogger()
    get_intergenic(sys.argv[1:])
