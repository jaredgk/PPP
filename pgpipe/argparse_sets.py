import sys
import argparse

def addRegionArgs(parser):
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
