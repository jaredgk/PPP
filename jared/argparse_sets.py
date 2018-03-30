import sys
import argparse

def addRegionArgs(parser):
    parser.add_argument('--regionf', dest="region_name",
                        help="Name of gene region file")
    parser.add_argument('--gr1', dest="region_idx", action="store_true",
                        help="Region list is 1 indexed, not 0")
    parser.add_argument('--halfopen', dest="halfopen", action="store_false",
                        help="If set, use closed coordinates instead of half-open")
    parser.add_argument('--regcol', dest='colstr', help= (
                        "Comma-separated list of columns for gene region "
                        " data, format is start/end if no chromosome "
                        " data, start/end/chrom if so"))
