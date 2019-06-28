import sys
import argparse
import logging
from pgpipe.genome_region import RegionList, subtractBed

def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stat-file",dest="stat_file",type=str,help="Name of statistic BED file to be filtered")
    parser.add_argument("--filter-file",dest="filter_file",type=str,help="Name of file with regions that indicate ranges in stat file to be thrown")
    parser.add_argument("--window",type=int,help=("Expand regions in filter file by this number of bases"),default=0)
    parser.add_argument("--zero-ho",action="store_true",dest="zeroho")
    parser.add_argument("--out",type=str)
    return parser

def filter_stat(sysargv):
    parser = createParser()
    args = parser.parse_args(sysargv)
    stat_list = RegionList(filename=args.stat_file,zeroho=args.zeroho,keep_full_line=True)
    filter_list = RegionList(filename=args.filter_file,zeroho=args.zeroho)
    original_len = len(stat_list.regions)
    if args.window != 0:
        filter_list.expandRegions(args.window)
    subtractBed(stat_list,filter_list)
    logging.warning("%d of %d regions selected as non-overlapping"%(len(stat_list.regions),original_len))
    stat_list.printList(file_name=args.out)


if __name__ == "__main__":
    filter_stat(sys.argv[1:])