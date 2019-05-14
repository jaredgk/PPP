import sys
import argparse
from gene_region import RegionList, subtractBed

def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stat-file",dest="stat_file",type=str,help="Name of statistic BED file to be filtered")
    parser.add_argument("--filter-file",dest="filter_file",type=str,help="Name of file with regions that indicate ranges in stat file to be thrown")
    parser.add_argument("--window",type=int,help=("Expand regions in filter file by this number of bases"),default=0)
    parser.add_argument("--zero-ho",action="store_true",dest="zeroho")
    return parser

def filterStat(sysargv):
    parser = createParser()
    args = parser.parse_args(sysargv)
    stat_list = RegionList(filename=args.stat_file,zeroho=args.zeroho,keep_full_line=True)
    filter_list = RegionList(filename=args.filter_file,zeroho=args.zeroho)
    if args.window != 0:
        filter_list.expandRegions(args.window)
    subtractBed(stat_list,filter_list)
    stat_list.printList()


if __name__ == "__main__":
    filterStat(sys.argv[1:])