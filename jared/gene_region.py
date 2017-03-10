import sys


class Region:
    def __init__(self, start, end, chrom):
        self.start = start
        self.end = end
        self.chrom = chrom

    def __cmp__(self, other):
        return self.__dict__ != other.__dict__


class RegionList:

    def __init__(self, filename, oneidx=True, collist=None,
                 defaultchrom=None):
        if collist is None:
            collist = [1, 2, 0]
        if filename is None:
            raise Exception("Filename for gene region list not provided")
        self.regions = []
        with open(filename, 'r') as regionfile:
            for line in regionfile:
                la = line.strip().split()
                try:
                    start = int(la[collist[0]])
                    end = int(la[collist[1]])
                    if len(collist) == 3:
                        chrom = la[collist[2]]
                    else:
                        chrom = defaultchrom
                except:
                    sys.stderr.write("Column is missing")
                if oneidx:
                    start -= 1
                    end -= 1
                self.regions.append(Region(start, end, chrom))
