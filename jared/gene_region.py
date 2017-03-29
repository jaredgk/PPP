import sys


class Region:
    def __init__(self, start, end, chrom):
        self.start = start
        self.end = end
        self.chrom = chrom

    def __eq__(self, other):
        return self.__dict__ == other.__dict__


class RegionList:

    def __init__(self, filename, oneidx=True, colstr=None,
                 defaultchrom=None, halfopen=True):
        if colstr is None:
            self.collist = [1, 2, 0]
        else:
            self.parseCols(colstr)
        if filename is None:
            raise Exception("Filename for gene region list not provided")
        self.regions = []
        self.oneidx = oneidx
        self.halfopen = halfopen
        with open(filename, 'r') as regionfile:
            for line in regionfile:
                la = line.strip().split()
                try:
                    start = int(la[self.collist[0]])
                    end = int(la[self.collist[1]])
                    if len(self.collist) == 3:
                        chrom = la[self.collist[2]]
                    else:
                        chrom = defaultchrom
                except:
                    sys.stderr.write("Column is missing")
                if oneidx:
                    start -= 1
                    end -= 1
                if not halfopen:
                    end += 1
                self.regions.append(Region(start, end, chrom))

    def parseCols(self,cols):
        col_list = map(int,cols.split(','))
        if len(col_list) < 2 or len(col_list) > 3:
            raise Exception(("Column string %s is either too short or "
                             "too long" % (cols)))
        self.collist = col_list
