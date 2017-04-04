import sys
import re
from functools import total_ordering

def matchChr(ref_chr,list_chr):
    """Checks ref_chr for whether "chr" is first chars, then modifies
    list_chr to match"""
    if ref_chr is None:
        return list_chr
    has_chr = (ref_chr[0:3]=="chr")
    list_has_chr = (list_chr[0:3]=="chr")
    if has_chr == list_has_chr:
        return list_chr
    elif has_chr:
        return "chr"+list_chr
    return list_chr[3:]




def sortChrom(c1,c2):
    l1 = parseChrom(c1)
    l2 = parseChrom(c2)

def getChromKey(chrom):
    c = chrom
    if chrom[0:3] == 'chr':
        c = chrom[3:]
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    k = [ convert(i) for i in re.split('([0-9]+)',c)]
    return k

@total_ordering
class Region:
    def __init__(self, start, end, chrom):
        self.start = start
        self.end = end
        self.chrom = chrom

    def chromMod(self):
        if self.chrom[0:3] == 'chr':
            c = self.chrom[3:]
        else:
            c = self.chrom
        return c

    def getChromKey(self):
        c = self.chromMod()
        convert = lambda text: int(text) if text.isdigit() else text.lower()
        k = [convert(i) for i in re.split('([0-9]+)',c)]
        return k


    def __eq__(self, other):
        c = self.chromMod()
        oc = other.chromMod()
        return ((self.start == other.start) and (self.end == other.end)
                and (c == oc))

    def __lt__(self, other):
        k1 = self.getChromKey()
        k2 = other.getChromKey()
        return k1 < k2 and self.start < other.start and self.end < other.end


class RegionList:

    def __init__(self, filename, oneidx=False, colstr=None,
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
                        c = la[self.collist[2]]
                        chrom = matchChr(defaultchrom,c)
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
        col_list = [int(i) for i in cols.split(',')]
        if len(col_list) < 2 or len(col_list) > 3:
            raise Exception(("Column string %s is either too short or "
                             "too long" % (cols)))
        self.collist = col_list
