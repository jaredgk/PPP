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


def getChromKey(chrom):
    """Get key from chromosome string so that it can be naturally sorted
    (int < string, numbers grouped together)
    """
    c = chrom
    if chrom[0:3] == 'chr':
        c = chrom[3:]
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    k = [ convert(i) for i in re.split('([0-9]+)',c)]
    return k

def keyComp(k1,k2):
    """Returns True if k1 < k2, False if not. int < str"""
    for i in range(min(len(k1),len(k2))):
        e1 = k1[i]
        e2 = k2[i]
        if isinstance(e1,int) != isinstance(e2,int):
            return isinstance(e1,int)
        if e1 != e2:
            return e1 < e2
    return len(k1) < len(k2)


@total_ordering
class Region:
    def __init__(self, start, end, chrom):
        """Zero-based, half open coordinates and chromosome info for
        a region in a genome. Coords will be formatted according to
        flags passed to RegionList, then stored in a single format.
        """
        self.start = start
        self.end = end
        self.chrom = chrom

    def chromMod(self):
        """Removes 'chr' from front of chromosome. Used for sorting"""
        if self.chrom[0:3] == 'chr':
            c = self.chrom[3:]
        else:
            c = self.chrom
        return c

    def getChromKey(self):
        """Splits chromosone name for natural sort. Ints and strs are
        grouped with adjacent elements of the same type

        Example: 'front88place' returns ['front',88,'place']
        """
        c = self.chromMod()
        convert = lambda text: int(text) if text.isdigit() else text.lower()
        k = [convert(i) for i in re.split('([0-9]+)',c) if len(i) != 0]
        return k


    def __eq__(self, other):
        c = self.chromMod()
        oc = other.chromMod()
        return ((self.start == other.start) and (self.end == other.end)
                and (c == oc))

    def __lt__(self, other):
        """Sort key order: chrom-key, start position, end position
        """
        k1 = self.getChromKey()
        k2 = other.getChromKey()
        if k1 != k2:
            return keyComp(k1,k2)
        if self.start != other.start:
            return self.start < other.start
        return self.end < other.end


class RegionList:

    def __init__(self, filename=None, genestr=None, oneidx=False,
                 colstr=None, defaultchrom=None, halfopen=True,
                 sortlist=True):
        """Class for storing gene region information

        Will either read a file or take a single gene region in a string
        as input, then create an object with some metadata and a list of
        gene regions.

        Parameters
        ----------
        filename : str (None)
            Name of file with gene coordinate data
        genestr : str (None)
            Semicolon-separated list of region data, either in "start:end"
            or "start:end:chrom" format.
        oneidx : bool (False)
            If true, region will be read in as a one-index based string.
            This results in one position being subtracted from the start
            and end coordinates
        colstr : str (None)
            If reading in from a file with many columns, will indicate
            what columns hold start/end data and chrom data if three
            items are given.
        defaultchrom : str (None)
            If no chromosome data is provided in the region file or string,
            will set the chromosome to this value
        halfopen : bool (True)
            If true, the region range will be [start,end), so the start
            coordinate will be included in the region but the end
            will not

        Exceptions
        ----------
        Both filename and genestr are None
        colstr has less than 2 or more than 3 values

        """
        if colstr is None:
            self.collist = [1, 2, 0]
        else:
            self.parseCols(colstr)
        if filename is None and genestr is None:
            raise Exception(("Either a filename or gene region must be "
                            "provided for creating a gene region list"))
        self.regions = []
        self.oneidx = oneidx
        self.halfopen = halfopen
        if filename is not None:
            self.initList(filename, oneidx, colstr, defaultchrom,
                          halfopen, sortlist)
        else:
            self.initStr(genestr, oneidx, colstr, defaultchrom, halfopen)

    def initList(self, filename, oneidx, colstr, defaultchrom, halfopen,
                 sortlist):
        """Initialize RegionList with a region file
        """
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
        if sortlist:
            self.regions.sort()

    def initStr(self, genestr, oneidx, colstr, defaultchrom, halfopen):
        la = genestr.split(':')
        start = int(la[0])
        end = int(la[1])
        if len(la) == 3:
            c = la[2]
            chrom = matchChr(defaultchrom,c)
        if oneidx:
            start += 1
            end += 1
        if not halfopen:
            end += 1
        self.regions.append(Region(start,end,chrom))

    def parseCols(self,cols):
        col_list = [int(i) for i in cols.split(',')]
        if len(col_list) < 2 or len(col_list) > 3:
            raise Exception(("Column string %s is either too short or "
                             "too long" % (cols)))
        self.collist = col_list

    #TO do: add check for intersecting regions
    #Add sort method for strictly text based sorting
