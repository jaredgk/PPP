import sys
import re
import logging
from functools import total_ordering


def getChromKey(chrom):
    """Get key from chromosome string so that it can be naturally sorted
    (int < string, numbers grouped together)
    """
    c = chrom
    if chrom[0:3] == 'chr':
        c = chrom[3:]
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    k = [ convert(i) for i in re.split('([0-9]+)',c) if len(i) != 0]
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
    def __init__(self, start, end, chrom, natsort=True):
        """Zero-based, half open coordinates and chromosome info for
        a region in a genome. Coords will be formatted according to
        flags passed to RegionList, then stored in a single format.
        """
        self.start = start
        self.end = end
        if chrom[0:3] == 'chr':
            self.chrom = chrom[3:]
        else:
            self.chrom = chrom

    def getChromKey(self):
        """Splits chromosone name for natural sort. Ints and strs are
        grouped with adjacent elements of the same type

        Example: 'front88place' returns ['front',88,'place']
        """
        return getChromKey(self.chrom)


    def __eq__(self, other):
        return ((self.start == other.start) and (self.end == other.end)
                and (self.chrom == other.chrom))

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

    def containsRecord(self, rec):
        k1 = self.getChromKey()
        kr = getChromKey(rec.chrom)
        if k1 != kr:
            if keyComp(k1,kr):
                return 'before'
            return 'after'
        if rec.pos < self.start:
            return 'before'
        if rec.pos >= self.end:
            return 'after'
        return 'in'

    def toStr(self, halfopen=True, oneidx=False, sep=':'):
        start = self.start
        end = self.end
        if not halfopen:
            end -= 1
        if oneidx:
            start += 1
            end += 1
        return str(start)+sep+str(end)+sep+str(self.chrom)


class RegionList:

    def __init__(self, filename=None, genestr=None, oneidx=False,
                 colstr=None, reglist=None, defaultchrom=None, halfopen=True,
                 sortlist=True, checkoverlap=True, natsort=True):
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
        if filename is None and genestr is None and reglist is None:
            raise Exception(("Either a filename or gene region must be "
                            "provided for creating a gene region list"))
        self.regions = []
        self.oneidx = oneidx
        self.halfopen = halfopen
        if filename is not None:
            self.initFile(filename, oneidx, colstr, defaultchrom,
                          halfopen, sortlist, checkoverlap)
        elif genestr is None:
            self.initStr(genestr, oneidx, colstr, defaultchrom, halfopen)
        else:
            self.initList(reglist, oneidx, colstr, defaultchrom,
                          halfopen, sortlist, checkoverlap)

    def initFile(self, filename, oneidx, colstr, defaultchrom, halfopen,
                 sortlist, checkoverlap):
        """Initialize RegionList with a region file
        """
        with open(filename, 'r') as regionfile:
            for line in regionfile:
                if line[0] == '#':
                    continue
                la = line.strip().split()
                start = int(la[self.collist[0]])
                end = int(la[self.collist[1]])
                if len(self.collist) == 3:
                    chrom = la[self.collist[2]]
                elif defaultchrom is not None:
                    chrom = defaultchrom
                else:
                    raise Exception("Chromosome for region is not specified")
                if oneidx:
                    start -= 1
                    end -= 1
                if not halfopen:
                    end += 1
                self.regions.append(Region(start, end, chrom))
        if sortlist:
            self.regions.sort()
        if checkoverlap:
            if self.hasOverlap():
                raise Exception("Region overlap detected")

    def initStr(self, genestr, oneidx, colstr, defaultchrom, halfopen):
        la = genestr.split(':')
        start = int(la[0])
        end = int(la[1])
        if len(la) == 3:
            chrom = la[2]
        elif defaultchrom is not None:
            chrom = defaultchrom
        else:
            raise Exception("Region has no chromosome provided")
        if oneidx:
            start -= 1
            end -= 1
        if not halfopen:
            end += 1
        self.regions.append(Region(start,end,chrom))

    def initList(self, reglist, oneidx, colstr, defaultchrom, halfopen,
                 sortlist, checkoverlap):
        for reg in reglist:
            start = int(reg[self.collist[0]])
            end = int(reg[self.collist[1]])
            if len(self.collist) == 3:
                chrom = reg[self.collist[2]]
            elif defaultchrom is not None:
                chrom = defaultchrom
            else:
                raise Exception("Chromosome for region is not specified")
            if oneidx:
                start -= 1
                end -= 1
            if not halfopen:
                end += 1
            self.regions.append(Region(start, end, chrom))
        if sortlist:
            self.regions.sort()
        if checkoverlap:
            if self.hasOverlap():
                raise Exception("Region overlap detected")

    def parseCols(self,cols):
        col_list = [int(i) for i in cols.split(',')]
        if len(col_list) < 2 or len(col_list) > 3:
            raise Exception(("Column string %s is either too short or "
                             "too long" % (cols)))
        self.collist = col_list

    def toStr(self, idx, sep=':'):
        return self.regions[idx].toStr(halfopen=self.halfopen,
                                       oneidx=self.oneidx,
                                       sep=sep)

    def hasOverlap(self):
        region_hold = sorted(self.regions)
        for i in range(len(region_hold)-1):
            if region_hold[i].chrom == region_hold[i+1].chrom:
                if (region_hold[i].start <= region_hold[i+1].start
                    < region_hold[i].end):
                    logging.warning("Regions %s and %s overlap" %
                                    (self.toStr(i),self.toStr(i+1)))
                    return True
        return False

    def printList(self, file_handle=None, file_name=None, delim="\t"):
        if file_handle is None and file_name is None:
            raise Exception("Either file_handle or file_name must be set")
        if file_name is not None:
            file_handle = open(file_name, 'w')
        for region in self.regions:
            start = region.start
            end = region.end
            if not self.halfopen:
                end -= 1
            if self.oneidx:
                start += 1
                end += 1
            reg_str = region.chrom+delim+str(start)+delim+str(end)+'\n'
            file_handle.write(reg_str)


    #TO do:
    #Add sort method for strictly text based sorting
