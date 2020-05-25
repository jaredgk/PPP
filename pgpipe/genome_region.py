import sys
import re
import logging
from functools import total_ordering
from random import shuffle


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

def setRegionSort(sorttype,sortlist=None):
    valid_types = ['natural','string','list']
    if sorttype not in valid_types:
        raise Exception("Sorting type %s is not valid" % (sorttype))
    Region.sort_method = sorttype
    if sorttype == 'list' and sortlist is not None:
        setSortOrder(sortlist)

def setSortOrder(sortlist,trimchrom=False):
    hold_list = []
    for c in sortlist:
        if trimchrom and c[0:3] == 'chr':
            hold_list.append(c[3:])
        else:
            hold_list.append(c)
    Region.sort_order = hold_list

@total_ordering
class Region:
    sort_method = "natural"
    sort_order = None

    def __init__(self, start, end, chrom, fullline = None):
        """Zero-based, half open coordinates and chromosome info for
        a region in a genome. Coords will be formatted according to
        flags passed to RegionList, then stored in a single format.
        """
        self.start = start
        self.end = end
        self.chrom = chrom
        self.fullline = fullline
        #if chrom[0:3] == 'chr':
        #    self.chrom = chrom[3:]
        #else:
        #    self.chrom = chrom
    def __deepcopy__(self):
        return Region(self.start,self.end,self.chrom)

    def cloneRegion(self):
        return Region(self.start,self.end,self.chrom)


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
        if Region.sort_method == 'natural':
            k1 = self.getChromKey()
            k2 = other.getChromKey()
            if k1 != k2:
                return keyComp(k1,k2)
        elif Region.sort_method == "string":
            #k1 = self.chrom
            #k2 = other.chrom
            k1 = (self.chrom if self.chrom[:3] != 'chr' else self.chrom[3:])
            k2 = (other.chrom if other.chrom[:3] != 'chr' else other.chrom[3:])
            if k1 != k2:
                return k1 < k2
        else:
            if Region.sort_order is None:
                raise Exception("Sorting by list requires calling setSortOrder() and providing a list for the order")
            k1 = self.chrom
            k2 = other.chrom
            if k1 not in Region.sort_order:
                raise Exception("Key %s is not in sort order list" % (k1))
            if k2 not in Region.sort_order:
                raise Exception("Key %s is not in sort order list" % (k2))
            if k1 != k2:
                return Region.sort_order.index(k1) < Region.sort_order.index(k2)
        if self.start != other.start:
            return self.start < other.start
        return self.end < other.end

    def __len__(self):
        return self.end-self.start

    def containsRecord(self, rec):
        k1 = self.getChromKey()
        kr = getChromKey(rec.chrom)
        if k1 != kr:
            if keyComp(k1,kr):
                return 'before'
            return 'after'
        comp_pos = rec.pos-1
##        print(comp_pos,self.start,self.end)
        if comp_pos < self.start:
            return 'before'
        if comp_pos >= self.end:
            return 'after'
        return 'in'

    def toStr(self, zeroho=False, zeroclosed=False, sep=':'):
        if self.fullline is not None:
            return self.fullline
        start = self.start
        end = self.end
        if not zeroho:
            start += 1
        elif zeroclosed:
            end -= 1
        return str(self.chrom)+sep+str(start)+sep+str(end)
        #return str(start)+sep+str(end)+sep+str(self.chrom)

    def getReference(self, refseq):
        return refseq.fetch(self.chrom,self.start,self.end)




class RegionList:

    def __init__(self, filename=None, genestr=None, reglist=None,
                 zeroclosed=False, zeroho=False,
                 colstr=None, sortlist=True, checkoverlap=None,
                 sortmethod=None, sortorder=None, chromfilter=None,
                 list_template=None, randomize=False,keep_full_line=False):
        """Class for storing gene region information

        Will create a list that stores genomic regions imported from
        a file, string, or list. List can be sorted in various ways,
        filtered for regions on one or multiple chromosomes, convert
        between coordinate systems, and integrated with other PPP
        functions.

        Parameters
        ----------
        filename : str (None)
            Name of file with gene coordinate data
        genestr : str (None)
            Semicolon-separated list of region data, either in "start:end"
            or "start:end:chrom" format.
        zeroho : bool (False)
            If true, input is in zero-based, half-open coordinates. Since
            pysam uses this format, coordinates will be unaltered
            (default is to use one-based closed coordinates)
        zeroclosed : bool (False)
            If true, uses zero-based, closed coordinates, which adds 1
            to the end value for a region. Rarely if ever used option.
        colstr : str (None)
            Three-element string separated by semicolons, with
            each element being an integer corresponding to the 0-based
            index of the column of the start, end, and chromosome values
            for input region(s).
        sortlist : bool (True)
            Will sort list once initialized using assigned sorting scheme.
        checkoverlap : str (None, 'error','fix')
            If None, will not check for overlaps. If 'error', will raise
            an exception if overlap is found. If 'fix', overlapping
            regions will be merged. Note zero-based, half-open intervals
            do not overlap when adjacent ones share end/start position,
            but one-based closed intervals do.
        sortmethod : str (None,'natural','string','list')
            Determines method for sorting. 'natural' is pre-set method,
            passing None to the constructor will ensure sort type is
            unchanged from other calls. This sorts akin to 'sort -V' unix
            command, which breaks a string into ints and substrings for
            sorting. 'string' will sort via regular string method. 'list'
            will use list passed to sortorder as ordering of chromosomes.
        sortorder : list (None)
            If 'list' sortmethod is used, regions will be sorted according
            to the order specified in this list.
        chromfilter : str (None)
            If set, will allow only regions from passed chromosome
            into the region list.
        list_template : RegionList (None)
            If creating new list with regions from a different list,
            will use options from given list instead of arguments passed
            for the following arguments: zeroho, zeroclosed, sortlist,
            checkoverlap, and sortorder.
        randomize : bool (False)
            If set, will randomly shuffle regions in list.

        Exceptions
        ----------
        Both filename and genestr are None
        colstr has less than 2 or more than 3 values

        """
        self.collist = None
        if list_template is not None:
            zeroho = list_template.zeroho
            zeroclosed = list_template.zeroclosed
            sortlist = list_template.sortlist
            checkoverlap = list_template.checkoverlap
            sortorder = list_template.sortorder

        if self.collist is None:
            if colstr is None:
                self.collist = [1, 2, 0]
            else:
                self.parseCols(colstr)
        if filename is None and genestr is None and reglist is None:
            raise Exception(("A filename, region string, or list of "
                             "region strings must be provided for " "creating a gene region list"))
        if [filename,genestr,reglist].count(None) < 2:
            raise Exception(("Only one of filename, region string, and "
                            "region list can be specified for init"))
        if zeroho and zeroclosed:
            raise Exception(("Zero-halfopen and zero-closed options "
                             "cannot be invoked at the same time"))

        if sortlist and randomize:
            raise Exception("Sorting and randomizing are not compatible options for RegionList")
        if checkoverlap not in [None,'fix','error']:
            raise Exception("Invalid checkoverlap value: %s" % (checkoverlap))
        self.regions = []
        self.zeroho = zeroho
        self.zeroclosed = zeroclosed
        self.chromfilter = chromfilter
        self.keep_full_line = keep_full_line
        self.header = None
        if sortmethod is not None:
            setRegionSort(sortmethod,sortlist=sortorder)
        if filename is not None:
            self.initFile(filename)
        elif genestr is not None:
            self.initStr(genestr)
        else:
            self.initList(reglist)

        if sortlist:
            self.regions.sort()
        if checkoverlap is not None:
            if self.hasOverlap():
                if checkoverlap == "fix":
                    self.fixOverlap()
                elif checkoverlap == "error":
                    raise Exception("Region overlap detected")

        if randomize:
            shuffle(self.regions)

    def initFile(self,filename):
        """Initialize RegionList with a region file
        """
        past_header = False
        with open(filename, 'r') as regionfile:
            for line in regionfile:
                if line[0] == '#':
                    self.header = line.strip()
                    past_header = True
                    continue
                la = line.strip().split()
                if not past_header:
                    past_header = True
                    try:
                        t = int(la[self.collist[0]])
                    except ValueError as e:
                        self.header = line.strip()
                        continue
                self.initRegion(la)


    def initStr(self,genestr):
        la = genestr.split(':')
        self.initRegion(la)


    def initList(self,reglist):
        for reg in reglist:
            self.initRegion(reg)



    def initRegion(self,la):
        start = int(la[self.collist[0]])
        end = int(la[self.collist[1]])
        chrom = la[self.collist[2]]
        if self.chromfilter is not None and chrom != self.chromfilter:
            return
        if not self.zeroho:
            start -= 1
        elif self.zeroclosed:
            end += 1
        fullline = (('\t'.join(la)) if self.keep_full_line else None)
        self.regions.append(Region(start,end,chrom,fullline))


    def parseCols(self,cols):
        col_list = [int(i) for i in cols.split(',')]
        if len(col_list) != 3:
            raise Exception(("Column string %s must have three "
                             "values" % (cols)))
        self.collist = col_list

    def toStr(self, idx, sep=':'):
        return self.regions[idx].toStr(zeroho=self.zeroho,
                                       zeroclosed=self.zeroclosed,
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

    def fixOverlap(self):
        region_hold = []
        self.regions.sort()
        i = 0
        temp_reg = None
        for i in range(len(self.regions)):
            if temp_reg is None:
                temp_reg = self.regions[i].cloneRegion()
                continue
            cur_reg = self.regions[i].cloneRegion()
            if cur_reg.start < temp_reg.end and cur_reg.chrom == temp_reg.chrom:
                temp_reg.end = max(temp_reg.end,cur_reg.end)
            else:
                region_hold.append(temp_reg)
                temp_reg = cur_reg
        region_hold.append(temp_reg)
        self.regions = region_hold

    def printList(self, file_handle=None, file_name=None,
                  return_str=False, delim="\t", add_chr=False,
                  remove_chr=False):
        if file_handle is None and file_name is None:
            file_handle = sys.stdout
        if file_name is not None:
            file_handle = open(file_name, 'w')
        if return_str:
            out_str = ''
        if self.header is not None:
            file_handle.write(self.header+'\n')
        for region in self.regions:
            if region.fullline is not None:
                file_handle.write(region.fullline+'\n')
                continue
            start = region.start
            end = region.end
            if not self.zeroho:
                start += 1
            elif self.zeroclosed:
                end -= 1
            chrom = region.chrom
            if add_chr and chrom[:3] != 'chr':
                chrom = "chr"+region.chrom
            if remove_chr and chrom[:3] == 'chr':
                chrom = region.chrom[3:]
            reg_str = chrom+delim+str(start)+delim+str(end)+'\n'
            if return_str:
                out_str += reg_str
            else:
                file_handle.write(reg_str)
        if return_str:
            return out_str

    def filterByChrom(self,chrom_list,include=False):
        if include:
            self.regions = [r for r in self.regions if r.chrom in chrom_list]
        else:
            self.regions = [r for r in self.regions if r.chrom not in chrom_list]

    def filterOutXY(self):
        self.filterByChrom(['X','Y','chrX','chrY'])

    def expandRegions(self,expand_size):
        #Expands regions, but does not merge overlaps 
        #that may occur
        for i in range(len(self.regions)):
            self.regions[i].start = max(self.regions[i].start-expand_size,0)
            self.regions[i].end += expand_size

    #TO do:
    #Add sort method for strictly text based sorting

def getIntervalsBetween(region_list, padding=0, firstline=True):
    region_hold = []
    if len(region_list.regions) < 2:
        raise Exception("Region list for complement requires at least two regions")
    if firstline:
        r1 = region_list.regions[0]
        if r1.start - padding > 0:
            region_hold.append([str(r1.chrom),0,r1.start-padding])
    for i in range(len(region_list.regions)-1):
        r1 = region_list.regions[i]
        rn = region_list.regions[i+1]
        if r1.chrom != rn.chrom:
            continue
        new_start = r1.end+padding
        new_end = rn.start-padding
        if new_start < new_end:
            region_hold.append([str(r1.chrom),new_start,new_end])
    out_list = RegionList(reglist=region_hold, zeroho=True)
    out_list.zeroho = region_list.zeroho
    out_list.zeroclosed = region_list.zeroclosed
    return out_list

def subtractBed(stat_list, filter_list):
    #Sorts regions and removes regions from stat list 
    #that overlap with any region in filter_list
    stat_idx = 0
    filter_idx = 0
    stat_list.regions.sort()
    filter_list.regions.sort()
    drop_list = [False for i in stat_list.regions]
    stat_region = stat_list.regions[stat_idx]
    filter_region = filter_list.regions[filter_idx]
    while stat_idx < len(stat_list.regions):
        #try:
        #    print (stat_list.regions[stat_idx].toStr(),filter_list.regions[filter_idx].toStr())
        #except:
        #    break
        #print (stat_list.regions[stat_idx].toStr())
        while filter_idx < len(filter_list.regions) and stat_list.regions[stat_idx].end > filter_list.regions[filter_idx].start:
            #print (filter_list.regions[filter_idx].toStr())
            if stat_list.regions[stat_idx].start < filter_list.regions[filter_idx].end:
                #print ("drop "+str(stat_idx))
                drop_list[stat_idx] = True
                break
            filter_idx+=1
        stat_idx += 1
    stat_list.regions = [stat_list.regions[i] for i in range(len(stat_list.regions)) if drop_list[i] is False]
    return 
