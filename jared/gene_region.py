import sys

class Region:
    def __init__(self,start,end,chrom):
        self.start = start
        self.end = end
        self.chrom = chrom


class RegionList:
    
    def __init__(self,filename,oneidx=False,collist=None):
        if collist is None:
            collist = [0,1,2]
        if filename is None:
            raise Exception("Filename for gene region list not provided")
        self.regions = []
        with open(filename,'r') as regionfile:
            for line in regionfile:
                la = line.strip().split()
                try:
                    start = int(la[collist[1]])
                    end = int(la[collist[2]])
                    chrom = la[collist[0]]
                except:
                    sys.stderr.write("Column is missing")
                if oneidx:
                    start -= 1
                    end -= 1
                self.regions.append(Region(start,end,chrom))
