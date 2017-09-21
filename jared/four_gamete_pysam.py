
"""
    samples an interval based on four-gamete tests from aligned sequences
"""
"""
    written for python 3.6
    2 possible kinds of 4 gamete-based intervals:
        HandK85, minimal parsimonious set of intervals, such that each contains
            one or more recombination events:
            each interval has endpoints that are informative sites that are
            incompatible on the basis of 4gamete criterion
            each is inferred to contain at least 1 recombination event
            a set of intervals form a  parsimonious non-overlapping set that
            are consistent with the data regions not in those intervals show no
            evidence of recombination
        CONTIG4GPASS, set of most inclusive intervals, such that all the sites
            in each are fully compatible with each other:
            each interval is a stretch of bases that contain zero incompatible
            sites by the 4 gamete criterion
            on each side,  the left and right,  the next polymorphic site that
            is not in the interval must be incompatible with one of the
            polymorphic sites that is in the interval.

    input filetypes:   FASTA, SITE, VCF
    kinds of intervals:   HandK85  CONTIG4GPASS (default)
    kinds of returns:
            all intervals
            single interval
    single interval return options:
            random interval - uniform over list of intervals
            randombase interval - prob sampling proportional to interval length
            leftmost
            rightmost
            interval with at least n informative sites:
                random_n
                randombase_n
                leftmost_n
                rightmost_n
    Maximum number of informative SNPs is 1000
        set in buildlistsofsites()

"""
import sys
import argparse
import logging
from logging_module import initLogger
import vcf_load_4g as vcfload
import random
from gene_region import Region
import pysam


class BaseData():

    def __init__(self, format, filename, region=None):
        self.poslist = []
        self.seqs = []
        self.format = format
        if format == 'VCF':
            self.onlysnps = True
            vcff = pysam.VariantFile(filename)
            record_list = getRecordList(vcff, region)
            self.ploidy, self.seqs, self.poslist =
            vcfload.checkVcfRegionPysam(record_list)
            #self.seqcount = len(self.seqs)
            #self.numbases = listcheckallsamelength(self.seqs)
        elif format == "FASTA":
            self.onlysnps = False
            f = open(filename,"r")
            seq = ""
            while True:
                s = f.readline().strip()
                if len(s) == 0:
                    if len(seq) > 0:
                        self.seqs.append(seq)
                    break
                else:
                    if s[0] == '>':
                        if len(seq) > 0:
                            self.seqs.append(seq)
                        seq = ""
                    else:
                        s = s.strip()
                        if len(s) > 0:
                            seq += s
        elif format == "SIT":
            self.onlysnps = False
            f = open(filename,"r")
            f.readline()
            while True:
                s = f.readline()
                if s[0] != '#':
                    break
            numseq = int(s.split()[0])
            seqlen = int(s.split()[1])
            s = f.readline().split()

            numnoncode = int(s[1])
            if numnoncode != seqlen:
                for i in range(numnoncode):
                    f.readline()
            s = f.readline()
            if s.find("All") < 0:
                numpops = int(s.strip())
                for i in range(numpops):
                    f.readline()
            ## now should be at beginning of data
            for i in range(numseq):
                self.seqs.append(f.readline().strip()[10:])
        self.numbases = listcheckallsamelength(self.seqs)
        self.seqcount = len(seqs)


def createParser():
    parser = argparse.ArgumentParser(description=("Given a file of aligned"
            "sequences or variable sites, will return intervals based"
            " on 4 gamete tests"))
    parser.add_argument("filename", help="Input filename")
    filetype_group = parser.add_mutually_exclusive_group(required=True)
    filetype_group.add_argument("--fasta", dest = "fasta",
            nargs="+", help = "fasta file(s) of"
            " multiple aligned sequences, all the same length")
    filetype_group.add_argument("--sit", dest = "sit",
            nargs="+", help = " SITES format file(s)")
    filetype_group.add_argument("--vcf", dest = "vcf", nargs=2
            help = "VCF file and region list, in that order")
    filetype_group.add_argument("--vcfs", dest="vcfs", nargs="+")
    output_group = parser.add_mutually_exclusive_group(required=True)
    output_group.add_argument("--out", dest="out")
    output_group.add_argument("--out-prefix", dest="out_prefix")
    output_group.add_argument("--out-reg", dest="out_reg")
    intervaltype_group = parser.add_mutually_exclusive_group(required=
            "--reti" in sys.argv)
    intervaltype_group.add_argument("--hk", dest="intervaltype",
            action='store_const', const="HandK85",help = "get a set of "
            "parsimonious intervals each containing at least 1 recombination"
            " event.  Follows appendix 2 of Hudson and Kaplan 1985")
    intervaltype_group.add_argument("--4gcompat", dest="intervaltype",
            action='store_const', const="CONTIG4GPASS",help = "get a set of "
            " intervals such that for each interval all the included bases are"
            " compatible with each other based on the 4 gamete test")
    returntype_group = parser.add_mutually_exclusive_group(required=True)
    returntype_group.add_argument("--reti", dest="returntype",
            action='store_const', const="returninterval",help = " return single "
            " interval")
    returntype_group.add_argument("--retl", dest="returntype",
            action='store_const', const="returnlist", help = " return a list of"
            " all intervals")
    returntype_group = parser.add_mutually_exclusive_group(required=False)
    returntype_group.add_argument("--rani", dest="intervalpicktype",
            action='store_const', const="randominterval",default =
            "randominterval", help = "return a randomly selected interval"
            " from the list of all intervals")
    returntype_group.add_argument("--ranb", dest="intervalpicktype",
            action='store_const', const="randombase", default="randominterval",
            help = "return a single randomly selected interval with probabiliy "
            "of selection proportional to interval length")
    returntype_group.add_argument("--left", dest="intervalpicktype",
            action='store_const', const="leftinterval", default="randominterval",
            help = "return the leftmost interval")
    returntype_group.add_argument("--right", dest="intervalpicktype",
            action='store_const', const="rightinterval",
            default="randominterval",help = "return the rightmost interval")
    parser.add_argument("--numinf",dest= "numinformativesites",
            action="store",nargs = 1,type = int, default = None,help =
            "limit the picking of an interval to only those with at least this"
            " many informative snps")
    parser.add_argument("--b2",dest="only2basesnps",action="store_true",
            default=False,help = "causes snps with more than 2 bases "
            "to be ignored, default is to include them")
    parser.add_argument("--incN", dest="includesnpmissingdata",
            action="store_true", default=False,help="causes snps with missing"
            " data to be included, default is to ignore them")
    parser.add_argument("--ranseed", dest="ranseed", type=int,
            action="store",nargs = 1,default = None, help="positive integer to"
            " serve as the random number seed "
            "for use with --ranb or --rani, not required but"
            " useful when debugging")
    return parser

def logArgs(args):
    logging.info('Arguments for vcf_region_write:')
    for k in vars(args):
        logging.info('Argument %s: %s' % (k, vars(args)[k]))

def compat(s1, s2):
    """
    s1 and s2 are lists of sets of line numbers, each representing a particular
    polymorphic site.
    Each list as 4 sets,  one for each base: A,C,G & T.
    Each set is all the line  numbers that share that particular base.
    So then a polymorphic site represented by s1 will have
    two or more nonempty sets in s1.
    Same for s2

    compat() counts the number of pairs of sets, betweeen s1 and s2,
    for all possible pairs, for which the size of the intersection is nonzero
    and less than the size of either set. If this is 4 or more then the
    two sites are not compatible

    this should work for 2,3 or 4 polymorphic bases in either site

    """
    count = 0
    for i1 in range(len(s1)):
        for i2 in range(len(s2)):
            if len(s1[i1]) > 1 and len(s2[i2]) > 1:
                intersect = s1[i1]  & s2[i2]
                lenintersect = len(intersect)
                if (lenintersect > 0 and lenintersect < len(s1[i1])
                    and lenintersect < len(s2[i2])):
                    count = count + 1
    if count >= 4 :
        return False
    else :
        return True

def listcheckallsamelength(lcheck):
    """
        takes a list of lists or strings and checks to see that
        they are all the same length returns the length that is shared by all
    """
    i = 0
    for x in lcheck:
        try:
            lval = len(x)
        except TypeError:
            print('position %d not compatible with len()'%i)
        if i>0:
            assert lval == lastlval
        lastlval = lval
        i+=1
    return lval

def getseqlist(filename,format):
    """
        format can be "FASTA" or "SIT" or "VCF"
        returns a list of strings,  all should be of same length
        ONLYSNPS :
            true if data is just snps (e.g. vcf files)
            false if data is full sequences
    """
    poslist = [] # only used for vcfs
    seqs = []
    if format.upper() == "FASTA":
        ONLYSNPS = False
        f = open(filename,"r")
        seq = ""
        while True:
            s = f.readline().strip()
            if len(s) == 0:
                if len(seq) > 0:
                    seqs.append(seq)
                break
            else:
                if s[0] == '>':
                    if len(seq) > 0:
                        seqs.append(seq)
                    seq = ""
                else:
                    s = s.strip()
                    if len(s) > 0:
                        seq += s
    if format.upper() == "SIT":
        ONLYSNPS = False
        f = open(filename,"r")
        f.readline()
        while True:
            s = f.readline()
            if s[0] != '#':
                break
        numseq = int(s.split()[0])
        seqlen = int(s.split()[1])
        s = f.readline().split()

        numnoncode = int(s[1])
        if numnoncode != seqlen:
            for i in range(numnoncode):
                f.readline()
        s = f.readline()
        if s.find("All") < 0:
            numpops = int(s.strip())
            for i in range(numpops):
                f.readline()
        ## now should be at beginning of data
        for i in range(numseq):
            seqs.append(f.readline().strip()[10:])

    if format.upper() == "VCF":
        #print (filename)
        ONLYSNPS = True
        rg = Region(0,6000000,"11")
        vcff = pysam.VariantFile(filename)
        ploidy, seqs, poslist = vcfload.checkVcfRegionPysam(vcff,rg)
        #print (len(seqs), len(seqs[0]))
        #phased,ploidy,seqs,poslist = vcfload.vcfToList(filename)
        #if not phased:
            #raise  Exception(("vcf file %s appears to be unphased.  Four gamete-based intervals require phased sequences"%filename))
    numbases = listcheckallsamelength(seqs)
    return len(seqs),numbases,seqs,poslist,ONLYSNPS

def getSeqListFasta()

def buildlistsofsites(filename,format,skipsiteswithmissingdata,skipsiteswith3or4bases) :
    """
        takes a sequence filename and format
        returns:
             the number of sequences
             the length of sequences
             a list of base positions of variable sites
             a list of base positions of informative sites
             a list of pairs of sites that are incompatible by the 4 gamete criterion
             an upper diagonal matrix for all pairs of informative sites, contains a 1 if they are compatible,  0 if not
    """
    MAXINFORMSNPS = 1000
    seqcount,numbases,seqs,list_of_positions,ONLYSNPS = getseqlist(filename,format)
    # Find polymorphic sites
    polysites = []
    informpolysites = []
    informpolyseqs = ['']*seqcount ## becomes a list of strings of informative sites,  same length as seqs,  not used as of 8/14/2017
    baseset = set(['a','c','g','t','A','C','G','T'])
    linesets = []
    for j in range(numbases):
        ok = True
        basecount = [0,0,0,0]
        for k in range(seqcount):
            if seqs[k][j].isalpha() == False:
                raise Exception ("non-alphabetic character %s"
                    " at position %d in sequence %d"%(seqs[k][j],k,j))
            if seqs[k][j] not in baseset: ## skip sites with missing data
                if skipsiteswithmissingdata == True:
                    ok = False
            if seqs[k][j].upper() == 'A' :
                basecount[0] = basecount[0] + 1
            if seqs[k][j].upper() == 'C' :
                basecount[1] = basecount[1] + 1
            if seqs[k][j].upper() == 'G' :
                basecount[2] = basecount[2] + 1
            if seqs[k][j].upper() == 'T' :
                basecount[3] = basecount[3] + 1
        if ok == False:
            pass
            # print ("missing data at position " + str(j))
        else:
            polybase = 0
            twobase = 0
            for bi in range(4):
                if basecount[bi] > 0 :
                    polybase += 1
                    if basecount[bi] > 1 :
                        twobase +=  + 1
            if  skipsiteswith3or4bases == True and polybase > 2 :
                pass
                # print ("more than two bases at position " + str(j))
            elif polybase > 1:
                polysites.append(j+1) ## use regular (start at 1) counting
                if (twobase > 1):
                    informpolysites.append(j+1)
                    for k in range(seqcount):
                        informpolyseqs[k] = informpolyseqs[k] + seqs[k][j]
                    lineset =[set([]),set([]),set([]),set([])]
                    for k in range(seqcount):
                        if seqs[k][j].upper() == 'A' :
                            lineset[0].add(k)
                        if seqs[k][j].upper() == 'C' :
                            lineset[1].add(k)
                        if seqs[k][j].upper() == 'G' :
                            lineset[2].add(k)
                        if seqs[k][j].upper() == 'T' :
                            lineset[3].add(k)
                    linesets.append(lineset)
    if len(informpolysites) > MAXINFORMSNPS:
        raise Exception("the number of informative snps: %d"
                " exceeds the maximum: %d"%(len(informpolysites,MAXINFORMSNPS)))
    c = 0
    sl = len(linesets)
    if sl > 1 :
        cmat = []
        cintervals = []
        for i in range(sl):
            cmat.append([])
            for j in range(sl):
                if i < j :
                    if compat(linesets[i],linesets[j]) == False :
                        cmat[i].append(0)
                        cintervals.append([informpolysites[i],informpolysites[j]])
                    else:
                        cmat[i].append(1)
                else:
                    if i==j:
                        cmat[i].append(1)
                    else:
                        cmat[i].append(-1)
    return (seqcount,numbases,polysites,informpolysites,cintervals,
        cmat,list_of_positions,ONLYSNPS)

def replace_positions(intervals,list_of_positions):
    """
        with vcf files the base number of snp locations is contained
        in list_of_positions and the value in intervals[] are simply
        the index positions of the snps in the vcf files.

        This function rewrites the intervals to reflect the original
        snp base positions
    """
    bintervals = []
    for iv in intervals:
        bintervals.append([list_of_positions[iv[0]],list_of_positions[iv[0]]])
    return bintervals

def hudsonkaplan85intervals(cintervals,numbases,list_of_positions) :
    """
        implements the HK algorithm from appendix 2 of Hudson & Kaplan 1985

        Rm is the number of recombination events that can be
        parsimoniously inferred from a sample of sequences
        D is a compatibility matrix d[i,j] == 1 if all 4 gametes are present,
             0 otherwise
        consider for convenience only i<j i.e. above the diagonal

        H&K algorithm prunes D to leave the minimum set of intervals
        necessary to explain the data

        pruning steps:
            for two intervals (i,j) and (m,n)  if m <= i< j <= n  remove (m,n)

            for remaining intervals do these steps to remove
            non disjoint intervals
                identify an interval that is not disjoint from all
                the others (if none can be found we are done)
                call this (i1,j1)

                now scan all intervals (m,n) and delete
                all those for which i1 < m < j2

                now search other intervals to find (i2,j2) such that i2 >= j1
                and this interval is not disjoint from all the
                remaining intervals
                now using (i2,j2) delete all intervals
                whose first component is <j2 and >i2
    """
    c = len(cintervals)
    removelist = []
    for ai in range(c-1) :
        for bi in range(ai+1,c) :
            assert (cintervals[ai][0] < cintervals[ai][1])
            assert (cintervals[bi][0] < cintervals[bi][1])
            i = cintervals[ai][0]
            j = cintervals[ai][1]
            m = cintervals[bi][0]
            n = cintervals[bi][1]
            if (m <= i)  and (j <= n):
                if cintervals[bi] not in removelist:
                    removelist.append(cintervals[bi])
            if (i <= m)  and (n <= j):
                if cintervals[ai] not in removelist:
                    removelist.append(cintervals[ai])
    for interval in removelist:
        cintervals.remove(interval)
    c = len(cintervals)
    ai = 0
    while ai < c-1:
        removelist = []
        for bi in range(ai+1,c):
            i1 = cintervals[ai][0]
            j1 = cintervals[ai][1]
            b = cintervals[bi][0]
            e = cintervals[bi][1]
            if  (i1 < b  < j1) or (b < i1 < e):
                removelist.append(cintervals[bi])
        for interval in removelist:
            cintervals.remove(interval)

        i1 = cintervals[ai][0]
        j1 = cintervals[ai][1]
        xi = ai + 1
        c = len(cintervals)
        while xi < c:
            if cintervals[xi][0] >= j1:
                removelist = []
                for bi in range(xi+1,c):
                    i2 = cintervals[xi][0]
                    j2 = cintervals[xi][1]
                    b = cintervals[bi][0]
                    e = cintervals[bi][1]
                    if  (i2 < b  < j2):
                        removelist.append(cintervals[bi])
            for interval in removelist:
                cintervals.remove(interval)
                c = len(cintervals)
            xi += 1
        c = len(cintervals)
        ai += 1
    # reset snp positions if needed
    if len(list_of_positions) > 0:
        cintervals = replace_positions(cintervals,list_of_positions)
    return cintervals

def compatibleintervals(cmat,informpolysites,numbases,
        list_of_positions,ONLYSNPS):

    """
     return a list of intervals
        CONTIG4GPASS, set of most inclusive intervals, such that all the sites
            in each are fully compatible with each other:
            each interval is a stretch of bases that contain zero incompatible
            sites by the 4 gamete criterion
            on each side,  the left and right,  the next polymorphic site that
            is not in the interval must be incompatible
            with one of the polymorphic sites that is in the interval.

            intervals are based only on informative sites
            it is possible for intervals could be extended to include flanking
            non-informative polymorphic sites
            probably would have little effect in most contexts
    """
    tempinformpolysites = list(informpolysites)
    cints = []
    nps = len(informpolysites)
    cintinformlistfull = []
    p = 0
    while p < nps:
        j = p
        donej = False
        cintinform = 1
        while True:
            j += 1
            if j >= nps or cmat[p][j] != 1:
                j -= 1
                donej = True
                break
            else:
                for i in range(p+1,j+1):
                    if cmat[i][j] != 1:
                        donej = True
                        j -= 1
                        break
            if donej:
                break
            else:
                cintinform += 1
        # if j > p:
        cints.append([tempinformpolysites[p],tempinformpolysites[j]])
        cintinformlistfull.append(cintinform)
        p += 1
    ## now remove all intervals that are overlapped by a larger one
    cintinformremovelist = []
    c = len(cints)
    removelist = []
    for ai in range(c-1) :
        for bi in range(ai+1,c) :
            assert (cints[ai][0] <= cints[ai][1])
            assert (cints[bi][0] <= cints[bi][1])
            i = cints[ai][0]
            j = cints[ai][1]
            m = cints[bi][0]
            n = cints[bi][1]
            if (m <= i)  and (j <= n):
                if cints[bi] not in removelist:
                    removelist.append(cints[ai])
                    cintinformremovelist.append(ai)
            if (i <= m)  and (n <= j):
                if cints[ai] not in removelist:
                    removelist.append(cints[bi])
                    cintinformremovelist.append(bi)
    for interval in removelist:
        cints.remove(interval)
    cintinformlist = []
    for i in range(c):
        if i not in cintinformremovelist:
            cintinformlist.append(cintinformlistfull[i])
    # now pad intervals with sequences to the flanking ones
    if ONLYSNPS:
        cints = replace_positions(cints,list_of_positions)
    cintspadded = []
    # if ONLYSNPS don't pad using 1 and numbases
    for ci in range(len(cints)):
        temp = [-1,-1]
        if ci == 0:
            temp[0] = cints[ci][0] if ONLYSNPS else 1
        else:
            temp[0] = cints[ci-1][1] + 1
        if ci == len(cints)-1:
            temp[1] = cints[ci][1] if ONLYSNPS else numbases
        else:
            temp[1] = cints[ci+1][0] - 1
        cintspadded.append(temp)
    return cintspadded,cintinformlist

def sampleinterval(picktype,numinf,intervals,infcounts):
    """
        picktype : randominterval, randombase,left,right
        numinf : integer or None
        intervals : list of intervals
        infcounts: list of informative site counts,  or None
    """
    if len(intervals) == 0:
        return None
    if picktype == "leftinterval":
        return intervals[0]
    if picktype == "rightinterval":
        return intervals[-1]
    if numinf != None:
        numinf = numinf[0]
    if numinf != None and infcounts != None:
        okints = [i for i in range(len(infcounts)) if infcounts[i] >= numinf ]
    else:
        okints = [i for i in range(len(intervals))]
    n = len(okints)
    if n == 0: # no intervals with that many informative sites
        return None
    if picktype == "randominterval":
        return  intervals[ okints[ random.randint(0,n-1)]]
    if picktype == "randombase":
        sumbase = 0
        for i in range(n):
            sumbase += intervals[okints[i]][1] -  intervals[okints[i]][0] + 1
        rbn = random.randint(1,sumbase)
        i = 0
        while True:
            if rbn > intervals[okints[i]][1] -  intervals[okints[i]][0] + 1:
                rbn -= intervals[okints[i]][1] -  intervals[okints[i]][0] + 1
                i += 1
            else:
                break
        return  intervals[ okints[i]]
    assert False  # should not get here
    return None

def getIntervalList(args, format, filename, region=None):
    basedata = BaseData(format, filename, region)
    (seqcount, numbases, polysites, informpolysites, cintervals, cmat,
        list_of_positions,ONLYSNPS) = buildlistsofsites(args.filename,
        args.inputfiletype,args.includesnpmissingdata,
        args.only2basesnps)

    if args.intervaltype == "HandK85":
        intervals=hudsonkaplan85intervals(cintervals,
            numbases,list_of_positions)
        #hk code doesn't return counts of the informatics sites in each interval
        if args.numinformativesites != None:
            raise Exception("--numinf not implemented when using --hk")
        intervalinformativecounts = None # not implemented  for hk
        numinf = None # set to None, hk use of # inform sites not implemented

    else: # argdic['intervaltype'] =='CONTIG4GPASS'
        intervals,intervalinformativecounts= compatibleintervals(cmat,
            informpolysites, numbases,list_of_positions,ONLYSNPS)

    # print (argdic)
    if args.returntype == 'returnlist':
        return intervals
    else:
        interval = sampleinterval(args.intervalpicktype,
                                  args.numinformativesites,
                                  intervals,intervalinformativecounts)
        return interval


def parse_input(args):
    in_type = ""


def sample_fourgametetest_intervals(sys_args):
    """Returns interval(s) from aligned sequences or snps.

    Given a set of aligned sequences or a vcf file, this function
    returns intervals that have been identified based on four-gamete
    incombpatibility.   Intervals can either flank a recombination event or
    be sure not to contain an incombatibility.

    usage:  [-h] (--fasta | --sit | --vcf) (--hk | --4gcompat) (--reti | --retl)
        [--rani | --ranb | --left | --right] [--numinf NUMINFORMATIVESITES]
        [--b2] [--incN]
        filename

    Parameters
    ----------

    positional arguments:
    filename              Input filename

    optional arguments:
    -h, --help            show this help message and exit
    --fasta               fasta file of multiple aligned sequences, all the same
                            length
    --sit                 SITES format file
    --vcf                 vcf format file
    --hk                  get a set of parsimonious intervals each containing at
                        least 1 recombination event. Follows appendix 2 of
                        Hudson and Kaplan 1985
    --4gcompat            get a set of intervals such that for each interval all
                        the included bases are compatible with each other
                        based on the 4 gamete test
    --reti                return a single interval
    --retl                return a list of all intervals
    --rani                return a randomly selected interval from the list of
                        all intervals
    --ranb                return a single randomly selected interval with
                        probabiliy of selection proportional to interval
                        length
    --left                return the leftmost interval
    --right               return the rightmost interval
    --numinf NUMINFORMATIVESITES
                        limit the picking of an interval to only those with at
                        least this many informative snps
    --b2                  causes snps with more than 2 bases to be ignored,
                        default is to include them
    --incN                causes snps with missing data to be included, default
                        is to ignore them
    --ranseed RANSEED   positive integer to serve as the random number seed
                        for use with --ranb or --rani, not required but useful
                        when debugging

    """
    parser = createParser()
    # parser.print_help()
    if len(sys_args) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(sys_args)
    logArgs(args)
    #argdic = vars(args)
    if args.ranseed != None:
        random.seed(int(args.ranseed[0]))
    if args.vcf is not None:
        vcfname = args.vcf[0]
        regname = args.vcf[1]
        region_list = RegionList(filename=regname)
        for region in region_list:

    intervals = getIntervalList(args)
    return intervals
        #


if __name__ == '__main__':
    initLogger()

    ### some testing command lines
    #sample_fourgametetest_intervals([])
    #exit()
    testlist = []
    vcf_name = "chr11subsamples4gtest.vcf.gz"
    testlist.append( ["test #" + str(len(testlist)),vcf_name, "--vcf", "--4gcompat", "--reti", "--ranseed", "123"])
    testlist.append( ["test #" + str(len(testlist)),vcf_name, "--vcf", "--hk", "--retl"])
    testlist.append( ["test #" + str(len(testlist)),vcf_name, "--vcf", "--hk",
            "--reti","--left"])
    testlist.append( ["test #" + str(len(testlist)),vcf_name, "--vcf", "--4gcompat", "--retl"])
    testlist.append( ["test #" + str(len(testlist)),vcf_name, "--vcf", "--4gcompat",
            "--reti","--ranseed","123","--ranb"])
    testlist.append( ["test #" + str(len(testlist)),"fourgfastatest.fasta", "--fasta", "--4gcompat",
            "--retl"])
    testlist.append( ["test #" + str(len(testlist)),"fourgfastatest.fasta", "--fasta", "--4gcompat",
            "--retl","--incN","--b2"])
    testlist.append( ["test #" + str(len(testlist)),"fourgfastatest.fasta", "--fasta", "--hk",
            "--retl"])
    testlist.append( ["test #" + str(len(testlist)),"fourgfastatest.fasta", "--fasta", "--hk",
            "--retl","--incN","--b2"])
    testlist.append( ["test #" + str(len(testlist)),"fourgfastatest.fasta", "--fasta", "--hk",
            "--reti","--incN","--left"])
    testlist.append( ["test #" + str(len(testlist)),"fourgfastatest.fasta", "--fasta", "--4gcompat",
            "--reti","--incN","--numinf","2","--ranb"])
    testlist.append( ["test #" + str(len(testlist)),"fourgfastatest.fasta", "--fasta", "--4gcompat",
            "--reti","--incN","--numinf","2","--rani"])
    testlist.append( ["test #" + str(len(testlist)),"fourgfastatest.fasta", "--fasta", "--4gcompat",
            "--reti","--incN","--numinf","4","--ranb"])
    testlist.append( ["test #" + str(len(testlist)),"fourgfastatest.fasta", "--fasta", "--4gcompat",
            "--reti","--incN","--numinf","100","--ranb"])
    print ("%d test jobs"%(len(testlist))                   )
    for tl in testlist:
        sys.argv = tl
        a = sample_fourgametetest_intervals(sys.argv[1:])
        print(sys.argv[0],a)
    exit()
