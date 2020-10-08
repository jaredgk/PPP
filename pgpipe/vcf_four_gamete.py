#!/usr/bin/env python
''' 
    The four-gamete test is a method for determining whether or not 
    there has been recombination between a pair of variants. To do this, 
    all individuals must have haplotypes defined as the variants at the 
    two sites.

    .. image:: ../../PPP_assets/PPP_FGT.png
        :width: 100 %
        :align: center

    In this illustration of four-gamete test, the haplotypes of the samples
    from 197337 to 199256 (highlighted in green) pass the four-gamete test. 
    In comparison, the haplotypes from 196944 to 197337 and from 199256 to 
    199492 (highlighted in red) both fail the four-gamete test as all possible 
    haplotypes are observed.

    Given phased input with individual variants over a region of the genome, 
    four_gamete generates an interval within those variants that passes the 
    four-gamete filtering criteria, then return either that interval or an 
    output file with variants in that interval.

    Common usage for this function is to input a VCF file that contains
    variants for individuals at a single locus, with output returned being
    a VCF that contains a subsample of these variants. A full VCF can be 
    used with --vcfreg, where the second argument is a BED file with one
    or more regions, output will be either a VCF for four-gamete passing 
    regions or a new BED file with the truncated regions.

    ###############
    Input Arguments
    ###############
    **--vcfs** *<input_vcf_1>*...*<input_vcf_n>*
        Input name of one or more VCF files, where each VCF represents
        a locus. 
    **--vcfreg** *<input_vcf>* *<BED file>*
        Input name of VCF file containing genome data and name of BED file
        with regions to be analyzed.

    ###############
    Output Aguments
    ###############
    **--out** *<output_filename>*
        Name for output file.
    **--out-prefix** *<ouput_prefix>*
        If multiple files are output, this option is required to set a 
        prefix for the output files.

    ##################
    Interval Arguments
    ##################
    **--numinf** *<minimum informative site count>*
        Region returned must have at least n informative sites, defaults to 1
    **--hk**
        If set, returns intervals with at least one recombination event 
        instead of regions with no recombination.
    **--reti**
        This script will generate a list of valid regions with no
        recombination. Selecting this option will return a single interval
        as specified by other arguments
    **--retl**
        Returns all valid intervals, either as a list of intervals or 
        multiple output files


    ################################
    Single Returned Region Arguments
    ################################
    Select one of:
    **--rani**
        Returns random interval (default)
    **--ranb**
        Returns random interval, with probability of interval 
        proportional to interval length
    **--left**
        Return first interval with enough informative sites
    **--right**
        Return last interval with enough informative sites
    **--maxlen**
        Return interval with most informative sites

    ###############
    Other Arguments
    ###############
    **--remove-multiallele**
        Removes multi-alleleic sites from analysis
    **--include-missing**
        Include sites with missing data in analysis
    **--ovlps**
        Extend region to include non-informative variants between
        an edge variant and a variant that breaks the four-gamete
        criteria
    **--ovlpi**
        Include informative variants from overlapping
        regions
'''

import sys
import argparse
import logging
import random
import pysam

from pgpipe.logging_module import initLogger
from pgpipe.genome_region import Region, RegionList
from pgpipe.vcf_reader_func import getRecordList, vcfRegionName, getRecordsInRegion, VcfReader
from pgpipe.misc import argprase_kwargs


class BaseData():

    def __init__(self, args, format, filename, region=None):
        self.poslist = []
        self.seqs = []
        self.format = format
        self.header = None
        self.records = []
        self.extend_noninf = args.extend_noninf
        self.extend_inf = args.extend_inf
        if format == 'VCF':
            self.onlysnps = True
            vcff = VcfReader(filename,index=args.tabix_index)
            self.allrecords = vcff.getRecordList(region=region, chrom=args.chrom)
            self.checkVcfRegion()
            self.header = vcff.reader.header
            vcff.close()
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
            f.close()
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
            f.close()
        self.numbases = listcheckallsamelength(self.seqs)
        logging.info('Sequence length: %d' % self.numbases)
        self.seqcount = len(self.seqs)
        logging.info('Sequence count: %d' % self.seqcount)
        self.buildlistofsites(args.includesnpmissingdata, args.only2basesnps)

    def buildlistofsites(self, includesnpmissingdata, only2basesnps):
        MAXINFORMSNPS = 1000
        # Find polymorphic sites
        self.polysites = []
        self.informpolysites = []
        informpolyseqs = ['']*self.seqcount ## becomes a list of strings of informative sites,  same length as seqs,  not used as of 8/14/2017
        baseset = set(['a','c','g','t','A','C','G','T'])
        linesets = []
        for j in range(self.numbases):
            ok = True
            basecount = [0,0,0,0]
            for k in range(self.seqcount):
                if self.seqs[k][j].isalpha() == False:
                    raise Exception ("non-alphabetic character %s"
                        " at position %d in sequence %d"%(self.seqs[k][j],k,j))
                if self.seqs[k][j] not in baseset: ## skip sites with missing data
                    if not includesnpmissingdata:
                        ok = False
                if self.seqs[k][j].upper() == 'A' :
                    basecount[0] = basecount[0] + 1
                if self.seqs[k][j].upper() == 'C' :
                    basecount[1] = basecount[1] + 1
                if self.seqs[k][j].upper() == 'G' :
                    basecount[2] = basecount[2] + 1
                if self.seqs[k][j].upper() == 'T' :
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
                if  only2basesnps == True and polybase > 2 :
                    pass
                    # print ("more than two bases at position " + str(j))
                elif polybase > 1:
                    self.polysites.append(j+1) ## use regular (start at 1) count
                    if (twobase > 1):
                        self.informpolysites.append(j+1)
                        lineset =[set([]),set([]),set([]),set([])]
                        for k in range(self.seqcount):
                            if self.seqs[k][j].upper() == 'A' :
                                lineset[0].add(k)
                            if self.seqs[k][j].upper() == 'C' :
                                lineset[1].add(k)
                            if self.seqs[k][j].upper() == 'G' :
                                lineset[2].add(k)
                            if self.seqs[k][j].upper() == 'T' :
                                lineset[3].add(k)
                        linesets.append(lineset)
                        if len(self.informpolysites) >= MAXINFORMSNPS:
                            logging.warning('Maximum SNP count of %d has been hit'%(MAXINFORMSNPS))
                            break
        logging.info('Informative site count: %d' % len(self.informpolysites))
        c = 0
        sl = len(linesets)
        if sl > 1 :
            self.cmat = []
            self.cintervals = []
            for i in range(sl):
                self.cmat.append([])
                for j in range(sl):
                    if i < j :
                        if self.records[i].chrom != self.records[j].chrom or compat(linesets[i],linesets[j]) == False :
                            self.cmat[i].append(0)
                            self.cintervals.append([self.informpolysites[i],self.informpolysites[j]])
                        else:
                            self.cmat[i].append(1)
                    else:
                        if i==j:
                            self.cmat[i].append(1)
                        else:
                            self.cmat[i].append(-1)


    def compatibleintervals(self):
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
        tempinformpolysites = list(self.informpolysites)
        cints = []
        nps = len(self.informpolysites)
        cintinformlistfull = []
        p = 0
        while p < nps:
            j = p
            donej = False
            cintinform = 1
            while True:
                j += 1
                if j >= nps or self.cmat[p][j] != 1:
                    j -= 1
                    donej = True
                    break
                else:
                    for i in range(p+1,j+1):
                        if self.cmat[i][j] != 1:
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
        cintspadded = []
        if False:
            cints = replace_positions(cints,self.poslist)
            return cints,cintinformlist
        if self.onlysnps:
        ## using a vcf file,  intervals are 1-based  but positions in list_of_positions are 0-based
            if not self.extend_inf:
                for ci in range(len(cints)):
                    tempcint = [-1,-1]
                    if ci == 0:
                        tempcint[0] = cints[ci][0]
                    else:
                        tempcint[0] = self.informpolysites[self.informpolysites.index(cints[ci-1][1])+1]
                    tempcint[1] = cints[ci][1]
                    cints[ci] = tempcint
                    cintinformlist[ci] = self.informpolysites.index(cints[ci][1]) - self.informpolysites.index(cints[ci][0])+1
            for ci in range(len(cints)):
                tempcint = [-1,-1]
                if not self.extend_noninf:
                    tempcint[0] = self.poslist[cints[ci][0]-1]
                    if ci == len(cints)-1:
                        tempcint[1] = self.poslist[self.informpolysites[-1]]-1
                    else:
                        tempcint[1] = self.poslist[self.informpolysites[self.informpolysites.index(cints[ci][1])+1]-1]-1
                else:
                    if ci == 0:
                        assert cints[ci][0] > 0
                        tempcint[0] = self.poslist[cints[ci][0] - 1]  ## -1 because list_of_positions is 0-based
                    else:
                        tempcint[0] = self.poslist[self.informpolysites[self.informpolysites.index(cints[ci][0]) -1] -1] + 1  ## position of previous informative snp + 1
                    if ci == len(cints)-1:
                        tempcint[1] = self.poslist[cints[ci][1] -1]
                    else:
                        tempcint[1] = self.poslist[self.informpolysites[self.informpolysites.index(cints[ci][1]) + 1] -1] -1   ## position of next informative snp -1
                cintspadded.append(tempcint)
        else:
            tempcint = [-1,-1]
            if not self.extend_inf:
                for ci in range(len(cints)):
                    if ci == 0:
                        tempcint[0] = cints[ci][0]
                    else:
                        tempcint[0] = self.informpolysites[self.informpolysites.index(cints[ci-1][1])+1]
                    tempcint[1] = cints[ci][1]
                    cints[ci] = tempcint
                    cintinformlist[ci] = self.informpolysites.index(cints[ci][1]) - self.informpolysites.index(cints[ci][0])+1
            else:
                for ci in range(len(cints)):
                    if ci == 0:
                        tempcint[0] = 1
                    else:
                        tempcint[0] = self.informpolysites[self.informpolysites.index(cints[ci][0]) -1] + 1   ## position of previous informative snp + 1
                    if ci == len(cints)-1:
                        tempcint[1] = self.numbases
                    else:
                        tempcint[1] = self.informpolysites[self.informpolysites.index(cints[ci][1]) + 1] -1    ## position of next informative snp -1
                    cintspadded.append(tempcint)
        return cintspadded,cintinformlist
        #return cints,cintinformlist

    #def hudsonkaplan85intervals(cintervals,numbases,list_of_positions) :
    def hk85intervals(self):
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
        c = len(self.cintervals)
        removelist = []
        for ai in range(c-1) :
            for bi in range(ai+1,c) :
                assert (self.cintervals[ai][0] < self.cintervals[ai][1])
                assert (self.cintervals[bi][0] < self.cintervals[bi][1])
                i = self.cintervals[ai][0]
                j = self.cintervals[ai][1]
                m = self.cintervals[bi][0]
                n = self.cintervals[bi][1]
                if (m <= i)  and (j <= n):
                    if self.cintervals[bi] not in removelist:
                        removelist.append(self.cintervals[bi])
                if (i <= m)  and (n <= j):
                    if self.cintervals[ai] not in removelist:
                        removelist.append(self.cintervals[ai])
        for interval in removelist:
            self.cintervals.remove(interval)
        c = len(self.cintervals)
        ai = 0
        while ai < c-1:
            removelist = []
            for bi in range(ai+1,c):
                i1 = self.cintervals[ai][0]
                j1 = self.cintervals[ai][1]
                b = self.cintervals[bi][0]
                e = self.cintervals[bi][1]
                if  (i1 < b  < j1) or (b < i1 < e):
                    removelist.append(self.cintervals[bi])
            for interval in removelist:
                self.cintervals.remove(interval)

            i1 = self.cintervals[ai][0]
            j1 = self.cintervals[ai][1]
            xi = ai + 1
            c = len(self.cintervals)
            while xi < c:
                if self.cintervals[xi][0] >= j1:
                    removelist = []
                    for bi in range(xi+1,c):
                        i2 = self.cintervals[xi][0]
                        j2 = self.cintervals[xi][1]
                        b = self.cintervals[bi][0]
                        e = self.cintervals[bi][1]
                        if  (i2 < b  < j2):
                            removelist.append(self.cintervals[bi])
                for interval in removelist:
                    self.cintervals.remove(interval)
                    c = len(self.cintervals)
                xi += 1
            c = len(self.cintervals)
            ai += 1
        # reset snp positions if needed
        if len(self.poslist) > 0:
            self.cintervals = replace_positions(self.cintervals,self.poslist)
        return self.cintervals

    def checkVcfRegion(self):
        """
        Checks that variants in region are SNPs and creates a list with valid sites. Also sets ploidy and position list of valid sites.

        """
        valid_bases = ['A','C','G','T']
        skiplines = []
        self.ploidy = -1
        self.poslist = []
        base_list = []
        self.records = []
        for record in self.allrecords:
            lineok=True
            for allele in record.alleles:
                if allele not in valid_bases:
                    lineok=False
                    break
            if not lineok:
                continue
            somedata = False
            #data_present = 0
            var_list = []
            #for indiv in record.samples:
            for i in range(len(record.samples)):
                indiv = record.samples[i]
                if self.ploidy == -1:
                    self.ploidy = len(indiv.alleles)
                if not somedata and indiv.alleles.count(None) == 0:
                    somedata = True
                for samp_allele in indiv.alleles:
                    if samp_allele is None:
                        var_list.append('N')
                    else:
                        var_list.append(samp_allele)
            if somedata:
                base_list.append(var_list)
                self.poslist.append(record.pos)
                self.records.append(record)
        self.seqs = [[base_list[row][col] for row in range(0, len(base_list))] for col in range(0, len(base_list[0]))]


#def createParser(passed_arguments = []):
def parseArguments(passed_arguments = []):
    parser = argparse.ArgumentParser(description=("Given a file of aligned"
            "sequences or variable sites, will return intervals based"
            " on 4 gamete tests"))
    #parser.add_argument("filename", help="Input filename")
    filetype_group = parser.add_mutually_exclusive_group(required=True)
    filetype_group.add_argument("--vcfreg", dest = "vcfreg", nargs=2,
            help = "VCF file and region list, in that order")
    filetype_group.add_argument("--vcfs", dest="vcfname", nargs="+",
                                help=("One or more input VCF files"))
    filetype_group.add_argument("--fasta", dest = "fasta",
            nargs="+", help = "fasta file(s) of"
            " multiple aligned sequences, all the same length")
    filetype_group.add_argument("--sit", dest = "sit",
            nargs="+", help = " SITES format file(s)")
    output_group = parser.add_mutually_exclusive_group()
    output_group.add_argument("--out", dest="out",help=("Output filename "
                            "for single input"))
    output_group.add_argument("--out-prefix", dest="out_prefix",help=("Prefix "
                            "for multi-file output"))
    intervaltype_group = parser.add_mutually_exclusive_group(required=
            "--reti" in sys.argv)
    intervaltype_group.add_argument("--hk", dest="intervaltype",
            action='store_const', const="HandK85",help = "get a set of "
            "parsimonious intervals each containing at least 1 recombination"
            " event.  Follows appendix 2 of Hudson and Kaplan 1985")
    intervaltype_group.add_argument("--fourgcompat", dest="intervaltype",
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
            "leftinterval", help = "return a randomly selected interval"
            " from the list of all intervals")
    returntype_group.add_argument("--ranb", dest="intervalpicktype",
            action='store_const', const="randombase", default="leftinterval",
            help = "return a single randomly selected interval with probabiliy "
            "of selection proportional to interval length")
    returntype_group.add_argument("--left", dest="intervalpicktype",
            action='store_const', const="leftinterval", default="leftinterval",
            help = "return the leftmost interval")
    returntype_group.add_argument("--right", dest="intervalpicktype",
            action='store_const', const="rightinterval",
            default="leftinterval",help = "return the rightmost interval")
    returntype_group.add_argument("--maxlen", dest="intervalpicktype",action="store_const", const="maxregion",default="leftinterval",help="return interval with most SNPs")
    parser.add_argument("--numinf",dest= "numinformativesites",
            action="store",nargs = 1,type = int, default = None,help =
            "limit the picking of an interval to only those with at least this"
            " many informative snps")
    parser.add_argument("--remove-multiallele",dest="only2basesnps",action="store_true",
            default=False,help = "causes snps with more than 2 bases "
            "to be ignored, default is to include them")
    parser.add_argument("--include-missing", dest="includesnpmissingdata",
            action="store_true", default=False,help="causes snps with missing"
            " data to be included, default is to ignore them")
    parser.add_argument("--ranseed", dest="ranseed", type=int,
            action="store",nargs = 1,default = None, help="positive integer to"
            " serve as the random number seed "
            "for use with --ranb or --rani, not required but"
            " useful when debugging")
    parser.add_argument("--chrom", dest="chrom", help="Select variants "
            "from a single specified chromosome")
    parser.add_argument("--ovlps",dest="extend_noninf",action="store_true",
                        help=("Extend region to overlapping non-informative sites"))
    parser.add_argument("--ovlpi",dest="extend_inf",action="store_true",
                        help="Extend region to overlapping informative sites")
    #add checks for correct input type
    parser.add_argument('--tbi', dest="tabix_index",help=("Filepath for tabix "
                        "index file if using a single bgzipped VCF"))
    parser.add_argument('--log', dest="log_name", help=("Filepath for log file"))
    #return parser
    if passed_arguments:
        return vars(parser.parse_args(passed_arguments))
    else:
        return vars(parser.parse_args())

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
        #logging.info('Index: %d %d' % (iv[0],iv[1]))
        #logging.info('Length: %d' % len(list_of_positions))
        bintervals.append([list_of_positions[iv[0]-1],list_of_positions[iv[1]-1]+1])
    return bintervals


def sampleinterval(picktype,numinf,intervals,infcounts):
    """
        picktype : randominterval, randombase,left,right
        numinf : integer or None
        intervals : list of intervals
        infcounts: list of informative site counts,  or None
    """
    if len(intervals) == 0:
        return None

    #if picktype == "leftinterval":
    #    return intervals[0]
    #if picktype == "rightinterval":
    #    return intervals[-1]
    if numinf != None:
        numinf = numinf[0]
    if numinf != None and infcounts != None:
        okints = [i for i in range(len(infcounts)) if infcounts[i] >= numinf ]
    else:
        okints = [i for i in range(len(intervals))]
    n = len(okints)
    if n == 0: # no intervals with that many informative sites
        return None
    if picktype == "leftinterval":
        return intervals[okints[0]]
    if picktype == "rightinterval":
        return intervals[okints[-1]]
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
    if picktype == "maxregion":
        #print (okints)
        #print (infcounts)
        #print (intervals)
        return intervals[infcounts.index(max(infcounts))]
    assert False  # should not get here
    return None

def getIntervalList(args, basedata):
    """
    Returns list of intervals from either hk85 or 4gcompat function.
    Also controls return based on args provided (list or single interval)
    """

    if args.intervaltype == "HandK85":
        #hk code doesn't return counts of the informatics sites in each interval
        intervals = basedata.hk85intervals()
        if args.numinformativesites != None:
            raise Exception("--numinf not implemented when using --hk")
        intervalinformativecounts = None # not implemented  for hk
        numinf = None # set to None, hk use of # inform sites not implemented

    else:
        intervals, intervalinformativecounts = basedata.compatibleintervals()

    if args.returntype == 'returnlist':
        return intervals
    else:
        interval = sampleinterval(args.intervalpicktype,
                                  args.numinformativesites,
                                  intervals,intervalinformativecounts)
        return interval


def outputSubregion(args, interval, basedata, region=None, filename=None):
    #create new region from interval
    #find records in range of new interval
    #open file with prefix and region name
    #write records, close file
    """
    Outputs records in identified sub-interval to desired output file
    """
    if interval is None:
        return
    if region is None:
        subregion = Region(interval[0]-1,interval[1],basedata.records[0].chrom)
    else:
        subregion = Region(interval[0]-1,interval[1],region.chrom)
    subrecords = getRecordsInRegion(subregion, basedata.records)
    if filename is None:
        subfn = vcfRegionName(args.out_prefix,subregion,"vcf.gz")
    else:
        subfn = filename
    #print (interval)
    #print (subfn)
    subf = pysam.VariantFile(subfn, 'w', header=basedata.header)
    for record in subrecords:
        subf.write(record)
    subf.close()





#def sample_fourgametetest_intervals(sys_args):
def sample_fourgametetest_intervals(**kwargs):
    """Returns interval(s) from aligned sequences or snps.

    Given a set of aligned sequences or a vcf file, this function
    returns intervals that have been identified based on four-gamete
    incombpatibility.   Intervals can either flank a recombination event or
    be sure not to contain an incombatibility.


    Parameters
    ----------

    Input (one required):
    --fasta : str
        Filename(s) of input FASTA file for testing
    --vcfreg : str, str
        Filename of input VCF file and BED file with list of
        regions to run tests on
    --vcfs : str
        Filename(s) of VCF files, each corresponding to a region, that will be subsampled for compatible intervals
    --sit : str
        Filename(s) of input SITES format files for testing

    Output types:
    --out : str
        Name of single output file.
    --out-prefix : str
        Prefix for name of output files
    --out-reg : str
        Name of file that stores passing interval per interval in input BED file (currently not working)

    Tests (one required)
    --hk : bool
        Will find intervals that contain at least one recombination event
    --fourgcompat : bool
        Will find intervals with zero recombination events in them

    Return options:
    --reti : bool
        Output will be of a single region, indicated by rani/ranb/left/right
    --retl : bool
        Output will be of all passing regions

    For single region:
    --left : bool (default)
        Single output region will be left-most region with enough informative sites
    --right : bool
        Output region will be right-most region
    --rani : bool
        Region will be randomly selected, with each region having equal chance of selection
    --ranb : bool
        Region will be randomly selected, with each region having selection odds proportional to size
    --maxlen : bool
        Region returned is that of one with most informative SNPs

    General options:
    --numinf : int
        Number of sites in an identified interval required for reporting (default is 1)
    --incN : bool
        Include SNPs with missing data (default is to ignore them)
    --ranseed : int
        Seed for picking random interval from list
    --tbi : str
        Path to tabix index for input VCF file if it's name doesn't match
    --ovlps : bool
        Extend region to overlapping non-informative sites
    --ovlpi : bool
        Extend region to overlapping informative sites

    """
    #parser = createParser()
    # parser.print_help()
    #if len(kwargs) == 0:
        #parser.print_help()
        #sys.exit(1)
    #args = parser.parse_args(sys_args)
    if __name__ != "__main__":
        kwargs = argprase_kwargs(kwargs, parseArguments)

    args = argparse.Namespace(**kwargs)
    if args.log_name is not None:
        initLogger(filename=log_name)
    logArgs(args)
    #argdic = vars(args)
    interval_list = []
    multiple_regions = False
    region_output = False
    if args.ranseed != None:
        random.seed(int(args.ranseed[0]))
    filetype = ''
    out_list = (args.out is None and args.out_prefix is None)

    if args.vcfreg is not None:
        filetype = "VCF"
        vcfname = args.vcfreg[0]
        regname = args.vcfreg[1]
        if regname != '-':
            region_list = RegionList(filename=regname)
            #region_output = True
            multiple_regions = (len(region_list.regions) > 1)
            if multiple_regions and not out_list and args.out is not None:
                raise Exception(("VCF with multiple regions cannot have "
                                "single output file (--out)"))
            for region in region_list.regions:
                basedata = BaseData(args, "VCF",vcfname, region)
                intervals = getIntervalList(args, basedata)
                interval_list.append(intervals)
                if args.out_prefix is not None:
                    if args.returntype == 'returnlist':
                        for iv in intervals:
                            outputSubregion(args, iv, basedata)
                    else:
                        outputSubregion(args, intervals, basedata)
                elif args.out is not None:
                    if args.returntype == 'returnlist':
                        raise Exception('Multiple outputs directed to single file')
                    outputSubregion(args, intervals, basedata, region=region, filename=args.out)

        else:
            basedata = BaseData(args, "VCF",vcfname)
            intervals = getIntervalList(args,basedata)
            interval_list.append(intervals)
            if args.out_prefix is not None:
                if args.returntype == 'returnlist':
                    for iv in intervals:
                        outputSubregion(args, iv, basedata)
                else:
                    outputSubregion(args, intervals, basedata)
            elif out_list:
                return interval_list
            elif args.out is not None:
                if args.returntype == 'returnlist':
                    raise Exception(('If --retl is specified, output must not '
                                    'be --out'))
                outputSubregion(args, intervals, basedata, filename=args.out)

    elif args.vcfname is not None:
        filetype="VCF"
        if len(args.vcfname) > 1 and args.out_prefix is not None:
            raise Exception(("Multiple VCFs require --out-prefix flag"))
        for vcfname in args.vcfname:
            basedata = BaseData(args, "VCF",vcfname)
            intervals = getIntervalList(args,basedata)
            #sys.stderr.write(str(intervals)+'\n')
            interval_list.append(intervals)
            if args.out_prefix is not None:
                if args.returntype == 'returnlist':
                    for iv in intervals:
                        outputSubregion(args, iv, basedata)
                else:
                    outputSubregion(args, intervals, basedata)
            elif out_list:
                return interval_list
            elif args.out is not None:
                if args.returntype == 'returnlist':
                    raise Exception(('If --retl is specified, output must not '
                                    'be --out'))
                outputSubregion(args, intervals, basedata, filename=args.out)
    elif args.fasta is not None:
        filetype="FASTA"
        if len(args.fasta) > 1 and args.out_prefix is not None:
            raise Exception(("Multiple FASTAs require --out-prefix flag"))
        for fastaname in args.fasta:
            basedata = BaseData(args, "FASTA", fastaname)
            intervals = getIntervalList(args, basedata)
            interval_list.append(intervals)
    elif args.sit is not None:
        filetype="SIT"
        if len(args.sit) > 1 and args.out_prefix is not None:
            raise Exception(("Multiple SITES files require --out-prefix flag"))
        if args.out_reg is not None:
            raise Exception(("SITES files cannot be output to region file"))
        for sitname in args.sit:
            intervals = getIntervalList(args, "SIT", sitname)
            interval_list.append(intervals)




    return interval_list
        #


if __name__ == '__main__':
    #initLogger()

    #sample_fourgametetest_intervals(sys.argv[1:])
    sample_fourgametetest_intervals(**parseArguments())

