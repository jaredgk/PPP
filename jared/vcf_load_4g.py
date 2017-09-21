import sys
import logging
from gene_region import Region
from vcf_reader_func import getRecordList

""" used with four_gamete_work.py, loads an uncompressed vcf
"""

def checkvcf(vcf_name):
    """
        do some basic checks of a vcf file
        assess ploidy (1 or 2)
        if ploidy is 2, assess whether or not file is phased
        count snps on lines not missing data
        build list of lines with missing data
        assess how much missing data
        return:
            phased
            ploidy
            nsnps
            skiplines
            propmissing
    """
    bases = ['A','C','G','T']
    skiplines = []
    ploidy = -1
    phased = False
    nsnps  = 0
    linenum = 0
    missing = 0
    notmissing = 0
    alleles = []
    f = open(vcf_name,'r')
    for line in f:
        if line[0] == '#':
            continue
        la = line.strip().split()
        lineok = True
        ref = la[3]
        if ref in bases:
            a = ['0']
            linealleles = [ref]
            alt = la[4].split(',')
            numalleles = len(alt) + 1
            for alta in alt:
                if alta.upper() not in bases:
                    skiplines.append(linenum)
                    lineok = False
                    break
                else:
                    a.append(str(len(a)))
                    linealleles.append(alta.upper())
                if lineok:
                    somedata = False
                    for i in range(9, len(la)):
                        g = la[i][0:3] # genotype part of the string
                        if g[0] in a:
                            notmissing += 1
                            somedata = True
                        if g[1] == '|':
                            phased = True
                            assert ploidy == 2 or ploidy == -1
                            ploidy = 2
                        elif g[1] != '/' and g[1] != ':':
                            # no '|' or '/'  so must be haploid
                            assert ploidy ==1 or ploidy == -1
                            ploidy = 1
                        elif g[1] == '/':
                            assert ploidy == 2 or ploidy == -1
                            ploidy = 2
                        if ploidy ==2 and g[2] in a:
                            notmissing += 1
                            somedata = True
                        missing += g.count('.')
                    if somedata == False:
                        skiplines.append(linenum)
                    else:
                        nsnps += 1
                        alleles.append(linealleles)
        else:
            skiplines.append(linenum)
        # print (linenum)
        linenum += 1
    propmissing =  missing/float(missing + notmissing)
    return phased,ploidy,nsnps,skiplines,propmissing,alleles

def getBaseList(ploidy, la):
    """
        The allele values are 0 for the reference
        allele (what is in the REF field),  1 for the first allele listed
        in ALT, 2 for the second allele list in ALT and
        so  on.
    """
    base_list = []
    base_list2 = []
    for i in range(9,len(la)):
        base_list.append(la[i][0])
        if ploidy==2:
            base_list2.append(la[i][2])
    return base_list,base_list2

def checkVcfRegionPysam(record_list):
    #record_list = getRecordList(vcf_reader, region)
    valid_bases = ['A','C','G','T']
    skiplines = []
    ploidy = -1
    pos_list = []
    base_list = []
    for record in record_list:
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
            if ploidy == -1:
                ploidy = len(indiv.alleles)
            if not somedata and indiv.alleles.count(None) == 0:
                somedata = True
            for samp_allele in indiv.alleles:
                if samp_allele is None:
                    var_list.append('.')
                else:
                    var_list.append(samp_allele)
        if somedata:
            base_list.append(var_list)
            pos_list.append(record.pos)
    transposed_bases = [[base_list[row][col] for row in range(0, len(base_list))] for col in range(0, len(base_list[0]))]
    return ploidy, transposed_bases, pos_list



def vcfToList(vcf_name):
    """Creates two lists: a 2d list with base information like a VCF file,
    and a list with position information from a given index.
    """
    phased, ploidy, nsnps, skiplines, propmissing,alleles = checkvcf(vcf_name)
    assert len(alleles) == nsnps
    f = open(vcf_name,'r')
    pos_list = []
    base_list = []
    li = 0
    for line in f:
        if line[0] == '#':
            continue
        if li not in skiplines:
            la = line.strip().split()
            pos = int(la[1])
            pos_list.append(pos)
            bl,bl2 = getBaseList(ploidy,la)
            if ploidy == 1:
                base_list.append(bl)
            else:
                assert len(bl) == len(bl2)
                mergedbl = []
                for i in range(len(bl)):
                    mergedbl.append(bl[i])
                    mergedbl.append(bl2[i])
                base_list.append(mergedbl)
        li += 1
    # transpose base_list
    # print (len(base_list),len(base_list[0]))
    # if sys.version_info > (2,):
    #     alignedsnparray = list(map(list, zip(*base_list)))
    # else:
    #     alignedsnparray = map(list, zip(*base_list))
    alignedsnparray = [[base_list[row][col] for row in range(0, len(base_list))] for col in range(0, len(base_list[0]))]
    numseqs = len(alignedsnparray)
    # print (numseqs, len(alignedsnparray[0]))
    # assert nsnps == len(alignedsnparray[0])
    for i in range(numseqs):
        for j in range(nsnps):
            if  alignedsnparray[i][j] == '.' :
                alignedsnparray[i][j] = 'N'
            else:
                alignedsnparray[i][j] = alleles[j][int(alignedsnparray[i][j])]
    return phased, ploidy, alignedsnparray, pos_list

if __name__ == '__main__':
    phased, ploidy, alignedsnparray, pos_list = vcfToList("merged_chr1_1000.vcf")
    numseqs = len(alignedsnparray)
    nsnps = len(alignedsnparray[0])
    print (phased,ploidy)

    print (numseqs,nsnps,len(pos_list))
    for i in range(numseqs):
        for j in range(nsnps):
            print (alignedsnparray[i][j],end = '')
        print ()
