import sys
import logging

def getBaseList(la):
    base_list = []
    for i in range(9,len(la)):
        base_list.append(la[i][0])
        base_list.append(la[i][2])
    return base_list

def vcfToList(vcf_name):
    """Creates two lists: a 2d list with base information like a VCF file,
    and a list with position information from a given index.


    """
    f = open(vcf_name,'r')
    pos_list = []
    base_list = []
    for line in f:
        if line[0] == '#':
            continue
        la = line.strip().split()
        pos = int(la[1])
        pos_list.append(pos)
        bl = getBaseList(la)
        base_list.append(bl)


    return base_list, pos_list
