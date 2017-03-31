import sys
import pysam
import logging
from random import sample



def getRecordList(vcf_reader, region, chrom):
    """Returns list for use in subsampling from input file"""
    var_sites = vcf_reader.fetch(chrom, region.start, region.end)
    lst = []
    for rec in var_sites:
        lst.append(rec)
    return lst


def getRecordListUnzipped(vcf_reader, region, chrom, prev_last_rec):
    lst = []
    if (prev_last_rec is not None and
        region.start <= prev_last_rec.pos < region.end):
        lst.append(prev_last_rec)
    elif prev_last_rec is not None and prev_last_rec.pos >= region.end:
        return []
    rec = next(vcf_reader,None)
    while rec is not None and rec.pos < region.end:
        if rec.pos >= region.start:
            lst.append(rec)
        rec = next(vcf_reader,None)
    prev_last_rec = rec
    return lst, prev_last_rec


def checkRecordIsSnp(rec):
    """Checks if this record is a single nucleotide variant, returns bool."""
    if len(rec.ref) != 1:
        return False
    for allele in rec.alts:
        if len(allele) != 1:
            return False
    return True


def getSubsampleList(vcfname, ss_count):
    vcf_o = pysam.VariantFile(vcfname)
    rec = next(vcf_o)
    vcf_o.close()
    lst = []
    for samp in rec.samples:
        lst.append(samp)
    return lst[:int(ss_count)]


def getVcfReader(args):
    subsamp_list = None
    if args.subsamp_num is not None:
        subsamp_list = getSubsampleList(args.vcfname, args.subsamp_num)
    elif args.subsamp_fn is not None:
        subsamp_file = open(args.subsamp_fn,'r')
        subsamp_list = [l.strip() for l in subsamp_file.readlines()]
        subsamp_file.close()
    vcf_reader = pysam.VariantFile(args.vcfname)
    if subsamp_list is not None:
        logging.debug('Subsampling %d individuals from VCF file' %
        (len(subsamp_list)))
        vcf_reader.subset_samples(subsamp_list)
    return vcf_reader
