import sys
import pysam
import logging
import struct
from random import sample
from collections import defaultdict
import os
import gzip

def checkIfGzip(filename):
    try:
        gf = gzip.open(filename)
        gl = gf.readline()
        gf.close()
        vcf_check = b'##fileformat=VCF'
        if gl[0:3] == b'BCF':
            return 'bcf'
        elif gl[:len(vcf_check)] == vcf_check:
            return checkHeader(filename)
        else:
            return 'other'
    except:
        return 'nozip'

def checkHeader(filename):
    f = open(filename,'rb')
    l = f.readline()
    f.close()
    BGZF_HEADER=b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00'
    #BGZF_HEADER=b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff'
    GZF_HEADER=b'\x1f\x8b'
    if l[:len(BGZF_HEADER)] == BGZF_HEADER:
        return 'bgzip'
    if l[:len(GZF_HEADER)] == GZF_HEADER:
        return 'gzip'
    return 'nozip'

def checkFormat(vcfname):
    """Checks header of given file for compression type


    Given a filename, opens file and reads first line to check if
    file has BGZF or GZIP header. May be extended to check for BCF format

    Parameters
    ----------
    filename : str
        Name of file to be checked

    Returns
    -------
    extension : str {'bgzip','gzip','vcf','other'}
        File extension as indicated by header

    """
    typ = checkIfGzip(vcfname)
    if typ != 'nozip':
        return typ
    f = open(vcfname)
    l = f.readline()
    f.close()
    VCF_TAG='##fileformat=VCF'
    if l[:len(VCF_TAG)] == VCF_TAG:
        return 'vcf'
    return 'other'

def checkIfCpG(record,fasta_ref,offset=0,add_chr=False):
    dr = None
    pos = record.pos
    c = record.chrom
    if record.alts is None:
        return False
    if add_chr:
        c = 'chr'+record.chrom
    if record.ref == 'C' and 'T' in record.alts:
        seq = fasta_ref.fetch(c,pos-1,pos+1)
        if seq[0].upper() != 'C':
            logging.warning('%s %d has bad base %s' % (record.chrom,record.pos,seq[0]))
            #raise Exception("checkIfCpG function not lining up properly")
        if seq[1].upper() == 'G':
            return True
        return False
    elif record.ref == 'G' and 'A' in record.alts:
        seq = fasta_ref.fetch(c,pos-2,pos)
        if seq[1].upper() != 'G':
            logging.warning('%s %d has bad base %s' % (record.chrom,record.pos,seq[1]))
            #raise Exception("checkIfCpg function not lining up on negative strand")
        if seq[0].upper() == 'C':
            return True
        return False
    return False

def checkForDuplicates(rec_list,pass_list):
    for i in range(len(rec_list)-1):
        if rec_list[i].pos == rec_list[i+1].pos:
            pass_list[i] = False
            pass_list[i+1] = False

def checkForMultiallele(rec_list,pass_list):
    for i in range(len(rec_list)):
        if i != len(rec_list)-1 and rec_list[i].pos == rec_list[i+1].pos:
            pass_list[i] = False
            pass_list[i+1] = False
        if len(rec_list[i].alleles) > 2:
            pass_list[i] = False

def flipChrom(chrom):
    if chrom[0:3] == 'chr':
        return chrom[0:3]
    return 'chr'+chrom

def getAlleleCountDict(rec):
    alleles = defaultdict(int)
    total_sites = 0
    missing_inds = 0
    for j in range(len(rec.samples)):
        samp = rec.samples[j]
        if None in samp.alleles:
            missing_inds += 1
        for k in range(len(samp.alleles)):
            b = samp.alleles[k]
            if b is not None:
                alleles[b] += 1
            total_sites+=1
    return alleles, total_sites, missing_inds

def isInformative(rec, mincount=2, alleles=None):
    count = 0
    if alleles is None:
        alleles, total_sites, missing_inds = getAlleleCountDict(rec)
    if len(alleles) != 2:
        return False
    i1,i2 = alleles.keys()
    return (alleles[i1] >= mincount and alleles[i2] >= mincount)

def getPassSites(record_list, remove_cpg=False, remove_indels=True,
                 remove_multiallele=True, remove_missing=0,
                 inform_level=2, fasta_ref=None):
    pass_list = [True for r in record_list]
    if remove_cpg == True and fasta_ref is None:
        raise Exception("CpG removal requires a reference")
    if inform_level > 2 or inform_level < 0:
        raise Exception("Inform level %d must be between 0 and 2" % inform_level)
    if remove_multiallele:
        checkForMultiallele(record_list,pass_list)
    for i in range(len(record_list)):
        rec = record_list[i]
        logging.info(rec.pos)
        if remove_indels and not checkRecordIsSnp(rec):
            pass_list[i] = False
        if remove_cpg and checkIfCpG(rec,fasta_ref):
            pass_list[i] = False
        alleles,total_sites,missing_inds = getAlleleCountDict(rec)
        if remove_missing != -1 and missing_inds > remove_missing:
            pass_list[i] = False
        if inform_level != 0 and not isInformative(rec,mincount=inform_level,alleles=alleles):
            pass_list[i] = False
    #pp = zip([r.pos for r in record_list],pass_list)
    #for ppp in pp:
    #    logging.info(ppp)
    return pass_list


def filterSites(record_list, remove_cpg=False, remove_indels=True,
                remove_multiallele=True, remove_missing=0, inform_level=2,
                fasta_ref=None):
    pass_list = getPassSites(record_list,remove_cpg,remove_indels,remove_multiallele,remove_missing,inform_level,fasta_ref)
    out_list = []
    for i in range(len(pass_list)):
        if pass_list[i]:
            out_list.append(record_list[i])
    return out_list



class VcfReader():
    def __init__(self, vcfname, compress_flag=False, subsamp_num=None,
                 subsamp_fn=None, subsamp_list=None, index=None, popmodel=None, use_allpop=False):

        ext = checkFormat(vcfname)
        if ext in ['gzip','other'] :
            raise Exception(('Input file %s is gzip-formatted, must be either '
                             'uncompressed or zipped with bgzip' % vcfname))
        self.file_uncompressed = (ext == 'vcf')
        self.reader_uncompressed = (self.file_uncompressed and not compress_flag)
        self.popmodel = None
        self.popkeys = None
        if popmodel is not None and use_allpop:
            raise Exception("Popmodel and allpop cannot both be specified")
        if compress_flag and file_uncompressed:
            vcfname = compressVcf(vcfname)
        if subsamp_num is not None:
            if subsamp_list is not None:
                raise Exception('Multiple subsampling options called in getVcfReader')
            subsamp_list = getSubsampleList(vcfname, subsamp_num)
        elif subsamp_fn is not None:
            if subsamp_list is not None:
                raise Exception('Multiple subsampling options called in getVcfReader')
            subsamp_file = open(subsamp_fn,'r')
            subsamp_list = [l.strip() for l in subsamp_file.readlines()]
            subsamp_file.close()

        if index is None:
            self.reader = pysam.VariantFile(vcfname)
        else:
            self.reader = pysam.VariantFile(vcfname, index_filename=index)
        if popmodel is not None:
            self.popmodel = popmodel
            popsamp_list = popmodel.inds
            self.reader.subset_samples(popsamp_list)
            self.setPopIdx()
        if use_allpop:
            self.setAllPop()
        if subsamp_list is not None:
            logging.debug('Subsampling %d individuals from VCF file' %
            (len(subsamp_list)))
            self.reader.subset_samples(subsamp_list)
        self.prev_last_rec = next(self.reader)
        self.chr_in_chrom = (self.prev_last_rec.chrom[0:3] == 'chr')

    def fetch(self, chrom=None, start=None, end=None):
        return self.reader.fetch(chrom, start, end)

    def getRecordList(self, region=None, chrom=None, start=None,
                      end=None):
        if self.reader_uncompressed:
            ret, self.prev_last_rec = getRecordListUnzipped(self.reader, self.prev_last_rec, region, add_chr=self.chr_in_chrom)
            return ret
        else:
            return getRecordList(self.reader, region, chrom, start, end, self.chr_in_chrom)

    def setPopIdx(self):
        self.popkeys = {}
        sample_names = [l for l in self.reader.header.samples]
        for p in self.popmodel.pop_list:
            self.popkeys[p] = []
            for ind in self.popmodel.ind_dict[p]:
                self.popkeys[p].append(sample_names.index(ind))

    def close(self):
        self.reader.close()

    def setAllPop(self):
        self.popkeys = {'ALL':[]}
        for i in range(len(self.reader.header.samples)):
            self.popkeys['ALL'].append(i)

def modChrom(c,vcf_chr):
    if c is None:
        return None
    if vcf_chr and c[:3] != 'chr':
        return 'chr'+c
    if not vcf_chr and c[:3] == 'chr':
        return c[3:]
    return c

def getRecordList(vcf_reader, region=None, chrom=None, start=None,
                  end=None, add_chr=False):
    """Returns list for use in subsampling from input file"""
    if region is not None:
        c = modChrom(region.chrom,add_chr)
        var_sites = vcf_reader.fetch(c, region.start, region.end)
    else:
        c = modChrom(chrom,add_chr)
        var_sites = vcf_reader.fetch(c, start, end)
    lst = []
    for rec in var_sites:
        lst.append(rec)
    return lst


def getRecordListUnzipped(vcf_reader, prev_last_rec, region=None, chrom=None,
                          start=None, end=None, add_chr=False):
    """Method for getting record list from unzipped VCF file.

    This method will sequentially look through a VCF file until it finds
    the given `start` position on `chrom`, then add all records to a list
    until the `end` position has been reached. Note that `prev_last_rec`
    must be kept track of externally to ensure that if consecutive regions
    are called, the record of the first variant outside the first region
    is not lost.

    Parameters
    ----------
    vcf_reader : pysam VariantFile object
        VCF reader initialized from other function
    region : region object
        Region object with start and end coordinates of region of interest.
    prev_last_rec : VariantRecord object
        Variable with last record read from VcfReader. Stored here so that
        if two adjacent regions are called, the overflow record from the
        first region is still included in the next region

    Returns
    -------
    lst : list
        List of records in given gene region
    prev_last_rec : VariantRecord object
        First record after target region, for use in next call

    """
    lst = []
    if region is None:
        lst.append(prev_last_rec)
        for rec in vcf_reader:
            lst.append(rec)
        return lst, lst[-1]

    if (prev_last_rec is not None and
        region.containsRecord(prev_last_rec) == 'in'):
        lst.append(prev_last_rec)
    elif (prev_last_rec is not None and
         region.containsRecord(prev_last_rec) == 'after'):
        return []
    rec = next(vcf_reader,None)
    if rec is None:
        return lst,None
    place = region.containsRecord(rec)
    while rec is not None and place != 'after':
        if place == 'in':
            lst.append(rec)
        rec = next(vcf_reader,None)
        if rec is None:
            break
        place = region.containsRecord(rec)
    prev_last_rec = rec
    return lst, prev_last_rec


def checkRecordIsSnp(rec):
    """Checks if this record is a single nucleotide variant, returns bool."""
    if len(rec.ref) != 1:
        return False
    if rec.alts is None:
        return False
    for allele in rec.alts:
        if len(allele) != 1:
            return False
    return True


def getSubsampleList(vcfname, ss_count):
    """Returns a list of the first `ss_count` individuals in `vcfname`

    """

    vcf_o = pysam.VariantFile(vcfname)
    rec = next(vcf_o)
    vcf_o.close()
    lst = []
    for samp in rec.samples:
        lst.append(samp)
    return lst[:int(ss_count)]


def compressVcf(vcfname,forceflag=False,remove=False):
    """Runs bgzip and tabix on input VCF file.

    Using the pysam library, this function runs the bgzip and tabix utilities
    on the given input file. By default, this will not overwrite an existing
    zipped file, but will overwrite an existing index. `remove` can be set to
    delete the unzipped file.

    Parameters
    ----------
    vcfname : str
        Name of uncompressed VCF file
    forceflag : bool (False)
        If true, will overwrite (vcfname).gz if it exists
    remove : bool (False)
        If true, will delete uncompressed source file

    Returns
    -------
    cvcfname : str
        Filepath to compressed VCF file
    """
    cvcfname = vcfname+".gz"
    pysam.tabix_compress(vcfname,cvcfname,force=forceflag)
    pysam.tabix_index(cvcfname,preset="vcf",force=True)
    if remove:
        os.remove(vcfname)
    return cvcfname

def vcfRegionName(prefix, region, ext, oneidx=False,
                  halfopen=True, sep='-'):
    chrom = region.toStr(halfopen, oneidx, sep)
    return prefix+'_'+chrom+'.'+ext

def getRecordsInRegion(region, record_list):
    sub_list = []
    for i in range(len(record_list)):
        loc = region.containsRecord(record_list[i])
        if loc == "in":
            sub_list.append(record_list[i])
        elif loc == "after":
            break
    return sub_list





#def getVcfReader(args):
def getVcfReader(vcfname, compress_flag=False, subsamp_num=None,
                 subsamp_fn=None, subsamp_list=None, index=None):
    """Returns a reader for a given input VCF file.

    Given a filename, filetype, compression option, and optional Subsampling
    options, will return a pysam.VariantFile object for iteration and
    a flag as to whether this file is compressed or uncompressed.

    Parameters
    ----------
    vcfname : str
        Filename for VCF file. The extension of this file will be used to
        determine whether it is compressed or not unless `var_ext` is set.
    var_ext : str (None)
        Extension for VCF file if it is not included in the filename.
    compress_flag : bool (False)
        If filetype is uncompressed and this is set to true, will run
        compressVcf function.
    subsamp_num : int (None)
        If set, will randomly select `subsamp_num` individuals (not
        genotypes) from the input VCF file and return a reader with
        only those data.
    subsamp_fn : str (None)
        If set, will return a reader with only data from the samples listed
        in the file provided. Cannot be used with other subsampling options.
    subsamp_list : list (None)
        If set, will return reader with records containing only
        individuals named in the list. Cannot be used with other subsampling
        options.

    Returns
    -------
    vcf_reader : pysam.VariantFile
        A reader that can be iterated through for variant records. If
        compressed, it will be able to use the pysam fetch method, otherwise
        it must be read through sequentially
    reader_uncompressed : bool
        If True, VCF reader is uncompressed. This means the fetch method
        cannot be used and region access must be done using the
        "getRecordListUnzipped" method.

    """
    ext = checkFormat(vcfname)
    if ext in ['gzip','other'] :
        raise Exception(('Input file %s is gzip-formatted, must be either '
                         'uncompressed or zipped with bgzip' % vcfname))
    file_uncompressed = (ext == 'vcf')
    reader_uncompressed = (file_uncompressed and not compress_flag)
    if compress_flag and file_uncompressed:
        vcfname = compressVcf(vcfname)
    #subsamp_list = None
    if subsamp_num is not None:
        if subsamp_list is not None:
            raise Exception('Multiple subsampling options called in getVcfReader')
        subsamp_list = getSubsampleList(vcfname, subsamp_num)
    elif subsamp_fn is not None:
        if subsamp_list is not None:
            raise Exception('Multiple subsampling options called in getVcfReader')
        subsamp_file = open(subsamp_fn,'r')
        subsamp_list = [l.strip() for l in subsamp_file.readlines()]
        subsamp_file.close()
    if index is None:
        vcf_reader = pysam.VariantFile(vcfname)
    else:
        vcf_reader = pysam.VariantFile(vcfname, index_filename=index)
    if subsamp_list is not None:
        logging.debug('Subsampling %d individuals from VCF file' %
        (len(subsamp_list)))
        vcf_reader.subset_samples(subsamp_list)
    return vcf_reader, reader_uncompressed
