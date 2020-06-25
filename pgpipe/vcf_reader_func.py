import sys
import pysam
import logging
import struct
from random import sample
from collections import defaultdict
import os
import gzip
import numpy as np


#sys.path.insert(0,os.path.abspath(os.path.join(os.pardir, 'andrew')))
from pgpipe.logging_module import initLogger
from pgpipe.model import read_model_file
#from pgpipe import test_cython

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
    except OSError:
        return 'nozip'
    except IOError:
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
    if vcfname == '-':
        return 'vcf' #may want to note its a stream
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



def checkIfCpG(record,fasta_ref,add_chr=False):
    #Note: will not detect multiallelic/indel combos
    dr = None
    pos = record.pos
    c = record.chrom
    if record.alts is None:
        return False
    if add_chr:
        c = 'chr'+record.chrom
    if record.ref == 'C' and 'T' in record.alts:
        try:
            seq = fasta_ref.fetch(c,pos-1,pos+1)
        except KeyError:
            seq = fasta_ref.fetch(flipChrom(c),pos-1,pos+1)
        if seq[0].upper() != 'C':
            logging.warning('%s %d has bad base %s' % (record.chrom,record.pos,seq[0]))
        if seq[1].upper() == 'G':
            return True
        return False
    elif record.ref == 'G' and 'A' in record.alts:
        try:
            seq = fasta_ref.fetch(c,pos-2,pos)
        except KeyError:
            seq = fasta_ref.fetch(flipChrom(c),pos-2,pos)
        if seq[1].upper() != 'G':
            logging.warning('%s %d has bad base %s' % (record.chrom,record.pos,seq[1]))
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

def checkPopWithoutMissing(rec_list,model,pop_keys,min_per_pop=5):
    for pop in pop_keys:
        present_count = 0
        for i in range(len(pop_keys[pop])):
            tidx = pop_keys[pop][i]
            if tidx != -1:
                all_data = True
                for rec in rec_list:
                    rec_alleles = rec.samples[tidx].alleles
                    if rec_alleles[0] in [None,'N']:
                        all_data = False
                        break
                if all_data:
                    present_count += 1
                if present_count >= min_per_pop:
                    break
        if present_count < min_per_pop:
            return False
    return True

def readSinglePopModel(popfilename,popname=None):
    popmodels = read_model_file(popfilename)
    if len(popmodels) != 1:
        if popname is None:
            raise Exception("Model file %s needs one model to be specified")
        return popmodels[popname]
    else:
        pp = list(popmodels.keys())
        return popmodels[pp[0]]


def flipChrom(chrom):
    if chrom[0:3] == 'chr':
        return chrom[0:3]
    return 'chr'+chrom

#def getAlleleCountCython(rec,idxlist=None):
#    alleles = defaultdict(int)
#    srec = str(rec)
#    if idxlist is None:
#        d_array = test_cython.getAlleleCountArray(srec.encode(),len(srec))
#    else:
#        d_array = test_cython.getAlleleCountArraySub(srec.encode(),len(srec),idxlist,len(idxlist))
#    for i,a in enumerate(rec.alleles):
#        alleles[a] = d_array[i+48]
#    if d_array[46] != 0:
#        alleles['N'] = d_array[46]
#    return alleles

def getAlleleCountDict(rec,idx_list=None):
    """
    Returns dict with allele counts for all
    """
    alleles = defaultdict(int)
    total_sites = 0
    missing_inds = 0
    if idx_list is None:
        idx_list = range(len(rec.samples))
    for j in idx_list:
        samp = rec.samples[j]
        if None in samp.alleles:
            alleles['N'] += len(samp.alleles)
            #missing_inds += 1
        for k in range(len(samp.alleles)):
            b = samp.alleles[k]
            if b is not None:
                alleles[b] += 1
            total_sites+=1
    return alleles

def getAlleleStats(rec):
    acd = getAlleleCountDict(rec)
    missing_inds = (acd['N'] if 'N' in acd.keys() else 0)
    total_sites = sum(acd.values())-missing_inds
    return acd,total_sites,missing_inds

#def parseRecString(recstr):



def getAlleleStatsAlt(rec):
    ala = str(rec).strip().split()[9:]
    l = []
    allele_idx = rec.alleles
    d = {}
    for a in allele_idx:
        d[a] = 0
    d['.'] = 0
    for j in range(len(ala)):
        for k in ala[j].split(':')[0].split(ala[j][1]):
            #d[allele_idx]
            if k == '.':
                d['N'] += 1
            else:
                d[allele_idx[int(k)]] += 1
    #print (str(d))
    return d
        #d[allele_idx[ii]] += 1 for ii in ala[j].split(':').split('|'))
        #l.extend(ala.split(':')[0].split('|'))

    #npa = np.array(l)

#def getAlleleStatsCython(rec,idx_list=None):
#    acd = getAlleleCountCython(rec,idx_list)
#    missing_inds = (acd['.'] if '.' in acd.keys() else 0)
#    total_sites = sum(acd.values())-missing_inds
#    return acd,total_sites,missing_inds



def isInvariant(rec):
    alleles, total_sites, missing_inds = getAlleleStats(rec)
    if len(alleles) == 1:
        return True
    return False

def isInformative(rec, mincount=2, alleles=None):
    count = 0
    if alleles is None:
        alleles, total_sites, missing_inds = getAlleleStats(rec)
    if 'N' in alleles.keys(): #jh edit 6/17/2020 alleles may have N's 
        del(alleles['N'])
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
        if remove_indels and not checkRecordIsSnp(rec):
            pass_list[i] = False
        if remove_cpg and checkIfCpG(rec,fasta_ref):
            pass_list[i] = False
        alleles,total_sites,missing_inds = getAlleleStats(rec)
        if remove_missing != -1 and missing_inds > int(remove_missing):
            pass_list[i] = False
        if inform_level != 0 and not isInformative(rec,mincount=inform_level,alleles=alleles):
            pass_list[i] = False
    #pp = zip([r.pos for r in record_list],pass_list)
    #for ppp in pp:
    #    logging.info(ppp)
    return pass_list

def checkRecordPass(rec, remove_cpg=False, remove_indels=True, 
                    remove_multiallele=True, remove_missing=0,
                    inform_level=2,fasta_ref=None):
    if remove_cpg and fasta_ref is None:
        raise Exception("CpG removal requires a reference")
    if inform_level > 2 or inform_level < 0:
        raise Exception("Inform level %d must be between 0 and 2" % inform_level)
    if remove_indels and not checkRecordIsSnp(rec):
        return False
    if remove_cpg and checkIfCpG(rec,fasta_ref):
        return False
    if remove_missing != -1 or inform_level != 0:
        alleles,total_sites,missing_inds = getAlleleStats(rec)
##        if len(alleles) > 1:
##            return True
        if remove_missing != -1 and missing_inds > int(remove_missing):
            return False
##  jh 6/17/2020  was calling isInformative which was removing variants,  replace with isInvariant
##        if inform_level != 0 and not isInformative(rec,mincount=inform_level,alleles=alleles):
        if inform_level != 0 and isInvariant(rec):
            return False
    return True


def filterSites(record_list, remove_cpg=False, remove_indels=True,
                remove_multiallele=True, remove_missing=0, inform_level=2,
                fasta_ref=None):
    pass_list = getPassSites(record_list,remove_cpg,remove_indels,remove_multiallele,remove_missing,inform_level,fasta_ref)
    out_list = []
    for i in range(len(pass_list)):
        if pass_list[i]:
            out_list.append(record_list[i])
    return out_list

def crossModelAndVcf(pop_list,vcf_samples,allow_missing_inds=True):
    missing_list = []
    present_list = []
    for p in pop_list:
        if p not in vcf_samples:
            missing_list.append(p)
        else:
            present_list.append(p)
    if len(missing_list) != 0:
        if allow_missing_inds:
            logging.warning("Samples %s are missing from VCF" % (','.join(missing_list)))
        else:
            raise Exception("Samples %s are missing from VCF" % (','.join(missing_list)))
    return present_list

def getIndsWithAllData(rec_list):
    idx_list = []
    for i in range(len(rec_list[0].samples)):
        has_data = True
        for rec in rec_list:
            if rec.samples[i].alleles[0] in [None,'N']:
                has_data = False
        if has_data:
            idx_list.append(i)
    return idx_list


class VcfReader():
    def __init__(self, vcfname, compress_flag=False, subsamp_num=None,
                 subsamp_fn=None, subsamp_list=None, index=None, 
                 popmodel=None, use_allpop=False, allow_missing_inds=True):

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
        if compress_flag and self.file_uncompressed:
            vcfname = compressVcf(vcfname)
        subsamp_list = None
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
        
        self.openSetInds(vcfname,index,popmodel,use_allpop,subsamp_list)
        self.info_rec = next(self.reader)
        self.prev_last_rec = None #next(self.reader)
        self.reader.close()
        self.openSetInds(vcfname,index,popmodel,use_allpop,subsamp_list)
        self.chr_in_chrom = (self.info_rec.chrom[0:3] == 'chr')

    def openSetInds(self, vcfname, index, popmodel, use_allpop, subsamp_list):
        if index is None:
            self.reader = pysam.VariantFile(vcfname)
        else:
            self.reader = pysam.VariantFile(vcfname, index_filename=index)
        if popmodel is not None:
            self.popmodel = popmodel
            popsamp_list = popmodel.inds
            vcf_names = [l for l in self.reader.header.samples]
            present_list = crossModelAndVcf(popsamp_list,vcf_names)
            #self.reader.subset_samples(popsamp_list)
            self.reader.subset_samples(present_list)
            self.setPopIdx(present_list)
        if use_allpop:
            self.setAllPop()
        if subsamp_list is not None:
            logging.debug('Subsampling %d individuals from VCF file' %
            (len(subsamp_list)))
            self.reader.subset_samples(subsamp_list)

    def fetch(self, chrom=None, start=None, end=None):
        return self.reader.fetch(chrom, start, end)

    def getRecordList(self, region=None, chrom=None, start=None,
                      end=None):
        if self.reader_uncompressed:
            ret, self.prev_last_rec = getRecordListUnzipped(self.reader, self.prev_last_rec, region, add_chr=self.chr_in_chrom)
            return ret
        else:
            return getRecordList(self.reader, region, chrom, start, end, self.chr_in_chrom)

    def setPopIdx(self,present_list):
        self.popkeys = {}
        sample_names = [l for l in self.reader.header.samples]
        for p in self.popmodel.pop_list:
            self.popkeys[p] = []
            for indiv in self.popmodel.ind_dict[p]:
                if indiv in present_list:
                    self.popkeys[p].append(sample_names.index(indiv))
                else:
                    self.popkeys[p].append(-1)

    def close(self):
        self.reader.close()

    def setAllPop(self):
        self.popkeys = {'ALL':[]}
        for i in range(len(self.reader.header.samples)):
            self.popkeys['ALL'].append(i)

    def returnNames(self,index_list):
        return [self.reader.header[i] for i in index_list]

    def getRegionIterUnzipped(self,region=None,add_chr=False):
        if region is None:
            a = self.getNext()
            while a is not None:
                #self.prev_last_rec = a
                yield a
                self.prev_last_rec = a
                a = self.getNext()
            self.prev_last_rec = None
            return
        trec = (self.prev_last_rec if self.prev_last_rec is not None else self.getNext())
        if trec is None:
            self.prev_last_rec = trec
            return

        self.prev_last_rec = trec
        stat = region.containsRecord(trec)
        if stat == 'after':
            return
        if stat == 'in':
            yield trec
            self.prev_last_rec = trec
            trec = self.getNext()
            if trec is None:
                self.prev_last_rec = None
                return
            stat = region.containsRecord(trec)
        elif stat == 'before':
            while stat == 'before':
                self.prev_last_rec = trec
                trec = self.getNext()
                if trec is None:
                    self.prev_last_rec = trec
                    return
                stat = region.containsRecord(trec)
        while stat != 'after':
            yield trec
            self.prev_last_rec = trec
            trec = self.getNext()
            if trec is None:
                self.prev_last_rec = trec
                return
            stat = region.containsRecord(trec)
        self.prev_last_rec = trec

        
    def getRegionIterZipped(self,region=None):
        if region is None:
            treader = self.reader.fetch()
        else:
            treader = self.reader.fetch(region.chrom,region.start,region.end)
        for rec in treader:
            yield rec
        return

    def getRegionIter(self,region=None):
        if self.reader_uncompressed:
            return self.getRegionIterUnzipped(region=region)
        else:
            return self.getRegionIterZipped(region=region)
        
    def getNext(self,set_prev=True):
        try:
            trec = next(self.reader)
            # jh 6/17/2020 added this crude trap for reading a line with only a newline symbol
            if trec.chrom == '':
                return None
            else:
                return trec
        except StopIteration as e:
            return None
            

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
        #lst.append(prev_last_rec)
        for rec in vcf_reader:
            lst.append(rec)
        return lst, lst[-1]

    if (prev_last_rec is not None and
        region.containsRecord(prev_last_rec) == 'in'):
        lst.append(prev_last_rec)
    elif (prev_last_rec is not None and
         region.containsRecord(prev_last_rec) == 'after'):
        return [],prev_last_rec  #jh added ',prev_last_rec'  6/5/2020
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
