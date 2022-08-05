'''
    The VcfReader class is a pysam-based wrapper for reading and writing 
    records from VCF files that integrates with other PPP functionality, 
    including the GenomeRegion class for handling BED regions/files, and the 
    Model class for specifying subpopulations and samples. It can read input 
    from unzipped VCF, bg-zipped VCF, and BCF files. It c
'''
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
from logging_module import initLogger
from model import read_model_file
#from pgpipe import test_cython

def checkIfGzip(filename):
    '''
    File Format Checker
    Identifies the file format of the given file.
    Parameters
    ----------
    filename : string
        Filepath of target file
    Returns
    -------
    type : string
        'BCF', 'bgzip', 'gzip', 'nozip', or 'other', where all -zip options 
        indicate a VCF file.
    '''
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
    '''
    VCF File Format Checker
    Identifies the format of a VCF file with given filename
    Parameters
    ----------
    filename : string
        Filepath of target file
    Returns 
    -------
    type : string
        'bgzip', 'gzip', or 'nozip'. Note that VcfReader cannot read gzipped VCF files.
    '''
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
    '''
    Checks if a given record in variant file is a CpG. Currently does not 
    support checking sites where there are two or more alternate alleles and
    one of them is an indel.
    Parameters
    ----------
    record : pysam.VariantRecord
        Record from VCF file to be checked
    fasta_ref : pysam.FastaFile
        Reference genome sequence, must be from start of chromosome
    add_chr : boolean (optional)
        If true, prefixes chromosome pulled from record with 'chr' to resolve
        situations where record doesn't have 'chr' in chromosome name but
        reference does.
        Default: False
    Return
    ------
    Boolean
        Returns true if record is located at a CpG, and false otherwise. 
    '''
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
    '''
    Finds if a list of VariantRecords has multiple records at same position
    Parameters
    ----------
    rec_list : list (pysam.VariantRecord)
        List of target records
    pass_list : list (Boolean)
        List with an entry for each record, either with each entry preset to 
        True, or with results from a different filtering function. Will be 
        modified to include False entries for duplicate records.
    '''
    for i in range(len(rec_list)-1):
        if rec_list[i].pos == rec_list[i+1].pos:
            pass_list[i] = False
            pass_list[i+1] = False

def checkForMultiallele(rec_list,pass_list):
    '''
    Finds if a list of VariantRecords has multi-allelic sites, and marks the 
    matching entry in a 'pass_list' of booleans to indicate a failure.
    Parameters
    ----------
    rec_list : list (pysam.VariantRecord)
        List of target records
    pass_list : list (Boolean)
        List with an entry for each record, either with each entry preset to 
        True, or with results from a different filtering function. Will be 
        modified to include False entries for multi-allelic sites.
    '''
    for i in range(len(rec_list)):
        if i != len(rec_list)-1 and rec_list[i].pos == rec_list[i+1].pos:
            pass_list[i] = False
            pass_list[i+1] = False
        if len(rec_list[i].alleles) > 2:
            pass_list[i] = False

def checkPopWithoutMissing(rec_list,model,pop_keys,min_per_pop=5):
    '''
    Checks that given populations in a model each have a minimum number of 
    individuals with data for the region in question. Returns true if 
    each population has at least min_per_pop individuals with no missing
    data in target region, false otherwise.
    Parameters
    ----------
    rec_list : list (pysam.VariantRecord)
        List of target records in a region
    model : model file
        Model file
    pop_keys : list (string)
        List of subpopulations in model file to examine
    min_per_pop : int (optional)
        Minimum number of individuals in a subpopulation that have no missing
        data.
        Default: 5
    Returns
    -------
    boolean
        True if enough individuals in each population have all their data,
        false otherwise.
    
    '''
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


def flipChrom(chrom):
    '''
    Flips whether a chomosome name is prefixed with 'chr'. Can be 
    used in try blocks as an easy check, yielding the expected
    exception if the modified chromosome name cannot be found. 
    '''
    if chrom[0:3] == 'chr':
        return chrom[0:3]
    return 'chr'+chrom

#def getAlleleCountCython(rec,idxlist=None):
#    #Removed until cython operation is re-enabled
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
    Returns dict with allele counts for all alleles with at least one 
    occurrence at the input record, plus an 'N' entry for number of alleles
    missing data. 
    Parameters
    ----------
    rec : pysam.VariantRecord
        Target record
    idx_list : list (integer) (optional)
        List of indexes of individuals to be looked at (defaults to looking
        at all individuals)
        Default: None
    Returns
    -------
    alleles : dict
        Dictionary with keys for all alleles present in record, plus an
        entry for missing data keyed 'N' if missing data is present.
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

def getAlleleStats(rec, idx_list = None):
    '''
    For a given input record, returns the allele count dictionary from 
    getAlleleCountDict, plus a count of the total number of genotypes with
    data and missing data (totaling to number of sites at record).
    Parameters
    ----------
    rec : pysam.VariantRecord
        Target record
    idx_list : list (integer) (optional)
        List of indexes of individuals to be looked at (defaults to looking
        at all individuals)
        Default: None
    Returns
    -------
    alleles : dict
        Dictionary with keys for all alleles present in record, plus an
        entry for missing data keyed 'N' if missing data is present.
    total_sites : int
        Count of genotypes that contain an allele (not missing data)
    missing_sites : int
        Count of genotypes that are missing data
    '''
    acd = getAlleleCountDict(rec, idx_list)
    missing_sites = (acd['N'] if 'N' in acd.keys() else 0)
    total_sites = sum(acd.values())-missing_sites
    return acd,total_sites,missing_sites


def getAlleleStatsAlt(rec):
    #development method for getting data from record string
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



def isInvariant(rec, alleles=None):
    '''
    Given a record, returns whether the record is invariant 
    (and has all genotypes with data)
    Parameters
    ----------
    rec : pysam.VariantRecord
        Record to be checked
    alleles : dictionary (optional, from getAlleleCountData)
        If provided, uses prior ACD instead of calculating a new one. 
    Returns
    -------
    bool
       True if record has a single allele and no missing data,
       false otherwise. 
    '''
    if alleles is None: 
        alleles, _, _ = getAlleleStats(rec)
    if len(alleles) == 1:
        return True
    return False

def isInformative(rec, mincount=2, alleles=None):
    '''
    Given a record, returns whether the site is considered informative in the
    four-gamete test criteria. Generally defined as whether a site is 
    bi-allelic where each allele has at least two copies. Missing data does
    not affect this calculation. 
    Parameters
    ----------
    rec : pysam.VariantRecord
        Record to be checked
    mincount : integer (optional)
        Minimum number of each allele required for record to be considered
        informative.
        Default: 2
    alleles : dictionary (optional, from getAlleleCountData)
        If provided, uses prior ACD instead of calculating a new one. 
    Returns
    -------
    boolean
        True if record is informative, false otherwise. 
        
    '''
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
    '''
    Given a list of VCF records, constructs a list of booleans corresponsing
    to the input list where a matched index is True if the record passes given
    criteria, and false otherwise. Can filter for CpGs, indels, multi-allelic 
    sites, sites with missing data, and non-informative sites. 
    Parameters
    ----------
    record_list : list (pysam VariantRecord)
        Input list of records to be checked
    Filter Options
    --------------
    remove_cpg : boolean
        If set to true, will fail records that are CpGs. Requires fasta_ref.
        Default: False
    remove_indels : boolean
        If set to true, will fail records that have at least one indel allele.
        Default: True
    remove_multiallele : boolean
        If set to true, will fail records that have more than two alleles.
        Default: True
    remove_missing : integer
        Will fail records with at least this many genotypes missing data. 
        0 will fail anything missing data, -1 will allow for any amount
        of missing data.
        Default: 0
    inform_level : integer
        Will set the informative site filter to require at least this many
        alleles from each record. Set to 0 to prevent informative site check.
        Default: 2
    fasta_ref : pysam FastaFile
        If CpG checking is set to true, use this object as the reference
        genome.
    Returns
    -------
    pass_list : list (boolean)
        For each record in the input record list, contains True if that
        record passes every specified filter and False if it fails at least
        one. 
    '''
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
    return pass_list

def checkRecordPass(rec, remove_cpg=False, remove_indels=True, 
                    remove_multiallele=True, remove_missing=0,
                    inform_level=2,fasta_ref=None):
    '''
    Returns whether a single record passes all of a given criteria, including
    options for a site not being one or more of the following: a CpG, an indel,
    multi-allelic, non-informative (four-gamete critera), or with some threshhold
    of missing data.
    Parameters
    ----------
    rec : pysam VariantRecord
        Input record to be checked
    Filter Options
    --------------
    remove_cpg : boolean
        If set to true, will fail records that are CpGs. Requires fasta_ref.
        Default: False
    remove_indels : boolean
        If set to true, will fail records that have at least one indel allele.
        Default: True
    remove_multiallele : boolean
        If set to true, will fail records that have more than two alleles.
        Default: True
    remove_missing : integer
        Will fail records with at least this many genotypes missing data. 
        0 will fail anything missing data, -1 will allow for any amount
        of missing data.
        Default: 0
    inform_level : integer
        Will set the informative site filter to require at least this many
        alleles from each record. Set to 0 to prevent informative site check.
        Default: 2
    fasta_ref : pysam FastaFile
        If CpG checking is set to true, use this object as the reference
        genome.
    Returns
    -------
    pass : boolean
        Returns true if records passes all filters, false otherwise. 
    '''
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
    '''
    Given a record list, returns a new list with records that pass specified
    filtering methods. Methods of filtering include removal of CpGs (with 
    reference genome), indels, multi-allelic sites, sites with missing data,
    and sites that do not pass the four-gamete test criterion (which requires
    each allele in a biallelic site have at least two copies at a site).
    Parameters
    ----------
    record_list : list (pysam VariantRecord)
        Input list of records to be checked
    Filter Options
    --------------
    remove_cpg : boolean
        If set to true, will fail records that are CpGs. Requires fasta_ref.
        Default: False
    remove_indels : boolean
        If set to true, will fail records that have at least one indel allele.
        Default: True
    remove_multiallele : boolean
        If set to true, will fail records that have more than two alleles.
        Default: True
    remove_missing : integer
        Will fail records with at least this many genotypes missing data. 
        0 will fail anything missing data, -1 will allow for any amount
        of missing data.
        Default: 0
    inform_level : integer
        Will set the informative site filter to require at least this many
        alleles from each record. Set to 0 to prevent informative site check.
        Default: 2
    fasta_ref : pysam FastaFile
        If CpG checking is set to true, use this object as the reference
        genome.
    Returns
    -------
    out_list : list (VariantRecord)
        Parsed list containing only records that pass all specified filtering
        criteria. 
    '''
    pass_list = getPassSites(record_list,remove_cpg,remove_indels,remove_multiallele,remove_missing,inform_level,fasta_ref)
    out_list = []
    for i in range(len(pass_list)):
        if pass_list[i]:
            out_list.append(record_list[i])
    return out_list

def crossModelAndVcf(pop_list,vcf_samples,allow_missing_inds=True):
    '''
    Given a list of individuals from a population or subsample, returns a list
    of all individuals from a given record that are present in the population 
    list. If allow_missing_inds is set to true, will warn if an individual in 
    the population list is missing from the record, otherwise will throw an
    exception. 
    Parameters
    ----------
    pop_list : list (string)
        Names of individuals to check for in the record
    
    vcf_samples : list (string)
        Names of all individuals from a target record
    allow_missing_inds : boolean (optional)
        If set to false, will throw exception if one or more individuals
        from the population list is not found in the record.
        Default: True
    Returns
    -------
    present_list : list (string)
        Names of individuals in both pop_list and vcf_samples. 
    '''
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
    #is this ever used? 
    idx_list = []
    for i in range(len(rec_list[0].samples)):
        has_data = True
        for rec in rec_list:
            if rec.samples[i].alleles[0] in [None,'N']:
                has_data = False
        if has_data:
            idx_list.append(i)
    return idx_list


def matchChrom(c,vcf_chr):
    '''
    Matches a chromosome for a region with a chromosome from a VCF, prefixing 'chr' if
    vcf_chr is True and 'chr' is not already a prefix, and removing 'chr' if vcf_chr is
    false and 'chr' is present. 
    Parameters
    ----------
    c : string
        Chromosome name, usually from reference or genome region
    vcf_chr : bool
        True if 'chr' should be added, false if it should be removed
    Returns
    -------
    new_c : string
        New chromosome name matching format specified by vcf_chr
    '''
    if c is None:
        return None
    if vcf_chr and c[:3] != 'chr':
        return 'chr'+c
    if not vcf_chr and c[:3] == 'chr':
        return c[3:]
    return c


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
    forceflag : bool (optional)
        If true, will overwrite (vcfname).gz if it exists
        Default: False
    remove : bool (optional)
        If true, will delete uncompressed source file
        Default: False
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
    '''
    Given an already output record list, output the records that fall in a
    given region.
    Parameters
    ----------
    region : Region object
        Target region for records.
    record_list : list (pysam VariantRecord)
        List of records to parse.
    Returns
    -------
    sub_list : list (pysam VariantRecord)
        List of records contained in the given region. 
        
    '''
    sub_list = []
    for i in range(len(record_list)):
        loc = region.containsRecord(record_list[i])
        if loc == "in":
            sub_list.append(record_list[i])
        elif loc == "after":
            break
    return sub_list



class VcfReader():
    def __init__(self, vcfname, compress_flag=False, subsamp_num=None,
                 subsamp_fn=None, subsamp_list=None, index=None, 
                 popmodel=None, use_allpop=False, allow_missing_inds=True):
        '''
        Creates a pysam VariantFile parser for bgzipped or uncompressed VCFs,
        as well as BCFs. Can be subsampled for individuals via a string list
        or a population model. 
        Parameters
        ----------
        vcfname : string
            Name of target input VCF/BCF file.
        Optional Parameters
        -------------------
        compress_flag : boolean
            If set to true, will compress the input VCF. Useful if trying to
            access specific genomic regions via index lookup. 
        subsamp_num : int
            If set, will subsample the first n genotypes from a given 
            VariantFile, in order from the first. 
        subsamp_fn : string
            Name of file that contains line-separated list of individuals
            to subsample from the input VariantFile
        subsamp_list : list (string)
            List of individuals to subsample from the input VariantFile
        index : string
            If input VariantFile is compressed but the index is not named
            (filename).tbi, use the provided index filename instead. 
        popmodel : PopModel
            Population model (single) to use for subsampling Variant File. 
            Will pull individuals from population specified, and store
            population data in field of VcfReader.
        use_allpop : boolean
            If set, will treat all individuals in the given VCF as a single 
            population. 
        allow_missing_ids : boolean
            If set, will allow subsampling methods that try to sample 
            individuals not present in the reader to proceed anyway.
        Returns
        -------
        VcfReader object, compressed appropriately and with only the 
        individuals specified being output to subsequent records. 
        '''

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
        '''
        Internal function in VcfReader initialization for opening
        a VCF file with only the given individuals, as specified by a 
        subsampling list or a population model. 
        '''
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
        '''
        Shortcut to access record from VcfReader pysam VariantFile object
        '''
        return self.reader.fetch(chrom, start, end)

    def getRecordList(self, region=None, chrom=None, start=None,
                      end=None):
        '''
        Given a genomeRegion or manually input start/end/chromosome values,
        will output all records from the given VCF in that region. Note that
        if the input VCF in uncompressed, regions must be accessed in sorted
        order. (genomeRegion's list will do this by default)
        Parameters
        ----------
        If a region is provided, it will take precedence over start/end/chrom
        arguments. 
        region : PPP region object
            Region for fetching records
        chrom : string
            If set, will fetch records from given chromosome. Will fetch
            all records if start and end are not specified
        start : int
            If set with chrom, will fetch records on given chromosome starting
            from this position (one-indexed)
        end : int
            If set with chrom and start, genome position to stop fetching of
            records at.
        Returns
        -------
        ret : list (pysam VariantRecord)
            List of all records in given region, or in entire VCF if no 
            region specified.  
        '''
        if self.reader_uncompressed:
            return self.getRecordListUnzipped(region=region, chrom=chrom,
                                              start=start, end=end)
        else:
            return self.getRecordListZipped(region=region, chrom=chrom,
                                            start=start, end=end)

    def getRecordListZipped(self, region=None, chrom=None, start=None,
                            end=None):
        """Returns list for use in subsampling from input file"""
        if region is not None:
            c = matchChrom(region.chrom,self.chr_in_chrom)
            var_sites = self.reader.fetch(c, region.start, region.end)
        else:
            c = matchChrom(chrom,self.chr_in_chrom)
            var_sites = self.reader.fetch(c, start, end)
        lst = []
        for rec in var_sites:
            lst.append(rec)
        return lst

    


    def getRecordListUnzipped(self, region=None, chrom=None,
                              start=None, end=None):
        """Method for getting record list from unzipped VCF file.
        This method will sequentially look through a VCF file until it finds
        the given `start` position on `chrom`, then add all records to a list
        until the `end` position has been reached. Note that `prev_last_rec`
        must be kept track of externally to ensure that if consecutive regions
        are called, the record of the first variant outside the first region
        is not lost.
        Parameters
        ----------
        If a region is provided, it will take precedence over start/end/chrom
        arguments. 
        
        region : region object
            Region object with start and end coordinates of region of interest.
        chrom : string
            If set, will fetch records from given chromosome. Will fetch
            all records if start and end are not specified
        start : int
            If set with chrom, will fetch records on given chromosome starting
            from this position (one-indexed)
        end : int
            If set with chrom and start, genome position to stop fetching of
            records at.
        Returns
        -------
        lst : list
            List of records in given gene region
        """
        lst = []
        if region is None:
            #lst.append(prev_last_rec)
            for rec in self.reader:
                lst.append(rec)
            self.prev_last_rec = lst[-1]
            return lst
            #return lst, lst[-1]

        if (self.prev_last_rec is not None and
            region.containsRecord(self.prev_last_rec) == 'in'):
            lst.append(self.prev_last_rec)
        elif (self.prev_last_rec is not None and
            region.containsRecord(self.prev_last_rec) == 'after'):
            return []
        rec = next(self.reader,None)
        if rec is None:
            self.prev_last_rec = None
            return lst
            #return lst,None
        place = region.containsRecord(rec)
        while rec is not None and place != 'after':
            if place == 'in':
                lst.append(rec)
            rec = next(self.reader,None)
            if rec is None:
                break
            place = region.containsRecord(rec)
        self.prev_last_rec = rec
        return lst
        #return lst, prev_last_rec

    def setPopIdx(self,present_list):
        '''
        Internal method for setting indices at which members of a population
        are accessible from a record.
        '''
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
        '''
        Sets all individuals in a VCF to be members of a single population. 
        Useful for some analyses where populations are required but a model
        file has not been created.
        '''
        self.popkeys = {'ALL':[]}
        for i in range(len(self.reader.header.samples)):
            self.popkeys['ALL'].append(i)

    def returnNames(self,index_list):
        '''
        Returns list of all individuals being output in VCF records.
        '''
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
        '''
        Method for returning iterator over records in a region. Note:
        currently untested.
        '''
        if self.reader_uncompressed:
            return self.getRegionIterUnzipped(region=region)
        else:
            return self.getRegionIterZipped(region=region)
        
    def getNext(self,set_prev=True):
        '''
        Helper method for preventing StopIteration exceptions when
        VcfReader is out of records. 
        '''
        try:
            trec = next(self.reader)
            # jh 6/17/2020 added this crude trap for reading a line with only a newline symbol
            if trec.chrom == '':
                return None
            else:
                return trec
        except StopIteration as e:
            return None

    def getSnpDensityFromHeader(self):
        for x in self.reader.header.records:
            if x.key == 'SNPDensity':
                return float(x.value)
        return None
