import sys
import pysam
import argparse
import os.path
import logging
from logging_module import individualFunctionLogger
from random import sample
from gene_region import Region, RegionList

#Input: VCF file, reference sequence, region list (possibly .bed file)
#Output: Sequences with reference genome overlayed with VCF SNP calls


def createParser():
    parser = argparse.ArgumentParser(description=("Generates sequences from"
                                     " samples from a VCF file, a reference"
                                     " genome, and a list of gene regions."))
    parser.add_argument("--vcf", dest="vcfname", help="Input VCF filename")
    parser.add_argument("--ref", dest="refname", help="Reference FASTA file")
    parser.add_argument("--gr", dest="genename",
                        help="Name of gene region file")
    parser.add_argument("--gr0", dest="gene_idx", action="store_false",
                        help="Gene Region list is 0 index based, not 1")
    parser.add_argument("--indels", dest="indel_flag", action="store_true",
                        help="Include indels when reporting sequences")
    parser.add_argument("--trim-to-ref-length", dest="trim_seq",
                        action="store_true",
                        help=("Trims sequences if indels cause them to be "
                              "longer than reference"))
    subsamp_group = parser.add_mutually_exclusive_group()
    subsamp_group.add_argument('--subsamp_list', dest="subsamp_fn",
                               help="List of sample names to be used")
    subsamp_group.add_argument('--subsamp_num', dest="subsamp_num",
                               help=("Number of individuals to be randomly "
                                     "subsampled from VCF file"))
    return parser


def logArgs(args):
    logging.info('Arguments for vcf_to_seq:')
    for k in vars(args):
        logging.info('Argument %s: %s' % (k, vars(args)[k]))

def validateFiles(args):
    """Validates that files provided to args all exist on users system"""
    for var in ['vcfname', 'refname', 'genename']:
        f = vars(args)[var]
        if not os.path.exists(f):
            raise ValueError('Filepath for %s not found at %s' %
                            (var, f))


def checkRecordIsSnp(rec):
    """Checks if this record is a single nucleotide variant, returns bool."""
    if len(rec.ref) != 1:
        return False
    for allele in rec.alts:
        if len(allele) != 1:
            return False
    return True


def getMaxAlleleLength(alleles):
    """If an indel, returns length of longest allele (returns 1 for snp)"""
    return max([len(r) for r in alleles])


def getNextIdx(rec, prev_indiv, prev_idx):
    """Using record sample array, find individual and haplotype indices for
    next sample. Will work for any ploidy. Returns -1's when all haplotypes
    have been iterated through"""
    if len(rec.samples[prev_indiv].alleles) > prev_idx + 1:
        return prev_indiv, prev_idx+1
    if len(rec.samples) > prev_indiv+1:
        return prev_indiv+1, 0
    return -1, -1


def checkRefAlign(vcf_r, fasta_ref, chrom):
    """Compares sequence from record to reference FASTA sequence"""
    vcf_seq = vcf_r.ref
    pos = vcf_r.pos-1
    fasta_seq = fasta_ref.fetch(chrom, pos, pos+len(vcf_seq))
    if vcf_seq != fasta_seq:
        raise Exception(("VCF bases and reference bases do not match.\n "
                        "VCF reference: %s\nFASTA reference: "
                        "%s") % (vcf_seq, fasta_seq))


def getRecordList(vcf_reader, region, chrom):
    """Returns list for use in subsampling from input file"""
    var_sites = vcf_reader.fetch(chrom, region.start, region.end)
    lst = []
    for rec in var_sites:
        lst.append(rec)
    return lst


def generateSequence(rec_list, ref_seq, fasta_ref,
                     region, chrom, indiv, idx, args):
    """Fetches variant sites from a given region, then outputs sequences
    from each individual with their correct variants. Will print sequence
    up to the variant site, then print correct variant. After last site,
    will output the rest of the reference sequence."""
    #var_sites = vcf_reader.fetch(chrom,region.start,region.end)
    fl = 0
    seq = ''
    prev_offset = 0

    for vcf_record in rec_list:
        issnp = checkRecordIsSnp(vcf_record)
        if not args.indel_flag and not issnp:
            continue

        pos_offset = vcf_record.pos - 1 - region.start
        for i in range(prev_offset, pos_offset-1):
            seq += ref_seq[i]
        checkRefAlign(vcf_record, fasta_ref, chrom)
        allele = vcf_record.samples[indiv].alleles[idx]
        if allele is None:
            raise Exception(("Individual %d at position %d is missing "
            "data") % (vcf_record.pos,indiv))
        if issnp:
            seq += vcf_record.samples[indiv].alleles[idx]
            prev_offset = pos_offset
        else:
            max_indel = getMaxAlleleLength(vcf_record.alleles)
            allele = vcf_record.samples[indiv].alleles[idx]
            for i in range(len(allele), max_indel):
                allele += '_'
            seq += allele
            indel_offset = len(vcf_record.ref)-1
            prev_offset = pos_offset+indel_offset

    for i in range(prev_offset, len(ref_seq)):
        seq += ref_seq[i]

    if args.trim_seq:
        return seq[:len(ref_seq)]
    return seq


def getHeader(record_count, chrom, region, oneidx=True):
    start = region.start
    end = region.end
    if oneidx:
        start += 1
        end += 1
    return '>'+str(record_count)+' '+chrom+' '+str(start)+':'+str(end)


def getFastaFilename(vcfname):
    for ext in ['.vcf.gz', '.vcf', '.bcf', 'vcf.bgz']:
        if ext in vcfname:
            if ext == '.vcf':
                raise Exception(('VCF filetype requires compression'
                                ' with bgzip and indexing with tabix'))
            return vcfname[:-1*len(ext)]+'.fasta', ext
    raise Exception('VCF filename %s has no valid extension' %
                    fasta_filename)


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


def vcf_to_seq(sys_args):
    parser = createParser()
    args = parser.parse_args(sys_args)
    logArgs(args)
    validateFiles(args)
    region_list = RegionList(args.genename, oneidx=args.gene_idx)
    fasta_filename, input_ext = getFastaFilename(args.vcfname)
    fasta_file = open(fasta_filename, 'w')

    vcf_reader = getVcfReader(args)
    first_el = next(vcf_reader)
    chrom = first_el.chrom

    fasta_ref = pysam.FastaFile(args.refname)
    record_count = 1
    logging.info('Total individuals: %d' % (len(first_el.samples)))
    logging.info('Total regions: %d' % (len(region_list.regions)))
    for region in region_list.regions:
        if region.chrom is not None:
            chrom = region.chrom
        rec_list = getRecordList(vcf_reader, region, chrom)
        logging.debug('Region %d to %d: %d variants' %
                      (region.start,region.end,len(rec_list)))
        ref_seq = fasta_ref.fetch(chrom, region.start, region.end)
        fasta_header = getHeader(record_count, chrom, region)
        fasta_file.write(fasta_header+'\n')

        indiv, idx = 0, 0
        while indiv != -1:
            seq = generateSequence(rec_list, ref_seq, fasta_ref,
                                   region, chrom, indiv, idx, args)
            fasta_file.write(seq+'\n')
            indiv, idx = getNextIdx(rec_list[0], indiv, idx)
        record_count += 1
    fasta_file.close()

if __name__ == "__main__":
    individualFunctionLogger()
    vcf_to_seq(sys.argv[1:])
