import sys
import pysam
import argparse
import os.path
import logging
from logging_module import initLogger
from random import sample
from gene_region import Region, RegionList
import vcf_reader_func as vf
#from tabix_wrapper import prepVcf

#Input: VCF file, reference sequence, region list (possibly .bed file)
#Output: Sequences with reference genome overlayed with VCF SNP calls


def createParser():
    parser = argparse.ArgumentParser(description=("Generates sequences from"
                                     " samples from a VCF file, a reference"
                                     " genome, and a list of gene regions."))
    parser.add_argument("--vcf", dest="vcfname", help="Input VCF filename")
    parser.add_argument("--ref", dest="refname", help="Reference FASTA file")
    parser.add_argument("--rl", dest="genename",
                        help="Name of gene region file")
    parser.add_argument("--gr1", dest="gene_idx", action="store_true",
                        help="Gene Region list is 1 index based, not 0")
    parser.add_argument("--indels", dest="indel_flag", action="store_true",
                        help="Include indels when reporting sequences")
    parser.add_argument("--trim-to-ref-length", dest="trim_seq",
                        action="store_true",
                        help=("Trims sequences if indels cause them to be "
                              "longer than reference"))
    parser.add_argument("--output", dest="output_name", help= (
                        "Optional name for output other than default"))
    parser.add_argument("--gene-col", dest="gene_col", help= (
                        "Comma-separated list of columns for gene region "
                        " data, format is start/end if no chromosome "
                        " data, start/end/chrom if so"))
    parser.add_argument("--ext", dest="var_ext", help=(
                        "Format for variant file if filename doesn't "
                        "contain extension"))
    parser.add_argument("--compress-vcf", dest="compress_flag",
                        action="store_true", help=("If input VCF is not "
                        "compressed, will compress and use zip search"))
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
        issnp = vf.checkRecordIsSnp(vcf_record)
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


def getHeader(record_count, chrom, region, oneidx=False, halfopen=True):
    start = region.start
    end = region.end
    if oneidx:
        start += 1
        end += 1
    if not halfopen:
        start -= 1
    return '>'+str(record_count)+' '+chrom+' '+str(start)+':'+str(end)


def getFastaFilename(args):
    vcfname = args.vcfname
    for ext in ['vcf.gz', 'vcf', 'bcf', 'vcf.bgz']:
        #if ext in vcfname:
        if ext == vcfname[-1*len(ext):]:
            if args.output_name is None:
                return vcfname[:-1*len(ext)]+'fasta', ext
            else:
                return args.output_name, ext
    if args.var_ext is not None:
        ext = args.var_ext
        if args.output_name is None:
            raise Exception(("If filetype is not implicit by filename, "
                    "an output filename must be provided"))
            #return vcfname[:-1*len(ext)]+'fasta', ext
        else:
            return args.output_name, ext
    raise Exception('VCF filename %s has no valid extension' %
                    vcfname)


def vcf_to_seq(sys_args):
    """Returns a FASTA file with seqs from individuals in given gene regions

    Given an input VCF file, a reference FASTA file, and a list of gene
    regions, will output a FASTA file with sequence data for all individuals
    in the regions given. The reference FASTA file must be a full file from
    one or multiple chromosomes, starting at the first base. The gene
    region file must have start and end coordinates (half-open), with an
    optional column for chromosome data if the VCF input has multiple
    chromosomes.

    Parameters
    ----------
    --vcf : str
        Filename for VCF input file. If it does not end with extension
        'vcf(.gz)', a value for --ext must be provided.
    --ref : str
        Filename for FASTA reference file. This file can contain multiple
        chromosomes but must start from the first base, as there is currently
        no way to offset the sequences when pulling from a Region
    --rl : str
        Filename for gene region file. Requires columns for start and end
        coordinates, with option for chromosome. Additional data may be
        included, the columns with relevant data can be specified with the
        --gene-col option
    --indels : bool, optional
        If set, indels will be included in the output sequences
    --output : str, optional
        If set, the default output name of (inputprefix).fasta will be
        replaced with the given string
    --gene-col : str, optional
        Comma-separated string with two or three elements (chromosome is
        optional). If length is 2, elements are the indices for columns
        in the input gene region file corresponding to the start/end
        coordinates of a region. If length 3, the third element
        specifies the index of the chromosome column. Default is "1,2,0",
        to match column order in a BED file.
    --ext : str ['vcf','vcf.gz'], optional
        Required if VCF filename does not end with the typical extension.



    Other Parameters
    ----------------
    --gr1 : bool, optional (False)
        If set, indicates that the genome coordinate data is in base 1
    --trim-to-ref-length: bool, optional (False)
        If set, the sequences output will always match the length of the
        region they are found in. For example, a sequence with an insertion
        will cause the sequence to be an additional length of n-1, with n
        being the length of the insertion.
    --compress-vcf : bool, optional
        If set, will use bgzip and tabix to compress and index given VCF
        file
    --subsamp_list : str, optional
        Name of single-column file with names of individuals to subsample
        from input VCF file.
    --subsamp_num : int, optional
        Number of individuals to be randomly subsampled from VCF file

    Output
    ------
    FASTA file
        Will be named either '--output' value or (vcfinput).fasta.
        Contains full sequence for given regions and individuals with
        appropriate SNP/indel data included.
    """
    parser = createParser()
    if len(sys_args) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(sys_args)
    logArgs(args)
    validateFiles(args)
    fasta_filename, input_ext = getFastaFilename(args)
    fasta_file = open(fasta_filename, 'w')

    vcf_reader, uncompressed = vf.getVcfReader(args)
    first_el = next(vcf_reader)
    chrom = first_el.chrom
    #compressed = (input_ext != 'vcf')

    region_list = RegionList(filename=args.genename, oneidx=args.gene_idx,
                            colstr=args.gene_col, defaultchrom=chrom)

    fasta_ref = pysam.FastaFile(args.refname)
    record_count = 1
    prev_last_rec = first_el
    logging.info('Total individuals: %d' % (len(first_el.samples)))
    logging.info('Total regions: %d' % (len(region_list.regions)))
    for region in region_list.regions:
        if region.chrom is not None:
            chrom = region.chrom
        if not uncompressed:
            rec_list = vf.getRecordList(vcf_reader, region, chrom)
        else:
            rec_list, prev_last_rec = vf.getRecordListUnzipped(vcf_reader,
                                      region, chrom, prev_last_rec)
        if len(rec_list) == 0:
            logging.warning(("Region from %d to %d has no variants "
                            "in VCF file") % (region.start,region.end))
        logging.debug('Region %d to %d: %d variants' %
                      (region.start,region.end,len(rec_list)))
        ref_seq = fasta_ref.fetch(chrom, region.start, region.end)
        fasta_header = getHeader(record_count, chrom, region,
                                 oneidx = args.gene_idx)
        fasta_file.write(fasta_header+'\n')

        indiv, idx = 0, 0
        while indiv != -1:
            seq = generateSequence(rec_list, ref_seq, fasta_ref,
                                   region, chrom, indiv, idx, args)
            fasta_file.write(seq+'\n')
            indiv, idx = getNextIdx(first_el, indiv, idx)
        record_count += 1
    fasta_file.close()

if __name__ == "__main__":
    initLogger(filename="/home/jared/workspace/ppp/galaxy.log")
    vcf_to_seq(sys.argv[1:])
