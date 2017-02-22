import Bio
import sys
import pysam
import argparse





#Input: VCF file, reference sequence, region list (possibly .bed file)
#Output: Sequences with reference genome overlayed with VCF SNP calls

def checkRecordIsSnp(rec):
    if len(rec.ref) != 1:
        return False
    for allele in rec.alts:
        if len(allele) != 1:
            return False
    return True

def replaceVar(ref_seq,rec,offset):
    ref_base = rec.ref
	
def indivIdx(indiv):
    return indiv/2,indiv%2

def generateSequence(vcf_reader,ref_seq,region,chrom,indiv):

    var_sites = vcf_reader.fetch(chrom,region[0],region[1])
    fl = 0
    seq = ''
    prev_offset = 0

    while fl == 0:
        try:
            vcf_record = next(var_sites)
        except:
            fl = 1
            break

        if not checkRecordIsSnp(vcf_record):
            continue

        pos_offset = vcf_record.pos - 1 - region[0]
        for i in xrange(prev_offset,pos_offset-1):
            seq += ref_seq[i]

        idv,idx = indivIdx(indiv)
        seq += vcf_record.samples[idv].alleles[idx]
        prev_offset = pos_offset

    for i in xrange(prev_offset,len(ref_seq)):
        seq += ref_seq[i]

    return seq
        
def readGeneRegionList


parser = argparse.ArgumentParser(description=("Generates sequences from samples"
                                 "from a VCF file, a reference genome, and a"
                                 "list of gene regions."))
parser.add_argument("--vcf", dest="vcfname",help="Input VCF filename")
parser.add_argument("--ref",dest="refname",help="Reference FASTA file")
args = parser.parse_args()

region_list = [[176200,176250]]

vcf_reader = pysam.VariantFile(args.vcfname)
first_el = next(vcf_reader)
chrom = first_el.chrom
sample_size = len(first_el.samples)*2

fasta_ref = pysam.FastaFile(args.refname)

for region in region_list:
    ref_seq = fasta_ref.fetch(chrom,region[0],region[1])
    print ">TestHeader"
    for i in xrange(sample_size):
        seq = generateSequence(vcf_reader,ref_seq,region,chrom,i)
        print seq
    #var_sites = vcf_reader.fetch(chrom,region[0],region[1],reopen=True)

