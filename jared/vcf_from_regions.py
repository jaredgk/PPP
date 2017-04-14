import sys
import pysam
import argparse
import logging
from logging_module import initLogger
from gene_region import Region, RegionList
import vcf_reader_func as vf



def createParser():
    parser= argparse.ArgumentParser(description=("Given a range or a file of "
                                "ranges and a VCF file, will generate one or "
                                "more VCF files with variants only from the "
                                "region(s) specified. For regions, start and "
                                "end coordinates required (default is zero-"
                                "based, half-open intervals, chromosome "
                                "optional unless VCF has more than one"))
    parser.add_argument("vcfname", help="Input VCF name")
    region_group = parser.add_mutually_exclusive_group()
    region_group.add_argument("--r", dest="gene_str", help=("Semicolon "
                              "separated string, formatted [chromosome:"
                              "start:end] or [start:end]"))
    region_group.add_argument("--rl", dest="genename", help=("Tab-delimited"
                              " file with two (start-end) or three "
                              "(chromosome) columns"))
    parser.add_argument("--gr1", dest="gene_idx", action="store_true",
                        help="Gene Region list is 1 index based, not 0")
    parser.add_argument("--noindels", dest="indel_flag", action="store_false",
                        help="Don't include indels in ouput VCF file(s)")
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
                        "compressed, will compress and use fetch search"))
    parser.add_argument("--multi-out", dest="multi_out", action="store_true",
                        help="Produces multiple output VCFs instead of one")
    subsamp_group = parser.add_mutually_exclusive_group()
    subsamp_group.add_argument('--subsamp_list', dest="subsamp_fn",
                               help="List of sample names to be used")
    subsamp_group.add_argument('--subsamp_num', dest="subsamp_num",
                               help=("Number of individuals to be randomly "
                                     "subsampled from VCF file"))
    return parser

def logArgs(args):
    logging.info('Arguments for vcf_region_write:')
    for k in vars(args):
        logging.info('Argument %s: %s' % (k, vars(args)[k]))

def getOutputName(args):
    vcfname = args.vcfname
    if args.output_name is not None:
        return args.output_name
    for ext in ['vcf.gz','vcf']:
        offset = -1*len(ext)
        if ext == vcfname[offset:]:
            return args.vcfname[:offset]+'subregions.'+ext
    raise Exception(("Either --output needs a value or the vcf "))

def vcf_region_write(sys_args):
    parser = createParser()
    if len(sys_args) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(sys_args)
    logArgs(args)
    vcf_reader, uncompressed = vf.getVcfReader(args)
    header = vcf_reader.header
    first_el = next(vcf_reader)
    chrom = first_el.chrom
    output_name = getOutputName(args)
    vcf_out = pysam.VariantFile(output_name,'w',header=header)

    if args.gene_str is not None:
        region_list = RegionList(genestr=args.gene_str, oneidx=args.gene_idx,
                                 defaultchrom=chrom)
    elif args.genename is not None:
        region_list = RegionList(filename=args.genename,oneidx=args.gene_idx,
                                colstr=args.gene_col, defaultchrom=chrom)
    else:
        raise Exception(("No value provided for region filename or "
                         "single region"))
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
        for rec in rec_list:
            issnp = vf.checkRecordIsSnp(rec)
            if not args.indel_flag and not issnp:
                continue
            vcf_out.write(rec)

if __name__ == '__main__':
    initLogger()
    vcf_region_write(sys.argv[1:])
