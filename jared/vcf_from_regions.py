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
    region_group.add_argument("--r", dest="gene_str", help=("Comma "
                              "separated string, formatted \"start,end,"
                              "chrom\" or \"start,end\""))
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
    subsamp_group.add_argument('--subsamp-list', dest="subsamp_fn",
                               help="List of sample names to be used")
    subsamp_group.add_argument('--subsamp-num', dest="subsamp_num",
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
    """Returns a VCF file with variants from regions in a list

    Given a VCF file and gene region information (from either a file with one
    or more regions or a string with a single one), will output a subset of
    variants that occur in the given regions.

    Parameters
    ----------
    --vcf : str
        Filename for VCF input file. If it does not end with extension
        'vcf(.gz)', a value for --ext must be provided.
    --r : str, optional
        Comma-separates string with either two or three elements, the first
        two corresponding to the start and end coordinates of the target
        region. If there is a third element, it is taken as chromosome data.
        Must use either this or --rl flag. (Not both)
    --rl : str, optional
        Filename for gene region file. Requires columns for start and end
        coordinates, with option for chromosome. Additional data may be
        included, the columns with relevant data can be specified with the
        --gene-col option. Must use either this or --r flag. (Not both)
    --noindels : bool, optional
        If set, indels will not be included in the output VCF
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
    --compress-vcf : bool, optional
        If set, will use bgzip and tabix to compress and index given VCF
        file
    --multi-out: bool, optional
        If set, each region will be written to a separate VCF file, tagged
        with (prefix).region[n].vcf(.gz)
    --subsamp-list : str, optional
        Name of single-column file with names of individuals to subsample
        from input VCF file.
    --subsamp-num : int, optional
        Number of individuals to be randomly subsampled from VCF file

    output
    ------
    VCF file(s)
        Will be named either '--output' value or (vcfinput).subregion.vcf(.gz)
        Contains variants from regions specified in region list from input
        VCF.


    """
    parser = createParser()
    if len(sys_args) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(sys_args)
    logArgs(args)
    vcf_reader, uncompressed = vf.getVcfReader(args)
    logging.info('VCF file read')
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
    logging.info('Region read')
    prev_last_rec = first_el
    logging.info('Total individuals: %d' % (len(first_el.samples)))
    logging.info('Total regions: %d' % (len(region_list.regions)))
    for region in region_list.regions:
        if not uncompressed:
            rec_list = vf.getRecordList(vcf_reader, region)
        else:
            rec_list, prev_last_rec = vf.getRecordListUnzipped(vcf_reader,
                            region, chrom, prev_last_rec)
        if len(rec_list) == 0:
            logging.warning(("Region from %d to %d has no variants "
                            "in VCF file") % (region.start,region.end))
        for rec in rec_list:
            issnp = vf.checkRecordIsSnp(rec)
            if not args.indel_flag and not issnp:
                continue
            vcf_out.write(rec)

if __name__ == '__main__':
    initLogger()
    vcf_region_write(sys.argv[1:])
