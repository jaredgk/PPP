import sys
import pysam
import argparse
import logging
import os
#sys.path.insert(0,os.path.abspath(os.path.join(os.pardir, 'andrew')))
from pgpipe.logging_module import initLogger
from pgpipe.genome_region import Region, RegionList
import pgpipe.vcf_reader_func as vf
from pgpipe.model import Model, read_model_file
from pgpipe.misc import argprase_kwargs

def parseArguments(passed_arguments = []):
    parser= argparse.ArgumentParser(description=("Given a range or a file of "
                                "ranges and a VCF file, will generate one or "
                                "more VCF files with variants only from the "
                                "region(s) specified. For regions, start and "
                                "end coordinates required (default is zero-"
                                "based, half-open intervals, chromosome "
                                "optional unless VCF has more than one"))
    parser.add_argument("--vcf",dest="vcfname", help="Input VCF name")
    region_group = parser.add_mutually_exclusive_group()
    region_group.add_argument("--r", dest="gene_str", help=("Comma "
                              "separated string, formatted \"start,end,"
                              "chrom\" or \"start,end\""))
    region_group.add_argument("--bed", dest="genename", help=("Tab-delimited"
                              " file with two (start-end) or three "
                              "(chromosome) columns"))
    parser.add_argument('--zero-ho', dest="zeroho", action="store_true")
    parser.add_argument('--zero-closed', dest="zeroclosed", action="store_true")
    parser.add_argument("--out", dest="output_name", help= (
                        "Optional name for output other than default"))
    parser.add_argument("--out-prefix",dest="out_prefix",help=("Specifies"
                        " output sent to multiple files starting with "
                        "given prefix"))
    parser.add_argument("--bed-column-index", dest="gene_col", help= (
                        "Comma-separated list of columns for gene region "
                        " data, format is start/end if no chromosome "
                        " data, start/end/chrom if so"))
    parser.add_argument("--compress-vcf", dest="compress_flag",
                        action="store_true", help=("If input VCF is not "
                        "compressed, will compress and use fetch search"))
    parser.add_argument("--parsecpg", dest="refname")
    parser.add_argument("--compress-out", dest="compress", action="store_true")
    parser.add_argument("--remove-indels",dest="remove_indels", 
                        action="store_true", 
                        help=("Removes indels from output VCF files"))
    parser.add_argument("--remove-multi", dest="remove_multiallele",
                        action="store_true")
    parser.add_argument("--remove-missing", dest="remove_missing", default=-1, 
                        help=("Will filter out site if more than the given number of "
                        "individuals (not genotypes) are missing data. 0 removes sites"
                        " with any missing data, -1 (default) removes nothing"))
    parser.add_argument("--informative-count", dest="informative_count",
                        type=int, default=0)
    parser.add_argument("--tbi", dest="tabix_index", help="Path to bgzipped "
                        "file's index if name doesn't match VCF file")
    parser.add_argument("--remove-missing-inds",dest="remove_missing_inds",
                        action="store_true", help=("Will remove individuals "
                        "with missing data from a loci's VCF file"))
    parser.add_argument("--multi-start-num",type=int,default=0,help=("If "
                        "multiple output files, start numbering with this number"))
    subsamp_group = parser.add_mutually_exclusive_group()
    subsamp_group.add_argument('--subsamp-list', dest="subsamp_fn",
                               help="List of sample names to be used")
    subsamp_group.add_argument('--subsamp-num', dest="subsamp_num",
                               help=("Number of individuals to be randomly "
                                     "subsampled from VCF file"))
    subsamp_group.add_argument('--model-file',dest="modelname",help="Model file for selecting individuals for writing")
    parser.add_argument('--model',dest="poptag",help="If model file is used, will use model with this name")
    parser.add_argument("--forceempty",dest="forceempty",action="store_true",help=("Will create empty VCF if a region is empty rather than throw an error"))
    if passed_arguments:
        return vars(parser.parse_args(passed_arguments))
    else:
        return vars(parser.parse_args())
    #return parser

def logArgs(args):
    logging.info('Arguments for vcf_region_write:')
    for k in vars(args):
        logging.info('Argument %s: %s' % (k, vars(args)[k]))

def getOutputName(args):
    vcfname = args.vcfname
    #nocompress = args.nocompress
    compress_output = args.compress
    if args.output_name is not None:
        if compress_output and args.output_name[-3:] != '.gz':
            return args.output_name+'.gz'
        return args.output_name
        #might need check for .gz if shouldnt be compressed
    end_ext = 'vcf.gz'

    if not compress_output:
        end_ext = 'vcf'
    for ext in ['vcf.gz','vcf']:
        offset = -1*len(ext)
        if ext == vcfname[offset:]:
            return args.vcfname[:offset]+'subregions.'+end_ext
    raise Exception(("Either --output needs a value or the vcf input needs "
                     "an extension of vcf or vcf.gz"))

def getOutputPrefix(args):
    if args.out_prefix is not None:
        return args.out_prefix
    if args.output_name is not None:
        return args.output_name
    for ext in ['vcf.gz','vcf']:
        offset = -1*len(ext)
        if ext == args.vcfname[offset:]:
            return args.vcfname[:offset]
    return args.vcfname

def getMultiFileName(pref, rc, compress):
    #oc = (not uncompressed and not nocompress_flag)
    ext = '.vcf'
    if compress:
        ext += '.gz'
    return pref+str(rc)+ext

def writeRegion(args, vcf_reader, region, rc, filter_sites, remove_cpg, 
                fasta_ref,header,vcf_out=None):
    rec_list = vcf_reader.getRecordList(region)
    out_p = getOutputPrefix(args)
    if len(rec_list) == 0:
        if args.forceempty:
            logging.warning(("Region %s has no variants "
                        "in VCF file") % (region.toStr()))
        else:
            raise Exception("Region %s has no variants" % (region.toStr()))
    outname = getMultiFileName(out_p,rc,args.compress)
    if vcf_out is None:

        vcf_t = pysam.VariantFile(outname,'w',header=header) #check
    else:
        vcf_t = vcf_out
    if filter_sites:
        pass_list = vf.getPassSites(rec_list, remove_cpg=remove_cpg,
                        remove_indels=args.remove_indels,
                        remove_multiallele=args.remove_multiallele,
                        remove_missing=args.remove_missing,
                        inform_level=args.informative_count,
                        fasta_ref=fasta_ref)
    for i in range(len(rec_list)):
        if (not filter_sites) or pass_list[i]:
            vcf_t.write(rec_list[i])
    if vcf_out is None:
        vcf_t.close()

def writeFile(args, vcf_reader, filter_sites, remove_cpg, fasta_ref, header):
    outname = getOutputName(args)
    outfile = pysam.VariantFile(outname,'w',header=header)
    for rec in vcf_reader.reader:
        if (not filter_sites) or vf.checkRecordPass(rec,
                                 remove_cpg=remove_cpg,
                                 remove_indels=args.remove_indels,
                                 remove_multiallele=args.remove_multiallele,
                                 remove_missing=args.remove_missing,
                                 inform_level=args.informative_count,
                                 fasta_ref=fasta_ref):
            outfile.write(rec)


def vcf_region_write(**kwargs):
    """Returns a VCF file with variants from regions in a list

    Given a VCF file and gene region information (from either a file with one
    or more regions or a string with a single one), will output a subset of
    variants that occur in the given regions.

    Parameters
    ----------
    vcfname : str
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
        If set, the default output name of (vcfinput).subregion.vcf(.gz) will
        be replaced with the given string. If writing multiple files, this
        value will be the file prefix, with outputs being of form
        (output).region(k).vcf[.gz]
    --nocompress
        If set, output VCF will be uncompressed. This will prevent adding
        the .gz extension to the output filename, which is how pysam
        determines what compression to use while writing
    --gene-col : str, optional
        Comma-separated string with two or three elements (chromosome is
        optional). If length is 2, elements are the indices for columns
        in the input gene region file corresponding to the start/end
        coordinates of a region. If length 3, the third element
        specifies the index of the chromosome column. Default is "1,2,0",
        to match column order in a BED file.



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
    --subsamp-fn : str, optional
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

    #parser = createParser()
    #if len(sys_args) == 0:
    #    parser.print_help()
    #    sys.exit(1)
    #args = parser.parse_args(sys_args)
    if __name__ != "__main__":
        kwargs = argprase_kwargs(kwargs,parseArguments)
    args = argparse.Namespace(**kwargs)
    logArgs(args)
    popmodel = None
    if args.modelname is not None:
        popmodels = read_model_file(args.modelname)
        if len(popmodels) != 1:
            popmodel = popmodels[args.poptag]
        else:
            pp = list(popmodels.keys())
            popmodel = popmodels[pp[0]]
    vcf_reader = vf.VcfReader(args.vcfname,
                              compress_flag=args.compress_flag,
                              subsamp_num=args.subsamp_num,
                              subsamp_fn=args.subsamp_fn,
                              popmodel=popmodel,
                              index=args.tabix_index)
    logging.info('VCF file read')
    header = vcf_reader.reader.header
    #first_el = vcf_reader.prev_last_rec
    first_el = vcf_reader.info_rec
    chrom = first_el.chrom
    full_fileread = False
    vcf_out = None
    if args.out_prefix is None:
        output_name = getOutputName(args)
        vcf_out = pysam.VariantFile(output_name, 'w', header=header)

    if args.gene_str is not None:
        region_list = RegionList(genestr=args.gene_str,zeroho=args.zeroho,
                                 zeroclosed=args.zeroclosed)
    elif args.genename is not None:
        region_list = RegionList(filename=args.genename,zeroho=args.zeroho,
                                 zeroclosed=args.zeroclosed,
                                 colstr=args.gene_col,sortlist=False)
    else:
        full_fileread=True
    logging.info('Region read')

    fasta_ref = None
    remove_cpg = (args.refname is not None)
    filter_sites = ((args.refname is not None)
                   or args.remove_indels
                   or args.remove_multiallele
                   or (args.remove_missing != -1)
                   or (args.informative_count != 0))
    if args.refname is not None:
        fasta_ref = pysam.FastaFile(args.refname)
        remove_cpg = True

    logging.info('Total individuals: %d' % (len(first_el.samples)))

    if full_fileread:
        writeFile(args,vcf_reader,filter_sites,remove_cpg,fasta_ref,header)
    else:
        logging.info('Total regions: %d' % (len(region_list.regions)))
        for rc,region in enumerate(region_list.regions,start=args.multi_start_num):
            #if args.multi_out:
            writeRegion(args,vcf_reader,region,rc,filter_sites,remove_cpg,
                        fasta_ref,header,vcf_out)
                #continue
        #rec_list = vcf_reader.getRecordList(region)
        #if len(rec_list) == 0:
        #    if args.forceempty:
        #        logging.warning(("Region %s has no variants "
        #                    "in VCF file") % (region.toStr()))
        #    else:
        #        raise Exception("Region %s has no variants" % (region.toStr()))
        #if args.multi_out:
        #    try:
        #        vcf_out.close()
        #    except:
        #        pass
        #    outname = getMultiFileName(out_p, rc, args.compress)
        #    vcf_out = pysam.VariantFile(outname, 'w', header=header)
        #if filter_sites:
        #    pass_list = vf.getPassSites(rec_list, remove_cpg=remove_cpg,
        #                  remove_indels=args.remove_indels, 
        #                  remove_multiallele=args.remove_multiallele,
        #                  remove_missing=args.remove_missing,
        #                  inform_level=args.informative_count,
        #                  fasta_ref=fasta_ref)
        #for i in range(len(rec_list)):
            #make filter_sites an option
        #    if (not filter_sites) or pass_list[i]:
        #        vcf_out.write(rec_list[i])


if __name__ == '__main__':
    #initLogger()
    vcf_region_write(**parseArguments())
