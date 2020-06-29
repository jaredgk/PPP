#!/usr/bin/env python
''' 

    This function checks to make sure that a locus has enough sites to be 
    considered informative in either the four-gamete test or an IM run. 
    Given a BED file and a VCF file, informative_loci_filter will find regions 
    in the VCF that have a specified number of variant sites.

    .. image:: ../../PPP_assets/PPP_Loci_Filter.png
        :width: 100 %
        :align: center

    In this illustration of the locus filtering process, locus_0004 is
    removed due to only having three variant sites (highlighted in green)
    when the threshold is set to four.

    Because many variants are not considered useful in these situations, filters 
    are provided for removing sites with missing data, non-biallelic sites, 
    indels, CpGs, and singletons from determining if there are a sufficient 
    number of sites in the region. Output is either a BED file of a set number 
    of random regions that pass the criteria, or a file with all regions that pass. 
    It can also remove regions that are below a minimum specified length. If a 
    model file is specified, only the individuals in the selected population will 
    be considered for singleton and missing data filters.

    ###############
    Input Arguments
    ###############
    **--vcf** *<vcf_name>*
        Name of input VCF file
    **--bed** *<bed_file>*
        Name of input BED file

    ################
    Output Arguments
    ################
    **--out** *<out_name>*
        Name of output BED file
    **--randcount** *<number_of_regions>*
        If set, will output set number of regions randomly selected
        from those that pass the criteria. Default behavior is to output
        all regions that match criteria.

    ###################
    Filtering Arguments
    ###################
    **--remove-multi**
        Do not count tri-allelic+ sites toward number of valid variants in
        a region
    **--remove-missing** *<max_misscount>*
        Do not count sites with more than max_misscount missing individuals.
        Default of -1 indicates all sites are included, 0 indicates sites
        with any missing data are not counted..
    **--remove-indels**
        Do not count indels toward number of valid variants in a region
     **--parsecpg** *<fasta_reference_filename>*
        Optional argument that if set, will detect whether or not variants
        are CpGs. A check is made to make sure the positions in the FASTA
        line up with the correct variant reference allele.   
    **--informative-count** *<minimum_allele_count>*
        Minimum number of haplotypes with both alleles at a site. Default
        is 2, meaning there must be two of each of reference and alternate
        allele in the target individuals. Can be set to 1 to filter 
        out invariant sites.
    **--minsites** *<min_sites>*
        Minimum number of variants required for a region to pass the filtering
        criteria. Variants that match specified arguments will not be counted
        towards this total.
    **--min-length** *<min_length>*
        Minimum base length of region for region to be considered.

    ################
    Region Arguments
    ################
    **--bed-column-index** *<start_idx>*,*<end_idx>*,*<chrom_idx>*
        Comma-separated string of the zero-based indices of the start, end,
        and chromosome columns in the input file, so the file doesn't need to
        be reformatted. Default for a regular BED file is 1,2,0.
    **--zero-ho**
        If set, indicates input BED regions are on zero-based, half-open
        coordinate system, as opposed to one-based, closed. For example, the
        first million bases on a chromosome would be:
            Zero-based, half-open: 0,1000000
            One-based, closed:     1,1000000
    **--pad** *<pad_length>*
        If set, regions in input file will be extended by pad_length bases
        on both sides.
    **--keep-full-line**
        If set, regions output will be the same line as was present in the 
        input file. Default behavior is to output start/end/chrom columns
        in that order, without any other data. 

    ###############
    Model Arguments
    ###############
    **--model-file** *<model_filename>*
        Name of model file that contains individuals to be considered
        for filtering. 
    **--model** *<model_name>*
        If model file contains more than one model, name of model
        to be used.

    ###############
    Other Arguments
    ###############
    **--no-sorting**
        Will output regions in order they were in input file.
        Default behavior is to sort regions before filtering.
    **--tbi** *<tabix_index>*
        If input VCF is compressed and tabix file is not default,
        provide the tabix filename here.
    **--no-xy**
        Removes regions on X/Y chromosomes from consideration.
'''

import sys
import pysam
import argparse
import os
from random import shuffle

from pgpipe.vcf_reader_func import getPassSites, VcfReader
from pgpipe.genome_region import Region, RegionList
from pgpipe.logging_module import initLogger
from pgpipe.model import Model, read_single_model

def createParser():
    parser = argparse.ArgumentParser(description=("Given a VCF file and a"
            " list of intervals (BED style), will output regions that have"
            " over a certain number of qualifying variants in the region"))
    parser.add_argument("--vcf", dest="vcfname", help=("Input VCF name"))
    parser.add_argument("--bed", dest="bedname", help=("BED filename "
                        "with regions for inspection"))
    parser.add_argument("--bed-column-index",dest="gene_col",help=("Three "
                        "comma-separated integers indicating columns for "
                        "start,end,and chromosme in BED file (default 1,2,0"))
    parser.add_argument("--out",dest="outname",help=("Output filename,"
                        "default is stdout"))
    parser.add_argument("--zero-ho", dest="zeroho", action="store_true",
                        help="BED input in zero-based, half-open format")
    parser.add_argument("--zero-closed", dest="zeroclosed",action="store_true")
    parser.add_argument("--parsecpg", dest="refname",help=("If filtering "
                        "for CpGs, provide reference filename here"))
    parser.add_argument("--remove-indels", dest="remove_indels", 
                        action="store_true", help=("Do not count indels "
                        "toward informative site count"))
    parser.add_argument("--remove-multi", dest="remove_multiallele",
                        action="store_true",help=("Do not count multiallelic"
                        "sites toward informative site count"))
    parser.add_argument("--remove-missing", dest="remove_missing", 
                        default=-1, type=int, help=("Will filter out site if"
                        " more than the given number of individuals "
                        "(not genotypes) are missing data. 0 removes sites "
                        "with any missing data, -1 (default) removes nothing"))
    parser.add_argument("--informative-count", dest="informative_count", 
                        type=int, default=2, help=("Minimum allele count of"
                        "sites counted as informative. Used so sites are "
                        "considered valid in four-gamete filtering."))
    parser.add_argument("--minsites", dest="minsites", default=3, 
                        help=("Regions with at least this many variants "
                        "passing filters will be output"))
    parser.add_argument("--tbi", dest="tabix_index", help=("Path to bgzipped "
                        "file's index if name doesn't match VCF file"))
    parser.add_argument("--randcount",dest="randcount",type=int,default=-1,
                        help=("If set, will randomly draw from input region "
                        "file until randcount # of passing BED regions "
                        "are found"))
    parser.add_argument("--no-xy",dest="filter_xy",action="store_true",
                        help="Remove X/Y chromosomes from valid regions")
    parser.add_argument("--min-length",dest="min_length",type=int,default=1000,
                        help="Minimum length of valid regions")
    parser.add_argument("--model-file",dest="modelname",help=("Model file for "
                        "selecting samples from VCF"))
    parser.add_argument("--model",dest="poptag",help="Name of pop if model "
                        "has more than one")
    parser.add_argument("--keep-full-line",dest="keep_full_line",
                        action="store_true",help=("Output BED line from "
                        "original input BED file"))
    parser.add_argument("--no-sorting",dest="sort_lists",action="store_false",
                        help=("Guarantees output regions are in same order"
                        "as input file"))
    return parser


def filter_bed_regions(sys_args):
    #parser = argparse.parse_args(sys_args)
    parser = createParser()
    args = parser.parse_args(sys_args)
    if args.outname is not None:
        outf = open(args.outname,'w')
    else:
        outf = sys.stdout

    popmodel = None
    if args.modelname is not None:
        popmodel = read_single_model(args.modelname,args.poptag)
        #popmodels = read_model_file(args.modelname)
        #if len(popmodels) != 1:
        #    popmodel = popmodels[args.poptag]
        #else:
        #    pp = list(popmodels.keys())
        #    popmodel = popmodels[pp[0]]
    
    vcf_reader = VcfReader(args.vcfname,index=args.tabix_index,popmodel=popmodel)
    fasta_seq = None
    if args.refname is not None:
        fasta_seq = pysam.FastaFile(args.refname)

    #regions = RegionList(filename=args.bedname,zeroho=args.zeroho,zeroclosed=args.zeroclosed,sortlist=(not args.randcoun))
    randomize = False
    if args.randcount != -1:
        randomize = True

    regions = RegionList(filename=args.bedname,zeroho=args.zeroho,
                         zeroclosed=args.zeroclosed,sortlist=args.sort_lists,
                         colstr=args.gene_col,
                         keep_full_line=args.keep_full_line)
    if args.filter_xy:
        regions.filterOutXY()
    remove_cpg = (True if args.refname is not None else False)
    idx_list = [i for i in range(len(regions.regions))]
    if args.randcount != -1:
        shuffle(idx_list)
    regions_output = 0
    idx_output = []
    #for region in regions.regions:
    for i in idx_list:
        region = regions.regions[i]
        if len(region) < args.min_length:
            continue
        rec_list = vcf_reader.getRecordList(region)
        pass_list = getPassSites(rec_list, remove_cpg=remove_cpg,
                    remove_indels=args.remove_indels,
                    remove_multiallele=args.remove_multiallele,
                    remove_missing=args.remove_missing,
                    inform_level=args.informative_count,
                    fasta_ref=fasta_seq)
        if pass_list.count(True) >= int(args.minsites):
            #print (region.toStr(sep='\t'))
            idx_output.append(i)
            regions_output += 1
        if args.randcount != -1 and regions_output == args.randcount:
            break

    if args.randcount != -1 and regions_output != args.randcount:
        sys.stderr.write("Only %d of %d regions found\n"%(regions_output,args.randcount))
        #exit(1)

    if randomize:
        idx_output.sort()
    
    if regions.header is not None:
        outf.write(regions.header+'\n')
    for i in idx_output:
        outf.write(regions.regions[i].toStr(sep='\t')+'\n')


if __name__ == '__main__':
    initLogger()
    filter_bed_regions(sys.argv[1:])
