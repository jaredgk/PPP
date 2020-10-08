#!/usr/bin/env python
'''
    Create IMa input file from four-gamete filtered VCF files.

    An IM analysis requires an IM-formatted input file that contains multiple 
    loci, where each locus has sequence information for variant sites as 
    determined by multiple individuals (which are grouped by population).
    To produce this file, this script takes either a VCF file with a
    BED file that indicates loci to be used or a VCF file per locus. A model
    file will be used to split the samples in the VCF(s) into populations.

    For each locus, a header line is created that contains several pieces of
    information, including number of individuals per population at this locus,
    number of sites in sequences provided, mutation model, inheritance scalar,
    and mutation rate per year at locus (over locus, not per base pair). 
    Population sizes and sequence ordering are handled internally, the 
    inheritance scalar and mutation rate can be set via commandline.

    Loci provided as input must pass the four-gamete filtering criteria. 
    In addition, if filtering has not been previously done options are
    available to filter out indels (default), multiallelic sites, and 
    CpGs (with reference genome). Sites with missing data can either be
    filtered by dropping individuals missing data from a locus from 
    analysis at that locus, or by replacing the missing site with a 
    reference allele. 

    ###############
    Input Arguments
    ###############
    **--vcf** *<input_filename>*
        Filename for input VCF if using BED file with locus information
    **--vcfs** *<vcf_filename_1>*...*<vcf_filename_n>*
        One or multiple VCF input filenames where each file contains
        sequences for a single locus. A file with lines corresponding
        to filenames can be provided with **--vcfs** *@<vcf_filelist>*
    **--model-file** *<model_filename>*
        Filename of model file.
    **--model** *<model name>*
        If model file contains multiple models, use this argument
        to specify name of population to use.
    **--reference-fasta** *<reference_filename>*
        Filename for reference FASTA file. File can be uncompressed or
        bg-zipped, but must be indexed with faidx. When option is specified,
        default options are to include sequence in output loci but not filter
        for CpGs (use --parse-cpg)
    **--bed** *<bed_filename>*
        Filename for BED file specifying loci if only one VCF is provided. 
        Can be used with multiple VCFs if line count aligns, used for getting
        correct locus length.
    
    ################
    Output Arguments
    ################
    **--out** *<out_filename>*
        Output filename.
    
    #############
    Model Options
    #############
    **--mutrate** *<mutation rate>*
        Set mutation rate per base pair (default is 1e-9). This value is
        multiplied by locus length to get mutation rate per locus.
    **--inheritance-scalar** *<scalar>*
        Sets inheritance scalar for all loci. Default behavior is to 
        set scalar to 1 for non-X/Y/MT chromosomes, .75 for 'X' and
        'chrX', and .25 for 'y', 'chrY', 'MT', and 'chrMT'.

    #################
    Filtering Options
    #################
    **--remove-multiallele**
        Set all multiallelic sites to be reference. 
    **--drop-missing-sites** *<individual_count>*
        Drops all sites where more than 'individual_count' individuals
        are missing data. Default is -1 (no dropping), and 0 will 
        drop all sites missing data and replace them with the reference
        allele. 
    **--drop-missing-inds**
        If set, if an individual is missing data at a locus, that 
        individual will not be included at that locus and population
        counts for that locus will be adjusted.
    **--remove-cpg**
        Requires --reference-fasta. If set, will replace CpG sites with
        reference allele at site, setting them as invariant.

    #############
    Other Options
    #############
    **--zero-ho**
        Use if input BED file uses zero-based, half-open coordinates instead
        of one-based, closed.
    **--bed-column-index** *<start_col,end_col,chrom_col>*
        Comma-separated list of zero-based indexes of start, end, and 
        chromosome name columns in input BED file. Default value for
        traditionally structured BED is 1,2,0
    **--noseq**
        If set, and --reference-fasta is provided, will not output 
        invariant sites to IM file.

'''
import sys
import pysam
import argparse
import os
import logging
from random import sample
from collections import OrderedDict

import pgpipe.vcf_reader_func as vf
from pgpipe.logging_module import initLogger, logArgs
from pgpipe.genome_region import Region, RegionList
from pgpipe.parse_functions import defaultsDictForFunction, getConfigFilename, makeRequiredList, getArgsWithConfig
from pgpipe.model import Model, read_single_model
from pgpipe.misc import argprase_kwargs
#from tabix_wrapper import prepVcf

#Input: VCF file, reference sequence, region list (possibly .bed file)
#Output: Sequences with reference genome overlayed with VCF SNP calls

#sys.path.insert(0,os.path.abspath(os.path.join(os.pardir, 'andrew')))


class locus():
    #Name, pops, store for seqs, length, mut model, scalar, mutrate
    def __init__(self, gener, rec_list, popmodel, args):
        if args.inhet_sc is None:
            if gener.chrom in ['X','chrX']:
                self.inhet_sc = 0.75
            elif gener.chrom in ['Y','chrY','MT','chrMT']:
                self.inhet_sc = 0.25
            else:
                self.inhet_sc = 1
        else:
            self.inhet_sc = args.inhet_sc
        self.name = gener.chrom+':'+str(gener.start)+':'+str(gener.end)
        self.gene_len = gener.end - gener.start
        for rec in rec_list:
            self.gene_len += (getMaxAlleleLength(rec.alleles) - len(rec.ref))
        self.popmodel = popmodel
        self.popkeys = []
        self.seqs = {}
        for p in popmodel.pop_list:
            self.popkeys.append(p)
            self.seqs[p] = []
        self.mut_model = 'I'
        #if args.mut_model is not None:
        #    if args.mut_model not in ['I','H','S','J','IS']:
        #        raise Exception("%s is an invalid mutation model (must select from I,H,S,J,IS)" % (args.mut_model))
        #    self.mut_model = args.mut_model
        if args.refname is not None and args.printseq:
            self.seq_len = self.gene_len
        else:
            self.seq_len = len(rec_list)
        self.mutrate = self.gene_len * args.mutrate

    def addSeq(self,seq,pop):
        try:
            self.seqs[pop].append(seq)
        except IndexError:
            raise Exception("Pop %s is not in pop list" % pop)

    def printToFile(self, outf = sys.stdout):
        out_str = ''
        out_str += (self.name)+' '
        for p in self.popkeys:
            out_str += (str(len(self.seqs[p]))+' ')
        out_str += str(self.seq_len)+' '
        out_str += self.mut_model + ' '
        out_str += str(self.inhet_sc)+' '
        out_str += str("%.14f" % (self.mutrate))+'\n'
        outf.write(out_str)
        for p in self.popkeys:
            for seq in self.seqs[p]:
                outf.write(seq+'\n')

class outputBuffer():

    def __init__(self, popmodel, outfile=sys.stdout):
        self.header = None
        self.hold = True
        self.outfile = outfile
        self.popmodel = popmodel
        self.poptree = None
        self.loci = []

    def createHeader(self):
        self.header = {}
        self.header['title'] = 'Test IMa input'
        self.header['pops'] = [p for p in self.popmodel.pop_list]

    def writeHeader(self):
        self.outfile.write(self.header['title']+'\n')
        self.outfile.write(str(len(self.header['pops']))+'\n')
        self.outfile.write(' '.join(self.header['pops'])+'\n')
        if self.poptree is not None:
            self.outfile.write(self.poptree+'\n')
        self.outfile.write(str(len(self.loci)+'\n'))



def parseArguments(passed_arguments = []):
    parser = argparse.ArgumentParser(description=("Generates an IMa input "
                                     "file from a VCF file, a reference"
                                     " genome, a list of gene regions, "
                                     "and a population info file."),
                                     fromfile_prefix_chars="@")
    parser.add_argument("--vcf", dest="vcfname", help="Input VCF filename")
    parser.add_argument("--vcfs", dest="vcflist", metavar="VCF",nargs="+",
                        help=("Input VCF per loci"))
    parser.add_argument("--reference-fasta", dest="refname", help="Reference FASTA file, use this plus --remove-cpg to remove CpGs from loci")
    parser.add_argument("--bed", dest="genename", help=("If using a single "
                        "VCF, list of regions to generate loci from"))
    parser.add_argument("--bed-column-index", dest="gene_col",
                        help=("Comma-separated list of length 3 with 0-based"
                        " indexes of start, end, and chromosome data in input"
                        " BED file. Default for normal BED is 1,2,0"))
    parser.add_argument("--model-file", dest="popname", help=("Filename of "
                        "pop model file"))
    parser.add_argument("--model",dest="poptag",help=("If model file has "
                        "multiple models, use model with this name"))
    parser.add_argument("--zero-ho", dest="zeroho", action="store_true",
                        help="Region coordinates are zero-based, half-open")
    parser.add_argument("--zero-closed", dest="zeroclosed",action="store_true",
                        help="Region coordinates are zero-based, closed")
    parser.add_argument("--keep-indels", dest="indel_flag", action="store_true",
                        help="Include indels when reporting sequences")
    parser.add_argument("--remove-multiallele",dest="remove_multiallele",
                        action="store_true", help=("Remove sites with more "
                        "than two alleles"))
    #Need to figure out exactly how this option is working
    parser.add_argument("--drop-missing-sites",dest="remove_missing",
                        action="store_const",const=0,default=-1,help=("Will "
                        "use reference allele for sites where more than n "
                        "individuals are missing data. Missing data will "
                        "by default cause an error."))
    parser.add_argument("--drop-missing-inds",dest="drop_inds",action="store_true",
                        help="Drop individuals from loci if they are missing any data")
    parser.add_argument("--remove-cpg",dest="parsecpg",action="store_true",
                        help=("Will replace CpG sites with reference allele"
                        ", requires --reference-fasta"))
    parser.add_argument("--noseq",dest="printseq",action="store_false",
                        help=("Will only print variants when reference is "
                        "provided. Used so CpG filtering can be done if "
                        "invariant sites aren't desired in output"))
    parser.add_argument("--trim-to-ref-length", dest="trim_seq",
                        action="store_true",
                        help=("Trims sequences if indels cause them to be "
                        "longer than reference"))
    parser.add_argument("--out", dest="output_name", default="input.ima.u",
                        help= ("Optional name for output other than default"))
    parser.add_argument("--compress-vcf", dest="compress_flag",
                        action="store_true", help=("If input VCF is not "
                        "compressed, will compress and use zip search"))
    parser.add_argument("--conf", dest="config_name", help= ("Name of "
                        "file with configuration options"))
    parser.add_argument("--no-ref-check", dest="ref_check",
                        action="store_false", help=("Prevents exception "
                        "generated by mismatched reference alleles from "
                        "VCF file compared to reference"))
    parser.add_argument("--mutrate",dest="mutrate",type=float,default=1e-9,
                        help="Mutation rate per base pair (default is 1e-9)")
    parser.add_argument("--inheritance-scalar",dest="inhet_sc",type=float,
                        default=None,help=("Sets inheritance scalar for all "
                        "chromosomes (default is 1 for autosomes, .75 for X,"
                        " .25 for Y/MT)"))
    parser.add_argument("--output-fasta",dest="fasta",action="store_true",
                        help=("Output fasta file(s) instead of IMa input. "
                        "Only vaguely supported."))
    parser.add_argument("--out-prefix",dest="multi_out",type=str,
                        help=("If output is fasta, generate one file"
                        "per loci."))
    if passed_arguments:
        return vars(parser.parse_args(passed_arguments))
    else:
        return vars(parser.parse_args())


def checkArgs(args):
    if args.vcfname is None and args.vcflist is None:
        raise Exception("Must provide at least one VCF file to either --vcf or --vcfs")
    if args.vcfname is not None and args.vcflist is not None:
        raise Exception("Cannot use both arguments --vcf and --vcfs")
    if args.vcfname is not None and args.genename is None:
        raise Exception("If using --vcf, must provide a BED file with --bed")
    if args.refname is None and args.parsecpg:
        raise Exception("CpG parsing requires reference genome (--reference-fasta)")
    if args.popname is None and not args.fasta:
        raise Exception("IMa file requires model input with --model-file [filename]")


def validateFiles(args):
    """Validates that files provided to args all exist on users system"""
    for var in ['vcfname', 'refname', 'genename','popname']:
        f = vars(args)[var]
        if f is not None and not os.path.exists(f):
            raise ValueError('Filepath for %s not found at %s' %
                            (var, f))
    if args.vcflist is not None:
        for f in args.vcflist:
            if not os.path.exists(f):
                raise ValueError(('File not found at %s') % f)


def getMaxAlleleLength(alleles):
    """If an indel, returns length of longest allele (returns 1 for snp)"""
    return max([len(r) for r in alleles])


def checkRefAlign(vcf_recs, fasta_ref, chrom, ref_check):
    """Compares sequence from record to reference FASTA sequence"""
    for vcf_r in vcf_recs:
        vcf_seq = vcf_r.ref
        pos = vcf_r.pos-1
        try:
            fasta_seq = fasta_ref.fetch(vcf_r.chrom, pos, pos+len(vcf_seq))
        except KeyError:
            fasta_seq = fasta_ref.fetch(vf.flipChrom(vcf_r.chrom),pos,pos+len(vcf_seq))
        if vcf_seq.upper() != fasta_seq.upper():
            if ref_check:
                raise Exception(("VCF bases and reference bases do not match."
                        "\nVCF reference: %s\nFASTA reference: "
                        "%s\nPosition: %s") % (vcf_seq, fasta_seq, str(pos)))
            else:
                logging.warning("Bases at position %s do not match" %
                            str(pos))
                break


def generateSequence(rec_list, ref_seq, region, indiv, idx, args):
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
        #If reference included, add all bases between last and current variant
        pos_offset = vcf_record.pos - 1 - region.start
        if ref_seq is not None:
            for i in range(prev_offset, pos_offset):
##                print(vcf_record.pos,region.start,prev_offset,pos_offset,i)
                seq += ref_seq[i]

        allele = vcf_record.samples[indiv].alleles[idx]
        if allele is None and not args.N_if_missing:
            raise Exception(("Individual {} at position {} is missing "
                "data").format(indiv,vcf_record.pos))
        if issnp:
            #Place allele, move reference offset
            if allele is None and args.N_if_missing:
                seq += 'N' 
            else:
                seq += vcf_record.samples[indiv].alleles[idx]
            prev_offset = pos_offset+1
        else:
            #Find longest allele, pad others to its length
            #Offset by length of reference allele
            max_indel = getMaxAlleleLength(vcf_record.alleles)
            allele = vcf_record.samples[indiv].alleles[idx]
            if allele is None:
                allele = 'N'
            for i in range(len(allele), max_indel):
                allele += '_'
            seq += allele
            indel_offset = len(vcf_record.ref)
            prev_offset = pos_offset+indel_offset
    #Output remaining sequence
    if ref_seq is not None:
        for i in range(prev_offset, len(ref_seq)):
            seq += ref_seq[i]
        if args.trim_seq:
            return seq[:len(ref_seq)]
    return seq



def writeHeader(popmodel, loci_count, out_f, mutrate, header="Test IMa input",
                pop_string=None):
    out_f.write(header+'\n')
    out_f.write('#Mutation rate: '+str(mutrate)+'\n')
    out_f.write(str(popmodel.npop)+'\n')
    pops = ''
    for p in popmodel.pop_list:
        pops += (p+' ')
    out_f.write(pops+'\n')
    #Add default treestring when found in IMa3 source
    if pop_string is not None:
        out_f.write(pop_string+'\n')
    out_f.write(str(loci_count)+'\n')

def getLocusHeader(gener, popmodel, rec_list, mut_model="I", inhet_sc=None, mut_rate=1e-9,include_seq=False):
    if inhet_sc is None:
        if gener.chrom in ['X','chrX']:
            inhet_sc = 0.75
        elif gener.chrom in ['Y','chrY','MT','chrMT']:
            inhet_sc = 0.25
        else:
            inhet_sc = 1
    name = gener.chrom+':'+str(gener.start)+':'+str(gener.end)
    gene_len = gener.end-gener.start
    for rec in rec_list:
        gene_len += (getMaxAlleleLength(rec.alleles) - len(rec.ref))
    lh = name
    #for i in range(len(pop_data)):
    #    lh += ' '+str(len(pop_data[i][1]))
    for p in popmodel.pop_list:
        lh += ' '+str(2*popmodel.nind[p])
    if include_seq:
        lh += ' '+str(gene_len)
    else:
        lh += ' '+str(len(rec_list))
    if mut_model not in ['I','H','S','J','IS']:
        raise Exception("%s is an invalid mutation model (must select from I,H,S,J,IS)" % (mut_model))
    lh += ' '+mut_model
    lh += ' '+str(inhet_sc)
    mutlocus = mut_rate * gene_len
    lh += ' '+str("%.14f" % (mutlocus))
    return lh

def getMultiName(args):
    suffix = '.fasta' if args.fasta else '.u'
    i = 1
    while True:
        yield args.multi_out+'_region'+str(i)+suffix
        i += 1

def getOutputFilename(args):
    if args.output_name is not None:
        return args.output_name
    suffix = '.fasta' if args.fasta else '.u'
    va = args.vcfname.split('.')
    if len(va) == 1:
        return va+suffix
    ext_cut = -1
    if va[-1] == 'gz':
        ext_cut = -2
    return '.'.join(va[:ext_cut])+suffix

def hasMissingData(rec_list, indiv_idx):
    for rec in rec_list:
        if rec.samples[indiv_idx].alleles[0] is None:
            return True
    return False


def vcf_to_ima(**kwargs):
    """Returns an IMa input file given four-gamete filtered loci in one or 
    multiple VCFs.

    An IM analysis requires a file with variants for a set of individuals at
    multiple loci. These loci must be checked so that there are no 
    recombination events between SNPs, as determined by a four-gamete test.
    Variant input can be either a single VCF with a BED file listing regions
    that correspond to loci, or a single VCF per locus. Also required is a
    model file that lists what individuals in the VCF header correspond to
    which population. Filtering can be done in this step, for indels,
    multiallelic variants, sites with missing data, and CpGs (reference
    required). Mutation rate and inheritance scalars can also be set via
    commandline. 

    Parameters
    ----------
    --vcf : str
        Filename for VCF input file. 
    --vcfs : str (1+)
        Filenames for VCF locus files. Can also call '--vcfs @[filename.txt]'
        where filename.txt has paths for all desired VCF files.
    --reference-fasta : str
        Filename for FASTA reference file. File must contain entirety of 
        chromosomes specified. By default this will put invaraint sites 
        in the output loci. Required for CpG filtering.
    --model-file : str
        Filename for model file. This file contains population designations
        for individuals to be used in IM analysis. 
    --model : str
        If model file has more than one population, use model with this name
    --bed : str
        Filename for loci file with single VCF input. Requires columns with
        start, end, and chromosome of each locus.
    --bed-column-index : ints
        Comma-separated list of length 3 with 0-based indexes of start, end,
        and chromosome data in input BED file. Default for normal BED is
        "1,2,0".
    --zero-ho : bool
        If set, treats regions in BED file as 0-based, (h)alf (o)pen
        coordinates rather than 1-based, closed. 
    --out : str
        Name for output file. Defaults to 'input.ima.u"
    --remove-multiallele : bool
        If set, will set all multiallelic sites found as reference allele
    --drop-missing-inds : bool
        If set, if an individual at a locus is missing data, that individual
        will not be included at the locus.
    --remove-cpg : bool
        If set, will replace CpG sites with reference allele. Requires
        --reference-fasta.
    --noseq : bool
        If set and --reference-fasta is called, will only use for CpG
        checking and not output invariant sites in output file.



    Other Parameters
    ----------------
    --remove-multiallele : bool
        If set, will set all multiallelic sites found as reference allele
    --drop-missing-inds : bool
        If set, if an individual at a locus is missing data, that individual
        will not be included at the locus.
    --remove-cpg : bool
        If set, will replace CpG sites with reference allele. Requires
        --reference-fasta.
    --noseq : bool
        If set and --reference-fasta is called, will only use for CpG
        checking and not output invariant sites in output file.
    --trim-to-ref-length: bool, optional (False)
        If set, the sequences output will always match the length of the
        region they are found in. For example, a sequence with an insertion
        will cause the sequence to be an additional length of n-1, with n
        being the length of the insertion.
    --mutrate : float
        Use custom mutation rate for file (default is 1e-9)
    --inheritance-scalar : float
        Use fixed inheritance scalar for every output loci. Default is
        1 for autosomal chromosomes, .75 for X and .25 for Y/MT.


    Output
    ------
    IMa file
        Will be named either '--out' value or (vcfinput).ima.u.
        Contains variants in designated loci for IMa run.
    """
    #parser = createParser()
    #if len(sys_args) == 0:
    #    parser.print_help()
    #    sys.exit(1)

    #required_args = ['vcfname','popname']
    #args = getArgsWithConfig(parser,sys_args,required_args,'vcf_to_ima')
    if __name__ != "__main__":
        kwargs = argprase_kwargs(kwargs,parseArguments)
    args = argparse.Namespace(**kwargs)
    checkArgs(args)
    logArgs(args)
    #validateFiles(args)
    if args.multi_out is None:
        output_filename = getOutputFilename(args)
        output_file = open(output_filename, 'w')
    else:
        outnamegen = getMultiName(args)
    popmodel = None
    use_allpop = False
    if args.popname is not None:
        popmodel = read_single_model(args.popname,args.poptag)
    else:
        use_allpop = True

    filter_recs = (args.remove_multiallele or (args.remove_missing != -1) or not args.indel_flag or args.parsecpg)

    if args.vcfname is not None:
        single_file = True
    elif args.vcflist is not None:
        single_file = False
    else:
        raise Exception("VCF tag must be specified")
    if single_file:
        vn = args.vcfname
    else:
        vn = args.vcflist[0]
    vcf_reader = vf.VcfReader(vn,
                              compress_flag=args.compress_flag,
                              popmodel=popmodel,
                              use_allpop=use_allpop)
    logging.info('VCF file read')

    regions_provided = False
    if args.genename is not None:
        regions_provided = True
        region_list = RegionList(filename=args.genename, zeroho=args.zeroho,
                            zeroclosed=args.zeroclosed,
                            colstr=args.gene_col)
        logging.info('Region list read')
    fasta_ref = None
    if args.refname is not None:
        fasta_ref = pysam.FastaFile(args.refname)
    record_count = 1
    first_el = vcf_reader.info_rec

    logging.info('Total individuals: %d' % (len(vcf_reader.info_rec.samples)))
    if regions_provided:
        logging.info('Total regions: %d' % (len(region_list.regions)))
    total_regions = (len(region_list.regions) if regions_provided else len(args.vcflist))
    if not args.fasta:
        writeHeader(popmodel, total_regions, output_file, args.mutrate)
    if not single_file:
        vcf_reader.reader.close()
    for i in range(total_regions):
        #if regions_provided:
        #region = region_list.regions[i]
        if args.multi_out is not None:
            try:
                output_file.close()
            except:
                pass
            output_filename = next(outnamegen)
            output_file = open(output_filename,'w')
        if single_file:
            region = region_list.regions[i]
            rec_list = vcf_reader.getRecordList(region)
        else:
            vcf_reader = vf.VcfReader(args.vcflist[i],
                                      compress_flag=args.compress_flag,
                                      popmodel=popmodel,
                                      use_allpop=use_allpop)
            rec_list = vcf_reader.getRecordList()
            vcf_reader.reader.close()
            if regions_provided:
                region = region_list.regions[i]
            else:
                region = Region(rec_list[0].pos-1,rec_list[-1].pos,rec_list[0].chrom)
        if filter_recs:
            #Make this modify rec list in place
            t = vf.filterSites(rec_list,remove_cpg=args.parsecpg,remove_indels=(not args.indel_flag),remove_multiallele=args.remove_multiallele,remove_missing=args.remove_missing,inform_level=1,fasta_ref=fasta_ref)
            rec_list = t
        if len(rec_list) == 0:
            logging.warning(("Region %s has no variants "
                            "in VCF file") % (region.toStr()))
        logging.debug('Region %d to %d: %d variants' %
                      (region.start,region.end,len(rec_list)))
        ref_seq = None
        if fasta_ref is not None:
            try:
                ref_seq = fasta_ref.fetch(region.chrom, region.start, region.end)
            except KeyError:
                ref_seq = fasta_ref.fetch(vf.flipChrom(region.chrom),region.start,region.end)
            checkRefAlign(rec_list, fasta_ref, region.chrom, args.ref_check)
        if not args.fasta:
            temp_locus = locus(region, rec_list, popmodel, args)
            #reg_header = getLocusHeader(region, popmodel, rec_list,mut_rate=args.mutrate,inhet_sc=args.inhet_sc,include_seq=(fasta_ref is not None))
            #output_file.write(reg_header+'\n')
        popnum = 0
        #for p in popmodel.pop_list:
        for p in vcf_reader.popkeys.keys():
            for indiv,indiv_idx in enumerate(vcf_reader.popkeys[p]):
                if indiv_idx == -1:
                    continue
                for hap in range(len(first_el.samples[indiv_idx].alleles)):
                    if args.fasta:
                        seq = generateSequence(rec_list, ref_seq,
                                   region, indiv_idx, hap, args)
                        output_file.write('>'+first_el.samples[indiv_idx].name+'_'+str(hap)+'\n')
                        output_file.write(seq+'\n')
                        continue
                    hmd = hasMissingData(rec_list,indiv_idx)
                    if hmd:
                        if not args.drop_inds:
                            raise Exception("Individual %d from pop %s at loci %d, site %d is missing data" % (indiv_idx,p,record_count,i))
                        else:
                            continue

                    seq = generateSequence(rec_list, ref_seq,
                                   region, indiv_idx, hap, args)
                    seq_name = str(popnum)+':'+str(indiv)+':'+str(hap)
                    seq_name += ''.join([' ' for i in range(len(seq_name),10)])
                    allseq = seq_name + seq
                    temp_locus.addSeq(allseq,p)
                    #output_file.write(seq_name)

                    #output_file.write(seq+'\n')
                #indiv += 1
            popnum += 1
            indiv = 0
        temp_locus.printToFile(outf = output_file)
        record_count += 1
    output_file.close()

if __name__ == "__main__":
    #initLogger()
    vcf_to_ima(**parseArguments())
