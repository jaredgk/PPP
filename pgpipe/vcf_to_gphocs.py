#!/usr/bin/env python
'''
Generates an input sequence file for the G-Phocs program from a
vcf file and a fasta reference file.

G-Phocs can estimate the phylogenetic and demographic history of a
set of genomes,  each sampled at a large number of genomic regions or loci.

Gronau I, Hubisz MJ, Gulko B, Danko CG, Siepel A.   Bayesian inference of ancient
human demography from individual genome sequences.  Nature Genetics 43 1031-1034.   2011

https://github.com/gphocs-dev/G-PhoCS/blob/master/GPhoCS_Manual.pdf

##################
Required Arguments
##################

**--vcf** *<input_vcf_filename>*
    The name of the vcf file.  This can be a bgzipped vcf file. . 

**--model-file** *<model_file_name>*
    The name of a PPP model file. 

**--model** *<model_name>*
    The name of a model in the model file.  The treemix file to be
    generated will contain the allele counts for each SNP in each of
    the populations.  The treemix run will estimate the phylogeny
    for the populations in the model.

    
**--bed-file** *<BED_file_name>*
    The Bed file specifies the regions of the vcf file to be sampled.
    Each row of the BED file (region) correspondes to one locus in the
    G-Phocs sequence file.
    
    The BED file is a sorted UCSC-style bedfile containing chromosome locations of
    the SNPs to be included in the output files. The BED file has no header.
    The first column is the chromosome name (this must match the chromosome
    name in the vcf file).
    The second column is start position (0-based, open interval)
    The third column is end position (closed interval).
    Any other columns are ignored.

**--out** *<output file name>*
    Specifies the complete output filename.

**--reference** *<reference fasta file>*
    The reference genome fasta file is required in order to generate full
    sequences from the SNP data in the vcf file. 
    
#################
Optional Aguments
#################

**--diploid** *<True (default)/False>*
    By default G-Phocs works with a single sequence for each individual, where
    heterozygous positions are shown using IUPAC ambiguity codes.
    If this option is False, then only the first sequence of each
    individual is returned and heterozygous positions are not shown.

**--nloci** *<number of loci>*
    By default the output file will contain as many loci as there are regions
    in the BED file.  With this option,  the first nloci regions will be used.


#############
Example usage
#############
Example command-lines:

.. code-block:: bash

    vcf_to_gphocs.py -h

   
.. code-block:: bash

    vcf_to_gphocs.py --vcf pan_example.vcf.gz --reference pan_example_ref.fa --model-file panmodels.model --modelname 4Pop" --bed-file pan_example_regions.bed --outvcf_gphocs_test.out
    

'''

import sys
import os
import subprocess
import logging
import argparse
from pgpipe.logging_module import initLogger, logArgs
import pgpipe.vcf_bed_to_seq as vBs
from pgpipe.model import Model, read_model_file
from pgpipe.genome_region import Region
import pgpipe.vcf_reader_func as vr
from pgpipe.misc import argprase_kwargs


def make_gphocs_sequence_file(vcf,reference,BEDfile,ids,filename=None,diploid = True,nloci = None):
    """
        Generates returns a G-phocs sequence file (Gronau et al., 2011) to stdout

        should look like
        1000
        locus1 3 10
        one CCAGAGAGCT
        two CCAGAGAGCT
        three CCAGAGAGCT
        locus2 3 1000
        ...

        -vcf :is a vcf file, bgzipped or not, with one or more chromosomes
        -reference :is a fasta file with references sequences for chromosomes
            in the same order as for vcf
            chromosome name(s)s in the fasta file (string after '>' symbol up to eol or first space)
                must exactly match the chromosome name(s) in the vcf file
        -BEDfile :is a sorted UCSC-style bedfile containing chromosome locations
            There is no header
            first column is chromosome name (must match vcf file and fasta file)
            second column is start position (0-based, open interval)
            third column is end position (closed interval)
            other columns are ignored
        -ids :is either a list of individuals, or a PPP model, that specifies which
            sequences to include
            -if a list, then the contents are strings of individuals that are in the
            vcf file and that should be included in the gphocs sequence file
            - if a popmodel,  then the individuals in the model should be in the vcf file
        - filename :is the name of the file to contain the sequences
            if None, results are sent to stdout
        -diploid :specifies whether the sample is diploid (default).
            -if True then IUPAC ambiguity codes are used for heterozygous positions
            -if False, then only the first sequence of each individual is returned
        -nloci : number of loci to return,  if None returns all regions in the BEDfile 
            
    """
    #gphocs use ambiguity codes for diploid individuals, acgti and hetzcode are used by makehetzseq()
    def makehetzseq(seqs):
        ns = []
        for i in range(len(seqs[0])):
            if seqs[0][i] in acgti and seqs[1][i] in acgti:
                ns.append(hetzcode[acgti[seqs[0][i]]][acgti[seqs[1][i]]])
            else:
                ns.append('N')
        return ''.join(ns)
    acgti = {'A':0,'C':1,'G':2,'T':3}
    hetzcode = [['A','M','R','W'],['M','C','S','Y'],['R','S','G','K'],['W','Y','K','T']]

    if filename is None:
        file_handle = sys.stdout
    else:
        file_handle = open(filename, 'w')
        
    if isinstance(ids,list): # make a model to pass to get_model_sequences()
        # order of sequences in output gphocs sequence file will be the same as in ids
        idlist = ids
        popmodel = Model('gphocsmodel')
        popmodel.assign_pop('gphocspop',ids)  
    else:
        assert isinstance(ids,Model)
        popmodel = ids
        idlist = []
        # this loop over pop and ind mirrors that in vcf_BED_to_seqs.py, should put ids same order as sequences that get returned 
        for pop in popmodel.pop_list:
            for ind in popmodel.ind_dict[pop]:
                idlist.append(ind)

    # get number of lines in BEDfile
    c = 0
    with open(BEDfile) as infp:
        for line in infp:
            if line.strip():
                c += 1
    if nloci != None:
        if  c < nloci:
            raise Exception("nloci (%d) is greater than the number of intervals (%d) in %s"%(nloci,c,BEDfile))
    else:
        nloci = c

    file_handle.write("%d\n"%nloci)


    a = vBs.get_model_sequences(vcf=vcf,
                        popmodel=popmodel,fasta_reference=reference,
                        BED_filename = BEDfile,return_single = (diploid==False))
    nr = 0
    while True:
        try:
            regionstr,s = next(a)# should raise StopIteration if generator is exhausted
            if diploid == False:
                assert len(s) == len(idlist)
                file_handle.write("locus%d_%s %d %d\n"%(nr,regionstr,len(s),len(s[0])))
                for i,id in enumerate(idlist):
                    file_handle.write("%s %s\n"%(id,s[i]))
            else:
                assert len(s) == 2*len(idlist)
                file_handle.write("locus%d_%s %d %d\n"%(nr,regionstr,len(s)//2,len(s[0])))
                for i,id in enumerate(idlist):
                    seqmerged = makehetzseq([s[2*i].upper(),s[2*i+1].upper()])
                    file_handle.write("%s %s\n"%(id,seqmerged))
                
            nr += 1
            if nr == nloci:
                break
        except StopIteration:
            break   # end loop
    if nr < nloci:
        raise Exception("Number of loci requested (%d) less than number retrieved (%d)"%(nloci,nr))
    file_handle.close()
    return 

def gphocs_parser(passed_arguments=[]):
    '''gphocs Argument Parser - Assigns arguments from command line'''

    def parser_confirm_file ():
        '''Custom action to confirm file exists'''
        class customAction(argparse.Action):
            def __call__(self, parser, args, value, option_string=None):
                if not os.path.isfile(value):
                    raise IOError('%s not found' % value)
                setattr(args, self.dest, value)
        return customAction

    gphocs_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    gphocs_parser.add_argument('--vcf', help = "Defines the filename of the VCF", type = str, required = True, action = parser_confirm_file())
    gphocs_parser.add_argument('--reference', help = 'Defines the fasta referemce filename',required = True, type = str, action = parser_confirm_file())
    gphocs_parser.add_argument('--model-file', help = 'Defines the model filename',required = True, type = str, action = parser_confirm_file())
    gphocs_parser.add_argument('--modelname', help = 'Defines the model and the individual(s)/population(s) to include', required = True,type = str)
    gphocs_parser.add_argument("--bed-file",help="Defines the BED filename", type = str, required = True, action = parser_confirm_file())
    gphocs_parser.add_argument('--out', help = 'Defines the complete output filename', type = str)
    gphocs_parser.add_argument('--diploid',default = True,
                                help="specifies whether the sample is diploid."
                                "If True then IUPAC ambiguity codes are used for heterozygous positions."
                                "If False, then only the first sequence of each individual is returned.")
    gphocs_parser.add_argument("--nloci",type=int,help="Number of 'loci' (BED file regions) to return."
                               "if None returns all regions in the BED file")
    if passed_arguments:
        return vars(gphocs_parser.parse_args(passed_arguments))
    else:
        return vars(gphocs_parser.parse_args())


def run (**kwargs):
    # Grab gphocs arguments from command line
##    gphocs_args = gphocs_parser(passed_arguments)
    # Update kwargs with defaults
    if __name__ != "__main__":
        kwargs = argprase_kwargs(kwargs, gphocs_parser)
    # Assign arguments
    
    gphocs_args = argparse.Namespace(**kwargs)


    # Adds the arguments (i.e. parameters) to the log file
    logArgs(gphocs_args, func_name = 'make_gphocs_sequence_file')

    popmodels = read_model_file(gphocs_args.model_file)
    popmodel = popmodels[gphocs_args.modelname]
    make_gphocs_sequence_file(gphocs_args.vcf,gphocs_args.reference,gphocs_args.bed_file,
                              popmodel,filename=gphocs_args.out,diploid = gphocs_args.diploid,
                              nloci=gphocs_args.nloci)


if __name__ == "__main__":
    initLogger()
    run(**gphocs_parser())
    exit()
    debugargs = ['--vcf','..//jhtests//pan_example.vcf.gz',
            '--reference',"..//jhtests//pan_example_ref.fa",
            '--model-file',"..//jhtests//panmodels.model",
            '--modelname',"4Pop",
            '--bed-file',"..//jhtests//pan_example_regions.bed",
            '--out','..//jhtests//results//vcf_gphocs_test.out']
    run(debugargs)
    debugargs = ['--vcf','..//jhtests//pan_example2.vcf.gz','--reference',
                 "..//jhtests//chr22_pan_example2_ref.fa",
            '--model-file',"..//jhtests//panmodels.model",'--modelname',"4Pop",
            '--bed-file',"..//jhtests//pan_test_regions.bed",
            '--out','..//jhtests//results//vcf_gphocs_test2.out','--diploid','False']
    run(debugargs)
