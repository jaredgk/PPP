import sys
import os
import subprocess
import logging
import pgpipe.vcf_BED_to_seqs as vBs
from pgpipe.model import Model, read_model_file
from pgpipe.genome_region import Region
import pgpipe.vcf_reader_func as vr


def make_treemix_file(vcffile,popmodel, filename=None):
    """
        Makes a SNP data input file from a vcf file for the treemix program
            Pickrell JK, Pritchard JK: Inference of Population Splits and
            Mixtures from Genome-Wide Allele Frequency Data. PLoS Genet 2012, 8(11):e1002967.

        treemix assumes biallelic states,  sites with more than two alleles, or missing data, are skipped
        
        vcffile is a vcf, or bgzipped vcf file
        
        filename is the name of the file to contain the sequences
            if None, results are sent to stdout
        
        popmodel is an instance of Model 
    """

    if filename is None:
        file_handle = sys.stdout
    else:
        file_handle = open(filename, 'w')


    a = vr.VcfReader(vcffile,popmodel=popmodel['popmod1'])

    for pop in popmodel.pop_list:
        file_handle.write("%s\t"%pop)
    file_handle.write("\n")
    numinds = []
    for pop in popmodel.pop_list:
        numinds.append(len(popmodel.ind_dict[pop]))
    ploidy = 2
    while True:
        x = a.getNext()
        if type(x) == type(None):
            break
        if vr.checkRecordPass(x,remove_indels=True,remove_multiallele=True):
            buildstr = ""
            checkcount = [0,0]
            spacer = [',',' ']
            for pop in popmodel.pop_list:
                adic = vr.getAlleleCountDict(x,popmodel.ind_dict[pop])
                
                for ai in range(2):
                    allele = x.alleles[ai]
                    try:
                        count = adic[allele]
                        checkcount[ai] += count
                        buildstr += "%d%s"%(count,spacer[ai])
                    except KeyError:
                        buildstr += "0%s"%(spacer[ai])
            if checkcount[0] > 0 and checkcount[1] > 0 and sum(checkcount) == ploidy*sum(numinds):
                file_handle.write("%s\n"%buildstr)
    file_handle.close()



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
##                if seqs[0][i] != seqs[1][i]:
##                    print(seqs[0][i],seqs[1][i])
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
        
    if isinstance(ids,list): # make a model to pass to get_model_sequences_from_multiple_regions()
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


    a = vBs.get_model_sequences_from_multiple_regions(vcf_filename=vcf,
                        popmodel=popmodel,fasta_reference=reference,
                        BED_filename = BEDfile,return_one_sequence = (diploid==False),
                        useNifmissingdata=True)
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
    return 
        

            
# testing calls to   vBs.get_model_sequences_from_region()


# in  Pan_all_hicov_chr22_decrun_missingasref.vcf chromosome names are 22
vcf_filename = "Pan_all_hicov_chr22_decrun_missingasref.vcf"
vcf_filename = "Pan_all_hicov_chr22_decrun_missingasref.vcf.gz"
# in chr22_hg18.fa changed chromosome name to '22'
fasta_reference = "chr22_hg18.fa"
model_file = "panmodels.model"
BED_file = "pan_test.bed"


vcf_filename = "Pan_chr_21_22_test.vcf"
vcf_filename = "Pan_chr_21_22_test.vcf.gz"
fasta_reference= "twochr_test_ref.fa"
BED_file = "twochr_test.bed"

popmodels = read_model_file(model_file)
popmodel = popmodels['3Pop']

##regionstr = '22:20000000-20001000'
##s = vBs.get_model_sequences_from_region(vcf = vcf_filename,popmodel=popmodel,
##        seq_reference=fasta_reference,region=regionstr)
##print('1',len(s),len(s[0]),s[0][0:20])
##seqref = s[0]
##s = vBs.get_model_sequences_from_region(vcf = vcf_filename,popmodel=popmodel,
##        seq_reference=seqref,region=regionstr)
##print('2',len(s),len(s[0]),s[0][0:20])
##
##chrname = regionstr[:regionstr.find(':')]
##[start,end] = map(int,regionstr[regionstr.find(':')+1:].split(sep='-'))
##print('z',chrname,start,end)
##region = Region(start-1,end,chrname)
##s = vBs.get_model_sequences_from_region(vcf = vcf_filename,popmodel=popmodel,
##        seq_reference=fasta_reference,region=region)
##print('3',len(s),len(s[0]),s[0][0:20])
##


# testing calls to   vBs.get_model_sequences_from_multiple_regions()
a = vBs.get_model_sequences_from_multiple_regions(vcf_filename=vcf_filename,
                        popmodel=popmodel,fasta_reference=fasta_reference,
                        BED_filename = BED_file,useNifmissingdata=True)
while True:
    try:
        s = next(a)# should raise StopIteration if generator is exhausted
        print(len(s),len(s[0]))
    except StopIteration:
        break   # end loop 


# testing calls to make_gphocs_sequence_file()
idlist =["Pan_troglodytes_schweinfurthii-A911_Kidongo","Pan_troglodytes_schweinfurthii-A912_Nakuu","Pan_troglodytes_troglodytes-A957_Vaillant","Pan_troglodytes_troglodytes-A958_Doris","Pan_troglodytes_troglodytes-A959_Julie","Pan_troglodytes_verus-9668_Bosco","Pan_troglodytes_verus-A956_Jimmie","Pan_troglodytes_verus-Clint","Pan_troglodytes_verus-X00100_Koby","Pan_paniscus-A917_Dzeeta","Pan_paniscus-A918_Hermien","Pan_paniscus-A919_Desmond"]

idlist = ["Pan_paniscus-A919_Desmond","Pan_troglodytes_schweinfurthii-A911_Kidongo","Pan_troglodytes_schweinfurthii-A912_Nakuu","Pan_troglodytes_troglodytes-A957_Vaillant","Pan_troglodytes_troglodytes-A958_Doris","Pan_troglodytes_troglodytes-A959_Julie","Pan_troglodytes_verus-9668_Bosco","Pan_troglodytes_verus-A956_Jimmie","Pan_troglodytes_verus-Clint","Pan_troglodytes_verus-X00100_Koby","Pan_paniscus-A917_Dzeeta","Pan_paniscus-A918_Hermien"]
##make_gphocs_sequence_file(vcf_filename,fasta_reference,BED_file,idlist,
##        filename="gphocstest.out",diploid = True)

make_gphocs_sequence_file(vcf_filename,fasta_reference,BED_file,popmodel,
        filename="gphocstest.out",diploid = True,nloci = 2)
