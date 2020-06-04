import sys
import os
import subprocess
import logging
import pgpipe.vcf_BED_to_seqs as vBs
from pgpipe.model import Model, read_model_file
from pgpipe.genome_region import Region

# in  Pan_all_hicov_chr22_decrun_missingasref.vcf chromosome names are 22
vcf_filename = "Pan_all_hicov_chr22_decrun_missingasref.vcf"
vcf_filename = "Pan_all_hicov_chr22_decrun_missingasref.vcf.gz"
# in chr22_hg18.fa changed chromosome name to '22'
fasta_reference = "chr22_hg18.fa"
model_file = "panmodels.model"
BED_file = "pan_test.bed"


vcf_filename = "Pan_chr_21_22_test.vcf"
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

a = vBs.get_model_sequences_from_multiple_regions(vcf_filename=vcf_filename,
                        popmodel=popmodel,fasta_reference=fasta_reference,
                        BED_filename = BED_file,useNifmissingdata=True)
while True:
    try:
        s = next(a)# should raise StopIteration if generator is exhausted
        print(len(s),len(s[0]))
    except StopIteration:
        break   # end loop 
