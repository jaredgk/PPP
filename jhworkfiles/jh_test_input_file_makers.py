import sys
import os
import subprocess
import logging
import pgpipe.vcf_BED_to_seqs as vBs
from pgpipe.model import Model, read_model_file
from pgpipe.genome_region import Region
import pgpipe.vcf_reader_func as vr

        

            
# testing calls to   vBs.get_model_sequences_from_region()


# in  Pan_all_hicov_chr22_decrun_missingasref.vcf chromosome names are 22
##vcf_filename = "Pan_all_hicov_chr22_decrun_missingasref.vcf"
##vcf_filename = "Pan_all_hicov_chr22_decrun_missingasref.vcf.gz"
### in chr22_hg18.fa changed chromosome name to '22'
##fasta_reference = "chr22_hg18.fa"
##model_file = "panmodels.model"
##BED_file = "pan_test.bed"
##
##
##vcf_filename = "Pan_chr_21_22_test.vcf"
##vcf_filename = "Pan_chr_21_22_test.vcf.gz"
##fasta_reference= "twochr_test_ref.fa"
##BED_file = "twochr_test.bed"
##
##popmodels = read_model_file(model_file)
##popmodel = popmodels['3Pop']

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
##a = vBs.get_model_sequences_from_multiple_regions(vcf_filename=vcf_filename,
##                        popmodel=popmodel,fasta_reference=fasta_reference,
##                        BED_filename = BED_file,useNifmissingdata=True)
##while True:
##    try:
##        s = next(a)# should raise StopIteration if generator is exhausted
##        print(len(s),len(s[0]))
##    except StopIteration:
##        break   # end loop 
##

# testing calls to make_gphocs_sequence_file()
##idlist =["Pan_troglodytes_schweinfurthii-A911_Kidongo","Pan_troglodytes_schweinfurthii-A912_Nakuu","Pan_troglodytes_troglodytes-A957_Vaillant","Pan_troglodytes_troglodytes-A958_Doris","Pan_troglodytes_troglodytes-A959_Julie","Pan_troglodytes_verus-9668_Bosco","Pan_troglodytes_verus-A956_Jimmie","Pan_troglodytes_verus-Clint","Pan_troglodytes_verus-X00100_Koby","Pan_paniscus-A917_Dzeeta","Pan_paniscus-A918_Hermien","Pan_paniscus-A919_Desmond"]
##
##idlist = ["Pan_paniscus-A919_Desmond","Pan_troglodytes_schweinfurthii-A911_Kidongo","Pan_troglodytes_schweinfurthii-A912_Nakuu","Pan_troglodytes_troglodytes-A957_Vaillant","Pan_troglodytes_troglodytes-A958_Doris","Pan_troglodytes_troglodytes-A959_Julie","Pan_troglodytes_verus-9668_Bosco","Pan_troglodytes_verus-A956_Jimmie","Pan_troglodytes_verus-Clint","Pan_troglodytes_verus-X00100_Koby","Pan_paniscus-A917_Dzeeta","Pan_paniscus-A918_Hermien"]
##make_gphocs_sequence_file(vcf_filename,fasta_reference,BED_file,idlist,
##        filename="gphocstest.out",diploid = True)
##import
##make_gphocs_sequence_file(vcf_filename,fasta_reference,BED_file,popmodel,
##        filename="gphocstest.out",diploid = True,nloci = 2)

##vcf_filename = "Pan_chr_21_22_test.vcf"
##vcf_filename = "Pan_chr_21_22_test.vcf.gz"
##
##model_file = "panmodels.model"
##popmodels = read_model_file(model_file)
##popmodel = popmodels['4Pop']
####make_treemix_file(vcf_filename,popmodel,"treemixtesta.out")
##
##make_treemix_file(vcf_filename,popmodel,"treemixtest.out")
##BED_file = "twochr_test.bed"
##make_treemix_file(vcf_filename,popmodel,"treemixtestbed.out",BEDfilename=BED_file)
##

import pgpipe.vcf_to_sfs as vcfsfs



vcf_filename = "Pan_chr_21_22_test.vcf"
vcf_filename = "Pan_chr_21_22_test.vcf.gz"
##vcf_filename = "temp.vcf" 



model_file = "panmodels.model"
popmodels = read_model_file(model_file)
popmodel = popmodels['4Pop']
##for pop in popmodel.pop_list:
##    print(2*len(popmodel.ind_dict[pop]))
BED_file = "twochr_test.bed"
BED_file = None

##downsamplesizes=[5,6,4,3]
##downsamplesizes=[5,2,4]
##downsamplesizes = None

##folded = True
##folded = False
##randomsnpprop = 0.1
randomsnpprop = None
seed = 127

##vcf_filename = "Pan_all_hicov_chr22_decrun_missingasref.vcf.gz"
##vcf_filename = "temp.vcf"
##altreference= "chr22_hg18.fa"
altreference = None

##vcf_filename = "sfs_test.vcf"
##model_file = "pantest.model"
##popmodels = read_model_file(model_file)
##popmodel = popmodels['2Pop']

downsamplesizes = None

folded = True
##folded = False


sfs = vcfsfs.build_sfs(vcf_filename,popmodel,BEDfilename=BED_file,altreference = altreference,folded = folded,downsamplesizes = downsamplesizes,randomsnpprop = randomsnpprop, seed = seed,makeint=True)
print(popmodel.pop_list)
print(sfs.shape,sfs.sum())
print(sfs)
##print(vcfsfs.reducesfsdims(sfs,popmodel,['Troglodytes','Schweinfurthii']))
rsfs = vcfsfs.reducesfsdims(sfs,popmodel,['Troglodytes','Schweinfurthii','verus'])
print(rsfs.shape)
print(rsfs)
    
