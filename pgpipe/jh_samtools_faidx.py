##print("hello world")
import sys
import argparse
import pysam
import gzip

import os
import sys
import json
import subprocess
import argparse
import logging
import itertools
import copy

##from pgpipe.vcf_reader_func import VcfReader,checkRecordPass,checkPopWithoutMissing,getAlleleCountDict
import pgpipe.vcf_reader_func as vr
from pgpipe.model import Model, read_model_file


##    --model : str
##        Name of the model
##    --model-tree : str, str
##        Assign newick-formatted population tree to model
##        Example: --model-tree model_name newick_tree
##    --model-tree-file : str, str
##        Assign newick-formatted population tree file to model
##        Example: --model-tree-file model_name newick_file
##    --model-pop : str, str
##        Assign single population name to a model
##        Example: --model-pop model_name pop_name
##    --model-pops : str, list
##        Assign multiple population names to a model
##        Example: --model-pops model_name pop_name1 pop_name2 etc.
##    --model-pop-file : str, str
##        Assign population names to a model using a file
##        Example: --model-pop-file model_name pop_file
##    --pop-ind : str, str
##        Assign single individual name to a population
##        Example: --pop-ind pop_name ind_name
##    --pop-inds : str, list
##        Assign multiple population names to a model
##        Example: --pop-inds pop_name ind_name1 ind_name2 etc.
##    --pop-ind-file : str, str
##        Assign population names to a model using a file
##        Example: --pop-ind-file pop_name ind_file
##    --out : str
##        Filename of the model output


model_file = r"/home/jodyhey/temp/PPP/jhfiles/panmodels.model"

popmodel = read_model_file(model_file)
exit()

vcffile = r"/home/jodyhey/temp/PPP/jhfiles/Pan_all_hicov_chr22_decrun_missingasref.vcf"

reffile = r"/home/jodyhey/temp/PPP/jhfiles/chr22_hg18.fa"

a = vr.VcfReader(vcffile,popmodel=popmodel['popmod1'])
##x = a.getNext()
##
##
##for pop in popmodel['popmod1'].pop_list:
##    for ind in popmodel['popmod1'].ind_dict[pop]:
##        print(pop,ind,x.samples[ind].alleles)

tfn = "treemixfile.out"
tf = open(tfn,'w')
for pop in popmodel['popmod1'].pop_list:
    tf.write("%s\t"%pop)
tf.write("\n")
numinds = []
for pop in popmodel['popmod1'].pop_list:
    numinds.append(len(popmodel['popmod1'].ind_dict[pop]))
ploidy = 2
while True:
    x = a.getNext()
    if type(x) == type(None):
        break
    if vr.checkRecordPass(x,remove_indels=True,remove_multiallele=True):
        buildstr = ""
        checkcount = [0,0]
        spacer = [',',' ']
        for pop in popmodel['popmod1'].pop_list:
            adic = vr.getAlleleCountDict(x,popmodel['popmod1'].ind_dict[pop])
            
            for ai in range(2):
                allele = x.alleles[ai]
                try:
                    count = adic[allele]
                    checkcount[ai] += count
                    buildstr += "%d%s"%(count,spacer[ai])
                except KeyError:
                    buildstr += "0%s"%(spacer[ai])
        if checkcount[0] > 0 and checkcount[1] > 0 and sum(checkcount) == ploidy*sum(numinds):
            tf.write("%s\n"%buildstr)
tf.close()


def pipe_bcftools (bcftools_call_args):
    '''
        Calls bcftools with pipe output

        The output of this function is the stdout and stderr of bcftools. This
        function should only be used if bcftools is being used as the stdin of
        another function. Please note that this function does not check the for
        errors in the bcftools call. Please check for errors after the call is
        closed using check_bcftools_for_errors.

        Parameters
        ----------
        bcftools_stderr : str
            bcftools stderr

        Returns
        -------
        bcftools_call : PIPE
            Pipe of subprocess call, including both stdout and stderr

    '''

    # Confirm where the specifed executable is located
    bcftools_path = confirm_executable('bcftools')

    # Check if the executable was found
    if not bcftools_path:
        raise IOError('bcftools not found. Please confirm the executable is installed')

    # bcftools subprocess call
    bcftools_call = subprocess.Popen([bcftools_path] + list(map(str, bcftools_call_args)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    return bcftools_call


                
            
        
