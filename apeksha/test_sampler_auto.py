import unittest
import sys
import os
import logging
import filecmp
import shutil

import vcf_calc
import vcf_filter
import vcf_sampler
import vcftools
import ruffus
import itertools
from itertools import combinations


class calcTest(unittest.TestCase):
 
#=================================function to make dir and sub -dir #================================#

    def test_calc_run_1(self):
        for i in range(1,80):
           os.makedirs(os.path.join('Run', 'test' + str(i))) 

#================================================# Generate all combinations  #=======================================================================#

 
        comb = []

        a=[['merged_chr1_10000.vcf.gz'],['--model-file','input.model'],['--model','2Pop'],['--statistic-file','merged_chr1_10000.windowed.weir.fst'],['--sample-file','sampled_data.tsv'],['--out-dir','Sample_Files'],['--out-prefix','Sample'],['--out-format','vcf.gz'],['--overwrite'],['--calc-statistic','windowed-weir-fst'], ['--statistic-window-size','10000'],['--sampling-scheme','random'],['--uniform-bins','10'],['--sample-size','10'],['--random-seed','456']]

        for i in range(len(a)):

            b = list(combinations(a, i+1))
            c = map(list, b)
            comb.append(c)

        whole_list = []

        for x in comb:
            for y in x:
                b = list(itertools.chain.from_iterable(y))
                whole_list.append(b)
                #print b
        m = whole_list
        #print m


#=============================================# Filter level 1 function #========================================================================#

        filList1 = []
        filList2 = []
        filList3 = []


        def filFunc(name, lol, listname):
            for x in lol:
                for y in x:
                    if y== name:
                       listname.append(x) 
#=====================================================================================================================#

        filFunc('merged_chr1_10000.vcf.gz', m, filList1)
        fnameFilter = filList1
        #for t in fnameFilter:
           # print t
            #print "--------------"

        filFunc('merged_chr1_10000.windowed.weir.fst', fnameFilter, filList2)
        statFilter = filList2
        print len(statFilter)
        #for t in statFilter:
            #print t
            #print "--------------"
        
        for i in range(len(statFilter)):
            print i, statFilter[i]
            vcf_sampler.run(statFilter[i])
            self.assertTrue(os.path.isfile('sampled_data.tsv'))
            self.assertTrue(os.path.isdir('Sample_Files'))
            os.remove('sampled_data.tsv')
            shutil.rmtree('Sample_Files')



if __name__ == '__main__':
    unittest.main()


