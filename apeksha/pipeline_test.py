import unittest
import sys
import os
import logging
import filecmp
import shutil
import vcftools
import ruffus
import itertools
import model
import model_creator
import vcf_filter
import vcf_calc
import vcf_sampler
from itertools import combinations
from ruffus import *

import ruffus.cmdline as cmdline


#=================================function to make dir and sub -dir #================================#

  
for i in range(1,61):
    os.makedirs(os.path.join('Run', 'test' + str(i))) 

#============================# Create path names to move the generated o/p files #===================#
newpath = []
filenamee = 'test'
counter1 = 1
counter2 = 2
list1 = []
list2 = range(61)    #random number, given len needed
for x in list2:
    counter1 = str(counter1)
    full_name = (filenamee+counter1)
    list1.append(full_name)
    counter1 = counter2
    counter2+=1

for x in list1:
    y = "/home/staff/asukeshkall/Downloads/Test-3/Run/"+ x
    newpath.append(y)        

#================================================# Generate all combinations  #=======================================================================#

filter_list=[['merged_chr1_10000.vcf.gz'],['--out-prefix','merged_chr1_10000'],['--out-format','removed_sites'],['--out-format','kept_sites'],['--out-format','vcf'],['--out-format','bcf'],['--filter-min-alleles', '2'],['--filter-max-alleles', '4']]

calc_list=[['--out-prefix','merged_chr1_10000'],['--pop-file', 'Paniscus.txt'],['--pop-file', 'Troglodytes.txt'],['--model-file','input.model'],['--model','2Pop'],['--calc-statistic', 'weir-fst'],['--calc-statistic', 'windowed-weir-fst'],['--calc-statistic', 'TajimaD'],['--calc-statistic', 'pi'],['--calc-statistic', 'freq'],['--calc-statistic', 'het'], ['--statistic-window-size', '2'],['--statistic-window-step', '2']]

sampler_list=[['merged_chr1_10000.vcf.gz'],['--model-file','input.model'],['--model','2Pop'],['--statistic-file','merged_chr1_10000.windowed.weir.fst'],['--sample-file','sampled_data.tsv'],['--out-dir','Sample_Files'],['--out-prefix','Sample'],['--out-format','vcf.gz'],['--overwrite'],['--calc-statistic','windowed-weir-fst'], ['--statistic-window-size','10000'],['--sampling-scheme','random'],['--uniform-bins','10'],['--sample-size','10'],['--random-seed','456']]


comb = []
whole_list = []
retList = []

def combFunc(a):
    for i in range(len(a)):
        b = list(combinations(a, i+1))
        c = map(list, b)
        comb.append(c)

    for x in comb:
        for y in x:
             b = list(itertools.chain.from_iterable(y))
             whole_list.append(b)
                #print b
    retList = whole_list
    return retList

 

#=============================================# Filter level 1 function #========================================================================#

filList1 = []
filList2 = []
filList3 = []
filList4 = []
filList5 = []
filList6 = []

def filFunc(name, lol, listname):
    for x in lol:
        for y in x:
            if y== name:
               listname.append(x) 

#============================= input =====================#

filter_comb = combFunc(filter_list)
filFunc('merged_chr1_10000.vcf.gz', filter_comb, filList1)
fnameFilter = filList1
#for t in fnameFilter:
    #print t
    #print "--------------"

filFunc('merged_chr1_10000', fnameFilter, filList2)
prefixFilter = filList2
#for t in prefixFilter:
    #print t
    #print "--------------"

#============================= input =====================#

calc_comb = combFunc(calc_list)

filFunc('Paniscus.txt', calc_comb, filList3)
panFilter = filList3
#for t in panFilter:
    #print t
    #print "--------------"

filFunc('Troglodytes.txt', panFilter, filList4)
trogFilter = filList4
#for t in trogFilter:
    #print t
    #print "--------------"

#=============================input=======================#

sampler_comb = combFunc(sampler_list)
filFunc('merged_chr1_10000.vcf.gz', sampler_comb, filList5)
fnameSampler = filList5
#for t in fnameSampler:
    #print t
    #print "--------------"

filFunc('merged_chr1_10000.windowed.weir.fst', fnameSampler, filList6)
statFilter = filList6
#for t in range(0,200):
    #print t, statFilter[t]
    #print "--------------"

#======================# Filter level 2 - to filter the redundant options #===================================#

app = []
for m in fnameFilter:  
    for x,y in itertools.combinations(m, 2):  
        if x == y:
           app.append(m) 
        
filter_final_list = [d for d in prefixFilter if d not in app]
#for t in filter_final_list:
    #print t

#================================================================================================#

list1 = []

for x in filter_final_list:
    for y in x:
        if y == 'bcf':
           list1.append(x)

#for t in list1:
    #print t


list2 = [['--out-prefix', 'merged_chr1_10000', '--pop-file', 'Paniscus.txt', '--pop-file', 'Troglodytes.txt', '--model-file', 'input.model', '--model', '2Pop', '--calc-statistic', 'windowed-weir-fst'],
       ['--out-prefix', 'merged_chr1_10000', '--pop-file', 'Paniscus.txt', '--pop-file', 'Troglodytes.txt', '--model-file', 'input.model', '--model', '2Pop', '--calc-statistic', 'windowed-weir-fst','--statistic-window-size','2'],
       ['--out-prefix', 'merged_chr1_10000', '--pop-file', 'Paniscus.txt', '--pop-file', 'Troglodytes.txt', '--model-file', 'input.model', '--model', '2Pop', '--calc-statistic', 'windowed-weir-fst','--statistic-window-step','2'],
       ['--out-prefix', 'merged_chr1_10000', '--pop-file', 'Paniscus.txt', '--pop-file', 'Troglodytes.txt', '--model-file', 'input.model', '--model', '2Pop', '--calc-statistic', 'windowed-weir-fst','--statistic-window-size','2','--statistic-window-step','2']]



list3 = [['merged_chr1_10000.vcf.gz',
                         '--statistic-file', 'merged_chr1_10000.windowed.weir.fst',
                         '--out-prefix', 'merged_chr1_10000',
                         '--sample-size', '20',
                         '--random-seed', '1000',
                         '--overwrite'],
         ['merged_chr1_10000.vcf.gz',
                         '--statistic-file', 'merged_chr1_10000.windowed.weir.fst',
                         '--out-prefix', 'merged_chr1_10000',
                         '--sample-size', '10',
                         '--random-seed', '1000',
                         '--uniform-bins','10',
                         '--overwrite'],
         ['merged_chr1_10000.vcf.gz',
                         '--statistic-file', 'merged_chr1_10000.windowed.weir.fst',
                         '--out-prefix', 'merged_chr1_10000',
                         '--sample-size', '20',
                         '--random-seed','254',
                         '--overwrite'],
         ['merged_chr1_10000.vcf.gz',
                         '--statistic-file', 'merged_chr1_10000.windowed.weir.fst',
                         '--out-prefix', 'merged_chr1_10000',
                         '--sample-size', '20',
                         '--random-seed', '1000',
                         '--statistic-window-size','10000',
                         '--overwrite']]

#======================# Filter level 2 - to filter the redundant options #===================================#

calc_unmatched = []
for e in trogFilter:
    if e.count('--calc-statistic')>1 or e.count('--out-prefix')>1:
       calc_unmatched.append(e)
        #print unmatched
               
calc_final_list = [d for d in trogFilter if d not in calc_unmatched]  

#for t in calc_final_list:
    #print t


#=============================================================================================================#

DEBUG_do_not_define_tail_task = False
DEBUG_do_not_define_head_task = False


starting_files = ['merged_chr1_10000.vcf.gz']

#file_names = ['merged_chr1_10000.removed.sites', 'merged_chr1_10000.kept.sites','merged_chr1_10000.recode.vcf','merged_chr1_10000.recode.bcf']

def exp_handler(type,val,tb):
    logging.exception("Error: %s" % (str(val)))

sys.excepthook = exp_handler

def start():
    file1 = ['merged_chr1_10000.vcf.gz']
    return file1


@transform(
           starting_files,
           filter = suffix(".vcf.gz"),
           output = ".recode.vcf")
def run_vcf_to_filter(vcf_in,seq_out):
    args = ['merged_chr1_10000.vcf.gz',
            '--out-format', 'vcf',
            '--out-prefix', 'merged_chr1_10000']
    vcf_filter.run(args)

@follows(run_vcf_to_filter)
@transform( run_vcf_to_filter,
           filter = suffix(".recode.vcf"),
           output = ".windowed.weir.fst")
def run_vcf_to_calc(f_in,f_out):
    args = [f_in,
            '--calc-statistic', 'windowed-weir-fst',
            '--out-prefix', 'merged_chr1_10000',
            '--model-file','input.model',
            '--model', '2Pop',
            '--overwrite']
   # print (f_in)
    #print (args)
    vcf_calc.run(args)

@follows(run_vcf_to_calc)
@transform( run_vcf_to_calc,
           filter = suffix(".windowed.weir.fst"),
           output = ".tsv")
def run_vcf_to_sampler(f_in,f_out):
    args = [starting_files[0],
                         '--statistic-file',f_in,
                         '--calc-statistic','windowed-weir-fst',
                         '--sample-size', '20',
                         '--random-seed', '1000',
                         '--sample-file','sampled_data_1.tsv',
                         '--out-dir','Sample_Files_1',
                         '--out-prefix', 'Sample_1',
                         '--overwrite']
    #print (args)
    vcf_sampler.run(args)


#====================================================================#
@follows(run_vcf_to_sampler)
@transform(
           starting_files,
           filter = suffix(".vcf.gz"),
           output = ".recode.bcf")
def run_vcf_to_filter1(vcf_in,seq_out):
    args = ['merged_chr1_10000.vcf.gz',
            '--out-format', 'bcf',
            '--out-prefix', 'merged_chr1_10000']
    vcf_filter.run(args)


@follows(run_vcf_to_filter1)
@transform( run_vcf_to_filter1,
           filter = suffix(".recode.bcf"),
           output = ".windowed.weir.fst")
def run_vcf_to_calc1(f_in,f_out):
    args = [f_in,
            '--calc-statistic', 'windowed-weir-fst',
            '--out-prefix', 'merged_chr1_10000.1',
            '--model-file','input.model',
            '--model', '2Pop',
            '--overwrite']
    #print (f_in)
    #print (args)
    vcf_calc.run(args)

@follows(run_vcf_to_calc1)
@transform( run_vcf_to_calc1,
           filter = suffix(".windowed.weir.fst"),
           output = ".tsv")
def run_vcf_to_sampler1(f_in,f_out):
    args = [starting_files[0],
                         '--statistic-file',f_in,
                         '--calc-statistic','windowed-weir-fst',
                         '--sample-size', '20',
                         '--random-seed', '1000',
                         '--sample-file','sampled_data_2.tsv',
                         '--out-dir','Sample_Files_2',
                         '--out-prefix', 'Sample_2',
                         '--overwrite']
    #print (args)
    vcf_sampler.run(args)


#=====================================================++#


@follows(run_vcf_to_sampler1)
@transform(
           starting_files,
           filter = suffix(".vcf.gz"),
           output = ".recode.vcf.gz")
def run_vcf_to_filter2(vcf_in,seq_out):
    args = ['merged_chr1_10000.vcf.gz',
            '--out-format', 'vcf.gz',
            '--out-prefix', 'merged_chr1_10000']
    vcf_filter.run(args)


@follows(run_vcf_to_filter2)
@transform( run_vcf_to_filter2,
           filter = suffix(".recode.vcf.gz"),
           output = ".windowed.weir.fst")
def run_vcf_to_calc2(f_in,f_out):
    args = [f_in,
            '--calc-statistic', 'windowed-weir-fst',
            '--out-prefix', 'merged_chr1_10000.2',
            '--model-file','input.model',
            '--model', '2Pop',
            '--overwrite']
    #print (f_in)
    #print (args)
    vcf_calc.run(args)

@follows(run_vcf_to_calc2)
@transform( run_vcf_to_calc2,
           filter = suffix(".windowed.weir.fst"),
           output = ".tsv")
def run_vcf_to_sampler2(f_in,f_out):
    args = [starting_files[0],
                         '--statistic-file',f_in,
                         '--calc-statistic','windowed-weir-fst',
                         '--sample-size', '20',
                         '--random-seed', '1000',
                         '--sample-file','sampled_data_3.tsv',
                         '--out-dir','Sample_Files_3',
                         '--out-prefix', 'Sample_3',
                         '--overwrite']
    #print (args)
    vcf_sampler.run(args)

#=======================================================#

@follows(run_vcf_to_sampler2)
@transform(
           starting_files,
           filter = suffix(".vcf.gz"),
           output = ".recode.vcf")
def run_vcf_to_filter3(vcf_in,seq_out):
    args = ['merged_chr1_10000.vcf.gz',
            '--out-format', 'vcf',
            '--out-prefix', 'merged_chr1_10000']
    vcf_filter.run(args)


@follows(run_vcf_to_filter3)
@transform( run_vcf_to_filter3,
           filter = suffix(".recode.vcf"),
           output = ".windowed.pi")
def run_vcf_to_calc3(f_in,f_out):
    args = [f_in,
            '--calc-statistic', 'windowed-pi',
            '--out-prefix', 'merged_chr1_10000',
            '--model-file','input.model',
            '--model', '2Pop',
            '--overwrite']
    #print (f_in)
    #print (args)
    vcf_calc.run(args)


@follows(run_vcf_to_calc3)
@transform( run_vcf_to_calc3,
           filter = suffix(".windowed.pi"),
           output = ".tsv")
def run_vcf_to_sampler3(f_in,f_out):
    args = [starting_files[0],
                         '--statistic-file',f_in,
                         '--calc-statistic','windowed-pi',
                         '--sample-size', '20',
                         '--random-seed', '1000',
                         '--sample-file','sampled_data_4.tsv',
                         '--out-dir','Sample_Files_4',
                         '--out-prefix', 'Sample_4',
                         '--overwrite']
    #print (args)
    vcf_sampler.run(args)

#=======================================================#


@follows(run_vcf_to_sampler3)
@transform(
           starting_files,
           filter = suffix(".vcf.gz"),
           output = ".recode.bcf")
def run_vcf_to_filter4(vcf_in,seq_out):
    args = ['merged_chr1_10000.vcf.gz',
            '--out-format', 'bcf',
            '--out-prefix', 'merged_chr1_10000']
    vcf_filter.run(args)


@follows(run_vcf_to_filter4)
@transform( run_vcf_to_filter4,
           filter = suffix(".recode.bcf"),
           output = ".windowed.pi")
def run_vcf_to_calc4(f_in,f_out):
    args = [f_in,
            '--calc-statistic', 'windowed-pi',
            '--out-prefix', 'merged_chr1_10000.1',
            '--model-file','input.model',
            '--model', '2Pop',
            '--overwrite']
    #print (f_in)
    #print (args)
    vcf_calc.run(args)


@follows(run_vcf_to_calc4)
@transform( run_vcf_to_calc4,
           filter = suffix(".windowed.pi"),
           output = ".tsv")
def run_vcf_to_sampler4(f_in,f_out):
    args = [starting_files[0],
                         '--statistic-file',f_in,
                         '--calc-statistic','windowed-pi',
                         '--sample-size', '20',
                         '--random-seed', '1000',
                         '--sample-file','sampled_data_5.tsv',
                         '--out-dir','Sample_Files_5',
                         '--out-prefix', 'Sample_5',
                         '--overwrite']
    #print (args)
    vcf_sampler.run(args)

#=======================================================#

@follows(run_vcf_to_sampler4)
@transform(
           starting_files,
           filter = suffix(".vcf.gz"),
           output = ".recode.vcf.gz")
def run_vcf_to_filter5(vcf_in,seq_out):
    args = ['merged_chr1_10000.vcf.gz',
            '--out-format', 'vcf.gz',
            '--out-prefix', 'merged_chr1_10000']
    vcf_filter.run(args)


@follows(run_vcf_to_filter5)
@transform( run_vcf_to_filter5,
           filter = suffix(".recode.vcf.gz"),
           output = ".windowed.pi")
def run_vcf_to_calc5(f_in,f_out):
    args = [f_in,
            '--calc-statistic', 'windowed-pi',
            '--out-prefix', 'merged_chr1_10000.2',
            '--model-file','input.model',
            '--model', '2Pop',
            '--overwrite']
    #print (f_in)
    #print (args)
    vcf_calc.run(args)


@follows(run_vcf_to_calc5)
@transform( run_vcf_to_calc5,
           filter = suffix(".windowed.pi"),
           output = ".tsv")
def run_vcf_to_sampler5(f_in,f_out):
    args = [starting_files[0],
                         '--statistic-file',f_in,
                         '--calc-statistic','windowed-pi',
                         '--sample-size', '20',
                         '--random-seed', '1000',
                         '--sample-file','sampled_data_6.tsv',
                         '--out-dir','Sample_Files_6',
                         '--out-prefix', 'Sample_6',
                         '--overwrite']
    #print (args)
    vcf_sampler.run(args)

#=======================================================#
@follows(run_vcf_to_sampler5)
@transform(
           starting_files,
           filter = suffix(".vcf.gz"),
           output = ".recode.vcf")
def run_vcf_to_filter6(vcf_in,seq_out):
    args = ['merged_chr1_10000.vcf.gz',
            '--out-format', 'vcf',
            '--out-prefix', 'merged_chr1_10000']
    vcf_filter.run(args)


@follows(run_vcf_to_filter6)
@transform( run_vcf_to_filter6,
           filter = suffix(".recode.vcf"),
           output = ".Tajima.D")
def run_vcf_to_calc6(f_in,f_out):
    args = [f_in,
            '--calc-statistic', 'TajimaD',
            '--out-prefix', 'merged_chr1_10000',
            '--model-file','input.model',
            '--model', '2Pop',
            '--overwrite']
    #print (f_in)
    #print (args)
    vcf_calc.run(args)


@follows(run_vcf_to_calc6)
@transform( run_vcf_to_calc6,
           filter = suffix(".Tajima.D"),
           output = ".tsv")
def run_vcf_to_sampler6(f_in,f_out):
    args = [starting_files[0],
                         '--statistic-file',f_in,
                         '--calc-statistic','TajimaD',
                         '--sample-size', '20',
                         '--random-seed', '1000',
                         '--sample-file','sampled_data_7.tsv',
                         '--out-dir','Sample_Files_7',
                         '--out-prefix', 'Sample_7',
                         '--overwrite']
    #print (args)
    vcf_sampler.run(args)

#=======================================================#

@follows(run_vcf_to_sampler6)
@transform(
           starting_files,
           filter = suffix(".vcf.gz"),
           output = ".recode.bcf")
def run_vcf_to_filter7(vcf_in,seq_out):
    args = ['merged_chr1_10000.vcf.gz',
            '--out-format', 'bcf',
            '--out-prefix', 'merged_chr1_10000']
    vcf_filter.run(args)


@follows(run_vcf_to_filter7)
@transform( run_vcf_to_filter7,
           filter = suffix(".recode.bcf"),
           output = ".Tajima.D")
def run_vcf_to_calc7(f_in,f_out):
    args = [f_in,
            '--calc-statistic', 'TajimaD',
            '--out-prefix', 'merged_chr1_10000.1',
            '--model-file','input.model',
            '--model', '2Pop',
            '--overwrite']
    #print (f_in)
    #print (args)
    vcf_calc.run(args)


@follows(run_vcf_to_calc7)
@transform( run_vcf_to_calc7,
           filter = suffix(".Tajima.D"),
           output = ".tsv")
def run_vcf_to_sampler7(f_in,f_out):
    args = [starting_files[0],
                         '--statistic-file',f_in,
                         '--calc-statistic','TajimaD',
                         '--sample-size', '20',
                         '--random-seed', '1000',
                         '--sample-file','sampled_data_8.tsv',
                         '--out-dir','Sample_Files_8',
                         '--out-prefix', 'Sample_8',
                         '--overwrite']
    #print (args)
    vcf_sampler.run(args)

#=======================================================#

@follows(run_vcf_to_sampler7)
@transform(
           starting_files,
           filter = suffix(".vcf.gz"),
           output = ".recode.vcf.gz")
def run_vcf_to_filter8(vcf_in,seq_out):
    args = ['merged_chr1_10000.vcf.gz',
            '--out-format', 'vcf.gz',
            '--out-prefix', 'merged_chr1_10000']
    vcf_filter.run(args)


@follows(run_vcf_to_filter8)
@transform( run_vcf_to_filter8,
           filter = suffix(".recode.vcf.gz"),
           output = ".Tajima.D")
def run_vcf_to_calc8(f_in,f_out):
    args = [f_in,
            '--calc-statistic', 'TajimaD',
            '--out-prefix', 'merged_chr1_10000.2',
            '--model-file','input.model',
            '--model', '2Pop',
            '--overwrite']
    #print (f_in)
    #print (args)
    vcf_calc.run(args)


@follows(run_vcf_to_calc8)
@transform( run_vcf_to_calc8,
           filter = suffix(".Tajima.D"),
           output = ".tsv")
def run_vcf_to_sampler8(f_in,f_out):
    args = [starting_files[0],
                         '--statistic-file',f_in,
                         '--calc-statistic','TajimaD',
                         '--sample-size', '20',
                         '--random-seed', '1000',
                         '--sample-file','sampled_data_9.tsv',
                         '--out-dir','Sample_Files_9',
                         '--out-prefix', 'Sample_9',
                         '--overwrite']
    #print (args)
    vcf_sampler.run(args)

#=======================================================#


class Test_ruffus(unittest.TestCase):


    def test_simpler (self):
        pipeline_run(pipeline= "main")

    #def func1(pipename, filterName, filterInput, filterOutput, calcName, calcInput, calcOutput):
       # test_pipeline = Pipeline("test")
        #test_pipeline.start()
      #  test_pipeline.transform(start, name = 'test:start', filter=suffix(".vcf.gz"),output = ".recode.vcf")
      #  test_pipeline.transform(run_vcf_to_filter, name = 'test:run_vcf_to_filter', filter = suffix(".recode.vcf"),output = ".windowed.weir.fst")
      #  test_pipeline.transform(run_vcf_to_calc, name = 'test:run_vcf_to_calc', filter = suffix(".windowed.weir.fst"),output = ".tsv")


    def test_newstyle_simpler1 (self):
        test_pipeline = Pipeline("test")
        #test_pipeline.start()
        test_pipeline.transform(start, name = 'test:start', filter=suffix(".vcf.gz"),output = ".recode.vcf")
        test_pipeline.transform(run_vcf_to_filter, name = 'test:run_vcf_to_filter', filter = suffix(".recode.vcf"),output = ".windowed.weir.fst")
        test_pipeline.transform(run_vcf_to_calc, name = 'test:run_vcf_to_calc', filter = suffix(".windowed.weir.fst"),output = ".tsv")


    def test_newstyle_simpler2 (self):
        test_pipeline = Pipeline("test1")
        #test_pipeline.start()
        test_pipeline.transform(start, name = 'test1:start', filter=suffix(".vcf.gz"),output = ".recode.bcf")
        test_pipeline.transform(run_vcf_to_filter, name = 'test1:run_vcf_to_filter1', filter = suffix(".recode.bcf"),output = ".windowed.weir.fst")
        test_pipeline.transform(run_vcf_to_calc, name = 'test1:run_vcf_to_calc1', filter = suffix(".windowed.weir.fst"),output = ".tsv")

    def test_newstyle_simpler3 (self):
        test_pipeline = Pipeline("test2")
        #test_pipeline.start()
        test_pipeline.transform(start, name = 'test2:start', filter=suffix(".vcf.gz"),output = ".recode.vcf.gz")
        test_pipeline.transform(run_vcf_to_filter, name = 'test2:run_vcf_to_filter2', filter = suffix(".recode.vcf.gz"),output = ".windowed.weir.fst")
        test_pipeline.transform(run_vcf_to_calc, name = 'test2:run_vcf_to_calc2', filter = suffix(".windowed.weir.fst"),output = ".tsv")

#==================================================================================================#

    def test_newstyle_simpler4 (self):
        test_pipeline = Pipeline("test3")
        #test_pipeline.start()
        test_pipeline.transform(start, name = 'test3:start', filter=suffix(".vcf.gz"),output = ".recode.vcf")
        test_pipeline.transform(run_vcf_to_filter, name = 'test3:run_vcf_to_filter3', filter = suffix(".recode.vcf"),output = ".windowed.pi")
        test_pipeline.transform(run_vcf_to_calc, name = 'test3:run_vcf_to_calc3', filter = suffix(".windowed.pi"),output = ".tsv")


    def test_newstyle_simpler5 (self):
        test_pipeline = Pipeline("test4")
        #test_pipeline.start()
        test_pipeline.transform(start, name = 'test4:start', filter=suffix(".vcf.gz"),output = ".recode.bcf")
        test_pipeline.transform(run_vcf_to_filter, name = 'test4:run_vcf_to_filter4', filter = suffix(".recode.bcf"),output = ".windowed.pi")
        test_pipeline.transform(run_vcf_to_calc, name = 'test4:run_vcf_to_calc4', filter = suffix(".windowed.pi"),output = ".tsv")

    def test_newstyle_simpler6 (self):
        test_pipeline = Pipeline("test5")
        #test_pipeline.start()
        test_pipeline.transform(start, name = 'test5:start', filter=suffix(".vcf.gz"),output = ".recode.vcf.gz")
        test_pipeline.transform(run_vcf_to_filter, name = 'test5:run_vcf_to_filter5', filter = suffix(".recode.vcf.gz"),output = ".windowed.pi")
        test_pipeline.transform(run_vcf_to_calc, name = 'test5:run_vcf_to_calc5', filter = suffix(".windowed.pi"),output = ".tsv")

#==================================================================================================#

    def test_newstyle_simpler7 (self):
        test_pipeline = Pipeline("test6")
        #test_pipeline.start()
        test_pipeline.transform(start, name = 'test6:start', filter=suffix(".vcf.gz"),output = ".recode.vcf")
        test_pipeline.transform(run_vcf_to_filter, name = 'test6:run_vcf_to_filter6', filter = suffix(".recode.vcf"),output = ".Tajima.D")
        test_pipeline.transform(run_vcf_to_calc, name = 'test6:run_vcf_to_calc6', filter = suffix(".Tajima.D"),output = ".tsv")


    def test_newstyle_simpler8 (self):
        test_pipeline = Pipeline("test7")
        #test_pipeline.start()
        test_pipeline.transform(start, name = 'test7:start', filter=suffix(".vcf.gz"),output = ".recode.bcf")
        test_pipeline.transform(run_vcf_to_filter, name = 'test7:run_vcf_to_filter7', filter = suffix(".recode.bcf"),output = ".Tajima.D")
        test_pipeline.transform(run_vcf_to_calc, name = 'test7:run_vcf_to_calc7', filter = suffix(".Tajima.D"),output = ".tsv")

    def test_newstyle_simpler9 (self):
        test_pipeline = Pipeline("test8")
        #test_pipeline.start()
        test_pipeline.transform(start, name = 'test8:start', filter=suffix(".vcf.gz"),output = ".recode.vcf.gz")
        test_pipeline.transform(run_vcf_to_filter, name = 'test8:run_vcf_to_filter8', filter = suffix(".recode.vcf.gz"),output = ".Tajima.D")
        test_pipeline.transform(run_vcf_to_calc, name = 'test8:run_vcf_to_calc8', filter = suffix(".Tajima.D"),output = ".tsv")

#==================================================================================================#

if __name__ == "__main__":
    unittest.main()
 


