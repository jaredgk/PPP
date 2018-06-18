import sys
from ruffus import *
import logging
import vcf_calc
import vcf_filter
import vcf_sampler


starting_files = ['merged_chr1_10000.vcf.gz']


def exp_handler(type,val,tb):
    logging.exception("Error: %s" % (str(val)))

sys.excepthook = exp_handler

@transform(starting_files,
           suffix(".vcf.gz"),
           ".recode.bcf")
def run_vcf_to_filter(vcf_in,seq_out):
    args = [vcf_in,
            '--out-format', 'bcf',
            '--out-prefix', 'merged_chr1_10000']
    vcf_filter.run(args)

@follows(run_vcf_to_filter)
@transform(run_vcf_to_filter,

           suffix(".recode.bcf"),
           ".windowed.weir.fst")
def run_vcf_to_calc(f_in,f_out):
    args = [f_in,
            '--calc-statistic', 'windowed-weir-fst',
            '--out-prefix', 'merged_chr1_10000',
            '--model-file','input.model',
            '--model', '2Pop',
            '--overwrite']

    vcf_calc.run(args)

@follows(run_vcf_to_calc)
@transform(run_vcf_to_calc,
           suffix(".windowed.weir.fst"),
           ".tsv")
def run_vcf_to_sampler(f_in,f_out):
    args = [starting_files[0],
                         '--statistic-file',f_in,
                         '--sample-size', '20',
                         '--random-seed', '1000',
                         '--overwrite']
    vcf_sampler.run(args)

@follows(run_vcf_to_sampler)
@transform(starting_files,
           suffix(".vcf.gz"),
           ".ima.u")
def vcf_to_ima_run(vcf_in,seq_out):
    args = ['--vcf', vcf_in,
            '--ref', 'human_g1k_chr11.fasta',
            '--bed', 'snp_region.txt',
            '--pop', 'testmodel.model']
    print (args)
    vcf_to_ima.vcf_to_ima(args)



def mainpipe():
    pipeline_run()


if __name__ == '__main__':
    mainpipe()
