import sys
import logging
import vcftools
import model
import model_creator
import vcf_filter
import vcf_calc
import vcf_sampler
import stat_sampler
import vcf_split
import vcf_phase
import os
import time
from ruffus import *

out_prefix = os.getcwd() + '/output/' + str(int(time.time())) + '/merged_chr1_10000'

input_dir = os.getcwd() + '/example/input/'


DEBUG_do_not_define_tail_task = False
DEBUG_do_not_define_head_task = False


def exp_handler(type, val, tb):
    logging.exception("Error: %s" % (str(val)))


sys.excepthook = exp_handler

starting_files = ['merged_chr1_10000.vcf.gz']

# ======================================================================================================================
# pipeline sequence: vcf_filter -> vcf_calc(windowed-weir-fst) -> stat_sampler -> vcf_split -> vcf_phase

subdir = str(int(time.time()))
os.mkdir('output/' + subdir)

@transform(starting_files, filter=suffix(".vcf.gz"), output=".recode.vcf", output_dir='output/' + subdir)
def run_vcf_to_filter(vcf_in, seq_out):
    args = ['--vcf', starting_files[0], '--out-format', 'vcf', '--out-prefix',
            'output/' + subdir + '/merged_chr1_10000', '--overwrite']
    vcf_filter.run(args)


@follows(run_vcf_to_filter)
@transform(run_vcf_to_filter, filter=suffix(".recode.vcf"), output=".windowed.weir.fst", output_dir='output/' + subdir)
def run_vcf_to_calc(f_in, f_out):
    args = ['--vcf', f_in, '--calc-statistic', 'windowed-weir-fst', '--out-prefix', 'output/' + subdir +'/merged_chr1_10000',
            '--statistic-window-size', '10000', '--statistic-window-step', '20000', '--model-file', 'input.model',
            '--model', '2Pop', '--overwrite']
    vcf_calc.run(args)


@follows(run_vcf_to_calc)
@transform(run_vcf_to_calc, filter=suffix(".windowed.weir.fst"), output=".sampled", output_dir='output/' + subdir)
def run_stat_to_sampler(f_in, f_out):
    args = ['--statistic-file', f_in, '--out-prefix', 'output/' + subdir +'/merged_chr1_10000',
            '--calc-statistic', 'windowed-weir-fst', '--sampling-scheme', 'uniform', '--uniform-bins', '5',
            '--sample-size', '5', '--random-seed', '100']
    stat_sampler.run(args)

@follows(run_vcf_to_calc)
@transform(run_vcf_to_calc, filter=suffix(".windowed.weir.fst"), output=".vcf.gz", output_dir='output/' + subdir)
def run_vcf_to_split(f_in, f_out):
    args = ['--vcf', starting_files[0], '--model-file', 'input.model', '--model', '2Pop', '--split-file',
            f_in, '--split-method', 'statistic-file', '--overwrite']
    vcf_split.run(args)

@follows(run_vcf_to_split)
@transform(run_vcf_to_split, filter=suffix(".vcf.gz"), output=".out", output_dir=os.getcwd())
def run_vcf_to_phase(f_in, f_out):
    args = ['--vcf', starting_files[0], '--beagle-path', '.', '--out-prefix', 'output/' + subdir +'/merged_chr1_10000',
            '--phase-algorithm', 'beagle', '--overwrite']
    vcf_phase.run(args)

# ======================================================================================================================
# pipeline sequence: vcf_filter -> vcf_calc(weir-fst) -> vcf_split -> vcf_phase

# subdir = time.time()
#
# @transform(starting_files, filter=suffix(".vcf.gz"), output=".recode.vcf")
# def run_vcf_to_filter(vcf_in, seq_out):
#     args = ['--vcf', starting_files[0], '--out-format', 'vcf', '--out-prefix',
#             'output/' + subdir + '/merged_chr1_10000', '--overwrite']
#     vcf_filter.run(args)
#
#
# @follows(run_vcf_to_filter)
# @transform(run_vcf_to_filter, filter=suffix(".recode.vcf"), output=".weir.fst")
# def run_vcf_to_calc(f_in, f_out):
#     args = ['--vcf', f_in, '--calc-statistic', 'weir-fst', '--out-prefix', 'output/' + subdir + '/merged_chr1_10000',
#             '--model-file', 'input.model', '--model', '2Pop', '--overwrite']
#     vcf_calc.run(args)
#
#
# @follows(run_vcf_to_calc)
# @transform(run_vcf_to_calc, filter=suffix(".weir.fst"), output=".vcf.gz")
# def run_vcf_to_split(f_in, f_out):
#     args = ['--vcf', starting_files[0], '--model-file', 'input.model', '--model', '2Pop', '--split-file', f_in,
#             '--split-method', 'statistic-file', '--overwrite']
#     vcf_split.run(args)
#
#
# @transform(run_vcf_to_filter, filter=suffix(".vcf.gz"), output=".out")
# def run_vcf_to_phase(f_in, f_out):
#     args = ['--vcf', starting_files[0], '--beagle-path', '.', '--out-prefix', 'output/' + subdir + '/merged_chr1_10000',
#             '--phase-algorithm', 'beagle', '--overwrite']
#     vcf_phase.run(args)

# ======================================================================================================================


pipeline_run(pipeline= "main")