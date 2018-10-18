import sys
import logging
import vcf_filter
import vcf_calc
import vcf_sampler
import stat_sampler
import vcf_split
import vcf_phase
import os
import shutil
from ruffus import *

root_dir = os.getcwd()
input_dir = os.getcwd() + '/example/input/'


DEBUG_do_not_define_tail_task = False
DEBUG_do_not_define_head_task = False

subdir = 1


def exp_handler(type, val, tb):
    logging.exception("Error: %s" % (str(val)))


sys.excepthook = exp_handler


for i in range(1, 10):
    os.makedirs(os.path.join('output', 'test' + str(i)))

starting_files = ['merged_chr1_10000.vcf.gz']

# ======================================================================================================================
# pipeline sequence: vcf_filter -> vcf_calc(windowed-weir-fst) -> stat_sampler(uniform) -> vcf_split -> vcf_phase

def setup1():
    global subdir
    shutil.copyfile(starting_files[0], root_dir + '/output/test' + str(subdir) + '/' + starting_files[0])
    os.chdir(root_dir + '/output/test' + str(subdir))


@follows(setup1)
@transform(starting_files[0], filter=suffix(".vcf.gz"), output=".recode.vcf")
def run_vcf_to_filter1(vcf_in, seq_out):
    args = ['--vcf', starting_files[0], '--out-format', 'vcf', '--out-prefix', 'merged_chr1_10000', '--overwrite']
    vcf_filter.run(args)


@follows(run_vcf_to_filter1)
@transform(run_vcf_to_filter1, filter=suffix(".recode.vcf"), output=".windowed.weir.fst")
def run_vcf_to_calc1(f_in, f_out):
    args = ['--vcf', f_in, '--calc-statistic', 'windowed-weir-fst', '--out-prefix', 'merged_chr1_10000',
            '--statistic-window-size', '10000', '--statistic-window-step', '20000', '--model-file',
            input_dir + 'input.model', '--model', '2Pop', '--overwrite']
    vcf_calc.run(args)


@follows(run_vcf_to_calc1)
@transform(run_vcf_to_calc1, filter=suffix(".windowed.weir.fst"), output=".sampled")
def run_stat_to_sampler1(f_in, f_out):
    args = ['--statistic-file', f_in, '--out-prefix', 'merged_chr1_10000', '--calc-statistic', 'windowed-weir-fst',
            '--sampling-scheme', 'uniform', '--uniform-bins', '5', '--sample-size', '5', '--random-seed', '100']
    stat_sampler.run(args)


@follows(run_vcf_to_calc1)
@transform(run_vcf_to_calc1, filter=suffix(".windowed.weir.fst"), output=".vcf.gz")
def run_vcf_to_split1(f_in, f_out):
    args = ['--vcf', starting_files[0], '--model-file', input_dir + 'input.model', '--model', '2Pop', '--split-file',
            f_in, '--split-method', 'statistic-file', '--statistic-window-size', '10000', '--overwrite']
    vcf_split.run(args)


# @follows(run_vcf_to_split1)
# @transform(run_vcf_to_split1, filter=suffix(".vcf.gz"), output=".out")
# def run_vcf_to_phase1(f_in, f_out):
#     args = ['--vcf', starting_files[0], '--beagle-path', input_dir, '--out-prefix', 'merged_chr1_10000',
#             '--phase-algorithm', 'beagle', '--overwrite']
#     vcf_phase.run(args)


# ======================================================================================================================
# pipeline sequence: vcf_filter -> vcf_calc(windowed-weir-fst) -> stat_sampler(random) -> vcf_split -> vcf_phase

@follows(run_vcf_to_split1)
def setup2():
    global subdir
    subdir = subdir + 1
    shutil.copyfile(starting_files[0], root_dir + '/output/test' + str(subdir) + '/' + starting_files[0])
    os.chdir(root_dir + '/output/test' + str(subdir))


@follows(setup2)
@transform(starting_files[0], filter=suffix(".vcf.gz"), output=".recode.vcf")
def run_vcf_to_filter2(vcf_in, seq_out):
    args = ['--vcf', starting_files[0], '--out-format', 'vcf', '--out-prefix', 'merged_chr1_10000', '--overwrite']
    vcf_filter.run(args)


@follows(run_vcf_to_filter2)
@transform(run_vcf_to_filter2, filter=suffix(".recode.vcf"), output=".windowed.weir.fst")
def run_vcf_to_calc2(f_in, f_out):
    args = ['--vcf', f_in, '--calc-statistic', 'windowed-weir-fst', '--out-prefix', 'merged_chr1_10000',
            '--statistic-window-size', '10000', '--statistic-window-step', '20000', '--model-file',
            input_dir + 'input.model', '--model', '2Pop', '--overwrite']
    vcf_calc.run(args)


@follows(run_vcf_to_calc2)
@transform(run_vcf_to_calc2, filter=suffix(".windowed.weir.fst"), output=".sampled")
def run_stat_to_sampler2(f_in, f_out):
    args = ['--statistic-file', f_in, '--out-prefix', 'merged_chr1_10000', '--calc-statistic', 'windowed-weir-fst',
            '--sampling-scheme', 'random', '--sample-size', '5', '--random-seed', '100', '--overwrite']
    stat_sampler.run(args)


@follows(run_vcf_to_calc2)
@transform(run_vcf_to_calc2, filter=suffix(".windowed.weir.fst"), output=".vcf.gz")
def run_vcf_to_split2(f_in, f_out):
    args = ['--vcf', starting_files[0], '--model-file', input_dir + 'input.model', '--model', '2Pop', '--split-file',
            f_in, '--split-method', 'statistic-file', '--statistic-window-size', '10000', '--overwrite']
    vcf_split.run(args)


# @follows(run_vcf_to_split2)
# @transform(run_vcf_to_split2, filter=suffix(".vcf.gz"), output=".out")
# def run_vcf_to_phase2(f_in, f_out):
#     args = ['--vcf', starting_files[0], '--beagle-path', input_dir, '--out-prefix', 'merged_chr1_10000',
#             '--phase-algorithm', 'beagle', '--overwrite']
#     vcf_phase.run(args)


# ======================================================================================================================
# pipeline sequence: vcf_filter -> vcf_calc(TajimaD) -> stat_sampler(uniform) -> vcf_split -> vcf_phase

@follows(run_vcf_to_split2)
def setup3():
    global subdir
    subdir = subdir + 1
    shutil.copyfile(starting_files[0], root_dir + '/output/test' + str(subdir) + '/' + starting_files[0])
    os.chdir(root_dir + '/output/test' + str(subdir))


@follows(setup3)
@transform(starting_files[0], filter=suffix(".vcf.gz"), output=".recode.vcf")
def run_vcf_to_filter3(vcf_in, seq_out):
    args = ['--vcf', starting_files[0], '--out-format', 'vcf', '--out-prefix', 'merged_chr1_10000', '--overwrite']
    vcf_filter.run(args)


@follows(run_vcf_to_filter3)
@transform(run_vcf_to_filter3, filter=suffix(".recode.vcf"), output=".Tajima.D")
def run_vcf_to_calc3(f_in, f_out):
    args = ['--vcf', f_in, '--calc-statistic', 'TajimaD', '--out-prefix', 'merged_chr1_10000',
            '--statistic-window-size', '10000', '--model-file', input_dir + 'input.model', '--model', '2Pop',
            '--overwrite']
    vcf_calc.run(args)


@follows(run_vcf_to_calc3)
@transform(run_vcf_to_calc3, filter=suffix(".Tajima.D"), output=".sampled")
def run_stat_to_sampler3(f_in, f_out):
    args = ['--statistic-file', f_in, '--out-prefix', 'merged_chr1_10000', '--calc-statistic',
            'TajimaD', '--sampling-scheme', 'uniform', '--uniform-bins', '5', '--sample-size', '5', '--random-seed',
            '100']
    stat_sampler.run(args)


@follows(run_vcf_to_calc3)
@transform(run_vcf_to_calc3, filter=suffix(".Tajima.D"), output=".vcf.gz")
def run_vcf_to_split3(f_in, f_out):
    args = ['--vcf', starting_files[0], '--model-file', input_dir + 'input.model', '--model', '2Pop', '--split-file',
            f_in, '--split-method', 'statistic-file', '--statistic-window-size', '10000', '--overwrite']
    vcf_split.run(args)


# @follows(run_vcf_to_split3)
# @transform(run_vcf_to_split3, filter=suffix(".vcf.gz"), output=".out")
# def run_vcf_to_phase3(f_in, f_out):
#     args = ['--vcf', starting_files[0], '--beagle-path', input_dir, '--out-prefix', 'merged_chr1_10000',
#             '--phase-algorithm', 'beagle', '--overwrite']
#     vcf_phase.run(args)


# ======================================================================================================================
# pipeline sequence: vcf_filter -> vcf_calc(TajimaD) -> stat_sampler(random) -> vcf_split -> vcf_phase

@follows(run_vcf_to_split3)
def setup4():
    global subdir
    subdir = subdir + 1
    shutil.copyfile(starting_files[0], root_dir + '/output/test' + str(subdir) + '/' + starting_files[0])
    os.chdir(root_dir + '/output/test' + str(subdir))


@follows(setup4)
@transform(starting_files[0], filter=suffix(".vcf.gz"), output=".recode.vcf")
def run_vcf_to_filter4(vcf_in, seq_out):
    args = ['--vcf', starting_files[0], '--out-format', 'vcf', '--out-prefix', 'merged_chr1_10000', '--overwrite']
    vcf_filter.run(args)


@follows(run_vcf_to_filter4)
@transform(run_vcf_to_filter4, filter=suffix(".recode.vcf"), output=".Tajima.D")
def run_vcf_to_calc4(f_in, f_out):
    args = ['--vcf', f_in, '--calc-statistic', 'TajimaD', '--out-prefix', 'merged_chr1_10000',
            '--statistic-window-size', '10000', '--model-file', input_dir + 'input.model', '--model', '2Pop',
            '--overwrite']
    vcf_calc.run(args)


@follows(run_vcf_to_calc4)
@transform(run_vcf_to_calc4, filter=suffix(".Tajima.D"), output=".sampled")
def run_stat_to_sampler4(f_in, f_out):
    args = ['--statistic-file', f_in, '--out-prefix', 'merged_chr1_10000', '--calc-statistic', 'TajimaD',
            '--sampling-scheme', 'random', '--sample-size', '5', '--random-seed', '100']
    stat_sampler.run(args)


@follows(run_vcf_to_calc4)
@transform(run_vcf_to_calc4, filter=suffix(".Tajima.D"), output=".vcf.gz")
def run_vcf_to_split4(f_in, f_out):
    args = ['--vcf', starting_files[0], '--model-file', input_dir + 'input.model', '--model', '2Pop', '--split-file',
            f_in, '--split-method', 'statistic-file', '--statistic-window-size', '10000', '--overwrite']
    vcf_split.run(args)


# @follows(run_vcf_to_split4)
# @transform(run_vcf_to_split4, filter=suffix(".vcf.gz"), output=".out")
# def run_vcf_to_phase4(f_in, f_out):
#     args = ['--vcf', starting_files[0], '--beagle-path', input_dir, '--out-prefix', 'merged_chr1_10000',
#             '--phase-algorithm', 'beagle', '--overwrite']
#     vcf_phase.run(args)


# ======================================================================================================================
# pipeline sequence: vcf_filter -> vcf_calc(window-pi) -> stat_sampler(uniform) -> vcf_split -> vcf_phase


@follows(run_vcf_to_split4)
def setup5():
    global subdir
    subdir = subdir + 1
    shutil.copyfile(starting_files[0], root_dir + '/output/test' + str(subdir) + '/' + starting_files[0])
    os.chdir(root_dir + '/output/test' + str(subdir))


@follows(setup5)
@transform(starting_files[0], filter=suffix(".vcf.gz"), output=".recode.vcf")
def run_vcf_to_filter5(vcf_in, seq_out):
    args = ['--vcf', starting_files[0], '--out-format', 'vcf', '--out-prefix', 'merged_chr1_10000', '--overwrite']
    vcf_filter.run(args)


@follows(run_vcf_to_filter5)
@transform(run_vcf_to_filter5, filter=suffix(".recode.vcf"), output=".windowed.pi")
def run_vcf_to_calc5(f_in, f_out):
    args = ['--vcf', f_in, '--calc-statistic', 'window-pi', '--out-prefix', 'merged_chr1_10000',
            '--statistic-window-size', '10000', '--statistic-window-step', '20000', '--overwrite']
    vcf_calc.run(args)


@follows(run_vcf_to_calc5)
@transform(run_vcf_to_calc5, filter=suffix(".windowed.pi"), output=".sampled")
def run_stat_to_sampler5(f_in, f_out):
    args = ['--statistic-file', f_in, '--out-prefix', 'merged_chr1_10000', '--calc-statistic', 'window-pi',
            '--sampling-scheme', 'uniform', '--uniform-bins', '5', '--sample-size', '5', '--random-seed', '100', '--overwrite']
    stat_sampler.run(args)


@follows(run_vcf_to_calc5)
@transform(run_vcf_to_calc5, filter=suffix(".windowed.pi"), output=".vcf.gz")
def run_vcf_to_split5(f_in, f_out):
    args = ['--vcf', starting_files[0], '--model-file', input_dir + 'input.model', '--model', '2Pop', '--split-file',
            f_in, '--split-method', 'statistic-file', '--statistic-window-size', '10000', '--overwrite']
    vcf_split.run(args)


# @follows(run_vcf_to_split5)
# @transform(run_vcf_to_split5, filter=suffix(".vcf.gz"), output=".out")
# def run_vcf_to_phase5(f_in, f_out):
#     args = ['--vcf', starting_files[0], '--beagle-path', input_dir, '--out-prefix', 'merged_chr1_10000',
#             '--phase-algorithm', 'beagle', '--overwrite']
#     vcf_phase.run(args)


# ======================================================================================================================
# pipeline sequence: vcf_filter -> vcf_calc(window-pi) -> stat_sampler(random) -> vcf_split -> vcf_phase


@follows(run_vcf_to_split5)
def setup6():
    global subdir
    subdir = subdir + 1
    shutil.copyfile(starting_files[0], root_dir + '/output/test' + str(subdir) + '/' + starting_files[0])
    os.chdir(root_dir + '/output/test' + str(subdir))


@follows(setup6)
@transform(starting_files[0], filter=suffix(".vcf.gz"), output=".recode.vcf")
def run_vcf_to_filter6(vcf_in, seq_out):
    args = ['--vcf', starting_files[0], '--out-format', 'vcf', '--out-prefix', 'merged_chr1_10000', '--overwrite']
    vcf_filter.run(args)


@follows(run_vcf_to_filter6)
@transform(run_vcf_to_filter6, filter=suffix(".recode.vcf"), output=".windowed.pi")
def run_vcf_to_calc6(f_in, f_out):
    args = ['--vcf', f_in, '--calc-statistic', 'window-pi', '--out-prefix', 'merged_chr1_10000',
            '--statistic-window-size', '10000', '--statistic-window-step', '20000', '--overwrite']
    vcf_calc.run(args)


@follows(run_vcf_to_calc6)
@transform(run_vcf_to_calc6, filter=suffix(".windowed.pi"), output=".sampled")
def run_stat_to_sampler6(f_in, f_out):
    args = ['--statistic-file', f_in, '--out-prefix', 'merged_chr1_10000', '--calc-statistic', 'window-pi',
            '--sampling-scheme', 'random', '--sample-size', '5', '--random-seed', '100', '--overwrite']
    stat_sampler.run(args)


@follows(run_vcf_to_calc6)
@transform(run_vcf_to_calc6, filter=suffix(".windowed.pi"), output=".vcf.gz")
def run_vcf_to_split6(f_in, f_out):
    args = ['--vcf', starting_files[0], '--model-file', input_dir + 'input.model', '--model', '2Pop', '--split-file',
            f_in, '--split-method', 'statistic-file', '--statistic-window-size', '10000', '--overwrite']
    vcf_split.run(args)


# @follows(run_vcf_to_split6)
# @transform(run_vcf_to_split6, filter=suffix(".vcf.gz"), output=".out")
# def run_vcf_to_phase6(f_in, f_out):
#     args = ['--vcf', starting_files[0], '--beagle-path', input_dir, '--out-prefix', 'merged_chr1_10000',
#             '--phase-algorithm', 'beagle', '--overwrite']
#     vcf_phase.run(args)


# ======================================================================================================================
pipeline_run(pipeline="main")
