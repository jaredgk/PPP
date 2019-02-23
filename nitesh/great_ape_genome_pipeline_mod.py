"""
Required input files:
1. merged.vcf.gz
2. input.model
"""

import sys
import os
# import graph_plotter
import admixture
import subprocess

#Insert ../jared into python path's first position so imports start from there
sys.path.insert(0,os.path.abspath(os.path.join(os.pardir,'pppipe')))

import four_gamete_pysam, vcf_to_ima
from logging_module import initLogger

#Importing Andrew's functions
sys.path.insert(0,os.path.abspath(os.path.join(os.pardir,'andrew')))

import vcf_filter, vcf_calc, vcf_sampler, vcf_phase, convert

#Importing Nitesh's functions
sys.path.insert(0,os.path.abspath(os.path.join(os.pardir,'nitesh')))

#Work for results, data_dir for storage
work_dir = '/home/students/saban001/PPP/nitesh/'
data_dir = '/home/students/saban001/PPP/nitesh/example/input/'
vcf_dir = work_dir+'great_ape_genome/'

#Different way of doing this, but your way is perfectly fine too, just personal preference
if not os.path.exists(vcf_dir):
    os.makedirs(vcf_dir)
    os.makedirs(vcf_dir+'four_gamete/')
    os.makedirs(vcf_dir + 'Sampled_nonmissing/')
    os.makedirs(vcf_dir + 'Phased/')

main_vcf_name = data_dir+'merged.vcf.gz'

filtered_vcf_pref = work_dir+'great_ape_genome/Pantrog_onlybiallelic_nomissing'

stat_file_pref = work_dir+'great_ape_genome/fst.calc'

#Uncomment for logging/debugging
#initLogger()

vcf_filter.run(['--vcf', main_vcf_name, '--filter-max-missing', '1.0', '--filter-include-indv-file', data_dir+'PaniscusTroglodytes.txt', '--filter-min-alleles', '2', '--filter-max-alleles', '2', '--out-format', 'vcf.gz', '--out-prefix', filtered_vcf_pref, '--filter-exclude-chr', 'chrX', 'chrY', '--overwrite'])

vcf_calc.run(['--vcf', filtered_vcf_pref + '.recode.vcf.gz', '--out-prefix', stat_file_pref, '--calc-statistic', 'windowed-weir-fst', '--model', '2Pop', '--statistic-window-size', '10000', '--statistic-window-step', '20000', '--model-file', data_dir + 'input.model', '--overwrite'])

vcf_sampler.run(['--vcf', filtered_vcf_pref + '.recode.vcf.gz', '--statistic-file', stat_file_pref + '.windowed.weir.fst', '--out-format', 'vcf.gz', '--calc-statistic', 'windowed-weir-fst', '--sampling-scheme', 'random', '--uniform-bins', '5', '--out-dir', work_dir + 'great_ape_genome/Sample_Files', '--overwrite'])

valid_files = []
phased_files = []

for i in range(200):
    try:
        sampfn = work_dir+'great_ape_genome/Sample_Files/Sample_'+str(i)+'.vcf.gz'
        nomiss_sampfn = work_dir+'great_ape_genome/Sampled_nonmissing/Sample_nomissing_'+str(i)
        phased_sampfn = work_dir+'great_ape_genome/Phased/phased_'+str(i)
        fgamfn = work_dir+'great_ape_genome/four_gamete/Sample_'+str(i)+'.vcf.gz'

        vcf_filter.run(['--vcf', sampfn, '--filter-max-missing', '1.0', '--out-format', 'vcf.gz', '--out-prefix', nomiss_sampfn, '--overwrite'])
        vcf_phase.run(['--vcf', nomiss_sampfn+'.recode.vcf.gz', '--beagle-path', data_dir, '--out-prefix', phased_sampfn, '--phase-algorithm', 'beagle', '--beagle-burn-iter', '12', '--beagle-iter', '24', '--beagle-states', '320', '--beagle-window', '20.0', '--beagle-overlap', '2.0', '--beagle-error', '0.0005', '--beagle-step', '0.05', '--beagle-nsteps', '5', '--random-seed', '100', '--overwrite'])
        os.system('tabix -p vcf ' + phased_sampfn + '.vcf.gz')
        four_gamete_pysam.sample_fourgametetest_intervals(['--vcfname', phased_sampfn+'.vcf.gz', '--out', fgamfn, '--4gcompat', '--reti', '--right', '--numinf', '2'])

    except IndexError:
        continue
    valid_files.append(fgamfn)
    phased_files.append(phased_sampfn+ ".vcf.gz")

phased_filenames = ' '.join([str(x) for x in phased_files])

subprocess.Popen('vcf-concat ' + phased_filenames + ' | bgzip -c > ' + work_dir + 'great_ape_genome/Phased/phased_merged.vcf.gz', shell=True, stdout=subprocess.PIPE)

convert.run(['--vcf', work_dir + 'great_ape_genome/Phased/phased_merged.vcf.gz', '--out-format', 'binary-ped', '--out-prefix', work_dir+'great_ape_genome/great_ape'])

admixture.run(['--file', work_dir+'great_ape_genome/great_ape.bed', '--pop', '2', ])

graph_plotter.bar_plot('great_ape.2.Q')

ima_filenames = ' '.join([str(y) for y in valid_files])

ima_args = ['--vcfs']
ima_args.extend(valid_files)
ima_args.extend(['--pop', data_dir + 'input.model', '--out', work_dir + 'ima_all_loci.ima.u'])

vcf_to_ima.vcf_to_ima(ima_args)

graph_plotter.do_plots(work_dir + 'ima_all_loci.ima.u')

