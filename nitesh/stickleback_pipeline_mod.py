"""
Required input files:
1. BEAGLE.Combined66.snpEff1Effect.SNPs.OnlyInCombined66.masked68.nocnvrs.recode.vcf.gz
2. stickleback_genes.bed
3. sticklebackmodel.txt
"""
import sys
import os
import subprocess
import fnmatch
# import graph_plotter
import admixture
#Insert ../jared into python path's first position so imports start from there
sys.path.insert(0,os.path.abspath(os.path.join(os.pardir,'pppipe')))

import four_gamete_pysam, vcf_to_ima
from logging_module import initLogger

#Importing Andrew's functions
sys.path.insert(0,os.path.abspath(os.path.join(os.pardir,'andrew')))

import vcf_filter, vcf_calc, vcf_sampler, convert

#Work for results, data_dir for storage
work_dir = '/home/students/saban001/PPP/nitesh/'
data_dir = '/home/students/saban001/PPP/nitesh/example/input/'
vcf_dir = work_dir+'stickleback_genome/'

#Different way of doing this, but your way is perfectly fine too, just personal preference
if not os.path.exists(vcf_dir):
    os.makedirs(vcf_dir)
    os.makedirs(vcf_dir+'four_gamete/')
#dir = 'stickleback_genome'
# if os.path.exists(dir):
#     shutil.rmtree(dir)
#os.makedirs(dir)


#Was unzipped on my system, so add .gz if compressed
main_vcf_name = data_dir+'BEAGLE.Combined66.snpEff1Effect.SNPs.OnlyInCombined66.masked68.nocnvrs.recode.vcf.gz'



filtered_vcf_pref = work_dir+'stickleback_genome/stickleback_biallelic_nomissing'

stat_file_pref = work_dir+'stickleback_genome/fst.calc'
#Uncomment for logging/debugging
#initLogger()

vcf_filter.run(['--vcf',main_vcf_name,'--filter-max-missing','1.0','--filter-min-alleles','2','--filter-max-alleles','2','--out-format','vcf.gz','--out-prefix',filtered_vcf_pref,'--filter-exclude-bed',data_dir+'stickleback_genes.bed','--overwrite'])

vcf_calc.run(['--vcf',filtered_vcf_pref+'.recode.vcf.gz','--out-prefix',stat_file_pref,'--calc-statistic','windowed-weir-fst','--model','2Pop','--statistic-window-size','10000','--statistic-window-step','20000','--model-file',data_dir+'sticklebackmodel.txt','--overwrite'])

vcf_sampler.run(['--vcf',filtered_vcf_pref+'.recode.vcf.gz','--statistic-file',stat_file_pref+'.windowed.weir.fst','--out-format','vcf.gz','--calc-statistic','windowed-weir-fst','--sampling-scheme','uniform','--uniform-bins','5','--out-dir',work_dir+'stickleback_genome/Sample_Files','--overwrite'])

valid_files = []

for i in range(166):
    try:
        sampfn = work_dir+'stickleback_genome/Sample_Files/Sample_'+str(i)+'.vcf.gz'
        fgamfn = work_dir+'stickleback_genome/four_gamete/Sample_'+str(i)+'.vcf.gz'
        four_gamete_pysam.sample_fourgametetest_intervals(['--vcfname',sampfn,'--out',fgamfn,'--4gcompat','--reti','--right','--numinf','2'])
        os.system('tabix -p vcf ' + fgamfn)
    except IndexError:
        continue
    valid_files.append(fgamfn)


valid_filenames = ' '.join([str(x) for x in valid_files])


subprocess.Popen('vcf-concat ' + valid_filenames + ' | bgzip -c > ' + work_dir + 'stickleback_genome/four_gamete/sample_merged.vcf.gz', shell=True, stdout=subprocess.PIPE)

print('vcf-merge ' + valid_filenames + ' | bgzip -c > ' + work_dir + 'stickleback_genome/four_gamete/sample_merged.vcf.gz')

convert.run(['--vcf', work_dir + 'stickleback_genome/four_gamete/sample_merged.vcf.gz', '--out-format', 'binary-ped', '--out-prefix', work_dir+'stickleback_genome/stickleback'])

admixture.run(['--file', work_dir+'stickleback_genome/stickleback.bed', '--pop', '2', ])

graph_plotter.bar_plot('stickleback.2.Q')

ima_filenames = ' '.join([str(x) for x in valid_files])

ima_args = ['--vcfs']
ima_args.extend(valid_files)
ima_args.extend(['--pop',work_dir+'sticklebackmodel.txt','--out',work_dir+'ima_all_loci.ima.u'])


vcf_to_ima.vcf_to_ima(ima_args)

