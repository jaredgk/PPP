import os
import shutil
import subprocess

dir = 'great_ape_genome'
if os.path.exists(dir):
    shutil.rmtree(dir)
os.makedirs(dir)

subdir1 = 'great_ape_genome/Sampled_nonmissing'
if os.path.exists(subdir1):
    shutil.rmtree(subdir1)
os.makedirs(subdir1)

subdir2 = 'great_ape_genome/Phased'
if os.path.exists(subdir2):
    shutil.rmtree(subdir2)
os.makedirs(subdir2)


process = subprocess.Popen("python vcf_filter.py --vcf merged.vcf.gz --filter-max-missing 1.0 --filter-include-indv-file example/input/PaniscusTroglodytes.txt --filter-min-alleles 2 --filter-max-alleles 2 --out-format vcf.gz --out-prefix great_ape_genome/Pantrog_onlybiallelic_nomissing  --filter-exclude-chr chrX chrY --overwrite", shell=True, stdout=subprocess.PIPE)
process.wait()
if process == 0:
    process2 = subprocess.Popen("python vcf_calc.py --vcf great_ape_genome/Pantrog_onlybiallelic_nomissing.recode.vcf.gz --out-prefix great_ape_genome/fst.calc --calc-statistic windowed-weir-fst --model-file example/input/input.model --model 2Pop --statistic-window-size 10000 --statistic-window-step 20000 --overwrite", shell=True, stdout=subprocess.PIPE)
    process2.wait()
    if process2.returncode == 0:
        print("FST calculation successful.\n")
        process3 = subprocess.Popen("python vcf_sampler.py --vcf great_ape_genome/merged.vcf.gz --statistic-file great_ape_genome/fst.calc.windowed.weir.fst --out-format vcf.gz --calc-statistic windowed-weir-fst --sampling-scheme random --uniform-bins 5 --out-prefix great_ape_genome/sampled --overwrite", shell=True, stdout=subprocess.PIPE)
        process3.wait()
        if process3.returncode == 0:
            print("Sampling successful.\n")
            for i in range(5):
                process4 = subprocess.Popen("python vcf_filter.py --vcf great_ape_genome/Sample_Files/sampled_" + str(i) + ".vcf.gz --filter-max-missing 1.0 --out-format vcf.gz --out-prefix great_ape_genome/Sampled_nonmissing/sampled_nomissing_" + str(i) + " --overwrite", shell=True, stdout=subprocess.PIPE)
                process4.wait()
                if process4.returncode == 0:
                    print("Filter successful.\n")
                    process5 = subprocess.Popen("python vcf_phase.py --vcf great_ape_genome/Sampled_nonmissing/sampled_nomissing_" + str(i) + ".recode.vcf.gz --beagle-path example/input/ --out-prefix great_ape_genome/Phased/phased_" + str(i) + " --phase-algorithm beagle --beagle-burn-iter 12 --beagle-iter 24 --beagle-states 320 --beagle-window 20.0 --beagle-overlap 2.0 --beagle-error 0.0005 --beagle-step 0.05 --beagle-nsteps 5 --random-seed 100 --overwrite",shell=True, stdout=subprocess.PIPE)
        else:
            print("Error while sampling.\n")
    else:
        print("Error while FST calculation.\n")
else:
    print("Error while filtering.\n")


