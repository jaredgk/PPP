"""
Required input files:
1. BEAGLE.Combined66.snpEff1Effect.SNPs.OnlyInCombined66.masked68.nocnvrs.recode.vcf.gz
2. stickleback_genes.bed
3. sticklebackmodel.txt
"""

import os
import subprocess
import fnmatch

dir = 'stickleback_genome'
# if os.path.exists(dir):
#     shutil.rmtree(dir)
os.makedirs(dir)

subdir1 = 'stickleback_genome/four_gamete'
os.makedirs(subdir1)


# Filtering out non biallelic and missing chromosomes
process = subprocess.Popen("python vcf_filter.py --vcf BEAGLE.Combined66.snpEff1Effect.SNPs.OnlyInCombined66.masked68.nocnvrs.recode.vcf.gz \
                            --filter-max-missing 1.0 --filter-min-alleles 2 --filter-max-alleles 2 --out-format vcf.gz \
                            --out-prefix stickleback_genome/stickleback_onlybiallelic_nomissing  --filter-exclude-bed stickleback_genes.bed --overwrite"
                           , shell=True, stdout=subprocess.PIPE)
process.wait()
if process.returncode == 0:

    # FST caculation
    process2 = subprocess.Popen("python vcf_calc.py --vcf stickleback_genome/stickleback_onlybiallelic_nomissing.recode.vcf.gz \
                                --out-prefix stickleback_genome/fst.calc --calc-statistic windowed-weir-fst --model 2Pop \
                                --statistic-window-size 10000 --statistic-window-step 20000 --model-file example/input/sticklebackmodel.txt --overwrite",
                                shell=True, stdout=subprocess.PIPE)
    process2.wait()
    if process2.returncode == 0:
        print("FST calculation successful.\n")

        # Sampling
        process3 = subprocess.Popen("python vcf_sampler.py --vcf BEAGLE.Combined66.snpEff1Effect.SNPs.OnlyInCombined66.masked68.nocnvrs.recode.vcf.gz \
                                    --statistic-file stickleback_genome/fst.calc.windowed.weir.fst --out-format vcf.gz \
                                    --calc-statistic windowed-weir-fst --sampling-scheme random --uniform-bins 5 --out-dir \
                                    stickleback_genome/Sample_Files --overwrite", shell=True, stdout=subprocess.PIPE)
        process3.wait()
        if process3.returncode == 0:
            print("Sampling successful.\n")
            for i in range(200):
                # four gamate test
                process4 = subprocess.Popen("python ../jared/four_gamete_pysam.py --vcfname stickleback_genome/Sample_Files/Sample_" + str(i) + ".vcf.gz \
                                            --out-prefix stickleback_genome/four_gamete/Sample_" + str(i) + " --4gcompat --reti --right --numinf 2",
                                            shell=True, stdout=subprocess.PIPE)

            loci_List = fnmatch.filter(os.listdir('stickleback_genome/Sample_Files'), '*.vcf.gz')  # gettibng list of sampled loci files
            new_loci_List = ["stickleback_genome/Sample_Files/" + s for s in loci_List]  # Adding relative path to each file

            process5 = subprocess.Popen("python ../jared/vcf_to_ima.py --vcfs " + (' '.join(str(x) for x in new_loci_List)) +
                                        " --ref ../nitesh/hg18.fasta --pop example/input/sticklebackmodel.txt --out \
                                        stickleback_genome/ima_all_loci.fasta --poptag 2Pop", shell=True, stdout=subprocess.PIPE)
            if process5.returncode == 0:
                # code for ima2p
                pass

        else:
            print("Error while sampling.\n")
    else:
        print("Error while FST calculation.\n")
else:
    print("Error while filtering.\n")