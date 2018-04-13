#Files downloaded from UCSC table browser

DIR="/home/jared/workspace/projects_ppp/african_simons/singlescript/"
SDIR="/home/jared/workspace/projects_ppp/african_simons/"

#Gene file, from UCSC
UCSC_GENE_FILE=$DIR"human_genes.bed"
#Repeat file, from UCSC
UCSC_REPEAT_FILE=$DIR"human_repeats.bed"

#Inverted files names
GENE_INVERT=$DIR"human_genes_inverted.bed"
REPEAT_INVERT=$DIR"human_repeats_inverted.bed"

#Files for use in multiple steps
VCF=$SDIR"jodyset2_minimerge.vcf.gz"
REF=$SDIR"human_g1k_v37.fasta"
MODEL=$DIR"human_pops.model"
LOCI=200

#Step 1: Invert files to get neutral regions
echo "Inverting files"
python find_intergenic_bed.py --regionf $UCSC_GENE_FILE --pad 10000 --out $GENE_INVERT --regcol 4,5,2  #--pad is used to get intergenic regions away from genes
python find_intergenic_bed.py --regionf $UCSC_REPEAT_FILE --out $REPEAT_INVERT --regcol 6,7,5

#Step 2: Use bedtools intersect to merge valid regions
INT_BED=$DIR"human_premissing.bed"
echo "Merging inverted files"
bedtools intersect -a $GENE_INVERT -b $REPEAT_INVERT | sort -V > $INT_BED

#3: Get areas without missing data from VCF input
echo "Generating regions without missing data"
MISSING_BED=$DIR"human_missing_regions.bed"
zcat $VCF | python get_nonmissing_chunks.py | awk '{ print $1,$3,$4 }' | tr " " "\t" > $MISSING_BED

#4: Merge missing data into bed file for check for informative sites
echo "Generating regions for informative check"
PREINFORM_BED=$DIR"human_preinformative.bed"
PREINFORM_SUBSAMP=$DIR"human_preinformative_subsamp.bed"
bedtools intersect -a $INT_BED -b $MISSING_BED | sort -V | grep -v "X" | grep -v "Y" > $PREINFORM_BED #Grep to remove X/Y chromosomes
GOOD_BED_REGIONS=$DIR"informative_regions.bed"

#To speed up informative region check, subsample 10x or more of desired final loci
shuf -n 5000 $PREINFORM_BED | sort -V > $PREINFORM_SUBSAMP

#5: Search subsampled intervals to make sure each has at least 3 informative sites
echo "Parsing for informative regions"
python informative_filter.py --vcf $VCF --bed $PREINFORM_SUBSAMP --parsecpg $REF --minsites 3 > $GOOD_BED_REGIONS

#6: Select 200 loci at random from list
SAMPLE_REGIONS=$DIR"final_sampled_regions.bed"
shuf -n $LOCI $GOOD_BED_REGIONS | sort -V > $SAMPLE_REGIONS


#7: From main VCF, pull regions and create a VCF file for each region with informative, biallelic, non-CpG sites
echo "Generating sub-VCFs"
mkdir $DIR"vcf/"
python vcf_from_regions.py $VCF --rl $SAMPLE_REGIONS --output $DIR"vcf/human_run_" --multi-out --parsecpg $REF --remove-indels --remove-multi

FILEPATH_F=$DIR"human_run_filepaths.txt"
REGION_F=$DIR"extended_human_regions.txt"
OUT_F=$DIR"human_pops.ima.u"

rm $FILEPATH_F
rm $REGION_F
#For each region:
for i in $(seq 1 $LOCI)
do
    START_F=$DIR"vcf/human_run_region"$i".vcf"
    PHASE_F=$DIR"vcf/human_run_phased_region"$i".vcf"
    GAM_F=$DIR"vcf/human_run_4g_region"$i".vcf"
    

    #Phase the VCF file
    cd ../andrew
    python vcf_phase.py --vcf $START_F --out $PHASE_F --phase-algorithm shapeit
    cd ../jared
    
    #Check for a subregion that passes the 4-gamete test
    python four_gamete_pysam.py --vcfname $PHASE_F --out $GAM_F --4gcompat --reti --right --numinf 2
    #Extend regions for more accurate lengths of loci
    python extend_with_vcf.py $VCF $GAM_F >> $REGION_F
    echo $GAM_F >> $FILEPATH_F #For input list passed to 
done

#Creates input file. Takes VCF list, extended VCF regions, and model file

python vcf_to_ima.py --vcfs @$FILEPATH_F --bed $REGION_F --pop $MODEL --output $OUT_F



