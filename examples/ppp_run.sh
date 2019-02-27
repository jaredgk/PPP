#DATA_DIR should contain source VCF and reference FASTA if applicable
DATA_DIR="/home/jared/"
#WORK_DIR contains model file and all intermediate/output files
WORK_DIR="/home/jared/"
TAG="tagset"

#must be bgzipped/tabixed
VCF=$DATA_DIR"vcfname.vcf.gz"
#must be faidx'd, can be bgzipped
REF=$DATA_DIR"ref.fasta"
SCRIPT_DIR="/home/jared/workspace/ppp/pgpipe/"


LOCI=200
MODEL=$WORK_DIR"model_file.model"
IMAOUT=$WORK_DIR$TAG".ima.u"

#Regions fit for analysis, may be filtered by genic regions, repeat regions,
#areas without missing data
PREINF_REGIONS=$WORK_DIR"preinformative_regiond.bed"
INFORMATIVE_REGIONS=$WORK_DIR$TAG"_informative_regions.bed"


FILTER_OPTS=" --model "$MODEL" --remove-indels --remove-multi --remove-missing 0"

VCF_DIR=$WORK_DIR"sub_vcfs/"
mkdir $VCF_DIR

#Finds regions that have enough informative sites to be likely to pass four-gamete test
python $SCRIPT_DIR"informative_filter.py" --vcf $VCF --bed $PREINF_REGIONS --minsites 3 --min-length 500 --randcount $LOCI --no-xy $FILTER_OPTS | sort -V > $INFORMATIVE_REGIONS

#Splits found regions into a VCF file per region
python $SCRIPT_DIR"vcf_from_regions.py" $VCF --bed $INFORMATIVE_REGIONS --output $VCF_DIR$TAG"_cpgfilter_" --multi-out --informative-count 1 $FILTER_OPTS

if [[ $? != 0 ]]
then
    exit
fi

EXTEND_REGIONS=$WORK_DIR$TAG"_final_regions.bed"
GOOD_FILES=$WORK_DIR$TAG"_final_filepaths.txt"

rm $EXTEND_REGIONS
rm $GOOD_FILES

for i in $(seq 1 $LOCI)
do
    echo $i
    FIRST_F=$VCF_DIR$TAG"_cpgfilter_region"$i".vcf"
    PHASE_F=$VCF_DIR$TAG"_phased_region"$i".vcf"
    FOURG_F=$VCF_DIR$TAG"_fourgamete_region"$i".vcf"
    
    #Phase region
    python $SCRIPT_DIR"vcf_phase.py" --vcf $FIRST_F --out $PHASE_F --phase-algorithm shapeit --overwrite --out-format vcf
    #Find subregion in VCF that has no recombination events amongst contained SNPs.
    #If no subregion is found, will not include locus in IMa input
    python $SCRIPT_DIR"four_gamete_pysam.py" --vcfs $PHASE_F --out $FOURG_F --4gcompat --reti --rani --numinf 2

    if [[ $? == 0 ]]
    then
        #Calculate halfway distance between SNPs flanking endpoints of subregion
        python $SCRIPT_DIR"extend_with_vcf.py" $VCF $FOURG_F >> $EXTEND_REGIONS
        echo $FOURG_F >> $GOOD_FILES
    fi



done
#Create IMa file
python $SCRIPT_DIR"vcf_to_ima.py" --vcfs @$GOOD_FILES --bed $EXTEND_REGIONS --pop $MODEL --out $IMAOUT