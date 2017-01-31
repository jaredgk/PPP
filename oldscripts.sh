GATK Notes

GenerateAltAlleleFasta
look this page up:
http://gatkforums.broadinstitute.org/discussion/1493/generatealtallelefasta

needs SVToolKit


java -jar GenomeAnalysisTK.jar -T SelectVariants --variant
stickleback.vcf -R stickleback.fasta -o BS2.vcf -se 'group*' -sn BS2

java -jar GenomeAnalysisTK.jar -T SelectVariants --variant
stickleback.vcf -R stickleback.fasta -o BS4.vcf -se 'group*' -sn BS4

java -jar GenomeAnalysisTK.jar -T SelectVariants --variant
stickleback.vcf -R stickleback.fasta -o BS6.vcf -se 'group*' -sn BS6

java -jar GenomeAnalysisTK.jar -T SelectVariants --variant
stickleback.vcf -R stickleback.fasta -o BS8.vcf -se 'group*' -sn BS8

java -jar GenomeAnalysisTK.jar -T SelectVariants --variant
stickleback.vcf -R stickleback.fasta -o BS10.vcf -se 'group*' -sn BS10

java -jar GenomeAnalysisTK.jar -T SelectVariants --variant
stickleback.vcf -R stickleback.fasta -o BS12.vcf -se 'group*' -sn BS12

java -jar GenomeAnalysisTK.jar -T SelectVariants --variant
stickleback.vcf -R stickleback.fasta -o BS1.vcf -se 'group*' -sn BS1

java -jar GenomeAnalysisTK.jar -T SelectVariants --variant
stickleback.vcf -R stickleback.fasta -o BS3.vcf -se 'group*' -sn BS3

java -jar GenomeAnalysisTK.jar -T SelectVariants --variant
stickleback.vcf -R stickleback.fasta -o BS5.vcf -se 'group*' -sn BS5

java -jar GenomeAnalysisTK.jar -T SelectVariants --variant
stickleback.vcf -R stickleback.fasta -o BS7.vcf -se 'group*' -sn BS7

java -jar GenomeAnalysisTK.jar -T SelectVariants --variant
stickleback.vcf -R stickleback.fasta -o BS9.vcf -se 'group*' -sn BS9

java -jar GenomeAnalysisTK.jar -T SelectVariants --variant
stickleback.vcf -R stickleback.fasta -o BS11.vcf -se 'group*' -sn BS11

--
This didn't work as I wanted:
Used GBSContextualSeq.py from Thomas Kono at UMN, with length of
flanking region 200 bases
This pulls out individual sequences around each SNP for each individual.
Thereon, I used this to pull out only fasta sequences with intersecting ID's

comm -12 <(comm -12 <(comm -12 <(comm -12 <(comm -12 <(comm -12 <(comm
-12 <(comm -12 <(comm -12 <(comm -12 <(comm -12 <(grep ">"
BS43_SNPs.fasta | sort) <(grep ">" BS44_SNPs.fasta | sort)) <(grep ">"
BS45_SNPs.fasta | sort)) <(grep ">" BS46_SNPs.fasta | sort)) <(grep
">" BS47_SNPs.fasta | sort)) <(grep ">" BS48_SNPs.fasta | sort))
<(grep ">" BS49_SNPs.fasta | sort)) <(grep ">" BS50b_SNPs.fasta |
sort)) <(grep ">" BS51_SNPs.fasta | sort)) <(grep ">" BS52b_SNPs.fasta
| sort)) <(grep ">" BS53_SNPs.fasta | sort)) <(grep ">"
BS54_SNPs.fasta | sort) > finalloci
--
Instead just settled for locus identification using
java -jar GenomeAnalysisTK.jar -R stickleback.fasta -T
FastaAlternateReferenceMaker -o bs43.fasta -L caregions.list --variant
bs43.vcf --use_IUPAC_sample BS43
repeated for all 10 loci

Then had to do some manipulations to the input fasta files:

1) added samplename_to every locus %s/>/>bssamplename_/g
2) sed '/bs43*/a_' bs43.fasta > blah
3) vim blah, %s/\n_/_/g
4) for i in {1..100}
do
for filename in bs43 bs44 bs45 bs46 bs47 bs48 bs49 bs50b bs51 bs52b bs53 bs54
do
grep -A 167 -h "$filename"_"$i"_ "$filename".fasta >> fasta/locus"$i".fasta
done
done

Now stickleback.msa contains 100 loci, multiple sequence alignments.
But I think this isn't needed for the python script from Jody



Now finalloci has all unique loci
Onto concatenating it all to a multiple sequence alignment to be phased, etc.
NONE OF THIS WORKED>>>> So!...
--
Nope - instead decided to work directly with VCF files. I created a
VCF to BED converter instead - small script that does something like
this:

sed -e 's/locus//' allca.vcf | awk '{OFS="\t"; if (!/^#/) {print
$1,$2-1,$2,$4,$5,$10.$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21}} >
allca.bed
Then series of manipulations:
%s/:[0-9]*//g
%s/,[0-9]*//g
%s/\.//g
%s/\/[0-9]*//g
%s/|/ /g
Also had to add 2 lines to the beginning.
allcaloci
chromosome  begin   end ref variantbase BS43    BS44    BS45    BS46
 BS47    BS48    BS49    BS50b   BS51    BS52b   BS53    BS54
Note that this file now has all sites - 10k for each locus. We only
want variant sites. So had to do some manipulation in R for this.
In R:
x=apply(allca[,6:17],1,unique)
y<-which(sapply(x,length)>1)
nodups<-allca[y,]
 write.table(nodups,"nodups.bed")

Then after that, I had to take the allca.bed file and sort it in Excel
- I can perhaps do this in R as well. Had to sort this based on the
chromosome name first, then followed by the starting site. Save this
as a separate BED file - call this allca_sorted.bed

Also, need to sort the variants file - I initially awk-ed the
samps.txt file (which was produced using vcftools while computing FST,
and then sampled from). Then sorted this also based on the same two
above fields.

This gives rise to allca_sorted.bed, and samps_vars (which I just
awked from the samps_sorted.txt).

More R:
samps_sorted<-read.table("samps_sorted",header=TRUE)
g<-data.frame(samps_sorted$CHROM,samps_sorted$BIN_START,samps_sorted$BIN_END,samps_sorted$BIN_END-samps_sorted$BIN_START,s)
 write.table(g,"g.txt")

Then awk only remaining columns to h.txt
Thereon, create a folder called "bedfiles". Then run these scripts to
create all BED files required by phase_four_gametes.py:

i=1
curr=0
for vals in $(cat samps_vars);
do
head -n 1 allca_sorted.bed > bedfiles/locus"$i".bed;
curr=$((curr+vals));
v=$((vals-1));
head -n $curr allca_sorted.bed | tail -n $v >> bedfiles/locus"$i".bed;
i=$((i+1)); done


for i in {1..100}
do
head -n $i h.txt | tail -n 1 | cat - bedfiles/locus"$i".bed > bedfiles/text.temp
mv bedfiles/text.temp bedfiles/locus"$i".bed
done

for i in {1..100}
do
head -n 1 locus"$i".bed > text.temp
sed '1d' locus"$i".bed | awk '{print
$1,$2,$3,$4,$5,$6,$8,$10,$12,$14,$16,$7,$9,$11,$13,$15,$17}' >>
text.temp
mv text.temp locus"$i".bed
done
dos2unix *.bed

Had to do some manipulations to the python script as well - please see
script for more details on this. Specifically:

1) change mutation rate
2) change the number of individuals per population
3) change population names
4) folders
5) file names
6)
And voila! Have .u files!

Stratified sampling - 5/28/2015
Had to create bins of Fst values first

cafst_cleaned<-read.table("Combined66_groups_sorted_cleaned_fst",header=FALSE)
cafst_cleaned$BINS<-c(cut(cafst_cleaned$V7,breaks=10))
Now BINS will contain 10 different bins.
I want 10 loci from each bin.
ids<-array(dim=c(1,100))
ids[1:10]<-sample(which(cafst_cleaned$BINS==1),10)
and so on...can run in loop
Now ids should have ID's of 10 sampled from each bin.
This was written as a new file
write.table(cafst_cleaned[ids,],file="ca_10fromeachbin_cleaned.txt",quote=FALSE,row.names=FALSE)

ca_10fromeachbin_cleaned.txt now has the correct 100 loci with
stratified sampling. Note however that the last bin only had 6 loci in
it.
So I sampled 14 from bin 9, and 6 from bin 10 instead.

cat ca_10fromeachbin_cleaned.txt | awk '{print $1":"$2"-"$3}' >
ca_10fromeachbin_cleaned_ranges.list
this file now contains all ranges

java -jar GenomeAnalysisTK.jar -T SelectVariants --variant
stickleback.vcf -R stickleback.fasta -o ca_10each_cleaned.vcf -se
'group*' -L ca_10fromeachbin_cleaned_ranges.list -sn BS44 -sn BS46 -sn
BS48 -sn BS50b -sn BS52b -sn BS54 -sn BS43 -sn BS45 -sn BS47 -sn BS49
-sn BS51 -sn BS53

Now ca_10each_cleaned.vcf contains all the sampled loci from
stratified sampling. Onto converting into BED files, etc.

Ran into troubles processing this single vcf/bed file - so just going
to create vcf files for all loci, then do the processing.

i=1
for vals in $(cat ca_10fromeachbin_cleaned_ranges.list);
do
echo $vals > temp.list
java -jar GenomeAnalysisTK.jar -T SelectVariants --variant
stickleback.vcf -R stickleback.fasta -o vcffiles_1/locus"$i".vcf -se
'group*' -L temp.list -sn BS44 -sn BS46 -sn BS48 -sn BS50b -sn BS52b
-sn BS54 -sn BS43 -sn BS45 -sn BS47 -sn BS49 -sn BS51 -sn BS53
i=$((i+1)); done
rm temp.list

okay now I have 100 vcf files, idx files, etc. Onto creating BED files
from these directly.


for i in {1..100}
do
sed -e 's/locus//' locus"$i".vcf | awk '{OFS="\t"; if (!/^#/) {print
$1,$2-1,$2,$4,$5,$11,$13,$15,$17,$19,$21,$10,$12,$14,$16,$18,$20}}' >
temp
sed -e 's/:[0-9]*//g' temp | sed -e 's/,[0-9]*//g' | sed -e 's/\.//g'
| sed -e 's/\/[0-9]*//g' | sed -e 's/|//g' >>
../bedfiles_1/locus"$i".bed
done

then loop inside R

for (i in 1:100) {
locusname<-paste("locus",i,".bed",sep="")
ca_10<-read.table(locusname,header=TRUE,skip=1)
x=apply(ca_10[,6:17],1,unique)
y<-which(sapply(x,length)>1)
nodups<-ca_10[y,]
write.table(nodups,locusname,col.names=FALSE,row.names=FALSE)
}

then in shell, remove all " values, also have to change all 0's to 00, 1's to 01

for i in {1..100}
do
sed -e 's/"//g' locus"$i".bed | sed -e 's/ 0 / 00 /g' | sed -e 's/ 1 /
01 /g' > temp
mv temp locus"$i".bed
done

for i in {1..100}
do
perl -p -e 's/ 0\n/ 00\n/' locus"$i".bed > temp
mv temp locus"$i".bed
done

for i in {1..100}
do
perl -p -e 's/ 1\n/ 01\n/' locus"$i".bed > temp
mv temp locus"$i".bed
done

for i in {1..100}
do
sed -e 's/ 0 / 00 /g' locus"$i".bed > temp
mv temp locus"$i".bed
done

for i in {1..100}
do
perl -p -e 's/ 1 / 01 /' locus"$i".bed > temp
mv temp locus"$i".bed
done



phew!!! done! now I need to write headers

i=1
for vals in $(cat ../ca_10fromeachbin_cleaned_ranges.list);
do
sed -i '1s/^/chromosome  begin   end ref variantbase  BS44 BS46 BS48
BS50b BS52b BS54 BS43 BS45 BS47 BS49 BS51 BS53\n/' locus"$i".bed
vars=$(wc -l locus"$i".bed | awk '{print $1}');
vars=$((vars-1));
sed -i "1s/^/$vals 9999 $vars\n/" locus"$i".bed
sed -i 's/:/ /g' locus"$i".bed
sed -i 's/group//g' locus"$i".bed
sed -i 's/-/ /g' locus"$i".bed
i=$((i+1));
done

DONE DONE DONE! NOW I have individual BED files to be run through the
fourgamete test, etc.
Once I run the python script:

for i in {1..100}
do
cat locus"$i".u >> ca_10each.u
done
then :s/ 6 6 / 12 12 /g
--


Now onto sampling 50 high Fst values, 50 low Fst values
ids2<-array(dim=c(1,100))
ids2[1:50]<-sample(which(cafst_cleaned$BINS==c(1,2,3,4,5)),50)
ids2[51:100]<-sample(which(cafst_cleaned$BINS==c(7,8,9,10)),50)
Now IDs

write.table(cafst_cleaned[ids2,],file="ca_50highlow_cleaned.txt",quote=FALSE,row.names=FALSE)

cat ca_50highlow_cleaned.txt | awk '{print $1":"$2"-"$3}' >
ca_50highlow_cleaned_ranges.list

i=1
for vals in $(cat ca_50highlow_cleaned_ranges.list);
do
echo $vals > temp.list
java -jar GenomeAnalysisTK.jar -T SelectVariants --variant
stickleback.vcf -R stickleback.fasta -o vcffiles_2/locus"$i".vcf -se
'group*' -L temp.list -sn BS44 -sn BS46 -sn BS48 -sn BS50b -sn BS52b
-sn BS54 -sn BS43 -sn BS45 -sn BS47 -sn BS49 -sn BS51 -sn BS53
i=$((i+1)); done
rm temp.list

for i in {1..100}
do
sed -e 's/locus//' locus"$i".vcf | awk '{OFS="\t"; if (!/^#/) {print
$1,$2-1,$2,$4,$5,$11,$13,$15,$17,$19,$21,$10,$12,$14,$16,$18,$20}}' >
temp
sed -e 's/:[0-9]*//g' temp | sed -e 's/,[0-9]*//g' | sed -e 's/\.//g'
| sed -e 's/\/[0-9]*//g' | sed -e 's/|//g' >>
../bedfiles_2/locus"$i".bed
done

then loop inside R

for (i in 1:100) {
locusname<-paste("locus",i,".bed",sep="")
ca_10<-read.table(locusname,header=TRUE,skip=1)
x=apply(ca_10[,6:17],1,unique)
y<-which(sapply(x,length)>1)
nodups<-ca_10[y,]
write.table(nodups,locusname,col.names=FALSE,row.names=FALSE)
}

then in shell, remove all " values, also have to change all 0's to 00, 1's to 01

for i in {1..100}
do
sed -e 's/"//g' locus"$i".bed | sed -e 's/ 0 / 00 /g' | sed -e 's/ 1 /
01 /g' > temp
mv temp locus"$i".bed
done

for i in {1..100}
do
perl -p -e 's/ 0\n/ 00\n/' locus"$i".bed > temp
mv temp locus"$i".bed
done

for i in {1..100}
do
perl -p -e 's/ 1\n/ 01\n/' locus"$i".bed > temp
mv temp locus"$i".bed
done

for i in {1..100}
do
sed -e 's/ 0 / 00 /g' locus"$i".bed > temp
mv temp locus"$i".bed
done

for i in {1..100}
do
perl -p -e 's/ 1 / 01 /' locus"$i".bed > temp
mv temp locus"$i".bed
done

for i in {1..100}
do
sed -e 's/ 1 / 01 /g' locus"$i".bed > temp
mv temp locus"$i".bed
done


phew!!! done! now I need to write headers

i=1
for vals in $(cat ../ca_50highlow_cleaned_ranges.list);
do
sed -i '1s/^/chromosome  begin   end ref variantbase  BS44 BS46 BS48
BS50b BS52b BS54 BS43 BS45 BS47 BS49 BS51 BS53\n/' locus"$i".bed
vars=$(wc -l locus"$i".bed | awk '{print $1}');
vars=$((vars-1));
sed -i "1s/^/$vals 9999 $vars\n/" locus"$i".bed
sed -i 's/:/ /g' locus"$i".bed
sed -i 's/group//g' locus"$i".bed
sed -i 's/-/ /g' locus"$i".bed
i=$((i+1));
done


--
samps1<-read.table("ca_10fromeachbin.txt",header=TRUE)
g<-data.frame(samps1$CHROM,samps1$BIN_START,samps1$BIN_END,samps1$BIN_END-samps1$BIN_START,samps1$N_VARIANTS)
write.table(g,"g1.txt",row.names=FALSE,quote=FALSE)
cat ca_10fromeachbin.txt | awk '{print $4}' > samps1_vars


java -jar GenomeAnalysisTK.jar -T SelectVariants --variant
stickleback.vcf -R stickleback.fasta -o ca_10each_cleaned.vcf -se
'group*' -L ca_10fromeachbin_regions_cleaned.txt -sn BS44 -sn BS46 -sn
BS48 -sn BS50b -sn BS52b -sn BS54 -sn BS43 -sn BS45 -sn BS47 -sn BS49
-sn BS51 -sn BS53


i=1
curr=0
for vals in $(cat samps1_vars);
do
head -n 1 allca_sorted.bed > bedfiles_1/locus"$i".bed;
curr=$((curr+vals));
v=$((vals-1));
head -n $curr allca_sorted.bed | tail -n $v >> bedfiles_1/locus"$i".bed;
i=$((i+1)); done


---

ARABIDOPSIS
To extract Fst:

First step was to create fasta files for all loci.
To do this, I had to manipulate all the *_4g.fs files:


for i in {0..28}
do
wordcount=$(wc -l disambig_aligned_Cluster_$i.fasta_4g.fs | awk '{print $1}');
wordcount=$((wordcount-6));
wordcount1=$((wordcount-1));
tail -$wordcount disambig_aligned_Cluster_"$i".fasta_4g.fs | head
-$wordcount1 > arabidopsis"$i".fasta
cat arabidopsis"$i".fasta | awk '{print ">"$1"\n"$2}' > temp
mv temp arabidopsis"$i".fasta
done


Inside R:

library(ape)
library(adegenet)
fsts<-rep(0,times=30)
for (x in 0:29) {
fname<-paste("arabidopsis",x,".fasta",sep="");
a0<-read.dna(fname,format="fasta")
a0_gi<-DNAbin2genind(a0)
#pop(a0_gi)<-rep("newpop",dim(a0)[1])
pops<-rep("newpop",times=dim(a0)[1])
for (i in 0:dim(a0)[1]) {
locname<-paste("Halleri",i,sep="")
j<-which(a0_gi@ind.names == locname)
if (length(j) > 0) {
pops[j] = "Halleri"
}
locname<-paste("Lyrata",i,sep="")
j<-which(a0_gi@ind.names == locname)
if (length(j) > 0) {
pops[j] = "Lyrata"
}
}
pop(a0_gi)=pops
fsts[x+1]=pairwise.fst(a0_gi)
}
--
ANOPHELES

library(ape)
library(adegenet)
fsts<-rep(0,times=37)
for (x in 1:36) {
fname<-paste("loc",x,"_pruned.fasta",sep="");
a0<-read.dna(fname,format="fasta")
a0_gi<-DNAbin2genind(a0)
pops<-rep("newpop",times=dim(a0)[1])
ms<-grep("M isolate",c(a0_gi@ind.names))
if (length(ms) == 0) {
ms<-grep("clone m",c(a0_gi@ind.names))
}
if (length(ms) == 0) {
ms<-grep("clone M",c(a0_gi@ind.names))
}
if (length(ms) == 0) {
ms<-grep("isolate M",c(a0_gi@ind.names))
}
if (length(ms) == 0) {
ms<-grep("M clone",c(a0_gi@ind.names))
}


ss<-grep("S isolate",c(a0_gi@ind.names))
if (length(ss) == 0) {
ss<-grep("clone s",c(a0_gi@ind.names))
}
if (length(ss) == 0) {
ss<-grep("clone S",c(a0_gi@ind.names))
}
if (length(ss) == 0) {
ss<-grep("isolate S",c(a0_gi@ind.names))
}
if (length(ss) == 0) {
ss<-grep("S clone",c(a0_gi@ind.names))
}


for (i in 1:length(ms)) {
pops[ms[i]]="M"
}
for (i in 1:length(ss)) {
pops[ss[i]]="S"
}
pop(a0_gi)=pops
fsts[x]=pairwise.fst(a0_gi)
}

--
STICKLEBACK

for i in {1..100}
do
wordcount=$(wc -l locus"$i"_4g.fs | awk '{print $1}');
wordcount=$((wordcount-6));
wordcount1=$((wordcount-1));
tail -$wordcount locus"$i"_4g.fs | head -$wordcount1 > stickleback"$i".fasta
cat stickleback"$i".fasta | awk '{print ">"$1"\n"$2}' > temp
mv temp stickleback"$i".fasta
done

in R

library(ape)
library(adegenet)
fsts<-rep(0,times=100)
for (x in 1:100) {
fname<-paste("stickleback",x,".fasta",sep="");
a0<-read.dna(fname,format="fasta")
a0_gi<-DNAbin2genind(a0)
#pop(a0_gi)<-rep("newpop",dim(a0)[1])
pops<-rep("newpop",times=dim(a0)[1])
for (i in 0:dim(a0)[1]) {
locname<-paste("L",i,sep="")
j<-which(a0_gi@ind.names == locname)
if (length(j) > 0) {
pops[j] = "Lake"
}
locname<-paste("R",i,sep="")
j<-which(a0_gi@ind.names == locname)
if (length(j) > 0) {
pops[j] = "River"
}
}
pop(a0_gi)=pops
fsts[x]=pairwise.fst(a0_gi)
}

--
MUSMUSCULUS dataset mining:

Got a hold of NEX files from C&H - same files as in the original paper.
Need FASTA files from this - so had to do some manipulation:

grep "AL[0-9]" XP_620246.nex | awk '{print ">"$1"\n"$2}' >> XP_620246.fasta
grep "D[0-9]" XP_620246.nex | awk '{print ">"$1"\n"$2}' >> XP_620246.fasta

etc

Then just ran my usual python scripts on this dataset.

To calculate Fsts:

library(ape)
library(adegenet)
fsts<-rep(0,times=14)
for (x in 1:14) {
fname<-paste("Locus",x,".fasta",sep="");
a0<-read.dna(fname,format="fasta")
a0_gi<-DNAbin2genind(a0)
#pop(a0_gi)<-rep("newpop",dim(a0)[1])
pops<-rep("newpop",times=dim(a0)[1])
for (i in 0:dim(a0)[1]) {
locname<-paste("ALtype",sep="")
j<-grep(locname,a0_gi@ind.names)
#j<-which(a0_gi@ind.names == locname)
if (length(j) > 0) {
pops[j] = "ALtype"
}
locname<-paste("Dtype",sep="")
j<-grep(locname,a0_gi@ind.names)
#j<-which(a0_gi@ind.names == locname)
if (length(j) > 0) {
pops[j] = "Dtype"
}
}
pop(a0_gi)=pops
fsts[x]=pairwise.fst(a0_gi)
}




--
HELICONIUS

Downloaded all BAM files from
https://usegalaxy.org/u/njnadeau/h/heliconius-sureselect-june-2011

Then had to merge BAM files using samtools:
had to download and concatenate both reference fastas - made heli.fasta
samtools mpileup -uf heli.fasta *.bam | bcftools view -bvcg - > var.raw.bcf
bcftools view var.raw.bcf > var.raw.vcf

--
This didn't work. So decided to find a different data set which
provided VCF files instead. This dataset comes from Martin et al.
(Genome-wide evidence for speciation with gene flow in Heliconius
butterflies)

Downloaded the VCF files in two parts from
datadryad.org/resource/doi:10.5061/dryad.dk712

Now just directly read the fst file into R for sampling according to bins.

fst<-read.csv("*.fst",header=TRUE)
fst$BINS<-c(cut(fst$am_timPFst,breaks=10))

ids<-array(dim=c(1,100))
ids[1:10]<-sample(which(fst$BINS==1),10)
ids[11:20]<-sample(which(fst$BINS==2),10)
ids[21:30]<-sample(which(fst$BINS==3),10)
ids[31:40]<-sample(which(fst$BINS==4),10)
ids[41:50]<-sample(which(fst$BINS==5),10)
ids[51:60]<-sample(which(fst$BINS==6),10)
ids[61:70]<-sample(which(fst$BINS==7),10)
ids[71:80]<-sample(which(fst$BINS==8),10)
ids[81:90]<-sample(which(fst$BINS==9),10)
ids[91:100]<-sample(which(fst$BINS==10),10)

Now ids should have ID's of 10 sampled from each bin.
This was written as a new file
write.table(fst[ids,],file="heliconius_10fromeachbin.txt",quote=FALSE,row.names=FALSE)

heliconius_10fromeachbin.txt now has the correct 100 loci with
stratified sampling.

cat heliconius_10fromeachbin.txt | awk '{print $1":"$3"-"$4}' >
heliconius_10fromeachbin_ranges.list
this file now contains all ranges

Since the two vcf files were just split in half, all I had to do was
cat the two into heliconius.vcf

Also downloaded the reference genome - H. melpomene from
ftp://ftp.ensemblgenomes.org/pub/metazoa/release-27/fasta/heliconius_melpomene/dna/
I downloaded the hard masked genome, and saved that as heliconius.fasta

Sigh - unfortunately not so easy...
java -jar /home/arun/stickleback/GenomeAnalysisTK.jar -T
SelectVariants --variant heliconius.vcf -R heliconius.fasta -o
heliconius_10each.vcf -se 'group*' -L
heliconius_10fromeachbin_ranges.list -sn tiP86 -sn tiP313 -sn tiP84
-sn tiP57 -sn am216 -sn am160 -sn am48 -sn am293

BLARG - at this point, I realized that the VCF file that was
downloaded from Dryad only had genotype calls, and not actual
location/genotype information. So had to give up on this and pick the
next dataset, which comes from Supple et al. (2013) -
http://dx.doi.org/10.1101/gr.150615.112

Data was downloaded - VCF (hopefully) from
http://datadryad.org/resource/doi:10.5061/dryad.rr65n/20

Also downloaded the reference FASTA from
ftp://ftp.ensemblgenomes.org/pub/metazoa/release-27/fasta/heliconius_melpomene/dna/Heliconius_melpomene.Hmel1.27.dna_rm.genome.fa.gz

This is hard-masked.

I am using two files - peru_aglaope*.vcf (rayed) and
peru_amaryllis*.vcf (postman) - contains data from Nadeau et al. 2012
basically

Had to do some gobbledygook here...
Had to add headers to both files first, then add a header note about
SB. Then remove the PL annotation using the below command:

/home/arun/Documents/stickleback/bcftools/bcftools annotate -x
FORMAT/PL,FORMAT/SB peru_aglaope_BD.vcf > pa.vcf


First step would be to merge these two vcfs.

I am going to do this using the vcf-merge script:

bgzip pam.vcf
bgzip pag.vcf
tabix -p vcf pam.vcf.gz
tabix -p vcf pag.vcf.gz
vcf-merge peru_am*.gz peru_ag*.gz > peruamag.vcf

Next step would be to compute Fsts:
vcftools --vcf merged.vcf.gz --weir-fst-pop Paniscus.txt
--weir-fst-pop Troglodytes.txt --out fst.txt --fst-window-size 1000000
--fst-window-step 200000

vcftools --vcf peruamag.vcf --weir-fst-pop Amaryllis.txt
--weir-fst-pop Aglaope.txt --out fst.txt --fst-window-size 15000
--fst-window-step 5000

This gave me only 121 windows - so might as well just keep them all
for further analyses.

cat fst.txt.windowed.weir.fst | awk '{print $1":"$3"-"$2}' > ranges.list

Since I'm using all the loci, I just have to go ahead and create BED
files as before...

 vcftools --vcf peruamag.vcf --out peruamagfiltered.vcf
--max-missing-count 0 --removeUnusedAlternates --recode

Had to index the fasta file first:


i=1
for vals in $(cat ranges.list);
do
echo $vals > temp.list
java -jar /home/arun/Documents/stickleback/GenomeAnalysisTK.jar -T
SelectVariants --variant peruamagfiltered.vcf.recode.vcf -ef -env -R
heliconius.fasta -o indivvcfs/locus"$i".vcf -L temp.list
i=$((i+1)); done
rm temp.list


Aglaope individuals first, then amaryllis

Now making BED files from all these VCF's (also note that only 118
loci were done - contigs were out of range for the rest of them)

for i in {1..118}
do
sed -e 's/locus//' locus"$i".vcf | awk '{OFS="\t"; if (!/^#/) {print
$1,$2-1,$2,$4,$5,$10,$11,$12,$13,$14,$15,$16,$17}}' > temp
sed -e 's/:[0-9]*//g' temp | sed -e 's/\.//g' | sed -e 's/|//g' | sed
-e 's/PASS//g' | sed -e 's/xqual//g' | sed -e 's/xlowcov//g' | sed -e
's/;//g' |  sed -e 's/,[0-9]*\t/\t/g' | sed -e 's/0,[0-9]*/0/g' | sed
-e 's/1,[0-9]*/1/g' | sed -e 's/2,[0-9]*/2/g' > locus"$i".bed
done

gen trying...

sed -e 's/locus//' locus1.vcf | awk '{OFS="\t"; if (!/^#/) {print
$1,$2-1,$2,$4,$5,$10,$11,$12,$13,$14,$15,$16,$17}}' > temp
sed -e 's/:[0-9]*//g' temp | sed -e 's/\.//g' | sed -e 's/|//g' | sed
-e 's/PASS//g' | sed -e 's/xqual//g' | sed -e 's/xlowcov//g' | sed -e
's/;//g' |  sed -e 's/,[0-9]*\t/\t/g' | sed -e 's/0,[0-9]*/0/g' | sed
-e 's/1,[0-9]*/1/g' | sed -e 's/2,[0-9]*/2/g' > blah


Then inside R - removing duplicate sites:

for (i in 71:118) {
locusname<-paste("locus",i,".bed",sep="")
heli<-read.table(locusname,header=TRUE,skip=1)
x=apply(heli[,6:13],1,unique)
y<-which(sapply(x,length)>1)
nodups<-heli[y,]
write.table(nodups,locusname,col.names=FALSE,row.names=FALSE)
}


In shell:

for i in {1..118}
do
sed -e 's/"//g' locus"$i".bed | sed -e 's/ 0 / 00 /g' | sed -e 's/ 1 /
01 /g' > temp
mv temp locus"$i".bed
done

for i in {1..118}
do
perl -p -e 's/ 0\n/ 00\n/' locus"$i".bed > temp
mv temp locus"$i".bed
done

for i in {1..118}
do
perl -p -e 's/ 1\n/ 01\n/' locus"$i".bed > temp
mv temp locus"$i".bed
done

for i in {1..118}
do
sed -e 's/ 0 / 00 /g' locus"$i".bed > temp
mv temp locus"$i".bed
done

for i in {1..118}
do
perl -p -e 's/ 1 / 01 /' locus"$i".bed > temp
mv temp locus"$i".bed
done

for i in {1..118}
do
sed -e 's/ 1 / 01 /g' locus"$i".bed > temp
mv temp locus"$i".bed
done

for i in {1..118}
do
perl -p -e 's/ 2 / 02 /' locus"$i".bed > temp
mv temp locus"$i".bed
done


sed -i "s/xhypercov//g" *.bed
sed -i "/,/d" *.bed
Now onto writing headers:

i=1
for vals in $(cat ../../ranges.list);
do
sed -i '1s/^/chromosome  begin   end ref variantbase
Melpomene_aglaope.09-246 Melpomene_aglaope.09-267
Melpomene_aglaope.09-268 Melpomene_aglaope.09-357
Melpomene_amaryllis.09-332 Melpomene_amaryllis.09-333
Melpomene_amaryllis.09-79 Melpomene_amaryllis.09-75\n/' locus"$i".bed
vars=$(wc -l locus"$i".bed | awk '{print $1}');
vars=$((vars-1));
sed -i "1s/^/$vals 1999 $vars\n/" locus"$i".bed
sed -i 's/:/ /g' locus"$i".bed
sed -i 's/group//g' locus"$i".bed
sed -i 's/-/ /g' locus"$i".bed
sed -i "s/\//g" locus"$i".bed
sed -i 's/\.09 /09_/g' locus"$i".bed
i=$((i+1));
done

sed -i 's/\.09 /09_/g' *.bed
sed -i "s/\//g" *.bed

Now that I have BED files, I can run my script on these...

Ran into issues running SITES on sites with > 2 alleles. So removed
all sites that had > 2 alleles.

for i in {1..118}
do
vcftools --vcf locus"$i".vcf --max-missing-count 0 --recode
--min-alleles 2 --max-alleles 2 --out temp
mv temp.recode.vcf locus"$i".vcf
done

Also had to shorten the windows:


i=1
for vals in $(cat ../ranges.list);
do
minwindow=$(echo $vals |  sed -e 's/:/ /g' | sed -e 's/-/ /g' | awk
'{print $2}');
maxwindow=$((minwindow+2000));
vcftools --vcf locus"$i".vcf --max-missing-count 0 --recode
--min-alleles 2 --max-alleles 2 --chr HE670865 --from-bp "$minwindow"
--to-bp "$maxwindow" --out temp
mv temp.recode.vcf 2kvcfs/locus"$i".vcf
i=$((i+1));
done

Again too big...so doing 2k window:

done!

now concatenating in order:

for i in {1..118}
do
cat locus"$i".u >> heliconius.u
done


geometric mean of theta = 4.288

So using this for priors:

theta = 21.44
m = 8.577
t = 0.466



--

CHIMP

Had to get a hold of the vcf files from
https://eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs for Pan
paniscus, and Pan troglodytes
There is no known migration between these two.

Then, had to first index the VCF files using tabix:

tabix -p vcf /home/arun/sharedCCGG/Pan_VCFs/Pan_paniscus.vcf.gz
tabix -p vcf /home/arun/sharedCCGG/Pan_VCFs/Pan_troglodytes.vcf.gz

Then the two indexed VCF's had to be merged using vcf-merge:

vcf-merge /home/arun/sharedCCGG/Pan_VCFs/Pan_paniscus.vcf.gz
/home/arun/sharedCCGG/Pan_VCFs/Pan_troglodytes.vcf.gz | bgzip -c >
/home/arun/sharedCCGG/Pan_VCFs/merged.vcf.gz

Now we can use the merged file for Fst calculations
vcftools --gzvcf merged.vcf.gz --weir-fst-pop Paniscus.txt
--weir-fst-pop Troglodytes.txt --out fst_trog.txt --fst-window-size
1000000 --fst-window-step 200000

Human genome assembly 37 was downloaded from UCSC:
ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
I downloaded the masked assembly.

Then now onto sampling 100 loci, like before:

fst<-read.table("*.fst",header=TRUE)
fst$BINS<-c(cut(fst$MEAN_FST,breaks=10))

ids<-array(dim=c(1,88))
ids[1:7]<-sample(which(fst$BINS==1),7)
ids[8:20]<-sample(which(fst$BINS==2),13)
ids[21:30]<-sample(which(fst$BINS==3),10)
ids[31:40]<-sample(which(fst$BINS==4),10)
ids[41:50]<-sample(which(fst$BINS==5),10)
ids[51:60]<-sample(which(fst$BINS==6),10)
ids[61:70]<-sample(which(fst$BINS==7),10)
ids[71:80]<-sample(which(fst$BINS==8),10)
ids[81:84]<-sample(which(fst$BINS==9),4)
ids[85:88]<-sample(which(fst$BINS==10),4)
write.table(fst[ids,],file="paniscustrog.txt",quote=FALSE,row.names=FALSE)

cat paniscustrog.txt | awk '{print $1":"$2"-"$3}' > paniscustrogranges.list


Bins 8 (3) ,9 (0) , and 10 (2) have very few loci - so need to sample them all.

ids[71:73]<-sample(which(fst$BINS==8),3)
ids[74:75]<-sample(which(fst$BINS==10),2)

Ended up with 88 loci in all.

Creating the reference FASTA:
java -jar /home/arun/Documents/stickleback/picard/dist/picard.jar
CreateSequenceDictionary R=hg38.fasta O=hg38.dict
samtools faidx hg38.fasta

Had to index the merged vcf file:

tabix -p vcf merged.vcf.gz

i=1
for vals in $(cat paniscustrogranges.list);
do
echo $vals > temp.list
java -jar /home/arun/Documents/stickleback/GenomeAnalysisTK.jar -T
SelectVariants --variant merged.vcf.gz -ef -env -R hg38.fasta -o
paniscustrogvcfs/locus"$i".vcf -L temp.list -sn
Pan_paniscus-9731_LB502 -sn Pan_paniscus-A914_Hortense -sn
Pan_paniscus-A915_Kosana -sn Pan_paniscus-A917_Dzeeta -sn
Pan_paniscus-A918_Hermien -sn Pan_paniscus-A919_Desmond -sn
Pan_paniscus-A922_Catherine -sn Pan_paniscus-A923_Kombote -sn
Pan_paniscus-A924_Chipita -sn Pan_paniscus-A925_Bono -sn
Pan_paniscus-A926_Natalie -sn Pan_paniscus-A927_Salonga -sn
Pan_paniscus-A928_Kumbuka -sn
Pan_troglodytes_troglodytes-A957_Vaillant -sn
Pan_troglodytes_troglodytes-A958_Doris -sn
Pan_troglodytes_troglodytes-A959_Julie -sn
Pan_troglodytes_troglodytes-A960_Clara
i=$((i+1)); done
rm temp.list


for i in {2..75}
do
vcftools --vcf locus"$i".vcf --max-missing-count 0 --recode --out temp
mv temp.recode.vcf locus"$i".vcf
done

for i in {1..75}
do
vcftools --vcf locus"$i".vcf --max-missing-count 0 --recode
--min-alleles 2 --max-alleles 2 --out temp
mv temp.recode.vcf locus"$i".vcf
done








for i in {1..75}
do
sed -e 's/locus//' locus"$i".vcf | awk '{OFS="\t"; if (!/^#/) {print
$1,$2-1,$2,$4,$5,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26}}'
> temp
sed -e 's/:[0-9]*//g' temp | sed -e 's/\.[0-9]*//g' | sed -e
's/,[0-9]*//g' > locus"$i".bed
done

gen trying:
sed -e 's/locus//' locus2.vcf | awk '{OFS="\t"; if (!/^#/) {print
$1,$2-1,$2,$4,$5,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45,$46,$47}}'
> temp

sed -e 's/:[0-9]*//g' temp | sed -e 's/\.[0-9]*//g' | sed -e
's/,[0-9]*//g' > locus2.bed







Then inside R - removing duplicate sites:

for (i in 12:75) {
locusname<-paste("locus",i,".bed",sep="")
pan<-read.table(locusname,header=TRUE,skip=1)
x=apply(pan[,6:22],1,unique)
y<-which(sapply(x,length)>1)
nodups<-pan[y,]
write.table(nodups,locusname,col.names=FALSE,row.names=FALSE)
}


In shell:

for i in {1..75}
do
sed -e 's/"//g' locus"$i".bed | sed -e 's/ 0 / 00 /g' | sed -e 's/ 1 /
01 /g' > temp
mv temp locus"$i".bed
done

for i in {1..75}
do
perl -p -e 's/ 0\n/ 00\n/' locus"$i".bed > temp
mv temp locus"$i".bed
done

for i in {1..75}
do
perl -p -e 's/ 1\n/ 01\n/' locus"$i".bed > temp
mv temp locus"$i".bed
done

for i in {1..75}
do
sed -e 's/ 0 / 00 /g' locus"$i".bed > temp
mv temp locus"$i".bed
done

for i in {1..75}
do
perl -p -e 's/ 1 / 01 /' locus"$i".bed > temp
mv temp locus"$i".bed
done

for i in {1..75}
do
sed -e 's/ 1 / 01 /g' locus"$i".bed > temp
mv temp locus"$i".bed
done

for i in {1..75}
do
perl -p -e 's/ 2 / 02 /' locus"$i".bed > temp
mv temp locus"$i".bed
done

WRiting headers:

i=1
for vals in $(cat ../ranges.list);
do
sed -i '1s/^/chromosome  begin   end ref variantbase
Pan_paniscus_9731_LB502 Pan_paniscus_A914_Hortense
Pan_paniscus_A915_Kosana Pan_paniscus_A917_Dzeeta
Pan_paniscus_A918_Hermien Pan_paniscus_A919_Desmond
Pan_paniscus_A922_Catherine Pan_paniscus_A923_Kombote
Pan_paniscus_A924_Chipita Pan_paniscus_A925_Bono
Pan_paniscus_A926_Natalie Pan_paniscus_A927_Salonga
Pan_paniscus_A928_Kumbuka Pan_troglodytes_ellioti_Akwaya_Jean
Pan_troglodytes_ellioti_Banyo Pan_troglodytes_ellioti_Basho
Pan_troglodytes_ellioti_Damian Pan_troglodytes_ellioti_Julie
Pan_troglodytes_ellioti_Kopongo Pan_troglodytes_ellioti_Koto
Pan_troglodytes_ellioti_Paquita Pan_troglodytes_ellioti_Taweh
Pan_troglodytes_ellioti_Tobi
Pan_troglodytes_schweinfurthii_100037_Vincent
Pan_troglodytes_schweinfurthii_100040_Andromeda
Pan_troglodytes_schweinfurthii_9729_Harriet
Pan_troglodytes_schweinfurthii_A910_Bwambale
Pan_troglodytes_schweinfurthii_A911_Kidongo
Pan_troglodytes_schweinfurthii_A912_Nakuu
Pan_troglodytes_troglodytes_A957_Vaillant
Pan_troglodytes_troglodytes_A958_Doris
Pan_troglodytes_troglodytes_A959_Julie
Pan_troglodytes_troglodytes_A960_Clara
Pan_troglodytes_verus_9668_Bosco Pan_troglodytes_verus_9730_Donald
Pan_troglodytes_verus_A956_Jimmie Pan_troglodytes_verus_Clint
Pan_troglodytes_verus_X00100_Koby\n/' locus"$i".bed
vars=0
vars=$(wc -l locus"$i".bed | awk '{print $1}');
vars=$((vars-1));
sed -i "1s/^/$vals 999999 $vars\n/" locus"$i".bed
sed -i 's/:/ /g' locus"$i".bed
sed -i 's/-/ /g' locus"$i".bed
i=$((i+1));
done

Note that some of the BED files have no sites - the filters removed them all.

--

Realized that I should only be using Troglodytes troglodytes!! :( So
going ahead and redoing these...Also had to shorten windows to say
10k:

i=1
for vals in $(cat ../paniscustrogranges.list);
do
chromosome=$(echo $vals |  sed -e 's/:/ /g' | sed -e 's/-/ /g' | awk
'{print $1}');
#minwindow=$(grep -A 1 "Clara" locus"$i".vcf | grep "chr" | awk '{print $2}');
#maxwindow=$((minwindow+10000));
vcftools --vcf locus"$i".vcf --max-missing-count 0 --recode
--min-alleles 2 --max-alleles 2 --out temp
mv temp.recode.vcf cleaned/locus"$i".vcf
#vcftools --vcf 10kvcfs/locus"$i".vcf --chr "$chromosome" --from-bp
"$minwindow" --to-bp "$maxwindow" --out temp
#mv temp.recode.vcf 10kvcfs/locus"$i".vcf
i=$((i+1));
done

i=1
for vals in $(cat ../paniscustrogranges.list);
do
chromosome=$(echo $vals |  sed -e 's/:/ /g' | sed -e 's/-/ /g' | awk
'{print $1}');
minwindow=$(grep -A 1 "Clara" cleaned/locus"$i".vcf | grep "chr" | awk
'{print $2}');
maxwindow=$((minwindow+10000));
vcftools --vcf cleaned/locus"$i".vcf --chr "$chromosome" --from-bp
"$minwindow" --to-bp "$maxwindow" --recode --out temp
mv temp.recode.vcf cleaned/shorts/locus"$i".vcf
i=$((i+1));
done

--

for i in {1..88}
do
sed -e 's/locus//' locus"$i".vcf | awk '{OFS="\t"; if (!/^#/) {print
$1,$2-1,$2,$4,$5,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26}}'
> temp
sed -e 's/:[0-9]*//g' temp | sed -e 's/\.[0-9]*//g' | sed -e
's/,[0-9]*//g' > locus"$i".bed
done




i=1
for vals in $(cat ../../../paniscustrogranges.list);
do
sed -i '1s/^/chromosome  begin   end ref variantbase
Pan_paniscus-9731_LB502 Pan_paniscus-A914_Hortense
Pan_paniscus-A915_Kosana Pan_paniscus-A917_Dzeeta
Pan_paniscus-A918_Hermien Pan_paniscus-A919_Desmond
Pan_paniscus-A922_Catherine Pan_paniscus-A923_Kombote
Pan_paniscus-A924_Chipita Pan_paniscus-A925_Bono
Pan_paniscus-A926_Natalie Pan_paniscus-A927_Salonga
Pan_paniscus-A928_Kumbuka Pan_troglodytes_troglodytes-A957_Vaillant
Pan_troglodytes_troglodytes-A958_Doris
Pan_troglodytes_troglodytes-A959_Julie
Pan_troglodytes_troglodytes-A960_Clara\n/' locus"$i".bed
vars=0
chromosome=$(echo $vals |  sed -e 's/:/ /g' | sed -e 's/-/ /g' | awk
'{print $1}');
minwindow=$(grep -A 1 "Clara" locus"$i".vcf | grep "chr" | awk '{print $2}');
maxwindow=$((minwindow+10000));
vars=$(wc -l locus"$i".bed | awk '{print $1}');
vars=$((vars-1));
sed -i "1s/^/$chromosome $minwindow $maxwindow 9999 $vars\n/" locus"$i".bed
sed -i 's/:/ /g' locus"$i".bed
sed -i 's/-/ /g' locus"$i".bed
i=$((i+1));
done

--

now can run our code!

-- done - see empties folder for loci that didn't cut it

for i in {1..87}
do
cat locus"$i".u >> paniscustrog.u
done

--
geometric mean = 9.63
q = 48.17
m = 19.26
t = 0.207

But...realized that the inheritance scalars should be according to
X/Y/autosome...gotta change that/...

Then ran IMa2p, and looks like there isn't enough data...so going to
go ahead and do some filtering steps BEFORE computing Fst.

This filtering keeps only the paniscus and troglodytes individuals,
with biallelic data, and with no missing data.

vcftools --gzvcf merged.vcf.gz --max-missing-count 0 --keep
PaniscusTroglodytes.txt --min-alleles 2 --max-alleles 2 --recode --out
Pantrog_onlybiallelic_nomissing.vcf


Fst:
vcftools --vcf Pantrog_onlybiallelic_nomissing.vcf.recode.vcf
--weir-fst-pop Paniscus.txt --weir-fst-pop Troglodytes.txt --out
fst_pantrog.txt --fst-window-size 10000 --fst-window-step 2000

fst_cut$BINS<-c(cut(fst_cut$MEAN_FST,breaks=10))

ids<-array(dim=c(1,200))
ids[1:20]<-sample(which(fst_cut$BINS==1),20)
ids[21:40]<-sample(which(fst_cut$BINS==2),20)
ids[41:60]<-sample(which(fst_cut$BINS==3),20)
ids[61:80]<-sample(which(fst_cut$BINS==4),20)
ids[81:100]<-sample(which(fst_cut$BINS==5),20)
ids[101:120]<-sample(which(fst_cut$BINS==6),20)
ids[121:140]<-sample(which(fst_cut$BINS==7),20)
ids[141:160]<-sample(which(fst_cut$BINS==8),20)
ids[161:180]<-sample(which(fst_cut$BINS==9),20)
ids[181:200]<-sample(which(fst_cut$BINS==10),20)
write.table(fst_cut[ids,],file="paniscustrog_cleaned.txt",quote=FALSE,row.names=FALSE)

cat paniscustrog_cleaned.txt | awk '{print $1":"$2"-"$3}' >
paniscustrogranges_cleaned.list

Got 200 loci!

gzip Pantrog_onlybiallelic_nomissing.vcf.recode.vcf
tabix -p Pantrog_onlybiallelic_nomissing.vcf.recode.vcf.gz


i=1
for vals in $(cat paniscustrogranges.list);
do
echo $vals > temp.list
java -jar /home/arun/Documents/stickleback/GenomeAnalysisTK.jar -T
SelectVariants --variant
Pantrog_onlybiallelic_nomissing.vcf.recode.vcf.gz -ef -env -R
hg38.fasta -o 10kvcfs/locus"$i".vcf -L temp.list
i=$((i+1)); done
rm temp.list
