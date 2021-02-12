# Useful info on Variant Calling: https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionVI/lessons/01_alignment.html
# http://people.duke.edu/~ccc14/duke-hts-2017/Statistics/08032017/GATK-pipeline-sample.html

source /home/jonbra/miniconda3/etc/profile.d/conda.sh
conda activate GATK
# Usage:
# bash gatk.sh 472Virus200604-1 

# Set variables
scriptdir="/home/jonbra/FHI/Prosjekter/HCV/"
humanrefs="/home/jonbra/FHI/Prosjekter/HCV/Data/Homo_sapiens_transcriptome/" # Folder need to contain reference fastas and indexes from bwa index
humanPrefix="Hsa_trans" # Prefix used in bwa mem
viralrefs="/home/jonbra/FHI/FHI-databaser_og_scripts/Referanser_HCV_ICTV_190508_clean/"
viralbase="HCVgenosubtypes_8.5.19_clean"
myrefdir="/home/jonbra/FHI/Prosjekter/HCV/Data/"$1
samplelabel=$1
minAgensRead=50000 # Uklart for meg hvordan dette tallet skal bestemmes. Og hvorfor man trenger det. Er det ikke en cutoff på coverage som brukes?
# NB! Reads need to be trimmed with trim-galore
R1=$(ls *val_1.fq)
R2=$(ls *val_2.fq)

# Info to Read Group in the bam file
myrgID=472Virus200604-1-LANE001 ### This needs to be a UNIQUE label!!
myrgPL=illumina ## Platform/technology
myrgPU=NA # Platform Unit
myrgSM=$1 # Sample
myrgLB=$1 # Lable

######## DEL 2 Mapping #### START ######

# Map to human transcriptome and keep unmapped reads
echo "Mapping to human transcriptome"
bwa mem -t 8 $humanrefs$humanPrefix $R1 $R2 > Hsa_transcriptome_map.sam 2> Hsa_map-bwasam.err 
samtools flagstat Hsa_transcriptome_map.sam | grep "itself" | cut -f 1 -d " "
echo "reads (not pairs) mapped"
echo "Done"

# Extract unmapped reads
echo "Extracting unmapped reads"
samtools view -bS -f 4 Hsa_transcriptome_map.sam > unmapped.bam 

#Separate into R1 and R2
samtools sort -n unmapped.bam -o unmapped_sorted.bam
bedtools bamtofastq -i unmapped_sorted.bam -fq newR1.fq -fq2 newR2.fq 
echo "Done"

# Remove unnecessary files
rm unmapped*.bam

# Map to the entire HCV reference database using BWA
echo "Mapping unmapped reads to entire HCV reference database"
bwa mem -M -t 8 ${viralrefs}${viralbase}.fa newR1.fq newR2.fq > mapped_to_HCV_database.sam 2> Map_to_database-bwasam.err # -M to maintain compatibility with Picard/GATK
echo "Done"
echo "Extracting major subtype"

# Get "best hit" (major subtype) and map to it 
newR4=mapped_to_HCV_database.sam
samtools view -bS ${newR4} | samtools sort -o ${newR4%.sam}_sorted.bam
rm ${newR4}
samtools index ${newR4%.sam}_sorted.bam

# Activate another environment to run weeSAM (due to pysam conflicts)
conda activate HCV
weeSAM --bam ${newR4%.sam}_sorted.bam --out ${newR4%.sam}_stats.txt 

# Use samtools idxstats? Need an indexed bam file. Columns out: reference_name reference_length #mapped_reads #unmapped_reads
# samtools idxstats mapped_to_HCV_database_sorted.bam | sort -k3 -n -r 
# Then use bedtools genomecov to calculate the percent coverage?

Rscript --vanilla ${scriptdir}Rscript_sumreads.R "${newR4%.sam}_stats.txt" "${newR4%.sam}_sumstats.txt" # Beregner også prosent av totalt antall agens read

# Re-activate GATK
conda deactivate

sort -t$'\t' -k3 -nr ${newR4%.sam}_stats.txt > ${newR4%.sam}_stats_sorted.txt #Ikke nødvendig, men gjør det lettere å gå tilbake å se på resultatene fra første mapping

#align vs. best hit
major=$(sed -n 2p mapped_to_HCV_database_sumstats.txt | cut -d " " -f1 | cut -d'"' -f2)
bestF1=$(sort -t$'\t' -k3 -nr ${newR4%.sam}_stats.txt | grep ${major} -m1 | cut -f1) #Finne første referanse i _stats.txt som inneholder "major" og bruke denne som referanse for mapping
bestF2="${samplelabel}_${bestF1%_*}" # brukes til navnsetting av outputfil
echo "Major subtype is ${major}, with ${bestF1} as the reference sequence"

echo "Mapping against major subtype with BWA"
# Husk å legge inn read group etc i BWA-mappingen
# Create bwa index of best reference
bwa index ${viralrefs}${bestF1}.fa
# Should I add option -aM to bwa mem? Check out
bwa mem -M -t 8 ${viralrefs}${bestF1}.fa newR1.fq newR2.fq > ${bestF2}_bwa_vbest.sam 2> ${bestF2}_bwa_vbest-bwasam.err
bestF3=${bestF2}_bwa_vbest.sam
#samtools view -bS ${bestF3} | samtools sort -o ${bestF3%.sam}_sorted.bam
#samtools index ${bestF3%.sam}_sorted.bam

# Add read groups to bam file, sort bam file by coordinate and index bam file
echo "Adding read group info to bam file, sorting and indexing"
gatk --java-options "-Xmx8g" AddOrReplaceReadGroups I=${bestF2}_bwa_vbest.sam O=${bestF3%.sam}_rg_sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true RGID=$myrgID RGSM=$myrgSM RGLB=$myrgLB RGPL=$myrgPL RGPU=$myrgPU > ${bestF2}_bwa_vbest-rgsort.out 2> ${bestF2}_bwa_vbest-rgsort.err
echo "Mapping done"

######## DEL 2 Mapping #### SLUTT ######

######## DEL 2b Mapping mot minority #### START ######

#align vs. next best genotype
sumAgensRead=$(awk 'FNR > 1 {print $2}' mapped_to_HCV_database_sumstats.txt | paste -sd+ | bc) # Summarise all mapped reads
minor=$(sed -n 3p mapped_to_HCV_database_sumstats.txt | cut -d " " -f1 | cut -d'"' -f2) # Find subtype with second most reads mapped
bestMinor=$(sort -t$'\t' -k3 -nr ${newR4%.sam}_stats.txt | grep ${minor}_ -m1 | cut -f1)     #Finne første referanse i _stats.txt som inneholder "minor" og bruke denne som referanse for mapping
bestMinor_percCov=$(sort -t$'\t' -k3 -nr ${newR4%.sam}_stats.txt | grep ${bestMinor} -m1 | cut -f5)     #Finner hvor godt dekt referansen var i første mapping
bestMinor_percCov2=${bestMinor_percCov/.*}          #Fjerner desimaler for at "if"-setningen skal gjenkjenne tallet
bestMinor2="${samplelabel}_${bestMinor%_*}" # Fjerner "_" og alt etter.


	if [ ${sumAgensRead} -gt ${minAgensRead} ] && [ ${bestMinor_percCov2} -gt 5 ]; then
    # Create bwa index of best reference
    bwa index ${viralrefs}${bestMinor}.fa
    bwa mem -M -t 8 ${viralrefs}${bestMinor}.fa newR1.fq newR2.fq > ${bestMinor2}_bwa_bestMinor.sam 2> ${bestMinor}_bwa_bestMinor-bwasam.err
    bestMinor3=$(ls *_bwa_bestMinor.sam)
    samtools view -bS ${bestMinor3} | samtools sort -o ${bestMinor3%.sam}_sorted.bam
    samtools index ${bestMinor3%.sam}_sorted.bam

    else
    echo "Møter ikke kriteriene for mapping mot minority"

	fi

echo "Mapping against minority done!"

######## DEL 2b Mapping mot minority #### SLUTT######


######## DEL 3 VariantCalling og Consensus #### START ######

# Lage konsensus for Main-genotype

# Mark duplicates in the bam file
# [] Inspect the duplication metrics. Plot histogram?
echo "Removing duplicates"
gatk --java-options "-Xmx8g" MarkDuplicates I=${bestF3%.sam}_rg_sorted.bam O=${bestF3%.sam}_rg_sorted.dedup.bam M=${bestF2}-duplication-metric.txt ASSUME_SORT_ORDER=coordinate CREATE_INDEX=true REMOVE_DUPLICATES=true > ${bestF2}-dedup.out 2> ${bestF2}-dedup.err
# Notice that duplicates are removed, not only marked

# Extract info on duplicates
grep -A 2 "## METRICS CLASS" ${bestF2}-duplication-metric.txt | cut -f 3 > tmp.txt
grep -A 2 "## METRICS CLASS" ${bestF2}-duplication-metric.txt | cut -f 7 >> tmp.txt
grep -A 2 "## METRICS CLASS" ${bestF2}-duplication-metric.txt | cut -f 9 >> tmp.txt
grep -A 2 "## METRICS CLASS" ${bestF2}-duplication-metric.txt | cut -f 10 >> tmp.txt

# Show duplicate info
cat tmp.txt

# Remove duplicate info
rm tmp.txt

readsmapped=$(samtools flagstat ${bestF3%.sam}_rg_sorted.dedup.bam | grep "itself" | cut -f 1 -d " ")
echo "${readsmapped} reads (not pairs) mapped after deduplication"

# Don't think I should do Base Score recalibration if I don't have a known sites file

# Run the HaplotypeCaller
echo "Running HaplotypeCaller"
gatk CreateSequenceDictionary R=${viralrefs}${bestF1}.fa O=${viralrefs}${bestF1}.dict # create reference dictionary
gatk --java-options "-Xmx8g" HaplotypeCaller -R ${viralrefs}${bestF1}.fa -I ${bestF3%.sam}_rg_sorted.dedup.bam -O ${bestF2}"-GATK-HC.vcf" -ploidy 1
echo "Done"

# Create consensus sequences with coverage less than 2 and 6 masked
echo "Creating consensus sequences"

# Get regions with low coverage
bedtools genomecov -bga -ibam ${bestF3%.sam}_rg_sorted.dedup.bam | grep -w '0$\|1$' > ${bestF2}_BWA-HC_regionswithlessthan2coverage.bed
bedtools genomecov -bga -ibam ${bestF3%.sam}_rg_sorted.dedup.bam | grep -w '0$\|1$\|2$\|3$\|4$\|5$' > ${bestF2}_BWA-HC_regionswithlessthan6coverage.bed

# Make consensus sequence from HC calls
gatk FastaAlternateReferenceMaker -R ${viralrefs}${bestF1}.fa -O tmp.fa -V ${bestF2}"-GATK-HC.vcf" 
sed -i "1s/.*/>${bestF1}/" tmp.fa # Replace fasta header with reference name

# Mask regions of low coverage with N's
bedtools maskfasta -fi tmp.fa -bed ${bestF2}_BWA-HC_regionswithlessthan2coverage.bed -fo ${bestF2}_BWA-HC_lessthan2masked_consensus.fa
sed -i "1s/.*/>less2mask_${bestF1}/" ${bestF2}_BWA-HC_lessthan2masked_consensus.fa # Change fasta header to reflect masking grade
bedtools maskfasta -fi tmp.fa -bed ${bestF2}_BWA-HC_regionswithlessthan6coverage.bed -fo ${bestF2}_BWA-HC_lessthan6masked_consensus.fa
sed -i "1s/.*/>less6mask_${bestF1}/" ${bestF2}_BWA-HC_lessthan6masked_consensus.fa # Change fasta header to reflect masking grade

rm tmp.*

echo "Consensus sequences created for bases with coverage less than 2 and 6"
conda deactivate
