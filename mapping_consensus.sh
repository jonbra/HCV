# Skript til lokalt bruk på Jons UiO Linux laptop
# Aktivere conda environment
conda activate HCV

basedir=$(pwd)
#runname=${basedir##*/}
runname=local-test

#husk å legge inn Rscript_sumreads.R
scriptdir=/home/jonbra/FHI/Prosjekter/HCV/
#tanotidir=/home/ngs2/Downloads/Tanoti-1.2-Linux/
#weesamdir=/home/ngs2/.fhiscripts/weeSAM/
script_name1=`basename $0`
#skille software fra rapportering
#VirusScriptDir=/home/ngs2/.fhiscripts/VirusScriptParts/


###### DATABASER/REFERANSESEKVENSER ########
HCV_RefDir=/home/jonbra/FHI/FHI-databaser_og_scripts/Referanser_HCV_ICTV_190508_clean
#HEV_RefDir=/media/data/Referanser_HEV
#Corona_RefDir=/media/data/Referanser_Corona
#Dengue_RefDir=/media/data/Referanser_Dengue
#Entero_RefDir=/media/data/Referanser_Entero
#TBEV_RefDir=/media/data/Referanser_TBEV

########## FYLL INN FOR AGENS ###################
# kan også legge til trimming-setinger her om man ønsker muligheten for at det skal være ulikt (phred-score og minimum lengde på read) 

#Skriv inn agens-navn (må være skrevet likt som i navnet på fasta-fil som inneholder databasen/referansesekvensene
Agens=HCV					#ingen mellomrom etter =

#husk å legge inn rett variabel for filbanen til databasen/referansesekvensene (se under "DATABASER/REFERANSESEKVENSER")
Refdir=${HCV_RefDir}	# f.eks. ${HCV_RefDir}

#presisere stringency for mapping, 1-100
String=85 				#Stringens i første mapping     ingen mellomrom etter =
String2=95                #Stringens i andre mapping (hoved og minor)

#Definere hvor mange read det må være mappet mot agens før det gjøres mapping mot minoritetsvariant, f.eks. 50000
minAgensRead=50000			#ingen mellomrom etter =





######## DEL 1 Trimming #### START ######

basedir=$(pwd)
runname=${basedir##*/} # ##*/ Fjerne alt foran basedir

R1=$(ls *_R1*.fastq.gz)
R2=$(ls *_R2*.fastq.gz)
#trim
trim_galore -q 30 --dont_gzip --length 50 --paired ${R1} ${R2}
    
echo "#"
echo "Read ferdig trimmet" 
echo "#"
echo "###################"


######## DEL 1 Trimming #### SLUTT ######


######## DEL 2 Mapping #### START ######
basedir=$(pwd)
runname=${basedir##*/}

R1=$(ls *_R1*.fastq.gz)
newR1=$(ls *val_1.fq)
newR2=$(ls *val_2.fq)

#align vs. entire db
tanoti -r ${Refdir}/${Agens}*.fa -i ${newR1} ${newR2} -o ${R1%%_*L001*}_tanoti.sam -p 1 -u 1 -m ${String} #dobbel % fjerner lengste mulige substring, enkelt % fjerner korteste mulige substring i ${variable%substring}
newR4=$(ls *_tanoti.sam)
samtools view -bS ${newR4} | samtools sort -o ${newR4%.sam}_sorted.bam
samtools index ${newR4%.sam}_sorted.bam
weeSAM --bam ${newR4%.sam}_sorted.bam --out ${newR4%.sam}_stats.txt 
Rscript --vanilla ${scriptdir}Rscript_sumreads.R "${newR4%.sam}_stats.txt" "${newR4%.sam}_sumstats.txt" # Beregner også prosent av totalt antall agens read
	
sort -t$'\t' -k3 -nr ${newR4%.sam}_stats.txt > ${newR4%.sam}_stats_sorted.txt #Ikke nødvendig, men gjør det lettere å gå tilbake å se på resultatene fra første mapping	
	
#align vs. best hit
major=$(sed -n 2p  *_tanoti_sumstats.txt | cut -d " " -f1 | cut -d'"' -f2)  
bestF1=$(sort -t$'\t' -k3 -nr ${newR4%.sam}_stats.txt | grep ${major} -m1 | cut -f1) #Finne første referanse i _stats.txt som inneholder "major" og bruke denne som referanse for mapping   
bestF2="${R1%%_*L001*}_${bestF1%_*}" # brukes til navnsetting av outputfil 
tanoti -r ${Refdir}/${bestF1}.fa -i ${newR1} ${newR2} -o ${bestF2}_tanoti_vbest.sam -p 1 -m ${String2}
bestF3=$(ls *_tanoti_vbest.sam)
samtools view -bS ${bestF3} | samtools sort -o ${bestF3%.sam}_sorted.bam
samtools index ${bestF3%.sam}_sorted.bam
   

cd "${basedir}"

echo "HEY HEY HEY, What's that sound?" 
echo "Mapping done!"



######## DEL 2 Mapping #### SLUTT ######


######## DEL 2b Mapping mot minority #### START ######

basedir=$(pwd)
runname=${basedir##*/}

R1=$(ls *_R1*.fastq.gz)
newR1=$(ls *val_1.fq)
newR2=$(ls *val_2.fq)


	
 #align vs. next best genotype 
newR4=$(ls *_tanoti.sam)
sumAgensRead=$(awk 'FNR > 1 {print $2}' *sumstats.txt| paste -sd+ | bc)
minor=$(sed -n 3p  *_tanoti_sumstats.txt | cut -d " " -f1 | cut -d'"' -f2)
bestMinor=$(sort -t$'\t' -k3 -nr ${newR4%.sam}_stats.txt | grep ${minor}_ -m1 | cut -f1)     #Finne første referanse i _stats.txt som inneholder "minor" og bruke denne som referanse for mapping 
bestMinor_percCov=$(sort -t$'\t' -k3 -nr ${newR4%.sam}_stats.txt | grep ${bestMinor} -m1 | cut -f5)     #Finner hvor godt dekt referansen var i første mapping
bestMinor_percCov2=${bestMinor_percCov/.*}          #Fjerner desimaler for at "if"-setningen skal gjenkjenne tallet
bestMinor2="${R1%%_*L001*}_${bestMinor%_*}"
    
	
	if [ ${sumAgensRead} -gt ${minAgensRead} ] && [ ${bestMinor_percCov2} -gt 5 ]; then      
   
    tanoti -r ${Refdir}/${bestMinor}.fa -i ${newR1} ${newR2} -o ${bestMinor2}_tanoti_bestMinor.sam -p 1 -m ${String2}
    bestMinor3=$(ls *_tanoti_bestMinor.sam)
    samtools view -bS ${bestMinor3} | samtools sort -o ${bestMinor3%.sam}_sorted.bam
    samtools index ${bestMinor3%.sam}_sorted.bam
    
    else
    echo "Møter ikke kriteriene for mapping mot minority"
    
	fi


    cd "${basedir}"


echo "HEY HEY HEY, What's that sound?" 
echo "Mapping against minority done!"

######## DEL 2b Mapping mot minority #### SLUTT######



######## DEL 3 VariantCalling og Consensus #### START ######

basedir=$(pwd)
runname=${basedir##*/}

# Lage konsensus for Main-genotype
	newR4=$(ls *_tanoti.sam) 
	major=$(sed -n 2p  *_tanoti_sumstats.txt | cut -d " " -f1 | cut -d'"' -f2)  
	bestF1=$(sort -t$'\t' -k3 -nr ${newR4%.sam}_stats.txt | grep ${major} -m1 | cut -f1)
	bestF3=$(ls *_tanoti_vbest.sam)

	samtools sort -n ${bestF3%.sam}_sorted.bam > ${bestF3%.sam}_sorted.byQuery.bam 
	samtools fixmate -m ${bestF3%.sam}_sorted.byQuery.bam ${bestF3%.sam}_sorted.fix.bam
	samtools sort ${bestF3%.sam}_sorted.fix.bam > ${bestF3%.sam}_sorted.fix_sorted.bam

	samtools markdup -r ${bestF3%.sam}_sorted.fix_sorted.bam ${bestF3%.sam}_sorted.marked.bam
	bcftools mpileup -Ou -f ${Refdir}/${bestF1}.fa ${bestF3%.sam}_sorted.marked.bam| bcftools call -mv -Ob -o calls.bcf.gz
	bcftools index calls.bcf.gz

#	bedtools genomecov -bga -ibam ${bestF3%.sam}_sorted.marked.bam| grep -w '0$' > regionswith0coverage.bed   # '0$\|1$\|2$\|3$\|4$\|5$' > regionswithlessthan6coverage
#	bcftools consensus -m regionswith0coverage.bed -f ${Refdir}${bestF1}.fa calls.vcf.gz -o cons.fa

    samtools index ${bestF3%.sam}_sorted.marked.bam

    bedtools genomecov -bga -ibam ${bestF3%.sam}_sorted.marked.bam| grep -w '0$\|1$\|2$\|3$\|4$\|5$' > regionswithlessthan6coverage.bed   
	bcftools consensus -m regionswithlessthan6coverage.bed -f ${Refdir}/${bestF1}.fa calls.bcf.gz -o cons.fa

	seqkit replace -p "(.+)" -r ${bestF3%%_*} cons.fa > ${bestF3%%_*}_consensus.fa #endrer navn fra referanse-navn til prøvenavn inne i fasta-fil
	
#sletter filer som ikke trengs videre: 
	rm *cons.fa 
	rm *calls*.bcf.gz
	rm *calls*.bcf.gz.csi 
	rm *regionswith*coverage.bed 
	rm *_sorted.byQuery.bam 
	rm *_sorted.fix.bam
	rm *_sorted.fix_sorted.bam



# Lage konsensus for minoritet-genotype

	sumAgensRead=$(awk 'FNR > 1 {print $2}' *sumstats.txt| paste -sd+ | bc)
    newR4=$(ls *_tanoti.sam)    
    minor=$(sed -n 3p  *_tanoti_sumstats.txt | cut -d " " -f1 | cut -d'"' -f2)    
    bestMinor=$(sort -t$'\t' -k3 -nr ${newR4%.sam}_stats.txt | grep ${minor}_ -m1 | cut -f1)
    bestMinor_percCov=$(sort -t$'\t' -k3 -nr ${newR4%.sam}_stats.txt | grep ${bestMinor} -m1 | cut -f5)    
    bestMinor_percCov2=${bestMinor_percCov/.*}  

	
	if [ ${sumAgensRead} -gt ${minAgensRead} ] && [ ${bestMinor_percCov2} -gt 5 ]; then 
		bestMinor=$(sort -t$'\t' -k3 -nr ${newR4%.sam}_stats.txt | grep ${minor}_ -m1 | cut -f1)
		bestMinor3=$(ls *_tanoti_bestMinor.sam)

		samtools sort -n ${bestMinor3%.sam}_sorted.bam > ${bestMinor3%.sam}_sorted.byQuery.bam 
		samtools fixmate -m ${bestMinor3%.sam}_sorted.byQuery.bam ${bestMinor3%.sam}_sorted.fix.bam
		samtools sort ${bestMinor3%.sam}_sorted.fix.bam > ${bestMinor3%.sam}_sorted.fix_sorted.bam

		samtools markdup -r ${bestMinor3%.sam}_sorted.fix_sorted.bam ${bestMinor3%.sam}_sorted.marked.bam
		bcftools mpileup -Ou -f ${Refdir}/${bestMinor}.fa ${bestMinor3%.sam}_sorted.marked.bam| bcftools call -mv -Ob -o calls.bcf.gz
		bcftools index calls.bcf.gz

		#bedtools genomecov -bga -ibam ${bestMinor3%.sam}_sorted.marked.bam| grep -w '0$' > regionswith0coverage.bed   # '0$\|1$\|2$\|3$\|4$\|5$' > regionswithlessthan6coverage
		#bcftools consensus -m regionswith0coverage.bed -f ${Refdir}${bestMinor}.fa calls.vcf.gz -o cons.fa
        
        samtools index ${bestMinor3%.sam}_sorted.marked.bam

        bedtools genomecov -bga -ibam ${bestMinor3%.sam}_sorted.marked.bam| grep -w '0$\|1$\|2$\|3$\|4$\|5$' > regionswithlessthan6coverage.bed   
		bcftools consensus -m regionswithlessthan6coverage.bed -f ${Refdir}/${bestMinor}.fa calls.bcf.gz -o cons.fa
   

		seqkit replace -p "(.+)" -r ${bestMinor3%%_*}_Minor cons.fa > ${bestMinor3%%_*}_Minor_consensus.fa #endrer navn fra referanse-navn til prøvenavn inne i fasta-fil
				
		
		#sletter filer som ikke trengs videre: 
		rm *cons.fa 
		rm *calls*.bcf.gz
		rm *calls*.bcf.gz.csi 
		rm *regionswith*coverage.bed 
		rm *_sorted.byQuery.bam 
		rm *_sorted.fix.bam
		rm *_sorted.fix_sorted.bam
	fi




cd "${basedir}"


echo "Consensus made" 

######## DEL 3 VariantCalling og Consensus #### SLUTT ######
