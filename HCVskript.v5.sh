### NY VERSJON v2 ######## I v2 er ICTV-databasen oppdatert og gaps og N'er er fjernet fra referansene.
						## I v2 er det lagt til opprettelse av consensus-sekvens

### Versjon 3 ble hoppet over ###

### NY VERSJON v4 ######## Duplikater fjernes før consensus lages


### NY VERSJON v5 ######## I v5 er det gjort store endringer på rekkefølgen av de ulike delene 
                        ## Stringens økt fra 85 til 95 i andre mapping
                        ## Lagt til coverage-plot og statistikk etter duplikater - info til pdf er nå uten duplikater
                        ## Consensus lages kun ved minimum 6 i dybde. Settes inn N'er ved 5 eller mindre. (etter duplikater er fjernet)
                        ## Endret å få ut % covered above 9 til above 5 og alle % covered er nå uten duplikater
                        ## Endret en del av parameterne som hentes ut til csv og pdf
                        ## GLUE-rapport lages fra bam-fil uten duplikater
                        ## For å kjøre igjennom minor-loop må den aktuelle referansen ha vært dekt med minst 5 % i føste runde med mapping

## Skript startes fra run-mappen (f.eks. Run443)

basedir=$(pwd)
runname=${basedir##*/}

#husk å legge inn Rscript_sumreads.R
scriptdir=/home/ngs2/.fhiscripts/
#tanotidir=/home/ngs2/Downloads/Tanoti-1.2-Linux/
#weesamdir=/home/ngs2/.fhiscripts/weeSAM/
script_name1=`basename $0`
#skille software fra rapportering
VirusScriptDir=/home/ngs2/.fhiscripts/VirusScriptParts/


###### DATABASER/REFERANSESEKVENSER ########
HCV_RefDir=/media/data/Referanser_HCV_ICTV_190508_clean
HEV_RefDir=/media/data/Referanser_HEV
Corona_RefDir=/media/data/Referanser_Corona
Dengue_RefDir=/media/data/Referanser_Dengue
Entero_RefDir=/media/data/Referanser_Entero
TBEV_RefDir=/media/data/Referanser_TBEV

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
runname=${basedir##*/}

for dir in $(ls -d Virus*/)
do
    cd ${dir}
    R1=$(ls *_R1*.fastq.gz)
    R2=$(ls *_R2*.fastq.gz)
#trim
    trim_galore -q 30 --dont_gzip --length 50 --paired ${R1} ${R2}
    
    cd "${basedir}"
done

echo "#"
echo "Read ferdig trimmet" 
echo "#"
echo "###################"


######## DEL 1 Trimming #### SLUTT ######


######## DEL 2 Mapping #### START ######
basedir=$(pwd)
runname=${basedir##*/}

for dir in $(ls -d Virus*/)
do
    cd ${dir}
	R1=$(ls *_R1*.fastq.gz)
    newR1=$(ls *val_1.fq)
    newR2=$(ls *val_2.fq)

#align vs. entire db
    tanoti -r ${Refdir}/${Agens}*.fa -i ${newR1} ${newR2} -o ${R1%%_*L001*}_tanoti.sam -p 1 -u 1 -m ${String} #dobbel % fjerner lengste mulige substring, enkelt % fjerner korteste mulige substring i ${variable%substring}
    newR4=$(ls *_tanoti.sam)
    samtools view -bS ${newR4} | samtools sort -o ${newR4%.sam}_sorted.bam
    samtools index ${newR4%.sam}_sorted.bam
    weeSAMv1.4 -b ${newR4%.sam}_sorted.bam -out ${newR4%.sam}_stats.txt 
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
done


echo "HEY HEY HEY, What's that sound?" 
echo "Mapping done!"



######## DEL 2 Mapping #### SLUTT ######


######## DEL 2b Mapping mot minority #### START ######

basedir=$(pwd)
runname=${basedir##*/}

for dir in $(ls -d Virus*/)
do
    cd ${dir}
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

done



echo "HEY HEY HEY, What's that sound?" 
echo "Mapping against minority done!"

######## DEL 2b Mapping mot minority #### SLUTT######



######## DEL 3 VariantCalling og Consensus #### START ######

basedir=$(pwd)
runname=${basedir##*/}

for dir in $(ls -d Virus*/)
do

cd ${dir}

# Lage konsensus for Main-genotype
	newR4=$(ls *_tanoti.sam) 
	major=$(sed -n 2p  *_tanoti_sumstats.txt | cut -d " " -f1 | cut -d'"' -f2)  
	bestF1=$(sort -t$'\t' -k3 -nr ${newR4%.sam}_stats.txt | grep ${major} -m1 | cut -f1)
	bestF3=$(ls *_tanoti_vbest.sam)

	samtools sort -n ${bestF3%.sam}_sorted.bam > ${bestF3%.sam}_sorted.byQuery.bam 
	samtools fixmate -m ${bestF3%.sam}_sorted.byQuery.bam ${bestF3%.sam}_sorted.fix.bam
	samtools sort ${bestF3%.sam}_sorted.fix.bam > ${bestF3%.sam}_sorted.fix_sorted.bam

	samtools markdup -r ${bestF3%.sam}_sorted.fix_sorted.bam ${bestF3%.sam}_sorted.marked.bam
	bcftools mpileup -f ${Refdir}/${bestF1}.fa ${bestF3%.sam}_sorted.marked.bam| bcftools call -mv -Ob -o calls.vcf.gz
	bcftools index calls.vcf.gz

#	bedtools genomecov -bga -ibam ${bestF3%.sam}_sorted.marked.bam| grep -w '0$' > regionswith0coverage.bed   # '0$\|1$\|2$\|3$\|4$\|5$' > regionswithlessthan6coverage
#	bcftools consensus -m regionswith0coverage.bed -f ${Refdir}${bestF1}.fa calls.vcf.gz -o cons.fa

    samtools index ${bestF3%.sam}_sorted.marked.bam

    bedtools genomecov -bga -ibam ${bestF3%.sam}_sorted.marked.bam| grep -w '0$\|1$\|2$\|3$\|4$\|5$' > regionswithlessthan6coverage.bed   
	bcftools consensus -m regionswithlessthan6coverage.bed -f ${Refdir}/${bestF1}.fa calls.vcf.gz -o cons.fa

	seqkit replace -p "(.+)" -r ${bestF3%%_*} cons.fa > ${bestF3%%_*}_consensus.fa #endrer navn fra referanse-navn til prøvenavn inne i fasta-fil
	
#sletter filer som ikke trengs videre: 
	rm *cons.fa 
	rm *calls*.vcf.gz
	rm *calls*.vcf.gz.csi 
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
		bcftools mpileup -f ${Refdir}/${bestMinor}.fa ${bestMinor3%.sam}_sorted.marked.bam| bcftools call -mv -Ob -o calls.vcf.gz
		bcftools index calls.vcf.gz

		#bedtools genomecov -bga -ibam ${bestMinor3%.sam}_sorted.marked.bam| grep -w '0$' > regionswith0coverage.bed   # '0$\|1$\|2$\|3$\|4$\|5$' > regionswithlessthan6coverage
		#bcftools consensus -m regionswith0coverage.bed -f ${Refdir}${bestMinor}.fa calls.vcf.gz -o cons.fa
        
        samtools index ${bestMinor3%.sam}_sorted.marked.bam

        bedtools genomecov -bga -ibam ${bestMinor3%.sam}_sorted.marked.bam| grep -w '0$\|1$\|2$\|3$\|4$\|5$' > regionswithlessthan6coverage.bed   
		bcftools consensus -m regionswithlessthan6coverage.bed -f ${Refdir}/${bestMinor}.fa calls.vcf.gz -o cons.fa
   

		seqkit replace -p "(.+)" -r ${bestMinor3%%_*}_Minor cons.fa > ${bestMinor3%%_*}_Minor_consensus.fa #endrer navn fra referanse-navn til prøvenavn inne i fasta-fil
				
		
		#sletter filer som ikke trengs videre: 
		rm *cons.fa 
		rm *calls*.vcf.gz
		rm *calls*.vcf.gz.csi 
		rm *regionswith*coverage.bed 
		rm *_sorted.byQuery.bam 
		rm *_sorted.fix.bam
		rm *_sorted.fix_sorted.bam
	fi




cd "${basedir}"


done

echo "Consensus made" 

######## DEL 3 VariantCalling og Consensus #### SLUTT ######


######## DEL 4 CoveragePlot og Statistikk #### START ######
basedir=$(pwd)
runname=${basedir##*/}

for dir in $(ls -d Virus*/)
do
    cd ${dir}
	
# Coverage plot og statistikkmed duplikater
	bestF3=$(ls *_tanoti_vbest.sam)
	weeSAMv1.4 -b ${bestF3%.sam}_sorted.bam -out ${bestF3%.sam}_stats.txt -plot ${bestF3%.sam}.pdf

# Coverage plot og statistikk uten duplikater	
	weeSAMv1.4 -b ${bestF3%.sam}_sorted.marked.bam -out ${bestF3%.sam}.marked_stats.txt -plot ${bestF3%.sam}_marked.pdf
	
	
    sumAgensRead=$(awk 'FNR > 1 {print $2}' *sumstats.txt| paste -sd+ | bc)
    newR4=$(ls *_tanoti.sam)    
    minor=$(sed -n 3p  *_tanoti_sumstats.txt | cut -d " " -f1 | cut -d'"' -f2)    
    bestMinor=$(sort -t$'\t' -k3 -nr ${newR4%.sam}_stats.txt | grep ${minor}_ -m1 | cut -f1)
    bestMinor_percCov=$(sort -t$'\t' -k3 -nr ${newR4%.sam}_stats.txt | grep ${bestMinor} -m1 | cut -f5)    
    bestMinor_percCov2=${bestMinor_percCov/.*}   

	if [ ${sumAgensRead} -gt ${minAgensRead} ] && [ ${bestMinor_percCov2} -gt 5 ]; then
	
# Coverage plot og statistikk med duplikater for minor
		bestMinor3=$(ls *_tanoti_bestMinor.sam)
		weeSAMv1.4 -b ${bestMinor3%.sam}_sorted.bam -out ${bestMinor3%.sam}_stats.txt -plot ${bestMinor3%.sam}.pdf

# Coverage plot og statistikk uten duplikater	for minor		
		weeSAMv1.4 -b ${bestMinor3%.sam}_sorted.marked.bam -out ${bestMinor3%.sam}.marked_stats.txt -plot ${bestMinor3%.sam}_marked.pdf
		
    fi


cd "${basedir}"


done

echo "Popped som plots"

######## DEL 4 CoveragePlot og Statistikk #### SLUTT ######


######## DEL 5 Identifisere parametere, lage summary og pdf-rapport for hver prøve #### START ######

basedir=$(pwd)
runname=${basedir##*/}

# Går inn i hver mappe og identifiserer ulike parametere og legger det inn i en csv fil 
for dir in $(ls -d Virus*/)
do
    cd ${dir}
#identify & log
    R1=$(ls *_R1*.fastq.gz)
    R2=$(ls *_R2*.fastq.gz)
    newR1=$(ls *val_1.fq)
    newR2=$(ls *val_2.fq)
    newR4=$(ls *_tanoti.sam)
    major=$(sed -n 2p  *_tanoti_sumstats.txt | cut -d " " -f1| cut -d'"' -f2)      
    bestF1=$(sort -t$'\t' -k3 -nr ${newR4%.sam}_stats.txt | grep ${major} -m1 | cut -f1)    
    bestF2="${R1%%_L001*}_v_${bestF1%_H*}"
    bestF3=$(ls *_tanoti_vbest.sam)
    bestF4=$(ls *_tanoti_vbest_sorted.bam) 
    readsb4=$(echo $(zcat ${R1}|wc -l)/2|bc)		#delt på 2 istedenfor 4 for å få reads for R1 og R2
    readsafter=$(echo $(cat ${newR1}|wc -l)/2|bc)
    readstrim=$(echo "scale=2 ; (($readsb4-$readsafter)/$readsb4)*100" | bc)
#   bpb4=$(zcat ${R1} | paste - - - - | cut -f2 | wc -c) 	 #gange 2 for å få bp for R1 og R2   
#   bpb4_2=$(echo "scale=2 ; $bpb4*2" | bc)
#   bpafter=$(cat ${newR1} | paste - - - - | cut -f2 | wc -c)
#   bpafter_2=$(echo "scale=2 ; $bpafter*2" | bc)
#   bptrim=$(echo "scale=2 ; (($bpb4-$bpafter) / $bpb4)*100" | bc)
    wee1113=$(sort -t$'\t' -k3 -nr *_tanoti_vbest_stats.txt | grep -m1 "" | cut -f3)
    mapreadsper=$(echo "scale=2 ; ($wee1113 / $readsafter) *100" | bc)
#   mapbp=$(awk '{s+=$4}END{print s}' ${bestF3%.sam}_sorted_aln.bam)
#   mapbpper=$(echo "scale=2 ; ($mapbp / $bpafter_2) *100" | bc)
    wee1114=$(sort -t$'\t' -k3 -nr *_tanoti_vbest_stats.txt | grep -m1 "" | cut -f5)
    wee1115=$(sort -t$'\t' -k3 -nr *_tanoti_vbest_stats.txt | grep -m1 "" | cut -f8)
   
    percmajor=$(sed -n 2p  *_tanoti_sumstats.txt | cut -d " " -f3)
    percmajor_2=$(echo "scale=2 ; $percmajor*100" | bc)
    sumAgensRead=$(awk 'FNR > 1 {print $2}' *sumstats.txt| paste -sd+ | bc)    

   # wee11=$(ls *_tanoti_vbest.pdf)
    # wee12=$(ls *_tanoti_bestMinor.pdf)    
     
   # bedtools genomecov -ibam ${bestF3%.sam}_sorted.bam -bga > ${bestF3%.sam}_sorted_aln.bam 
   # LengthBelowDepth6=$(awk '$4 <6' *vbest_sorted_aln.bam | awk '{a=$3-$2;print $0,a;}' | awk '{print $5}' | paste -sd+ | bc)
   # LengthBelowDepth30=$(awk '$4 <30' *vbest_sorted_aln.bam | awk '{a=$3-$2;print $0,a;}' | awk '{print $5}' | paste -sd+ | bc)
    RefLength=$(awk 'FNR == 2 {print $2}' *vbest_stats.txt)
   # PercCovAboveDepth5=$(echo "scale=5;(($RefLength-$LengthBelowDepth6)/$RefLength)*100" |bc)    
   # PercCovAboveDepth29=$(echo "scale=5;(($RefLength-$LengthBelowDepth30)/$RefLength)*100" |bc)


    # After removal of duplicates
    wee1120=$(sort -t$'\t' -k3 -nr *_tanoti_vbest.marked_stats.txt | grep -m1 "" | cut -f3)
    wee1121=$(sort -t$'\t' -k3 -nr *_tanoti_vbest.marked_stats.txt | grep -m1 "" | cut -f8)
    wee13=$(ls *_tanoti_vbest_marked.pdf)
    wee14=$(ls *_tanoti_bestMinor_marked.pdf)

    bedtools genomecov -ibam ${bestF3%.sam}_sorted.marked.bam -bga > ${bestF3%.sam}_sorted.marked_aln.bam 
    W_LengthBelowDepth6=$(awk '$4 <6' *vbest_sorted.marked_aln.bam | awk '{a=$3-$2;print $0,a;}' | awk '{print $5}' | paste -sd+ | bc)
    W_LengthBelowDepth30=$(awk '$4 <30' *vbest_sorted.marked_aln.bam | awk '{a=$3-$2;print $0,a;}' | awk '{print $5}' | paste -sd+ | bc)
    W_PercCovAboveDepth5=$(echo "scale=5;(($RefLength-$W_LengthBelowDepth6)/$RefLength)*100" |bc)    
    W_PercCovAboveDepth29=$(echo "scale=5;(($RefLength-$W_LengthBelowDepth30)/$RefLength)*100" |bc)


 #minor genotype
    minor=$(sed -n 3p  *_tanoti_sumstats.txt | cut -d " " -f1 | cut -d'"' -f2)
    percminor=$(sed -n 3p  *_tanoti_sumstats.txt | cut -d " " -f3)
    percminor_2=$(echo "scale=2 ; $percminor*100" | bc)
    bestMinor=$(sort -t$'\t' -k3 -nr ${newR4%.sam}_stats.txt | grep ${minor}_ -m1 | cut -f1)
    bestMinor_percCov=$(sort -t$'\t' -k3 -nr ${newR4%.sam}_stats.txt | grep ${bestMinor} -m1 | cut -f5)    
    bestMinor_percCov2=${bestMinor_percCov/.*}  

if [ ${sumAgensRead} -gt ${minAgensRead} ] && [ ${bestMinor_percCov2} -gt 5 ]; then     
    bestMinor3=$(ls *_tanoti_bestMinor.sam)
    minor2=$(sed -n 3p  *_tanoti_sumstats.txt | cut -d " " -f1 | cut -d'"' -f2)
	wee1116=$(sort -t$'\t' -k3 -nr *_tanoti_bestMinor_stats.txt | grep -m1 "" | cut -f3)    
    wee1117=$(sort -t$'\t' -k3 -nr *_tanoti_bestMinor_stats.txt | grep -m1 "" | cut -f5)
    wee1118=$(sort -t$'\t' -k3 -nr *_tanoti_bestMinor_stats.txt | grep -m1 "" | cut -f8)
   # bedtools genomecov -ibam ${bestMinor3%.sam}_sorted.bam -bga > ${bestMinor3%.sam}_sorted_aln.bam
   # M_LengthBelowDepth6=$(awk '$4 <6' *Minor_sorted_aln.bam | awk '{a=$3-$2;print $0,a;}' | awk '{print $5}' | paste -sd+ | bc)
   # M_LengthBelowDepth30=$(awk '$4 <30' *Minor_sorted_aln.bam | awk '{a=$3-$2;print $0,a;}' | awk '{print $5}' | paste -sd+ | bc)
    M_RefLength=$(awk 'FNR == 2 {print $2}' *bestMinor_stats.txt)
   # M_PercCovAboveDepth5=$(echo "scale=5;(($M_RefLength-$M_LengthBelowDepth6)/$M_RefLength)*100" |bc)    
   # M_PercCovAboveDepth29=$(echo "scale=5;(($M_RefLength-$M_LengthBelowDepth30)/$M_RefLength)*100" |bc)

    # After removal of duplicates
    wee1122=$(sort -t$'\t' -k3 -nr *_tanoti_bestMinor.marked_stats.txt | grep -m1 "" | cut -f3)    
    wee1123=$(sort -t$'\t' -k3 -nr *_tanoti_bestMinor.marked_stats.txt | grep -m1 "" | cut -f5)
    wee1124=$(sort -t$'\t' -k3 -nr *_tanoti_bestMinor.marked_stats.txt | grep -m1 "" | cut -f8)
 
    bedtools genomecov -ibam ${bestMinor3%.sam}_sorted.marked.bam -bga > ${bestMinor3%.sam}_sorted.marked_aln.bam
    WM_LengthBelowDepth6=$(awk '$4 <6' *Minor_sorted.marked_aln.bam | awk '{a=$3-$2;print $0,a;}' | awk '{print $5}' | paste -sd+ | bc)
    WM_LengthBelowDepth30=$(awk '$4 <30' *Minor_sorted.marked_aln.bam | awk '{a=$3-$2;print $0,a;}' | awk '{print $5}' | paste -sd+ | bc)
    WM_PercCovAboveDepth5=$(echo "scale=5;(($M_RefLength-$WM_LengthBelowDepth6)/$M_RefLength)*100" |bc)    
    WM_PercCovAboveDepth29=$(echo "scale=5;(($M_RefLength-$WM_LengthBelowDepth30)/$M_RefLength)*100" |bc)
else
    minor2=$('N/A')
    wee1122=$('N/A')
    wee1117=$('N/A')
    WM_PercCovAboveDepth5=$('N/A')
    wee1124=$('N/A')
    
fi


#write bit
echo "Parameters, ${dir%/}" >> ${dir%/}_summary.csv
#echo "Total_number_of_reads_before_trim:, ${readsb4}"  >> ${dir%/}_summary.csv
#echo "Total_number_of_reads_after_trim:, ${readsafter}" >> ${dir%/}_summary.csv
#echo "Percent_reads_trimmed_removed:, ${readstrim}" >> ${dir%/}_summary.csv
#echo "Total_mapped_${Agens}_reads:, ${sumAgensRead}" >> ${dir%/}_summary.csv
echo "Percent_mapped_reads_of_trimmed:, ${mapreadsper}" >> ${dir%/}_summary.csv # mot den enkelte referansen etter andre runde mapping
#echo "Total_bp_before_trim:, ${bpb4_2}" >> ${dir%/}_summary.csv
#echo "Total_bp_after_trim:, ${bpafter_2}" >> ${dir%/}_summary.csv
#echo "Percent_bp_trimmed_removed:, ${bptrim}" >> ${dir%/}_summary.csv
echo "Majority_genotype:, ${major}" >> ${dir%/}_summary.csv
#echo "Best_hit_from_database:, ${bestF1}" >> ${dir%/}_summary.csv
#echo "Percent_majority_genotype:, ${percmajor_2}" >> ${dir%/}_summary.csv
echo "Number_of_mapped_reads:, ${wee1113}" >> ${dir%/}_summary.csv
#echo "Mapped_bp:, ${mapbp}" >> ${dir%/}_summary.csv
#echo "Percent_mapped_bp_of_trimmed:, ${mapbpper}" >> ${dir%/}_summary.csv
echo "Percent_covered:, ${wee1114}" >> ${dir%/}_summary.csv
#echo "Average_depth:, ${wee1115}" >> ${dir%/}_summary.csv
#echo "Percent_covered_above_depth=5:, ${PercCovAboveDepth5}" >> ${dir%/}_summary.csv
#echo "Percent_covered_above_depth=29:, ${PercCovAboveDepth29}" >> ${dir%/}_summary.csv


# After removal of duplicates
echo "Number_of_mapped_reads_without_duplicates:, ${wee1120}" >> ${dir%/}_summary.csv
echo "Average_depth_without_duplicates:, ${wee1121}" >> ${dir%/}_summary.csv
echo "Percent_covered_above_depth=5_without_duplicates:, ${W_PercCovAboveDepth5}" >> ${dir%/}_summary.csv
echo "Percent_covered_above_depth=29_without_duplicates:, ${W_PercCovAboveDepth29}" >> ${dir%/}_summary.csv


echo "Most_abundant_minority_genotype:, ${minor}" >> ${dir%/}_summary.csv

if [ ${sumAgensRead} -gt ${minAgensRead} ] && [ ${bestMinor_percCov2} -gt 5 ]; then  
#echo "Best_hit_from_database_minor:, ${bestMinor}" >> ${dir%/}_summary.csv
echo "Percent_most_abundant_minority_genotype:, ${percminor_2}" >> ${dir%/}_summary.csv
echo "Number_of_mapped_reads_minor:, ${wee1116}" >> ${dir%/}_summary.csv
echo "Percent_covered_minor:, ${wee1117}" >> ${dir%/}_summary.csv
#echo "Average_depth_minor:, ${wee1118}" >> ${dir%/}_summary.csv
#echo "Percent_covered_above_depth=5_minor:, ${M_PercCovAboveDepth5}" >> ${dir%/}_summary.csv
#echo "Percent_covered_above_depth=29_minor:, ${M_PercCovAboveDepth29}" >> ${dir%/}_summary.csv
echo "Number_of_mapped_reads_minor_without_duplicates:, ${wee1122}" >> ${dir%/}_summary.csv
echo "Average_depth_minor_without_duplicates:, ${wee1124}" >> ${dir%/}_summary.csv
echo "Percent_covered_above_depth=5_minor_without_duplicates:, ${WM_PercCovAboveDepth5}" >> ${dir%/}_summary.csv
echo "Percent_covered_above_depth=29_minor_without_duplicates:, ${WM_PercCovAboveDepth29}" >> ${dir%/}_summary.csv

else
#echo "Best_hit_from_database_minor:, N/A" >> ${dir%/}_summary.csv 
echo "Percent_most_abundant_minority_genotype:, ${percminor_2}" >> ${dir%/}_summary.csv
echo "Number_of_mapped_reads_minor:, N/A" >> ${dir%/}_summary.csv
echo "Percent_covered_minor:, N/A" >> ${dir%/}_summary.csv
#echo "Average_depth_minor:, N/A" >> ${dir%/}_summary.csv
#echo "Percent_covered_above_depth=5_minor:, N/A" >> ${dir%/}_summary.csv
#echo "Percent_covered_above_depth=29_minor:, N/A" >> ${dir%/}_summary.csv
echo "Number_of_mapped_reads_minor_without_duplicates:, N/A" >> ${dir%/}_summary.csv
echo "Average_depth_minor_without_duplicates:,  N/A" >> ${dir%/}_summary.csv
echo "Percent_covered_above_depth=5_minor_without_duplicates:, N/A" >> ${dir%/}_summary.csv
echo "Percent_covered_above_depth=29_minor_without_duplicates:, N/A" >> ${dir%/}_summary.csv


fi


echo "Script_name:, ${script_name1}" >> ${dir%/}_summary.csv
 
#pdf rapport ved bruk av latex
 exec 3<> ${bestF3%%_*}_analysis.tex
        echo "\documentclass[a4paper,15pt]{article}" >&3
        echo "\usepackage[margin=3cm]{geometry}" >&3
        echo "\usepackage{fancyhdr}" >&3
        echo "\fancyhead{}" >&3         # Clears all page headers and footersecho 
        echo "\pagestyle{fancy}" >&3
        echo "\rhead{\today}" >&3   
        echo "\chead{${script_name1}}" >&3       
        echo "\lhead{${bestF3%%_*}}" >&3
        echo "\newcommand\tab[1][1cm]{\hspace*{#1}}" >&3       
        echo "\usepackage{graphicx}" >&3
echo "\makeatletter" >&3     
echo "\def\@seccntformat#1{%
  \expandafter\ifx\csname c@#1\endcsname\c@section\else
  \csname the#1\endcsname\quad
  \fi}" >&3
echo "\makeatother" >&3

echo "\begin{document}" >&3


echo "\section{Trimming and total ${Agens} reads}" >&3          
    echo "\textbf{Percent reads removed with trimming:} \tab" >&3
    echo "${readstrim} \newline" >&3
    echo "\textbf{Total mapped ${Agens} reads:} \tab[3.02cm]" >&3
    echo "${sumAgensRead} \tab" >&3
echo "\section{Main genotype (data without duplicates)}" >&3
        echo "\textbf{Main genotype:} \tab[2.83cm]" >&3
        echo "${major} \newline" >&3
        echo "\textbf{Mapped reads:} \tab[2.93cm]" >&3
        echo "${wee1120} \newline" >&3
        echo "\textbf{Percent covered min. 1 depth:} \tab[0.15cm]" >&3
        echo "${wee1114} \newline" >&3
        echo "\textbf{Percent covered min. 6 depth:} \tab[0.15cm]" >&3
        echo "${W_PercCovAboveDepth5} \newline" >&3
        echo "\textbf{Average depth:} \tab[2.85cm]" >&3
        echo "${wee1121} \newline \newline" >&3
        echo  "\IfFileExists{${wee13}}{\includegraphics[scale=0.65]{${wee13}}} \tab" >&3

echo "\section{Minority genotype (data without duplicates)}" >&3       
        echo "\textbf{Minority genotype:} \tab[2.13cm]" >&3
        echo "${minor2} \newline" >&3
        echo "\textbf{Mapped reads:} \tab[2.93cm]" >&3
        echo "${wee1122} \newline" >&3
        echo "\textbf{Percent covered min. 1 depth:} \tab[0.15cm]" >&3
        echo "${wee1117} \newline" >&3
        echo "\textbf{Percent covered min. 6 depth:} \tab[0.15cm]" >&3
        echo "${WM_PercCovAboveDepth5}  \newline" >&3
        echo "\textbf{Average depth:} \tab[2.85cm]" >&3
        echo "${wee1124} \newline \newline" >&3
        echo  "\IfFileExists{${wee14}}{\includegraphics[scale=0.65]{${wee14}}} \tab" >&3
           

echo "\end{document}" >&3
    exec 3>&-
    pdflatex ${bestF3%%_*}_analysis.tex

#sletter midlertidige filer
rm *analysis.aux
rm *analysis.log
#rm *analysis.tex


    cd "${basedir}"
done

######## DEL 5 Identifisere parametere, lage summary og pdf-rapport for hver prøve #### SLUTT ######


######## GLUE #### START ######
## bruker bam-filene uten duplikater


basedir=$(pwd)
runname=${basedir##*/}

for dir in $(ls -d Virus*/)
do
   cd ${dir}

   docker start gluetools-mysql #starter først gluetools-mysql docker (lagt inn fordi docker stopper å kjøre ved restart av pc)
   newR5=$(ls *_tanoti_vbest_sorted.marked.bam)
   pwd=$(pwd)
  docker run --rm --name gluetools -v ${pwd}:/opt/bams -w /opt/bams --link gluetools-mysql cvrbioinformatics/gluetools:latest gluetools.sh --console-option log-level:FINEST --inline-cmd project hcv module phdrReportingController invoke-function reportBamAsHtml ${newR5} 15.0 ${newR5%.bam}.html

    minor=$(sed -n 3p  *_tanoti_sumstats.txt | cut -d " " -f1 | cut -d'"' -f2)
    percminor=$(sed -n 3p  *_tanoti_sumstats.txt | cut -d " " -f3)
    percminor_2=$(echo "scale=2 ; $percminor*100" | bc)
    sumAgensRead=$(awk 'FNR > 1 {print $2}' *sumstats.txt| paste -sd+ | bc)
    newR4=$(ls *_tanoti.sam)    
    minor=$(sed -n 3p  *_tanoti_sumstats.txt | cut -d " " -f1 | cut -d'"' -f2)    
    bestMinor=$(sort -t$'\t' -k3 -nr ${newR4%.sam}_stats.txt | grep ${minor}_ -m1 | cut -f1)
    bestMinor_percCov=$(sort -t$'\t' -k3 -nr ${newR4%.sam}_stats.txt | grep ${bestMinor} -m1 | cut -f5)    
    bestMinor_percCov2=${bestMinor_percCov/.*} 


    if [ ${sumAgensRead} -gt ${minAgensRead} ] && [ ${bestMinor_percCov2} -gt 5 ]; then   
        M_newR5=$(ls *_tanoti_bestMinor_sorted.marked.bam)
        pwd=$(pwd)
        docker run --rm --name gluetools -v ${pwd}:/opt/bams -w /opt/bams --link gluetools-mysql cvrbioinformatics/gluetools:latest gluetools.sh --console-option log-level:FINEST --inline-cmd project hcv module phdrReportingController invoke-function reportBamAsHtml ${M_newR5} 15.0 ${M_newR5%.bam}.html

fi

    cd "${basedir}"
done

echo "DAA-resistant polymorphisms identified" 

######## GLUE #### STOPP ######




######## DEL 6 Sammenfatte resultater #### START ######


basedir=$(pwd)
runname=${basedir##*/}


mkdir "./${runname}_summaries"
#mkdir "./${runname}_IGV_bam_filer"
#mkdir "./${runname}_summaries/Konsensus-sekvenser"


for dir in $(ls -d Virus*/)
do

	cp ${dir}/*_summary.csv "./${runname}_summaries/"
	cp ${dir}/*.html "./${runname}_summaries/"
	cp ${dir}/*analysis.pdf "./${runname}_summaries/"   

#	cp ${dir}/*_consensus.fa "./${runname}_summaries/Konsensus-sekvenser/"  
	cp ${dir}/*_consensus.fa "./${runname}_summaries/"
	
	
	#cp ${dir}/*tanoti_*_sorted.bam "./${runname}_IGV_bam_filer"
    #cp ${dir}/*tanoti_*_sorted.bam.bai "./${runname}_IGV_bam_filer"
	#cp ${dir}/*sorted.marked.bam "./${runname}_IGV_bam_filer" 
    #cp ${dir}/*sorted.marked.bam.bai "./${runname}_IGV_bam_filer"  
	
done


#lager en fil .tmp for hver prøve hvor alle verdiene legges inn i (uten overskriftene) 
	cd "./${runname}_summaries"

	for f in $(ls *.csv) 
	do
		sed 's/\./,/g' $f | awk 'BEGIN {OFS=","} {print $2}' > $f-5.tmp
    done

echo "Parameters:" >> parameters                            # Lager en fil parameteres hvor alle oversikriftene legges 
#echo "Total number of reads before trim:"  >> parameters
#echo "Total number of reads after trim:" >> parameters
#echo "Percent reads removed with trimming:" >> parameters
#echo "Total mapped ${Agens} reads:" >> parameters
echo "Percent mapped reads of trimmed:" >> parameters
#echo "Total bp before trim:" >> parameters
#echo "Total bp after trim:" >> parameters
#echo "Percent bp trimmed/removed:" >> parameters
echo "Majority genotype:" >> parameters
#echo "Genotype /best hit in database:" >> parameters 
#echo "Percent majority genotype:" >> parameters
echo "Number of mapped reads:" >> parameters
#echo "Mapped bp:" >> parameters
#echo "Percent mapped bp of trimmed:" >> parameters
echo "Percent covered:" >> parameters
#echo "Average depth:" >> parameters
#echo "Percent covered above depth=5:" >> parameters
#echo "Percent covered above depth=29:" >> parameters

echo "Number of mapped reads without duplicates:" >> parameters
echo "Average depth without duplicates:" >> parameters
echo "Percent covered above depth=5 without duplicates:" >> parameters
echo "Percent covered above depth=29 without duplicates:" >> parameters

echo "Most abundant minority genotype:" >> parameters
#echo "Best hit for minor genotype:" >>parameters
echo "Percent most abundant minority genotype:" >> parameters
echo "Number of mapped reads minor:" >>parameters
echo "Percent covered minor:" >>parameters
#echo "Average depth minor:" >>parameters
#echo "Percent covered above depth=5 minor:" >> parameters
#echo "Percent covered above depth=29 minor:" >> parameters
echo "Number of mapped reads minor without duplicates:" >>parameters
echo "Average depth minor without duplicates:" >>parameters
echo "Percent covered above depth=5 minor without duplicates:" >> parameters
echo "Percent covered above depth=29 minor without duplicates:" >> parameters

echo "Script name:" >> parameters

paste parameters *.tmp >> ${runname}_summaries.csv    # verdiene og overskriftene limes inn i en og samme fil

find . -type f -name "*.tmp" -exec rm -f {} \;
find . -type f -name "parameters" -exec rm -f {} \; # sletter de midlertidige filene
rm *summary.csv

cd "${basedir}"

######## DEL 6 Sammenfatte resultater #### SLUTT ######


echo "Takk for at du brukte dette skriptet til å hente ut resultater for ${Agens}"
