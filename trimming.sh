source /home/jonbra/miniconda3/etc/profile.d/conda.sh
conda activate TrimGalore
# Example usage
# bash trimming.sh 472Virus200604-1_S7_L001_R1_001.fastq.gz 472Virus200604-1_S7_L001_R2_001.fastq.gz

######## DEL 1 Trimming #### START ######

echo "Starter trimming"
trim_galore -q 30 --dont_gzip --length 50 --paired $1 $2
echo "Read ferdig trimmet" 

######## DEL 1 Trimming #### SLUTT ######
conda deactivate
