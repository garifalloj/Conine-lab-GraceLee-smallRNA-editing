#!/bin/sh
#SBATCH -t 15:00
#SBATCH --mem=2G

#generate STAR mapping scripts for both Grace and Lexi RNAseq datasets
#must use MD tag in samtools view to get mismatches

fastq_dir=/miR-200_and_miR-465_112025/fastq/
cd $fastq_dir

for fastq_file in *_R1_001.fastq.gz; do
    (
    #cut "R1.fastq.gz" or "R2.fastq.gz" from the file name
    sample_name=$(echo $fastq_file | sed 's/_R1_001\.fastq\.gz//')
    echo -e "#!/bin/sh
#SBATCH -t 8:00:00
#SBATCH --mem=30G
#SBATCH -c 1

/STAR --runThreadN 4 --genomeDir /GRCm38_gencodeVM25.50bpRL --readFilesCommand zcat --outFilterMultimapNmax 10 --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.3 --outSAMattributes NH HI AS nM MD --readFilesIn ${fastq_dir}${fastq_file} ${fastq_dir}${sample_name}_R2_001.fastq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /miR-200_and_miR-465_112025/STAR/${sample_name}_" > /miR-200_and_miR-465_112025/STAR/STAR_${sample_name}.sh
    
    )
done 