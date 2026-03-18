#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH --mem=32G
#SBATCH -c 1


#script to use picard to mark and remove duplicates from BAM files
#only include mapped and properly paired reads
#then filter for mapping quality > 20


module load picard/2.26.10-Java-15.lua
module load SAMtools/1.16.1-GCC-11.3.0.lua



bam_file_directory=/miR-200_and_miR-465_112025/STAR


cd $bam_file_directory

for bam_file in *Aligned.sortedByCoord.out.bam; do 
    
    sample_name=$(echo $bam_file | sed 's/_S[0-9]\+_L001_Aligned.sortedByCoord.out.bam//')
    
    
    #use picard to remove PCR duplicates from coordinate sorted bams, output to Picard2
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT=$bam_file OUTPUT=/miR-200_and_miR-465_112025/picard_removeDups/${sample_name}_remDups.bam VALIDATION_STRINGENCY=LENIENT TMP_DIR=/miR-200_and_miR-465_112025/picard_removeDups/tmp METRICS_FILE=${sample_name}_dup.txt REMOVE_DUPLICATES=true

    #samtools filter to only include properly paired reads
    samtools view -f 3 -b /miR-200_and_miR-465_112025/picard_removeDups/${sample_name}_remDups.bam > /miR-200_and_miR-465_112025/picard_removeDups/${sample_name}_remDups_properPairs.bam

    #samtools filter to only include reads with mapping quality > 20
    samtools view -q 20 -b /miR-200_and_miR-465_112025/picard_removeDups/${sample_name}_remDups_properPairs.bam > /miR-200_and_miR-465_112025/picard_removeDups/${sample_name}_remDups_properPairs_mapq20.bam

    #also need to index the bam files
    #sailor pipeline will index if NOT present
    samtools index /miR-200_and_miR-465_112025/picard_removeDups/${sample_name}_remDups_properPairs_mapq20.bam

done