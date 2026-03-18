#!/bin/sh
#script will generate the per library sailor .json and sailor.sh files for workflow_sailor
#this uses the updated duplicates removed, properly paired, and mapq20 filtered BAM files


BAM_directory=/miR-200_and_miR-465_112025/picard_removeDups


cd $BAM_directory
for bam_file in *_remDups_properPairs_mapq20.bam; do
    (
				
                #remove "_remDups_properPairs_mapq20" from sample name
                sample_name=$(echo $bam_file | sed 's/_remDups_properPairs_mapq20.bam//')
            
				echo -e "#!/bin/sh
#SBATCH -t 4:00:00
#SBATCH --mem=10G
#SBATCH -c 2

#might need full path to the miniconda environment
source /home/user/miniconda3/etc/profile.d/conda.sh
conda activate marine_environment

python /MARINE/marine.py \
--bam_filepath ${BAM_directory}/${bam_file} \
--output_folder /miR-200_and_miR-465_112025/MARINE/${sample_name} \
--cores 2 \
--sailor \
--bedgraphs "AG" \
--annotation_bedfile_path /gencode.vM25.annotation.genes.bed \
--contigs "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrM","chrX","chrY" \
--strandedness 0 \
--paired_end" > /miR-200_and_miR-465_112025/MARINE/"${sample_name}".MARINE.sh;	 

    )
done

