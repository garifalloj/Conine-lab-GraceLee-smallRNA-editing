#filter out any SNPs from dbSNP from all potential MARINE outputs.
#Note: this removes any edit sites that are known SNP locations.


library(GenomicRanges)
library(rtracklayer)
library(dplyr)

sample_name_dirs <- list.dirs(".", recursive = F, full.names = F)


#read in each as granges objects, and examine
gr_dbSNP <- rtracklayer::import("/references/mm10_dbsnp_combined_reformat.noZeroLen.bed")

#edit filtered to remove these overlapping

#for loop 
for (sampleName_rna_edit_dir in sample_name_dirs){
  final_filtered_edits_annotated <- read.delim(paste0("./", sampleName_rna_edit_dir, "/final_filtered_site_info_annotated.tsv"), header = T, sep = '\t')
  
  gr_edits <- GRanges(
    seqnames = final_filtered_edits_annotated$contig,
    ranges   = IRanges(
      start = final_filtered_edits_annotated$position,
      end   = final_filtered_edits_annotated$position
    )
  )
  
  
  keep <- !overlapsAny(gr_edits, gr_dbSNP)
  
  
  
  edits_filtered <- final_filtered_edits_annotated[keep, ]
  write.table(edits_filtered, file = paste0(sampleName_rna_edit_dir, "/", sampleName_rna_edit_dir, "_filtered_site_annotated_dbSNPfilt.tsv"), row.names = FALSE, sep = '\t', quote = FALSE)
  
}




