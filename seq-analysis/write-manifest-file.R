# Read in SRA sample IDs and list of raw fastq files to generate a manifest file ready for qiime2 import

# Required libraries
library(tidyverse)

# Import all files in currect directory
paths <- as.data.frame(list.files(pattern = ".fastq.gz", full.names = TRUE))
colnames(paths)[1]<-"FASTQ"

# Get local path & add to dataframe
path_to_files <- getwd()
paths$PATH <- path_to_files

# Extract sample ID
paths_run <- paths %>% 
## Use for sampleid_L001_Rx_001.fastq.gz
#	mutate(SAMPLEID = str_replace(FASTQ, "(\\w*?)_L001_R((\\w+))_001.fastq.gz","\\1"))
## Use for _R1.fastq.gz
	mutate(SAMPLEID = str_replace(FASTQ, "(\\w*?)_R((\\d)).fastq.gz","\\1"))
#
## ^See wildcard options on this line to modify how R script pulls out your sample IDs from fastq files

paths_run$SAMPLEID <- gsub("./","", paths_run$SAMPLEID)
paths_run$FASTQ <- gsub("./","", paths_run$FASTQ)

# Write full path
paths_run$FULL_PATH <- paste(paths_run$PATH, paths_run$FASTQ, sep="/")

forward <- paths_run %>%
	filter(grepl("R1_001.fastq.gz|_1.fastq.gz|R1.fastq.gz", FASTQ)) %>%
	select(SAMPLEID, `forward-absolute-filepath` = FULL_PATH)

# Write manifest file
manifest <- paths_run %>%
	filter(grepl("R2_001.fastq.gz|_2.fastq.gz|R2.fastq.gz", FASTQ)) %>%
	select(SAMPLEID, `reverse-absolute-filepath` = FULL_PATH) %>%
	right_join(forward) %>%
	select('sample-id' = SAMPLEID, `forward-absolute-filepath`, `reverse-absolute-filepath`)

# Write output as a manifest file
write.table(manifest, "manifest",quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t")


