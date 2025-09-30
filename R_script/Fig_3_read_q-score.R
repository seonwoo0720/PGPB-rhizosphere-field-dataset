######################################################################################################################
####### Field-scale 16S rRNA dataset on rhizosphere microbial community under Rhodobacter capsulatus treatment #######
######################################################################################################################
# Date of creation: 2025-07-15

# ================================================================================================
# Data availability:
# - Raw sequencing data are not included in this repository due to large file sizes.  
# - Data can be downloaded from the NCBI Sequence Read Archive (SRA) under accession number [PRJNA1332407].
# ================================================================================================

# ================================================================================================
# Working directory:
# - Ensure that the "Q-mean" folder exists in your current working directory.
# - All FASTQ/FASTQ.GZ files for quality assessment should be placed inside "Q-mean" or its subfolders.
# - For more details, please refer to the README file.
# ================================================================================================

# Fig 3. Distribution of Mean Q-scores
# ================================================================================================
# This script calculates and visualizes the distribution of mean Q-scores from raw FASTQ files
# generated during 16S rRNA sequencing.
#
# Workflow:
#   1. Recursively search the base directory ("Q-mean") for FASTQ/FASTQ.GZ files.
#   2. For each file, extract quality scores using ShortRead.
#   3. Compute mean Q-scores per read and organize results into a combined dataset.
#   4. Generate density plots of mean Q-score distributions for each crop group.
#   5. Save the plot as a high-resolution PNG file for publication use.
#
# Key Notes:
#   - Input: FASTQ files (compressed or uncompressed) stored in subfolders by crop.
#   - Output: "mean_qscore_all_reads.png" (28 x 20.2 cm, 600 dpi).
#   - X-axis: Mean Q-score per read (range 1–40).
#   - Y-axis: Density (proportion of reads).
#   - Fill color: Crop-specific grouping (folder name).
#   - Interpretation: Provides an overview of sequencing read quality across experimental groups.
# ================================================================================================

# Load required libraries
library(ShortRead)
library(ggplot2)
library(data.table)

# Set your base directory here
base_dir = "Q-mean"

# Recursively find all FASTQ or FASTQ.GZ files
fastq_files = list.files(base_dir, pattern = "\\.fastq(\\.gz)?$", recursive = TRUE, full.names = TRUE)

# Function to compute mean Q-scores
get_mean_qscores = function(file) {
  fq = readFastq(file)
  quals = as(quality(fq), "matrix")
  means = rowMeans(quals)
  return(data.frame(
    file = file,
    folder = basename(dirname(file)),
    mean_qscore = means
  ))
}

#Process all files and combine

all_data = rbindlist(lapply(fastq_files, get_mean_qscores))

# Plot using ggplot2
all = ggplot(all_data, aes(x = mean_qscore, fill = folder)) +
  geom_density(alpha = 0.5) +
  theme_bw(base_size = 16) +
  scale_x_continuous(
    breaks = seq(0,40, by =5),
    limits = c(1,40)
  )+
  labs(
    title = "Distribution of Mean Q-scores (Overall)",
    x = "Mean Q-score of reads",
    y = "Proportion",
    fill = "Crops"
  )
all

ggsave("mean_qscore_all_reads.png", all, width = 28, height = 20.2, units = "cm", dpi = 600)

