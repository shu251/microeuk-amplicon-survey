# Survey of microbial eukaryotes at hydrothermal vents
## Details on sequencing methods
Main data analysis location <https://shu251.github.io/microeuk-amplicon-survey/>

Sequence survey spanning three deep-sea hydrothermal vent habitats.

*Data sets*
-   Axial Seamount
-   Gorda Ridge
-   Mid-Cayman Rise

## 1. QIIME2 analysis for individual sequence runs

Analysis run with [QIIME2](https://docs.qiime2.org/2021.11/) v2021.4.

### 1.1 Process individual sequence runs

Code available in this repo includes bash script submitted to HPC at WHOI (via SLURM).
Ahead of time, all sequences were trimmed and assessed for quality with fastqc (and multiqc) using [this snakemake workflow](https://github.com/shu251/qc-trim). Manifest files were created using the Rscript ```write-manifest-file.R```, when R script is run in the same directory as all input sequences, it will create a QIIME2 formatted manifest file.

#### Axial Seamount & Gorda Ridge

Process 18S rRNA gene sequences from Axial Seamount and Gorda Ridge.
Raw sequences can be found _here_.
```
# Grab trim stats from cutadapt
qiime demux summarize \
        --i-data gr-axial-pe-trimmed.qza \
        --o-visualization gr-axial-pe-trimmed.qzv

# Run dada2
qiime dada2 denoise-paired \
        --i-demultiplexed-seqs gr-axial-pe-trimmed.qza \
        --p-trunc-len-f 260 \
        --p-trunc-len-r 225 \
        --p-max-ee-f 2 \
        --p-max-ee-r 2 \
        --p-min-overlap 10 \
        --p-pooling-method independent \
        --p-n-reads-learn 100000 \
        --p-n-threads 10 \
        --p-chimera-method pooled \
        --o-table /vortexfs1/scratch/sarahhu/gr-axial-asv-table.qza \
        --o-representative-sequences /vortexfs1/scratch/sarahhu/gr-axial-ref-seqs.qza \
        --o-denoising-stats /vortexfs1/scratch/sarahhu/gr-axial-dada2-stats.qza

# Convert to TSV ASV table
qiime tools export \
        --input-path /vortexfs1/scratch/sarahhu/gr-axial-asv-table.qza \
        --output-path /vortexfs1/scratch/sarahhu/gr-axial-output/

biom convert -i /vortexfs1/scratch/sarahhu/gr-axial-output/feature-table.biom \
        -o /vortexfs1/scratch/sarahhu/gr-axial-output/gr-axial-asv-table.tsv \
        --to-tsv

# Get dada2 stats
qiime metadata tabulate \
       --m-input-file /vortexfs1/scratch/sarahhu/gr-axial-dada2-stats.qza \
       --o-visualization /vortexfs1/scratch/sarahhu/gr-axial-dada2-stats.qzv

```

#### Mid-Cayman Rise

```
# Import data
qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path manifest-files/manifest-mcr-samples \
        --output-path mcr-samples-pe.qza \
        --input-format PairedEndFastqManifestPhred33V2

# Get starting stats on input sequences
qiime demux summarize \
        --i-data mcr-samples-pe.qza \
        --o-visualization mcr-samples-pe.qzv

# Remove 18S V4 primers
qiime cutadapt trim-paired \
        --i-demultiplexed-sequences mcr-samples-pe.qza \
        --p-cores 4 \
        --p-front-f CCAGCASCYGCGGTAATTCC \
        --p-front-r ACTTTCGTTCTTGATYRA \
        --p-error-rate 0.1 \
        --p-overlap 3 \
        --p-match-adapter-wildcards \
        --o-trimmed-sequences mcr-samples-pe-trimmed.qza 

# Grab trim stats from cutadapt
qiime demux summarize \
        --i-data mcr-samples-pe-trimmed.qza \
        --o-visualization mcr-samples-pe-trimmed.qzv

# Run dada2
qiime dada2 denoise-paired \
        --i-demultiplexed-seqs /vortexfs1/omics/huber/shu/microeuk-survey-18S/tag-seq-analysis/mcr-samples/mcr-samples-pe-trimmed.qza \
        --p-trunc-len-f 260 \
        --p-trunc-len-r 225 \
        --p-max-ee-f 2 \
        --p-max-ee-r 2 \
        --p-min-overlap 10 \
        --p-pooling-method independent \
        --p-n-reads-learn 100000 \
        --p-n-threads 10 \
        --p-chimera-method pooled \
        --o-table /vortexfs1/scratch/sarahhu/mcr-samples-asv-table.qza \
        --o-representative-sequences /vortexfs1/scratch/sarahhu/mcr-samples-ref-seqs.qza \
        --o-denoising-stats /vortexfs1/scratch/sarahhu/mcr-samples-dada2-stats.qza

# Convert to TSV ASV table
qiime tools export \
        --input-path /vortexfs1/scratch/sarahhu/mcr-samples-asv-table.qza \
        --output-path /vortexfs1/scratch/sarahhu/mcr-samples-output/
        
biom convert -i /vortexfs1/scratch/sarahhu/mcr-samples-output/feature-table.biom \
        -o /vortexfs1/scratch/sarahhu/mcr-samples-output/mcr-samples-asv-table.tsv \
        --to-tsv

# Get dada2 stats
qiime metadata tabulate \
       --m-input-file /vortexfs1/scratch/sarahhu/mcr-samples-dada2-stats.qza \
       --o-visualization /vortexfs1/scratch/sarahhu/mcr-samples-dada2-stats.qzv

```
#### Additional controls

```
# Run dada2
qiime dada2 denoise-paired \
        --i-demultiplexed-seqs mcr-ctrl-pe-trimmed.qza \
        --p-trunc-len-f 260 \
        --p-trunc-len-r 225 \
        --p-max-ee-f 2 \
        --p-max-ee-r 2 \
        --p-min-overlap 10 \
        --p-pooling-method independent \
        --p-n-reads-learn 100000 \
        --p-n-threads 8 \
        --p-chimera-method pooled \
        --o-table /vortexfs1/scratch/sarahhu/mcr-ctrl-asv-table.qza \
        --o-representative-sequences /vortexfs1/scratch/sarahhu/mcr-ctrl-ref-seqs.qza \
        --o-denoising-stats /vortexfs1/scratch/sarahhu/mcr-ctrl-dada2-stats.qza

qiime tools export \
        --input-path /vortexfs1/scratch/sarahhu/mcr-ctrl-asv-table.qza \
        --output-path /vortexfs1/scratch/sarahhu/mcr-ctrl-output/

biom convert -i /vortexfs1/scratch/sarahhu/mcr-ctrl-output/feature-table.biom \
        -o /vortexfs1/scratch/sarahhu/mcr-ctrl-output/mcr-ctrl-asv-table.tsv \
        --to-tsv

# Get dada2 stats
qiime metadata tabulate \
       --m-input-file /vortexfs1/scratch/sarahhu/mcr-ctrl-dada2-stats.qza \
       --o-visualization /vortexfs1/scratch/sarahhu/mcr-ctrl-dada2-stats.qzv
```


## 2. Merge ASV tables & assign taxonomy

After merging ASV tables, assign taxonomy using the [PR2](https://pr2-database.org) database v.4.14 with vsearch.

```
qiime feature-table merge \
        --i-tables /vortexfs1/scratch/sarahhu/mcr-samples-asv-table.qza \
        --i-tables /vortexfs1/scratch/sarahhu/mcr-ctrl-asv-table.qza \
        --i-tables /vortexfs1/scratch/sarahhu/gr-axial-asv-table.qza \
        --o-merged-table /vortexfs1/scratch/sarahhu/microeuk-merged-asv-table.qza


qiime feature-table merge-seqs \
        --i-data /vortexfs1/scratch/sarahhu/mcr-samples-ref-seqs.qza \
        --i-data /vortexfs1/scratch/sarahhu/mcr-ctrl-ref-seqs.qza \
        --i-data /vortexfs1/scratch/sarahhu/gr-axial-ref-seqs.qza \
        --o-merged-data /vortexfs1/scratch/sarahhu/microeuk-merged-ref-seqs.qza

# Assign taxonomy
qiime feature-classifier classify-consensus-vsearch \
        --i-query /vortexfs1/scratch/sarahhu/microeuk-merged-ref-seqs.qza \
        --i-reference-reads /vortexfs1/omics/huber/db/pr2-db/pr2_version_4.14_seqs.qza \
        --i-reference-taxonomy /vortexfs1/omics/huber/db/pr2-db/pr2_version_4.14_tax.qza \
        --o-classification /vortexfs1/scratch/sarahhu/microeuk-merged-taxa.qza \
        --p-threads 8 \
        --p-maxaccepts 10 \
        --p-perc-identity 0.8 \
        --p-min-consensus 0.70

# Convert to TSV ASV table
qiime tools export \
        --input-path /vortexfs1/scratch/sarahhu/microeuk-merged-asv-table.qza \
        --output-path /vortexfs1/scratch/sarahhu/microeuk-merged-output/

biom convert -i /vortexfs1/scratch/sarahhu/microeuk-merged-output/feature-table.biom \
        -o /vortexfs1/scratch/sarahhu/microeuk-merged-output/microeuk-merged-asv-table.tsv \
        --to-tsv

qiime tools export \
        --input-path /vortexfs1/scratch/sarahhu/microeuk-merged-taxa.qza \
        --output-path /vortexfs1/scratch/sarahhu/microeuk-merged-output/

```

## 3. Import into R, compile

Move to R to compile, use phyloseq, and additional QC steps. All code is in ```asv-data-processing.Rmd```

## 4. Remove contaminate sequences