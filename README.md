# my_CUTnRUN

*A modular, reproducible pipeline for CUT&RUN processing, alignment, QC, peak calling, and downstream analysis.*

This repository provides an end-to-end workflow for analyzing CUT&RUN sequencing data—from raw FASTQ files to normalized signal tracks, peak calls, and annotated differential enrichment summaries. The pipeline is designed for HPC environments, supports spike-in (E. coli) normalization, and is fully modular so each step can be run independently.

---

## Features

* **Complete CUT&RUN workflow**

  * Adapter trimming with Trim Galore
  * Combined host + spike-in alignment (Bowtie2)
  * Host-only + spike-in–only BAM separation
  * MAPQ filtering and optional blacklist removal
* **QC & spike-in assessment**

  * Per-sample host/spike-in read counts
  * Useful for assessing enzyme activity or digestion consistency
* **Signal tracks**

  * `bamCoverage` RPGC-normalized bigWig files
  * Optional bigWig → bedGraph conversion
* **Peak calling**

  * SEACR (stringent or relaxed) with IgG controls
  * Per-sample peak sets + merged consensus sets
* **Downstream R utilities**

  * Gene annotation with ChIPseeker + TxDb
  * Unified peak tables across replicates
  * Basic differential CUT&RUN summaries (e.g., WT vs KO)

