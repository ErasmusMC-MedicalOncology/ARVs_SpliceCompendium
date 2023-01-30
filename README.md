# Quantification of AR splice-variants (ARVs) incl. UTR variants in mCRPC (CPCT-02).

This repository contains the custom code (R) used to perform all analysis of the study as further detailed by Isebia and van Riet et al. in <JOURNAL> (2023): [TITLE](https://www.google.com/).

All processed (and cleaned) data as-used in the presented analysis and figures has been published alongside the manuscript (Suppl. Table 1).

The raw and pre-processed data of the CPCT-02 WGS and RNA-Seq cohort can be requested under the data-request **DR-071** from the Hartwig Medical Foundation (HMF): https://www.hartwigmedicalfoundation.nl/en/applying-for-data/.

The commands and tools used in analyzing the RNA-Seq data (starting from paired-end .fastq) are described within the manuscript. In addition, the workflow can be found in this repository under *Misc/*.

All provided code in the R/ folder was performed on Debian 11 with the following *R* version:
```R
R version 4.2.2 (2022-10-31) -- "Innocent and Trusting"
Copyright (C) 2022 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)
```

Required R dependencies are listed on the top of each script.