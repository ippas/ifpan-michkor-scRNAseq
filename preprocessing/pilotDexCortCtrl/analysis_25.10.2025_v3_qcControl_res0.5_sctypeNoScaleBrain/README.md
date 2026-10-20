# scRNA-seq Analysis v3 — qcControl_res0.5_sctypeNoScaleBrain (25.10.2025)

This folder contains scripts and results for a **partial and exploratory re-analysis** of the scRNA-seq dataset.

## Overview
- **Version:** v3  
- **Resolution:** 0.5  
- **Annotation:** ScType (applied to existing clustering, no scaling)  
- **Dataset:** QC-filtered brain subset  

## Purpose
This version is **not a final analysis**.  
It was prepared **for presentation purposes** and to **define the goals and direction** of a more comprehensive scRNA-seq study.  

Specifically, it aims to:
- Validate and refine preliminary ScType cell-type annotations  
- Examine cluster-level QC and marker distributions  
- Identify analytical steps that will require improvement or re-run in the final version  

## Notes
- The starting point is a preprocessed Seurat object generated in `analysis_12.08.2025_v2`.  
- No new normalization, scaling, or reclustering is performed here.  
- Scripts assume that embeddings and metadata are already available.  

## Folder structure
