# PTCL ctDNA : Shallow WGS Analysis

**Thesis title:** Shallow whole genome sequencing of ctDNA in 
T-cell lymphoma focusing on copy number variations and fragmentomics

**Author:** Noor Shaban 
**Institution:** University of Leicester  
**Year:** 2026  

## Project Overview

Analysis pipeline for plasma sWGS ctDNA profiling in peripheral T-cell lymphoma (PTCL), across two analytical arms:
1. CNA profiling — tumour fraction estimation (ichorCNA), manual curation, absolute copy number scaling (RASCAL), recurrent CNA detection (GISTIC2), and arm-level aneuploidy scoring ( FGA and AS ). 
2. Fragmentomics — fragment-length metrics, tumour fraction correlation
   
A sub-analysis compares plasma ctDNA against 16 matched FFPE tissue sWGS samples.

Raw BAM files are not included (patient-level data).

# Cohort

Total plasma samples	34 (Batch 1: n=5, Batch 2: n=16, Batch 3: n=13)
Sequencing	illumina Novaseq, GRCh38
CNA cohort	n=20 (ichorCNA TF less than 3% excluded 13 and one was above threshold but excluded due to incorrect ploidy )
Fragmentomics cohort	n=33 (B2_S11 excluded "degraded sample")
Matched FFPE pairs	16 from Batch 2 only 


# Dependencies
Tool	                      Version
R	                          4.3.1
QDNAseq	                    1.36.0 (Bioconductor 3.17)
ichorCNA	                 broadinstitute/ichorCNA (GitHub)
RASCAL	                   0.7.0 (crukci-bioinformatics/rascal)
GISTIC2                     Apptainer container
SAMtools	                  1.18
Picard MarkDuplicates    	  2.18.23
HMMcopy readCounter       	Built from source

All scripts were run on the University of Leicester ALICE HPC cluster (Rocky Linux 9, SLURM).
