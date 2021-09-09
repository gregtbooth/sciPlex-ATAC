# sciPlex-ATAC-seq analysis

Gregory T. Booth


## Overview 

This repository contains scripts used to perform analyses and generate figures in the manuscript, "High Capacity Sample Multiplexing for Single Cell Chromatin Accessibility Profiling" by Booth et. al. 2021. 

Raw and processed data files used in the manuscript can be accessed through GEO (GSE178953). The aligned fragments files contained in GEO for each experiment can easily be used to create ArchR project objects, which are relied on by most scripts in this repository. The manual for ArchR can be found [here](https://www.archrproject.com/).


## Pre-processing raw sciATAC-seq data 

Because we have made use of the ArchR scATAC-seq analysis framework, preprocessing only needs to create aligned .bam file with cell barcodes contained in the headers for each read pair alignment. The bam file can then be used by the ArchR package to generate Arrow files containing information for all downstream analyses. Hashing data can be appended to the ArchR project to give sample information.

### Two-level sciATAC

For two-level sciPlex-ATAC-seq data, such as that of the chemical screen experiment, preprocessing instructions was performed as previously described by Cusanovich et. al. Nature 2018 and scripts can be found [here](https://github.com/shendurelab/fly-atac).

### Three-level sciATAC

For three-level sciPlex-ATAC-seq data, such as that of the mixed lymphocyte reaction experiment,  preprocessing instructions was performed as previously described by Domcke et. al. Science 2020 and scripts can be found [here](https://github.com/shendurelab/human-atac).


## Pre-processing raw hash-sequencing data

sciPlex-ATAC-seq libraries contain hash reads which identify sample origins. We have included scripts for extracting these hash reads from fastq files and organizing them with respect to cell barcodes in the "Preprocess" directory. Because of the differences in library construction, hash reads from two- and three-level sciPlex-ATAC-seq experiments require different handling.

### Two-level hash reads
In two-level experiments, hash reads are sequenced separately from ATAC reads. Fastq files are first created containing the all cell barcodes in the header, then these hash reads are organized and summarized into a table with hash ID counts per cell barcode. 

### Three-level hash reads
For three-level experiments, hash reads are sequenced along with ATAC fragments. Therefore, hash reads are simply extracted from fastq files, by looking for reads which match expected hash ID sequences. These hash reads are organized and summarized into a table with hash ID counts per cell barcode. 


## Downstream Analysis 

In general, the analyses performed took advantage of the ArchR R package (Granja et. al. 2021) to process and store much of the data. Monocle3 (Qui et. al. 2018) and Cicero (Pliner et. al. 2018) were also heavily used for select analyses.

Separate directories are included for each of the four experiments presented in the manuscript. Each directory contains "NoteBook" scripts with numbers corresponding to a general order in which to perform each analysis. For example, in many cases, later notebooks (e.g. NB3) will rely on the output of earlier notebooks (e.g. NB1). 

For the chemical screen and mixed lymphocyte experiments, The analysis starts with the following: 
- assigning nuclear labels to each cell barcode
- Creating an Arrow project from bam or fragment files
- Creating an ArchR project and filtering cells

The ArchR project containing hash information for each cell in the filtered set is then used as the input for most downstream analyses. 
