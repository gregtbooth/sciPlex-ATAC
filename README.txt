sciPlex-ATAC-seq analysis

Gregory T. Booth

This repository contains scripts and processed files used to generate analyses and figures in the manuscript, "High Capacity Sample Multiplexing for Chromatin Accessibility Profiling Within Single Cells" by Booth et. al. 2021. 

In general, the analyses performed took advantage of the ArchR R package (Granja et. al. 2021) to process and store much of the data. Monocle3 (Qui et. al. 2018) and Cicero (Pliner et. al. 2018) were also heavily used for select analyses.

Separate directories are included for each of the four experiments presented in the manuscript. Each directory contains "NoteBook" scripts with prefix numbers corresponding to a general order in which to perform each analysis. For example, in many cases, later notebooks (e.g. NB3) will rely on the output of earlier notebooks (e.g. NB1). 

For the chemical screen and mixed lymphocyte experiments, The analysis starts with the following: 
- assigning nuclear labels to each cell barcode
- Creating an Arrow project from bam or fragment files
- Creating an ArchR project and filtering cells

The ArchR project containing hash information for each cell in the filtered set is then used as the input for most downstream analyses. 

