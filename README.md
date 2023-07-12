# scViralQuant
The goal of this python package is to accurately align and quantify viral reads for 10x single cell data.  This software integrates viral UMI counts and viral UMI gene counts into 10x files (features.tsv.gz and matrix.mtx.gz in the filtered_feature_bc_matrix and  raw_feature_bc_matrix folders) so that softwares such as seurat can be used to analyze data. The software uses outside software samtools, bowtie2, and htseq to quantify viral reads.  Default parameters are used when aligning reads to the viral genome and and non default htseq-count include using the --intersection-nonempty parameter.

### Dependencies
----------------
scViralQuant is tested to work in python 3.6.8

* argsparse
* htseq
* pandas 
* scipy
* seaborn

Currently only set to run on linux & mac but has only been tested on a linux platform

### Installation
----------------
```bash
python setup.py install
```

### Running scViralQuant 
------------------------
Set the following parameters 
 
* -10x = Path/to/10x/sample
* -op = path/for/results 
* -p = number of processors (defualt = 1)
* -p2genome = path/to/viral/fastafilefolder - in this folder should be at most 2 files 1) the fasta file with the viral genome sequence and 2) viral gtf file (not mandatory this can just run with fasta file with the viral genome).  Once bowtie2 indexes are made folder can be reused with out having to remake bowtie indexes.  Note: In the fasta file the header will be used to quantify the number of viral copies, it is recommended that if the fasta header is a complicated name it be simplified (i.e. > HIV_virus)
* -align = alignment tool bbmap (default) or bowtie2
* -overwrite = will overwrite the files directly in the 10x filtered_feature_bc_matrix folder.  If not specified a copy of this folder will be made and then information about viral count information will be added to the files in this copied folder 
* --bbmap_params = parameters specific to the bbmap alignment tool  
* --bowtie2_params = parameters specific to the bowtie2 alignment tool
* --tmp_removal = if specified will remove the temporary directory (_tmp/) of files used by scviralquant (these files can be large so if space is an issue these should be deleted)

Example runs:
```bash 
 scviralquant -10x Path/to/10x/sample -op path/for/results -p 8 -p2genome path/to/viral/fastafilefolder
```
```bash 
 scviralquant -10x Path/to/10x/sample -op path/for/results -p 8 -p2genome path/to/viral/fastafilefolder --tmp_removal
```
```bash 
 scviralquant -10x Path/to/10x/sample -op path/for/results -p 8 -p2genome path/to/viral/fastafilefolder --bbmap_params "--semiperfectmode"
```
```bash 
 scviralquant -10x Path/to/10x/sample -op path/for/results -p 8 -p2genome path/to/viral/fastafilefolder  -align bowtie2 --bowtie2_params "--very-sensitive  --non-deterministic"
```

### Output files 
----------------
Output files by scViralQuant

* viral_copy.png - violin plot showing number of cells with virus
* viral_genes.png - violin plot showing cells with viral gene expression 
* virus_al_counts.csv - total number of reads mapping to the virus in each cell 
* virus_al_gene_counts.csv - number of reads mapping to viral genes in each cell 
* virus_al.bam - reads mapped to virus (no unmapped reads)
* virus_al_mapped.sam - reads mapped to only virus (no unmapped reads)
* virus_al_sort.bam - sorted reads mapped to virus 
* virus_al_sort.bam.bai - index file to virus_al_sort.bam
* virus_al_sort_counts.sam - htseq output reads mapping to virus
* virus_genes_al_sort_counts.sam - htseq output reads mapping to individual virus genes (will not be produced if viral gtf is not provided)
* Overwrites original 10x data provided to include viral counts and viral gene counts (if gtf file is provided)

### LICENSE
----------------
BSD - 3-Clause Copyright 2023 University of Washington