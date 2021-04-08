# Benchmarking-of-small-RNA-seq
Analysis workflow for benchmarking of small RNA-seq library preparation protocols

**Samples**   
Raw fastq files are available for download at Gene Expression Omnibus with accession [GSE149513](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149513). 

**Workflow**   
Initial trimming and mapping was done with bash scripts. After installation of packages in environment.yml file, whole bash part should be executable with main.sh script with files from GEO. Counting and analysis were done in R. Some of R scripts expect files which generation is not included in the workflow. Such files are uploaded here in "src" folder.

IsomiRROR package for isomiR mapping was downloaded from [gitlab](https://gitlab.lrz.de/Physio/isomiRROR).


