isomiRROR
--------------------------------
isomiR – references, overviews and readcounts

Description
------------------------------------------------

isomiRROR is a snakemake-based pipeline to analyze isoforms of miRNA from small RNA Seq data. It generates extensive isomiR references from fasta files of mature miRNA sequences and aligns it with trimmed small RNA NGS reads to create fully comprehensive readcount files of isomiRs. Informations are stored in easily accessible keys within the isomiR identifier and include (in order from left to right):

    hsa-let-7i-5p_Can5_0_X_template_Add3_2_GT_nontemplate_unique_Poly_1_13:T>C_Seed_GAGGTA_M8_G_End3_CTGTTGT_length_24

NNN-TGAGGTAGTAGTCTGTGCTGTT-GT

GGC-TGAGGTAGTAGTTTGTGCTGTT-GGU
 
mature hsa-let-7i-5p sequence (+ precursor sequence)

* name of the canonical mature miRNA, the isomiR originates from
* type of modification on 5’ end compared to the canonical miRNA
    * Add5:  one or more nucleotides added
    * Trim5: one or more nucleotides removed
    * Can5: no modification
* number of nucleotides on 5’ end added (1-3), removed ((-1)-(-6)) or no modification (0)
* sequence of nucleotides on 5’ end that were either added or trimmed in 5’ to 3’ orientation (X for no modification)
* modified 5’ sequence present on precursor sequence (template) or not (nontemplate)
* type of modification on 3’ end compared to the canonical miRNA
    * Add3:  one or more nucleotides added
    * Trim3: one or more nucleotides removed
    * Can3: no modification
* number of nucleotides on 3’ end added (1-3), removed ((-1)-(-6)) or no modification (0)
* sequence of nucleotides on 3’ end that were either added or trimmed in 5’ to 3’ orientation (X for no modification)
* modified 3’ sequence present on precursor sequence (template) or not (nontemplate)
* presence of more isomiRs with the same sequence originating from different canonical miRNA sequences (multi) or not (unique)
    * the most likely isomiR is displayed in the readcount files but other possibilities are stored in multireads.txt
* presence of polymorphic substitutions/ mismatches (Poly) or not (0mm)
* number of polymorphic substitutions/ mismatches
* position and type of nucleotide substituted
* seed sequence (nucleotide from position 2-7 of the isomiR sequence)
* m8 sequence (nucleotide at position 8 of the isomiR sequence)
* 3’ end sequence (last 7 nucleotide of the the isomiR sequence)
* length of the isomiR sequence

Advantages of isomiRROR over other isomiR analysis tools:
    
* works with any mature miRNA sequences independent of species or annotation status (i.e. novel miRNAs)
* availability of precursor sequences is optional and only necessary for defining template status of isomiR modifications (coming soon)
* user’s choice of alignment tool (coming soon)
* complete information on sequence modifications


Installing and running isomiRROR
---------------------------------------------

The sole requirement for running is conda, everything else will be delivered with isomiRROR.

### Step 1: Download and install miniconda3
First you need to download and install miniconda3:
for linux
    
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh

for mac os

    curl https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o Miniconda3-latest-MacOSX-x86_64.sh
    bash Miniconda3-latest-MacOSX-x86_64.sh

### Step 2: Clone the workflow
Clone the worflow to the folder where you run your pipelines

    git clone https://gitlab.lrz.de/Physio/isomiRROR.git

### Step 3: Create an environment based on environment.yaml

    cd isomiRROR
    conda update conda
    conda env create --file environment.yaml

### Step 4: Supply input files
isomiRROR needs 2 types of files to work on: fasta files of mature miRNA and optionally precursor hairpin sequences (with T(hymine) instead of U(racil)) and fastq files of sequencing reads of interest. Fasta files can be found on a number of online databases, most prominently miRBase, but can also be custom compositions of miRNAs of user’s choice.  Fasta files can be stored in the /refs subfolder or in a location of user’s choice, while sample fastq files should go into the /data subfolder.

### Step 5: Complete config.yaml and samples.csv with the missing information
In order to run the pipeline you will need to complete the config.yaml file and the samples.csv file. Both are located in the templates folder and should be moved to the root folder of the experiment.
1. config.yaml – reference files and isomiR length limits
The config.yaml stores the information of where to find the fasta files for the isomiR reference generation as well the minimum and maximum length of isomiRs to consider for alignment.

References
* reference_folder: refs or PATH/TO/LOCAL/FOLDER
* mature_db_file: miRNA
* hairpin_ref_file: hairpin

By default isomiRROR will look  for the reference sequences in the /refs subfolder but users can also specify a different folder on their workstation. In both cases the name of the mature miRNA fasta file (mature_db_file) and their corresponding precursor hairpin fasta file (hairpin_ref_file) need to be given without the .fa ending.

Input variables
* min_length: 16
* max_length: 30

Since sequences below a certain length are very prone to align to multiple locations and represent most likely degraded RNA from other RNA subspecies, it is recommended to limit the range of isomiRs before alignment. Likewise the number of isomiRs above a certain lenths is rapidly declining and computational resources and time can be saved . Out of personal experience I recommend using a minimum length (min_length) of 16 or higher  and a macimum length (max_length) of around 30, but users are encouraged to experiment and find a fitting length range for their own experiments. Since isomiR lengths are stored in the ID found in the final readcount file as well, filtering of certain lengths post isomiRROR is possible as well. 

2. samples.csv
The samples.csv simply holds a list of the sample names that should be alignend to the isomiR reference with each sample on a new line. Sample names should be given without the .fastq ending.

    samples
    sample_name1
    sample_name2

### Step 6: Run isomiRROR!
You are now ready to analyze isomiRs but I highly recommend to take a look at the options that are available since I won't cover everything here.
For now isomiRROR is only able to be run as a complete pipeline but running of single steps (isomiR reference generation, mapping with specific alignment tools etc.) will be implented soon.

Simply run :

    snakemake --use-conda

together with your choice of parametery (i.e --cores) in the root folder of your isomiRROR folder and it will be running from beginning to end.

### Output files

isomiRROR will generate a new isomiR reference fasta file, an alignment index corresponding to the mapping tool of your choice splitted by lengths, a text file containing the possible multireads (isomiRs with identical sequences but originating from different mature miRNA sequences) and of course a readcount file containing a table with isomiR ID in the first column and corresponding readcounts for each sample in the following.