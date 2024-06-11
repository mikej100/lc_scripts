Bash utilities for processing chromatin assay results such as ChIP-seq data.
Written by Mike Jennings, based on scripts from \<to be named>.

Set up
======

~/.bashrc
----------
In your `~/.bashrc file`, put the following line. This creates the SCRIPTS environment variable which is used by scripts.

    SCRIPTS=<path to directory with these scripts>
For ease of use, you may want put this directory in your PATH, either in your session, or in your .bashrc file.Then you can type the name of the script and it will run, regardless of your
current work directory. 

    PATH=$SCRIPTS:$PATH

Project directory
-----------------
Subdirectory for each track of data (model)
The name of the subdirectory will be used for track name of aggregated data 
such as peak call

`models.cfg` - file containing the names of each model directory in the project

Model subdirectory
---------------------
subdirectory named `slurm`. This is where the job output will be written.
If this directory does not exist, batch job will fail immediately with exit code 53.


`sra_ids.cfg` - file containing the list of SRAs be downloaded.

[to be developed]`sample_names.cfg` file containing the sample names to be integrate
will be used to name the individual .sam and .bam files
Prior to developing this feature, the sra_id is used

Analysis subdirectory
---------------------
Create subdirectory in the project with a name for the analysis to perform

Add a `models.cfg` with the names of the models to include in the analysis. 
This can be a subset of the models in the project.

Example directory structure
----------------------------
Here is an example project structure part-way through running the scripts.

        .
        ├── ana_1
        ├── ATAC
        │   ├── ATAC.bam
        │   ├── ATAC.bam.bai
        │   ├── ATAC.bw
        │   ├── ATAC.rpkm.bw
        │   ├── Lanceotron_output
        │   │   ├── ATAC_L-tron.bed
        │   │   ├── ATAC_L-tron_toppeaks.bed
        │   │   ├── ATA_L-tron_toppeaks.bed
        │   │   └── _L-tron_toppeaks.bed
        │   ├── slurm
        │   ├── SRR15970891_ter119_APH_spleen_rep1_R1.fastq.gz
        :│   └── SRR15970891_ter119_APH_spleen_rep1_R2.fastq.gz
        ├── enhancer_mcs_mm10_conservation.npz
        ├── H3K4me1
        │   ├── slurm
        │   ├── sra_ids.cfg
        │   ├── SRR14573546.bai
        │   ├── SRR14573546.bam
        │   ├── SRR14573546.bw
        │   ├── SRR14573546_R1.fastq.gz
        │   ├── SRR14573546_R2.fastq.gz
        │   ├── SRR14573547.bai
        │   ├── SRR14573547.bam
        │   ├── SRR14573547.bw
        │   ├── SRR14573547_R1.fastq.gz
        │   ├── SRR14573547_R2.fastq.gz
        │   ├── SRR14573548.bai
        │   ├── SRR14573548.bam
        │   ├── SRR14573548.bw
        │   ├── SRR14573548_R1.fastq.gz
        │   └── SRR14573548_R2.fastq.gz
        ├── H3K4me3
        │   ├── slurm
        │   ├── sra_ids.cfg
        │   ├── SRR5453541.bai
        │   ├── SRR5453541.bam
        │   ├── SRR5453541.bw
        │   ├── SRR5453541_R1.fastq.gz
        │   ├── SRR5453541_R2.fastq.gz
        │   ├── SRR5453542.bai
        │   ├── SRR5453542.bam
        │   ├── SRR5453542_R1.fastq.gz
        │   └── SRR5453542_R2.fastq.gz
        ├── Med1
        │   ├── Med1.bam
        │   ├── Med1.bam.bai
        │   ├── Med1.bw
        │   ├── R1.fastq.gz
        │   ├── R2.fastq.gz
        │   ├── slurm
        │   ├── sra_ids.cfg
        │   ├── SRR3199708.bai
        │   ├── SRR3199708.bam
        │   ├── SRR3199708.bw
        │   ├── SRR3199708_R1.fastq.gz
        │   └── SRR3199708_R2.fastq.gz
        ├── mini-mouse
        │   ├── ATAC
        │   │   ├── ATAC.bam
        │   │   └── slurm
        │   ├── H3K4me1
        │   │   ├── H3K4me1.bam
        │   │   └── slurm
        │   ├── H3K4me3
        │   │   ├── H3K4me3.bam
        │   │   └── slurm
        │   ├── Med1
        │   │   ├── Med1.bam
        │   │   └── slurm
        │   └── models.cfg
        ├── models.cfg

Run the scripts
===============

To run the upstream pipeline scripts (which contain Q, B, P in their name) `cd` to the subdirectory
 of the model. Then type the script. For example

        sbatch $SCRIPTS/Step1_B_fastq2bw.sh
If you have added $SCRIPTS to your $PATH, then just

        sbatch Step1_B_fastqbw.sh

Slurm logs will be written to the subfolder `/slurm` and the output file created in the current directory or subdirectories. You can use the utilities `sjobs` to view the jobs, `stree` to view data in the directory tree and `stidy-up` to tidy away slurm log files to the slurm subdirectory if you have from running other scripts.

The slurm log contains information about the precise version of the script which is run. This information can be used to inpsect and  a copy of that version of the script. Note that this relies on having committed any changes to git.


Utilities
=========
Slurm utilities
---------------
`sjobs` - list your jobs from today, includes exit codes

`stidy-up` - move things that look like slurm job output into subfolder named slurm. Create slurm folder if it does not exist.

`stree` - show tree of directories and files excluding things that look like slurm output

Other utitilities
-----------------

`make-mini-project.sh` - create a subproject with tiny subset of the data so that you can test out analysis and scripts. Currently, this is set up to select .bam data for the mouse alpha globin TAD.

`zip-fastq.sh` - gz zip any .fastq files in current directory.