Bash utilities for processing chromatin assay results such as ChIP-seq data.
Written by Mike Jennings, based on scripts from \<to be named>.

Set up
======

~/.bashrc
----------
    SCRIPTS=<path to directory with these scripts>
For ease of use, you may want put this directory in your PATH.
Then you can type the name of the script and it will run, regardless of your
current work directory. 
$ PATH=$SCRIPTS:$PATH

Project directory
-----------------
Subdirectory for each track of data (model)
The name of the subdirectory will be used for track name of aggregated data 
such as peak call

Model subdirectory
---------------------
subdirectory named slurm. This is where the job output will be written.
If this directory does not exist, batch job will fail immediately with exit code 53

sra_ids.cfg - file containing the list of SRAs be downloaded
[to be developed]sample_names.cfg file containing the sample names to be integrate
will be used to name the individual .sam and .bam files
Prior to developing this feature, the sra_id is used 

Utilities
=========
sjobs - list your jobs from today, includes exit codes

stree - show tree of directories and files excluding things that look like slurm output

make-mini-project - create a subproject with tiny subset of the data so that you can test out analysis and scripts. Currently, this is set up to select .bam data for the mouse alpha globin TAD \[REQUIRES FIXING FOR MM39]