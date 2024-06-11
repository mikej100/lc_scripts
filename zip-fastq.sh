#!/bin/bash
#
# Utility to zip fastq files.
for file in *.fastq ; do
    echo "G-zipping $file" ;
    gzip $file ;
done
