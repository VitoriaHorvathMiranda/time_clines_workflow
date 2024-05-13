#!/bin/bash

find . -type f -name "*.fastq.gz" | while read -r file
do
    md5sum $file >> MD5sum.txt
    
done