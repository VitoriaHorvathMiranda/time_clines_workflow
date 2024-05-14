#!/bin/bash

# Check if at least two arguments were provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <sufix> <MD5_txt>"
    exit 1
fi

sufix="$1"
output="$2"

find . -type f -name "*$sufix" | while read -r file
do
    md5sum $file >> $output
    
done