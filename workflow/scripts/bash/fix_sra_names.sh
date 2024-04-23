#!/bin/bash
# This script takes a directory with sra.fastq.gz files and sra info csv file (/time_clines_workflow/resources/SRA_samples_info.csv)
# and changes the name of the downloaded samples to a more intuitive one

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <whole path input_file> <directory with sra.fastq.gz>" #precisa ser o 'whole path'
    exit 1
fi

input_file="$1"

# Directory where your sequence files are located
sequence_dir="$2"

# Change to the specified directory
cd "$sequence_dir" || { echo "Directory not found: $sequence_dir"; exit 1; }
echo "Current directory: $PWD" 


# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: File '$input_file' not found."
    exit 1
fi

# Read each 11th column of the input file and store it in the sra array
while read line; do
    sample+=($line)
done < <(tail -n +3 "$input_file" | cut -d ',' -f 11 )


# Read each 1st column of the input file and store it in the sample array
while read line; do
    sra+=($line)
done < <(tail -n +3 "$input_file" | cut -d ',' -f 1 )


echo "sra codes: '${sra[@]}'"
echo "sample names: '${sample[@]}'"

# Iterate over indices of arrays
for ((i = 0; i < ${#sra[@]}; i++)); do
    # Construct file paths
    old_path="${sra[i]}.fastq.gz"
    new_path="${sample[i]}.fastq.gz"
    
    # Check if the file exists before renaming
    if [ -f "$old_path" ]; then
        #mv "$old_path" "$new_path"
        echo "Renamed: $old_path -> $new_path"
    else
        echo "Warning: File not found - $old_path"
    fi
done