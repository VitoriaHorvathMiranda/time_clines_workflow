#!/bin/bash

# Store arguments into an array
args=("$@")

echo "zip files: '${args[@]}'"

# Iterate over the array elements
for arg in "${args[@]}"; do
  echo "Argument: $arg"
  files=$(less $arg | awk '{print $8}' | grep -E "(1997|2009|2010|2017|2022)" | sed ':a;N;$!ba;s/\n/ /g')
  echo "these are teh files to unzip: $files"
  unzip -j $arg $files -d ../resources/climate_data

done