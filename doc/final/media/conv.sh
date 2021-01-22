#!/bin/bash

for eps_file in "$1"/*.eps; do
  epspdf "$eps_file"
done

for pdf_file in "$1"/*.pdf; do
  filename="${pdf_file%.*}"
  pdf2svg "$pdf_file" "$filename.svg" 
done
