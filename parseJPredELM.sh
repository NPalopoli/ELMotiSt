#!/bin/bash

dir='/home/npalopoli/20150924_ELM-Struct/JPred/elm_instances.fasta_dir_output'
ls $dir/*.fasta | xargs -n 1 basename | cut -d'.' -f 1 >files.lst

while read line
do
  accession=`head -1 "$dir"/"$line".fasta | cut -d'|' -f 2`
  ./parseJPredELM.py "$dir"/"$line".fasta "$dir"/"$line".jnet /home/npalopoli/SLiMBench/ELMmap/elm_instances.2015-08-27.tsv "$accession" >"$line".tsv
done <files.lst

rm files.lst
