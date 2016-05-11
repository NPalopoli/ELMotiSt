#!/bin/bash

###############################################################################
#
# Generate parsed table of ELM structural information
# File name:ELMotiSt_PDBssSEQATOM.sh
# Author: Nicolas Palopoli
# Date created: 2015/10/05
# Date last modified: 2015/11/13
# Bash version: 4.3.11(1)-release 
#
###############################################################################
#
# RUN: ./ELMotiSt_PDBssSEQATOM.sh >ELMotiSt_PDBssSEQATOM.out 2>&1 &
#
###############################################################################



### FUNCTIONS ###

# Test file/dir exists
function testexist {
  if [ ! -f $1 ] 
  then
    if [ ! -e $1 ]
    then
      echo "ERROR. Missing $1. Exit."
      exit
    fi    
  fi
}

# Remove file if exists
function testrm {
  if [ -f $1 ]
  then
    rm "$1"
  fi
}

# Print help
function usage {
  description="$(basename "$0")\nProgram to generate parsed table of ELM structural information.\nArguments:\n"
  arguments="-h|--help\tShow this help\n
             -e|--inelm\tPath to elm_instances.[date].tsv file\n
             -j|--dirjpred\tPath to directory with JPred output files\n
             -s|--insifts\tPath to parsed SIFTS file\n"
  echo -e -n $description $arguments
}

### SETUP ###

# Set default paths
dirjpred='/home/npalopoli/20150924_ELM-Struct/JPred/elm_instances.fasta_ALL'
#inelm='/home/npalopoli/SLiMBench/ELMmap/elm_instances.2015-08-27.tsv'
inelm='/home/npalopoli/SLiMBench/ELMmap/elm_instances.2015-08-27_addPDB.tsv'
#insifts='/home/npalopoli/SIFTSParse/parseSIFTS_withELMID.out'
insifts='/home/npalopoli/SIFTSParse/parseSIFTS_PDBssSEQATOM_withELMID.out'

# Parse arguments to override defaults
while [[ $# > 0 ]]  # number of positional parameters
do
  argument="$1"
  case $argument in
    -e|--inelm)
      inelm="$2"
      testexist $inelm
      shift  # past argument
      ;;
    -j|--dirjpred)
      dirjpred="$2"
      testexist $dirjpred
      shift
      ;;
    -s|--insifts)
      insifts="$2"
      testexist $insifts
      shift
      ;;
    -h|--help)
      usage
      exit
      ;;
    *)
      shift  # unknown option
      ;;
  esac
  shift # past argument or value
done

# Check required files
testexist ELMotiSt_PDBssSEQATOM.py

# Write list of files to process
if [ ! -f ./files.lst ]
then
#  ls "$dirjpred"/*.fasta | xargs -n 1 basename | cut -d'.' -f 1 >files.lst
# Sample subsets for testing
# ls /home/npalopoli/20150924_ELM-Struct/JPred/elm_instances.fasta_dir_output/sp-P43489-TNR4_HUMAN.jnet | xargs -n 1 basename | cut -d'.' -f 1 >files.lst
# ls /home/npalopoli/20150924_ELM-Struct/JPred/elm_instances.fasta_dir_output/sp-P03070-LT_SV40.jnet | xargs -n 1 basename | cut -d'.' -f 1 >files.lst
# ls $dir/*DROME.fasta | xargs -n 1 basename | cut -d'.' -f 1 >files.lst
# ls /home/npalopoli/20150924_ELM-Struct/JPred/elm_instances.fasta_dir_output/sp-P53141-MLC1_YEAST.fasta | xargs -n 1 basename | cut -d'.' -f 1 >files.lst
#  ls /home/npalopoli/20150924_ELM-Struct/JPred/elm_instances.fasta_dir_output/sp-P39748-FEN1_HUMAN.fasta | xargs -n 1 basename | cut -d'.' -f 1 >files.lst
  ls /home/npalopoli/20150924_ELM-Struct/JPred/elm_instances.fasta_dir_output/sp-P38936-CDN1A_HUMAN.fasta | xargs -n 1 basename | cut -d'.' -f 1 >files.lst
fi

# Check/create output directory
if [ ! -d ./ELMotiStELM_PDBssSEQATOM_output ]
then
  mkdir ./ELMotiStELM_PDBssSEQATOM_output 
fi

# Remove list of result files if existing
testrm files_tsv.lst

### START ###

# Run ELMotiSt_PDBssSEQATOM.py for each in files.lst
while read line
do
  accession=`head -1 "$dirjpred"/"$line".fasta | cut -d'|' -f 2`
  echo -e "RUN\t$line\t$accession"
  ./ELMotiSt_PDBssSEQATOM.py "$dirjpred"/"$line".fasta "$dirjpred"/"$line".jnet "$inelm" "$insifts" "$accession" >ELMotiStELM_PDBssSEQATOM_output/"$line".tsv
  echo "$line".tsv >>files_tsv.lst
done <files.lst

# Collect results in single output
sample=`head -1 files.lst`
head -1 ELMotiStELM_PDBssSEQATOM_output/"$sample".tsv>ELMotiStELM_PDBssSEQATOM_all.tsv
while read line
do
  tail -n +2 ELMotiStELM_PDBssSEQATOM_output/"$line" >>ELMotiStELM_PDBssSEQATOM_all.tsv
done<files_tsv.lst

# Cleanup
testrm files.lst
testrm files_tsv.lst
