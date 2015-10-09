#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
  Parse *.jnet files from JPred output
  File name: parseJPred.py
  Author: Nicolas Palopoli
  Date created: 2015/10/05
  Date last modified: 2015/10/05
  Python Version: 2.7
'''

import sys
from collections import OrderedDict
import csv

# Read input files
# TO-DO: replace argument parsing with argparse
try:
  infasta = open(sys.argv[1])
  injnet = open(sys.argv[2])
  inELMinstances = open(sys.argv[3])
  primaryAcc = sys.argv[4]
except IndexError:
  print("Input file(s) not specified. Format: ./parseJPred.py <in.fasta> <in.jnet> <elm_instances[.date].tsv> <primaryAcc>")
  exit()
except IOError:
  print("Input file(s) not found. Format: ./parseJPred.py <in.fasta> <in.jnet> <elm_instances[.date].tsv> <primaryAcc>")
  exit()

def readFasta(infasta):
  '''Store fasta sequences from file.'''
  seqs = OrderedDict()
#  seqs={}  # dict for raw seqs
  readFirstSeq = False
  for line in infasta:
    if line.strip() or line not in ['\n', '\r\n']:  # avoid empty or only whitespace lines
      line=line.rstrip()  # discard newline at the end (if any)
      if line[0]=='>':  # or line.startswith('>'); distinguish header
        if readFirstSeq:  # exit if more than 2 sequences
          break
        readFirstSeq = True
        words=line.split()
        name=words[0][1:]
        seqs['res']=''
      else :  # sequence, not header, possibly multi-line
        seqs['res'] = seqs['res'] + line
  seqs['res'] = list(seqs['res'])
  seqs['position'] = range(1,len(seqs['res'])+1)
  seqs['name'] = [name] * len(seqs['res'])
  return seqs

def readJPred(injnet):
  '''Store JPred predictions by program from file.'''
  entries={}  # dict for raw entries
  entries = OrderedDict()
  for line in injnet:
    if line.strip() or line not in ['\n', '\r\n']:  # avoid empty or only whitespace lines
      line = line.rstrip()  # discard newline at the end (if any)
      words = line.split(':')
      program = words[0]  # programs are keys
      entries[program] = ''
      values = words[1].split(',')  # predictions are values
      entries[program] = values[0:-1]
  tempValue = []
  for value in entries['JNETJURY']:  # replace empty values with dashes
    if value == '*':
      tempValue.append('*')
    else:
      tempValue.append('-')
  entries['JNETJURY'] = tempValue
  return entries

def readELMinstances(infile):
  '''Store ELM instances information as list of dicts'''
  elm = csv.DictReader(filter(lambda row: row[0]!='#', infile),delimiter='\t', quotechar='"')
  return elm

def mapELMpositions(parsedELM,primaryAcc):
  '''Make dict with [ELMAccession:(Start,End)]'''
  ELMpos = {}
  for row in parsedELM:
    if primaryAcc == row['Primary_Acc']:
      ELMpos[row['Accession']] = [row['Start'],row['End']]
#      ELMpos[row['Accession']] = (row['Start'],row['End'])
  return ELMpos

def placeELM(seq,ELMpos):
  '''Map ELM to fasta sequence'''
  seq['ELMpos'] = list('-' * len(seq['res']))
  seq['ELMacc'] = list('-' * len(seq['res']))
  for accession, limits in ELMpos.iteritems():
    for pos in range(int(limits[0])-1,int(limits[1])):
      seq['ELMpos'][pos] = seq['res'][pos]
      if '-' in seq['ELMacc'][pos]: 
        seq['ELMacc'][pos] = accession
      else:
        seq['ELMacc'][pos] = seq['ELMacc'][pos] + accession
  return seq

def printTable(results):
  '''Print parsing results as table'''
  for row in zip(*([key] + value for key, value in results.items())):
    print '\t'.join(map(str, row))
  return None

# Make dict of input files
seq = readFasta(infasta)
infasta.close()

parsedELM = readELMinstances(inELMinstances)
inELMinstances.close()
ELMpos = mapELMpositions(parsedELM,primaryAcc)
seq = placeELM(seq,ELMpos)

predictions = readJPred(injnet)
injnet.close()

results = seq.copy()  # merge tables
results.update(predictions)

printTable(results)

