#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
  Generate parsed table of ELM structural information with 2nd structure and SEQATOM/SEQRES mapping
  File name: ELMotiSt_PDBssSEQATOM.py
  Author: Nicolas Palopoli
  Date created: 2015/10/05
  Date last modified: 2016/03/03
  Python Version: 2.7
'''

import sys
from collections import OrderedDict
import csv
from Bio import SeqIO

# Read input files
try:
  infasta = open(sys.argv[1])
  injnet = open(sys.argv[2])
  inELMinstances = open(sys.argv[3])
  inSIFTSparse = sys.argv[4]
  primaryAcc = sys.argv[5]
except IndexError:
  print("Input file(s) not specified. Format: ./ELMotiSt_PDBssSEQATOM.py <in.fasta> <in.jnet> <elm_instances[.date].tsv> <parseSIFTS.out> <primaryAcc>")
  exit()
except IOError:
  print("Input file(s) not found. Format: ./ELMotiSt_PDBssSEQATOM.py <in.fasta> <in.jnet> <elm_instances[.date].tsv> <parseSIFTS.out> <primaryAcc>")
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
#  entries={}  # dict for raw entries
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
  '''Make dict with [ELMAccession:[Start,End,ELMType,ELMIdentifier]]'''
  ELMpos = {}
  for row in parsedELM:
    if primaryAcc == row['Primary_Acc']:  # Uniprot accession match
#      ELMpos[row['Accession']] = [row['Start'],row['End']]
      ELMpos[row['Accession']] = [row['Start'],row['End'],row['ELMType'],row['ELMIdentifier']]
#      ELMpos[row['Accession']] = (row['Start'],row['End'])
  return ELMpos

def placeELM(seq,ELMpos):
  '''Map ELM to fasta sequence'''
  seq['ELMpos'] = list('-' * len(seq['res']))
  seq['ELMacc'] = list('-' * len(seq['res']))
  seq['ELMType'] = list('-' * len(seq['res']))
  seq['ELMIdentifier'] = list('-' * len(seq['res']))
  for accession, vals in ELMpos.iteritems():
    for pos in range(int(vals[0])-1,int(vals[1])):  # iterate over Primary_Acc sequence from 'Start' to 'End' positions of ELM
      seq['ELMpos'][pos] = seq['res'][pos]  # if aa at position matches:
      if '-' in seq['ELMacc'][pos]:  # create list of ELMs at position
        seq['ELMacc'][pos] = accession
        seq['ELMType'][pos] = vals[2]
        seq['ELMIdentifier'][pos] = vals[3]
      else:  # append to list of ELMs at position if already exists
        seq['ELMacc'][pos] = seq['ELMacc'][pos] + ',' + accession
        seq['ELMType'][pos] = seq['ELMType'][pos] + ',' + vals[2]
        seq['ELMIdentifier'][pos] = seq['ELMIdentifier'][pos] + ',' + vals[3]
  return seq

'''
inSIFTSparse

>ELMI000068:P38936:4RJF:F:sequence
------------------------------------------------------------------------------------------------------------------------------------------GRKRRQTSMTDFFHSKRRLIFS
>ELMI000068:P38936:4RJF:F:secstr
-----------------------------------------------------------------------------------------------------------------------------------------------B--GGGTS-EEEEE---
>ELMI000068:P38936:4RJF:F:disorder
------------------------------------------------------------------------------------------------------------------------------------------XXXX------------------
>ELMI000068:P38936:4RJF:F:SEQRES
------------------------------------------------------------------------------------------------------------------------------------------GRKRRQTSMTDFFHSKRRLIFS
>ELMI000068:P38936:4RJF:F:SEQATOM
----------------------------------------------------------------------------------------------------------------------------------------------RQTSMTDFFHSKRRLIFS

>ELMI000068:P38936:1AXC:F:sequence
------------------------------------------------------------------------------------------------------------------------------------------GRKRRQTSMTDFYHSKRRLIFS
>ELMI000068:P38936:1AXC:F:secstr
-----------------------------------------------------------------------------------------------------------------------------------------------B--GGGTSEEEEEEE--
>ELMI000068:P38936:1AXC:F:disorder
----------------------------------------------------------------------------------------------------------------------------------------------------------------
>ELMI000068:P38936:1AXC:F:SEQRES
------------------------------------------------------------------------------------------------------------------------------------------GRKRRQTSMTDFYHSKRRLIFS
>ELMI000068:P38936:1AXC:F:SEQATOM
----------------------------------------------------------------------------------------------------------------------------------------------RQTSMTDFYHSKRRLIFS
'''

def readSIFTSparse(inSIFTSparse,ELMpos,seq):
  '''Read 2nd struct from PDB following SIFTS parsing'''
#  PDBss = OrderedDict()
  PDBss = {}
  PDBss['PDBid'] = list('.' * len(seq['res']))
  PDBss['PDBchain'] = list('.' * len(seq['res']))
  PDBss['PDBseq'] = list('.' * len(seq['res']))
  PDBss['PDBss'] = list('.' * len(seq['res']))
  PDBss['PDBseqres'] = list('.' * len(seq['res']))
  PDBss['PDBseqatom'] = list('.' * len(seq['res']))
  PDBss['PDBdis'] = list('.' * len(seq['res']))
  fastaseqs = SeqIO.parse(open(inSIFTSparse),'fasta')
  for fasta in fastaseqs:
#chk    if fasta.id[0:10] == ELMpos.keys()[0]:
    for keys in ELMpos:
      if fasta.id[0:10] != keys:
        continue
#      if 'sequence' in fasta.id:
#        for pos in range(0,len(seq['res'])):
#          if seq['ELMacc'][pos] == '-':
#            continue
#          if PDBss['PDBid'][pos] != '-':
#            PDBss['PDBid'][pos] = '%s,%s' % (PDBss['PDBid'][pos],fasta.id[18:22])
#            PDBss['PDBchain'][pos] = '%s,%s' % (PDBss['PDBchain'][pos],fasta.id[23])
#          else:
#            PDBss['PDBid'][pos] = fasta.id[18:22]
#            PDBss['PDBchain'][pos] = fasta.id[23]
#        PDBss['PDBid'][pos] = fasta.id[18:22]
#        PDBss['PDBchain'][pos] = fasta.id[23]
      if 'sequence' in fasta.id:
        for pos in range(0,len(fasta.seq)):
          if PDBss['PDBseq'][pos] != '.':
            PDBss['PDBseq'][pos] = '%s,%s' % (PDBss['PDBseq'][pos],fasta.seq[pos])
          elif fasta.seq[pos] != '-':
            PDBss['PDBseq'][pos] = fasta.seq[pos]
#chk      if PDBss['PDBseq'][pos] != '.':
          if fasta.seq[pos] != '-':
            if PDBss['PDBid'][pos] != '.':  # and fasta.id[18:22] not in PDBss['PDBid'][pos]:
              PDBss['PDBid'][pos] = '%s,%s' % (PDBss['PDBid'][pos],fasta.id[18:22])
              PDBss['PDBchain'][pos] = '%s,%s' % (PDBss['PDBchain'][pos],fasta.id[23])
#            elif fasta.seq[pos] != '-':
            else:
              PDBss['PDBid'][pos] = fasta.id[18:22]
              PDBss['PDBchain'][pos] = fasta.id[23]
      elif 'SEQRES' in fasta.id:
        for pos in range(0,len(fasta.seq)):
          if PDBss['PDBseqres'][pos] != '.':
            PDBss['PDBseqres'][pos] = '%s,%s' % (PDBss['PDBseqres'][pos],fasta.seq[pos])
          elif PDBss['PDBseq'][pos] != '.':
#chk          elif fasta.seq[pos] != '-':
            PDBss['PDBseqres'][pos] = fasta.seq[pos]
      elif 'SEQATOM' in fasta.id:
        for pos in range(0,len(fasta.seq)):
          if PDBss['PDBseqatom'][pos] != '.':
            PDBss['PDBseqatom'][pos] = '%s,%s' % (PDBss['PDBseqatom'][pos],fasta.seq[pos])
          elif PDBss['PDBseq'][pos] != '.':
#chk          elif fasta.seq[pos] != '-':
            PDBss['PDBseqatom'][pos] = fasta.seq[pos]
      elif 'secstr' in fasta.id:
        for pos in range(0,len(fasta.seq)):
          if PDBss['PDBss'][pos] != '.':
            PDBss['PDBss'][pos] = '%s,%s' % (PDBss['PDBss'][pos],fasta.seq[pos])
          elif PDBss['PDBseq'][pos] != '.':
            PDBss['PDBss'][pos] = fasta.seq[pos]
      elif 'disorder' in fasta.id:
        for pos in range(0,len(fasta.seq)):
          if PDBss['PDBdis'][pos] != '.':
            PDBss['PDBdis'][pos] = '%s,%s' % (PDBss['PDBdis'][pos],fasta.seq[pos])
          elif PDBss['PDBseq'][pos] != '.':
            PDBss['PDBdis'][pos] = fasta.seq[pos]
      
#  return PDBss
  for pos in range(0,len(seq['res'])):  # restore temporary initial '.' as '-'
    if PDBss['PDBid'][pos] == '.':
      PDBss['PDBid'][pos] = '-'
    if PDBss['PDBchain'][pos] == '.':
      PDBss['PDBchain'][pos] = '-'
    if PDBss['PDBseq'][pos] == '.':
      PDBss['PDBseq'][pos] = '-'
    if PDBss['PDBss'][pos] == '.':
      PDBss['PDBss'][pos] = '-'
    if PDBss['PDBseqres'][pos] == '.':
      PDBss['PDBseqres'][pos] = '-'
    if PDBss['PDBseqatom'][pos] == '.':
      PDBss['PDBseqatom'][pos] = '-'
    if PDBss['PDBdis'][pos] == '.':
      PDBss['PDBdis'][pos] = '-'
  return PDBss

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

PDBss = readSIFTSparse(inSIFTSparse,ELMpos,seq)
#print type(PDBss)

results = seq.copy()  # merge tables
results.update(PDBss)
results.update(predictions)

printTable(results)

