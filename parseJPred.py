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

# Read input files
try:
  infasta = open(sys.argv[1])
  injnet = open(sys.argv[2])
except IndexError:
  print("Input file(s) not specified. Format: ./parseJPred.py <in.fasta> <in.jnet>")
  exit()
except IOError:
  print("Input file(s) not found. Format: ./parseJPred.py <in.fasta> <in.jnet>")
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
        seqs[name]=''
      else :  # sequence, not header, possibly multi-line
        seqs[name] = seqs[name] + line
  seqs[name] = list(seqs[name])
  seqs['position'] = range(1,len(seqs[name])+1)
  seqs['name'] = [name] * len(seqs[name])
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

# Make dict of input files
seq = readFasta(infasta)
infasta.close()
predictions = readJPred(injnet)
injnet.close()
results = seq.copy()
results.update(predictions)

# Print table with results
for row in zip(*([key] + value for key, value in results.items())):
  print '\t'.join(map(str, row))


'''
print seq
for name in seq:
  print name,len(seq[name]),list(seq[name])

for program in entries:
  print program,len(entries[program]),entries[program]
'''

'''
jnetpred
JNETCONF
JNETSOL25
JNETSOL5
JNETSOL0
JNETHMM
JNETPSSM
JNETJURY
JNETPROPE
JNETPROPH
JNETPROPC
'''

'''
# Split sequences by length
seqsSplit = {}  # dict for split seqs
for i in seqs.keys():
  chunk_size = 800  # max sequence length
  if len(seqs[i]) >= chunk_size:
    seqsSplit[i] = [seqs[i][pos:pos+chunk_size] for pos in range(0, len(seqs[i]), chunk_size)]
  else:
    seqsSplit[i] = list(seqs[i].split())

# Save fragments as separate files
for entry in seqsSplit.keys():
  count = 1
  for i in seqsSplit[entry]:
    outfile = entry+'_'+str(count)+'.fasta'  # output filename: 'fastaID_fragmentNumber.fasta'
    with open(outfile, 'a') as f:
      f.write('>{0}\n{1}\n'.format(entry,i))
    count += 1
'''
