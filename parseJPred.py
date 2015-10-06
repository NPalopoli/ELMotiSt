#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
  Parse JPred output
  File name: parseJPred.py
  Author: Nicolas Palopoli
  Date created: 2015/10/05
  Date last modified: 2015/10/05
  Python Version: 2.7
'''

import sys

# Read input file
try:
  infile = open(sys.argv[1])
except IndexError:
  print("Input file name not specified. Exit.")
  exit()
except IOError:
  print("Input file {} does not exist. Exit.".format(sys.argv[1]))
  exit()

def readJPred(infile):
  '''Store JPred predictions by program from file.'''
  entries={}  # dict for raw entries
  for line in infile:
    if line.strip() or line not in ['\n', '\r\n']:  # avoid empty or only whitespace lines
      line = line.rstrip()  # discard newline at the end (if any)
      words = line.split(':')
      program = words[0]  # programs are keys
      entries[program] = ''
      values = words[1].split(',')  # predictions are values
      entries[program] = values[0:-1]
  return entries

entries = readJPred(infile)

infile.close()

for program in entries:
  print program,len(entries[program]),entries[program]

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
