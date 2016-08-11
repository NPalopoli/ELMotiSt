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
from Bio.Alphabet import IUPAC
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

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

def OLD_readSIFTSparse(inSIFTSparse,ELMpos,seq):
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
#  PDBss['PDBidchain'] = list('-')
# pdbidchain = list('-') * len(seq['res'])
  readflag = 0
  pdbseq = [[] for _ in seq['res']]
  pdbidchain = [[] for _ in seq['res']]
  pdbseqres = [[] for _ in seq['res']]
  pdbseqatom = [[] for _ in seq['res']]
  fastaseqs = SeqIO.parse(open(inSIFTSparse),'fasta',IUPAC.extended_protein)
  for fasta in fastaseqs:
#    if fasta.id[0:10] == ELMpos.keys()[0]:
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
#chk        for pos in range(0,len(seq['res'])):
          try:
            if fasta.seq[pos] != '-' and ((fasta.id[18:22],fasta.id[23]) not in pdbseq[pos]):
              pdbseq[pos].append((fasta.id[18:22],fasta.id[23]))
          except IndexError:
            print fasta.seq, seq['res']
            break
          else:
#            readflag = 1
#          if readflag == 1:
            if PDBss['PDBseq'][pos] != '.':
              PDBss['PDBseq'][pos] = '%s,%s' % (PDBss['PDBseq'][pos],fasta.seq[pos])
            elif fasta.seq[pos] != '-':
              PDBss['PDBseq'][pos] = fasta.seq[pos]
#      if PDBss['PDBseq'][pos] != '.':
#          if fasta.seq[pos] != '-':
          if fasta.seq[pos] != '-' and ((fasta.id[18:22],fasta.id[23]) not in pdbidchain[pos]):
            pdbidchain[pos].append((fasta.id[18:22],fasta.id[23]))
            readflag = 1
            if PDBss['PDBid'][pos] != '.':
              PDBss['PDBid'][pos] = '%s,%s' % (PDBss['PDBid'][pos],fasta.id[18:22])
              PDBss['PDBchain'][pos] = '%s,%s' % (PDBss['PDBchain'][pos],fasta.id[23])
#            elif fasta.seq[pos] != '-':
            else:
              PDBss['PDBid'][pos] = fasta.id[18:22]
              PDBss['PDBchain'][pos] = fasta.id[23]
          else:
            readflag = 0
      elif 'SEQRES' in fasta.id:
        for pos in range(0,len(fasta.seq)):
#chk          if PDBss['PDBseqres'][pos] != '.':
          if fasta.seq[pos] != '-' and ((fasta.id[18:22],fasta.id[23]) not in pdbseqres[pos]):
            pdbseqres[pos].append((fasta.id[18:22],fasta.id[23]))
            if PDBss['PDBseqres'][pos] != '.':
              PDBss['PDBseqres'][pos] = '%s,%s' % (PDBss['PDBseqres'][pos],fasta.seq[pos])
            elif PDBss['PDBseq'][pos] != '.':
  #          elif fasta.seq[pos] != '-':
              PDBss['PDBseqres'][pos] = fasta.seq[pos]
#          if PDBss['PDBseqres'][pos] != '.' and ((fasta.id[18:22],fasta.id[23]) not in pdbidchain[pos]):
#            PDBss['PDBseqres'][pos] = '%s,%s' % (PDBss['PDBseqres'][pos],fasta.seq[pos])
#          elif PDBss['PDBseq'][pos] != '.':
##          elif fasta.seq[pos] != '-':
#            PDBss['PDBseqres'][pos] = fasta.seq[pos]
      elif 'SEQATOM' in fasta.id:
        for pos in range(0,len(fasta.seq)):
          if readflag == 1:
            if PDBss['PDBseqatom'][pos] != '.':
              PDBss['PDBseqatom'][pos] = '%s,%s' % (PDBss['PDBseqatom'][pos],fasta.seq[pos])
            elif PDBss['PDBseq'][pos] != '.':
#            elif fasta.seq[pos] != '-':
              PDBss['PDBseqatom'][pos] = fasta.seq[pos]
      elif 'secstr' in fasta.id:
        for pos in range(0,len(fasta.seq)):
          if readflag == 1:
            if PDBss['PDBss'][pos] != '.':
              PDBss['PDBss'][pos] = '%s,%s' % (PDBss['PDBss'][pos],fasta.seq[pos])
            elif PDBss['PDBseq'][pos] != '.':
             PDBss['PDBss'][pos] = fasta.seq[pos]
      elif 'disorder' in fasta.id:
        for pos in range(0,len(fasta.seq)):
          if readflag == 1:
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


def testreadSIFTSparse(inSIFTSparse,ELMpos,seq):
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
#  PDBss['PDBidchain'] = list('-')
# pdbidchain = list('-') * len(seq['res'])
  readflag = 0
  pdbseq = [[] for _ in seq['res']]
  pdbidchain = [[] for _ in seq['res']]
  pdbseqres = [[] for _ in seq['res']]
  pdbseqatom = [[] for _ in seq['res']]
  fastaseqs = SeqIO.parse(open(inSIFTSparse),'fasta')
  names = ''
  for fasta in fastaseqs:
    minindex = 0
    maxindex = 0
#    if fasta.id[0:10] == ELMpos.keys()[0]:
    for keys in ELMpos:
      if fasta.id[0:10] != keys:
        break
      names = names + '\n' + fasta.id
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
#      continue
#      listseqres = str(fasta.seq) + '\n' + ''.join(seq['res'])
#      listseqres.extend(seq['res'])
#      break
      fastaseq = fasta.seq.strip("-")
      seqres = ''.join(seq['res']).strip("-")
      alignments = pairwise2.align.globalms(fastaseq,seqres,10,0,-10,-1)
#      alignments = pairwise2.align.globalms(str(fasta.seq),''.join(seq['res']),10,0,-10,-1)
#      alignments = pairwise2.align.globalxx(str(fasta.seq),''.join(seq['res']))
#      listfastaseq = alignments[0][0]
#      listfastaseq = list(alignments[0][0])
#      listseqres = alignments[0][1]
#      listseqres = list(alignments[0][1])
#      break
#      listindex = []
#      for index,pos in enumerate(listseqres):
#        if listseqres[index] != '-':
#          listindex.extend(index)
#      minindex = min(listindex)
#      maxindex = max(listindex)
#      return minindex, maxindex
  listfastaseq = list(alignments[0][0])
  listseqres = list(alignments[0][1])
  listindex = []
  for index,pos in enumerate(listseqres):  # o es listseqres?
    if listfastaseq[index] != '-':
#      listindex.extend(str(index))
      listindex.append(index)
#  minindex = min(map(int,listindex))
#  maxindex = max(map(int,listindex))
  minindex = min(listindex)
  maxindex = max(listindex)
  return 'FS\t' + fastaseq + '\nSR\t' + seqres + '\nA0\t' + alignments[0][0] + '\nA1\t' + alignments[0][1]  + '\nMI\t' + str(minindex) + '\nMX\t' + str(maxindex)
#  return 'FS\t' + fastaseq + '\nSR\t' + seqres + '\nA0\t' + alignments[0][0] + '\nA1\t' + alignments[0][1] + '\nLI\t' + ''.join(listindex) + '\nMI\t' + str(minindex) + '\nMX\t' + str(maxindex)
#  return 'FS\t' + fastaseq + '\nSR\t' + seqres + '\nA0\t' + alignments[0][0] + '\nA1\t' + alignments[0][1] + '\nLI\t' + ''.join(listindex) + '\nMI\t' + minindex + '\nMX\t' + maxindex
#  return 'FS\t' + fastaseq + '\nSR\t' + seqres + '\nA0\t' + alignments[0][0] + '\nA1\t' + alignments[0][1] + '\nLI\t' + listindex + '\nMI\t' + minindex + '\nMX\t' + maxindex
#  return str(names)
#  return fastaseq + '\n' + seqres + '\n' + ''.join(listfastaseq).replace('-','') + '\n' + ''.join(listseqres)


def getaliindex(fastaseq,fastaid,seqres):
  '''Align Uniprot sequence to SIFTS parsed sequence'''
  minindex = 0
  maxindex = 0
#  if len(fastaseq) == len(seqres):
#    
#  fastaseq = fastaseq.ungap('-')  # SIFTS sequence  # gap strip no esta funcionando
  fastaseq2 = str(fastaseq).strip('-')
  seqres2 = ''.join(seqres).strip('-')  # Uniprot sequence
  if len(fastaseq) == 0 or len(fastaseq2) == 0:  # cases when no SIFTS sequence: XXXX | ____ ; XXXX | ---- 
#    return (0, len(seqres), str(fastaseq), ''.join(seqres), str(fastaseq2), ''.join(seqres2))  # minindex, maxindex, fastaseq, seqres, ali1, ali2
    return (0, len(seqres) - 1, str(fastaseq), ''.join(seqres), '-' * len(seqres), ''.join(seqres))  # minindex, maxindex, fastaseq, seqres, string of '-', seqres
  elif seqres2 in fastaseq2:  # case: SIFTS includes Uniprot plus other residues
    fastaseqindex = str(fastaseq).index(fastaseq2)  # index of first position not gap in fastaseq
    fastaseq2index = fastaseq2.index(seqres2)  # index of first position of seqres2 (ungapped seqres) in fastaseq2 (ungapped fastaseq)
#    minindex = fastaseqindex + fastaseq2index  # minindex starts at position of fastaseq (gapped) where match with seqres (ungapped)
    return (0, len(seqres) - 1, str(fastaseq), ''.join(seqres), fastaseq[(fastaseqindex + fastaseq2index):(fastaseqindex + fastaseq2index + len(seqres))], ''.join(seqres))  # minindex, maxindex, fastaseq, seqres, sequence of fastaseq matching seqres2 (ungapped seqres), seqres
  elif fastaseq2 in seqres:  # case: Uniprot includes SIFTS plus other residues
    fastaseqindex = str(fastaseq).index(fastaseq2)  # index of first position not gap in fastaseq
    seqres2index = seqres2.index(fastaseq2)  # index of first position of fastaseq2 (ungapped fastaseq) in seqres2 (ungapped seqres)
    return (0, len(seqres) - 1, str(fastaseq), ''.join(seqres), ''.join('-' * seqres2index, str(fastaseq), '-' * (len(seqres2) - seqres2index - len(fastaseq))), ''.join(seqres))  # minindex, maxindex, fastaseq, seqres, sequence of fastaseq matching seqres2 completed with gaps, seqres
  else:  # cases: sequences difer in length or aa sequence and none is empty
    fastaseqindex = str(fastaseq).index(fastaseq2)  # index of first position not gap in fastaseq
    seqresindex = ''.join(seqres).index(seqres2)  # index of first position not gap in seqres
#  try:
    matrix = matlist.ident  # use identity matrix: http://biopython.org/DIST/docs/api/Bio.SubsMat.MatrixInfo-module.html#ident
    matrix.update(((a,'X'),6) for (a,b) in matrix.keys() if (a,'X') not in matrix)  # add tuples with 'X' in matrix
    matrix.update(((b,a),val) for (a,b),val in matrix.items())  # make matrix squared
    alignments = pairwise2.align.globalds(fastaseq2,seqres2,matrix,-10,-1, one_alignment_only=1)
#    print 'ALIGNED'
#    print fastaid
#    print 'FASTASEQ2:\n' + fastaseq2
#    print 'SEQRES2:\n' + seqres2
#    print type(alignments)
#    print alignments
#    print type(alignments[0])
#    print alignments[0]
#    print 'END'
#      listfastaseq = list(align1)
#      listseqres = list(align2)
  try:
    print 'TRYING'
  except:
    print 'ERROR_START:' + fastaid
    print 'ERROR_FASTASEQ2:\n' + fastaseq2
#    print 'ERROR_START:' + fastaid
    print 'ERROR_SEQRES2:\n' + seqres2
#    print 'ALIGNMENTS:\n' + alignments
    print 'ERROR_END\n'
    return (0, len(seqres), str(fastaseq), ''.join(seqres), str(fastaseq), ''.join(seqres))
  else:
    print 'ALIGNED_START:' + fastaid
    print 'ALIGNED_FASTASEQ2:\n' + fastaseq2
    print 'ALIGNED_SEQRES2:\n' + seqres2
    print type(alignments)
    print alignments
    print type(alignments[0])
    print alignments[0]
    print 'ALIGNED_END\n'
#    listfastaseq = list(alignments[0][0])
#    listseqres = list(alignments[0][1])
#    topaln = alignments[0]
#    listfastaseq, listseqres, score, begin, end = topaln
#    listfastaseq = list(listfastaseq)
#    listseqres = list(listseqres)
    alifastaseq, aliseqres, score, begin, end = alignments[0]
    listalifastaseq = list(alifastaseq)
    listaliseqres = list(aliseqres)
    listindex = []
    for index,pos in enumerate(listalifastaseq):  # o es listseqres?
      if listaliseqres[index] != '-':
        listindex.append(index)
#    minindex = min(listindex)
#    maxindex = max(listindex)
    minindex = min(listindex) + seqresindex
    maxindex = max(listindex) + seqresindex
    return (minindex, maxindex+1, fastaseq, seqres, alifastaseq, aliseqres)  # maxindex+1 deberia ser correcto en todos los casos
#  except:
#    print 'ERROR'
#    return ('0', '1', str(fastaseq), ''.join(seqres), str(fastaseq), ''.join(seqres))
#    return ('0', '1', str(fastaseq).strip('-'), ''.join(seqres).strip('-'), alignments[0][0], alignments[0][1])
#  return 'FS\t' + fastaseq + '\nSR\t' + seqres + '\nA0\t' + alignments[0][0] + '\nA1\t' + alignments[0][1]  + '\nMI\t' + str(minindex) + '\nMX\t' + str(maxindex)

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
#  PDBss['PDBidchain'] = list('-')
# pdbidchain = list('-') * len(seq['res'])
  readflag = 0
  pdbseq = [[] for _ in seq['res']]
  pdbidchain = [[] for _ in seq['res']]
  pdbseqres = [[] for _ in seq['res']]
  pdbseqatom = [[] for _ in seq['res']]
  fastaseqs = SeqIO.parse(open(inSIFTSparse),'fasta')
  for fasta in fastaseqs:
#    if fasta.id[0:10] == ELMpos.keys()[0]:
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

#      alignments = pairwise2.align.globalms(fasta.seq,seq['res'],100,-10,-100,-100)
#      listfastaseq = list(alignments[0][0])
#      listseqres = list(alignments[0][1])
#      listindex = []
#      for index,pos in enumerate(listseqres):
#        if listseqres[index] != '-':
#          listindex.extend(index)
#      minindex = min(listindex)
#      maxindex = max(listindex)

#      minindex, maxindex, fastaseq, seqres, ali0, ali1  = getaliindex(fasta.seq,fasta.id,seq['res'])
      
#      if 'sequence' not in fasta.id:
#        break
#      elif 'sequence' in fasta.id:

#      if len(fasta.seq) == len(seq['res']):
#        minindex = 0
#        maxindex = len(seq['res']) + 1
#      else:
#        if len(fasta.seq) < len(seq['res']):
#          diflen = len(seq['res']) - len(fasta.seq)
#          for addgappos in range(0,diflen):
#            fasta.seq[addgappos] = '-'
#        minindex, maxindex, fastaseq, seqres, ali1, ali2  = getaliindex(fasta.seq,fasta.id,seq['res'])

      if 'sequence' in fasta.id:
        fastaseqtmp = str(fasta.seq)
        seqrestmp = ''.join(seq['res'])
#        if ( '-' not in fastaseqtmp or '-' not in seqrestmp ) and len(fastaseqtmp) == len(seqrestmp):
        if fastaseqtmp == seqrestmp:  # case: XXXX | XXXX
          minindex = 0
          maxindex = len(seq['res']) - 1
          alifastaseq = str(fasta.seq)
          aliseqres = ''.join(seq['res'])
        else:  # cases: sequences difer in length or aa sequence
          minindex, maxindex, fastaseq, seqres, alifastaseq, aliseqres = getaliindex(fasta.seq,fasta.id,seq['res'])
#chk        for pos in range(0,len(fasta.seq)):
        minindex = int(minindex)
        maxindex = int(maxindex)
        for pos in range(minindex,maxindex):
#X STARTEST
          if pos < len(alifastaseq):
#X            if fasta.seq[pos] != '-' and ((fasta.id[18:22],fasta.id[23]) not in pdbidchain[pos]):
            if alifastaseq[pos] != '-' and ((fasta.id[18:22],fasta.id[23]) not in pdbidchain[pos]):
              pdbidchain[pos].append((fasta.id[18:22],fasta.id[23]))
              readflag = 1
              if PDBss['PDBid'][pos] != '.':
                PDBss['PDBid'][pos] = '%s,%s' % (PDBss['PDBid'][pos],fasta.id[18:22])
                PDBss['PDBchain'][pos] = '%s,%s' % (PDBss['PDBchain'][pos],fasta.id[23])
#              elif fasta.seq[pos] != '-':
              else:
                PDBss['PDBid'][pos] = fasta.id[18:22]
                PDBss['PDBchain'][pos] = fasta.id[23]
            else:
              readflag = 0  # '-' in position or PDBID:CHAIN already read
#X ENDTEST
          try:
#chk            if fasta.seq[pos] != '-' and ((fasta.id[18:22],fasta.id[23]) not in pdbseq[pos]):
#X            if seq['res'][pos] != '-' and ((fasta.id[18:22],fasta.id[23]) not in pdbseq[pos]):  # 
            if aliseqres[pos] != '-' and ((fasta.id[18:22],fasta.id[23]) not in pdbseq[pos]):  # 
              pdbseq[pos].append((fasta.id[18:22],fasta.id[23]))  # 18-22:PDBID; 23:PDBCHAIN 
          except IndexError:
            print 'INDEX_ERROR:', fasta.id, '\n', fasta.seq, '\n', ''.join(seq['res']), '\n', minindex, '\n', maxindex, '\n'
#            print fasta.seq, '\n', seq['res'], '\n', ''.join(seq['res']).strip("-"), '\n', minindex, '\n', maxindex, '\n', fastaseq, '\n', seqres, '\n', ali1, '\n', ali2
            break
          else:
#            readflag = 1
#          if readflag == 1:
#X            if pos >= len(fasta.seq):
            if pos >= len(alifastaseq):
              if PDBss['PDBseq'][pos] != '.':
                PDBss['PDBseq'][pos] = '%s,%s' % (PDBss['PDBseq'][pos],'-')
              else:
                PDBss['PDBseq'][pos] = '-'
##X            else:
            elif readflag == 1:
              if PDBss['PDBseq'][pos] != '.':
#X                PDBss['PDBseq'][pos] = '%s,%s' % (PDBss['PDBseq'][pos],fasta.seq[pos])
#X              elif fasta.seq[pos] != '-':
#X                PDBss['PDBseq'][pos] = fasta.seq[pos]
                PDBss['PDBseq'][pos] = '%s,%s' % (PDBss['PDBseq'][pos],alifastaseq[pos])
              elif alifastaseq[pos] != '-':
                PDBss['PDBseq'][pos] = alifastaseq[pos]
#      if PDBss['PDBseq'][pos] != '.':
#          if fasta.seq[pos] != '-':
#X          if pos < len(fasta.seq):
##X          if pos < len(alifastaseq):
#X            if fasta.seq[pos] != '-' and ((fasta.id[18:22],fasta.id[23]) not in pdbidchain[pos]):
##X            if alifastaseq[pos] != '-' and ((fasta.id[18:22],fasta.id[23]) not in pdbidchain[pos]):
##X              pdbidchain[pos].append((fasta.id[18:22],fasta.id[23]))
##X              readflag = 1
##X              if PDBss['PDBid'][pos] != '.':
##X                PDBss['PDBid'][pos] = '%s,%s' % (PDBss['PDBid'][pos],fasta.id[18:22])
##X                PDBss['PDBchain'][pos] = '%s,%s' % (PDBss['PDBchain'][pos],fasta.id[23])
#              elif fasta.seq[pos] != '-':
##X              else:
##X                PDBss['PDBid'][pos] = fasta.id[18:22]
##X                PDBss['PDBchain'][pos] = fasta.id[23]
##X            else:
##X              readflag = 0  # '-' in position or PDBID:CHAIN already read
      elif 'SEQRES' in fasta.id:
#        continue
#chk        for pos in range(0,len(fasta.seq)):
        for pos in range(int(minindex),int(maxindex)):
##X          if readflag == 1:  # testing readflag, now also in SEQRES; if not working, move all below one indentation left
#          if fasta.seq[pos] != '-' and ((fasta.id[18:22],fasta.id[23]) not in pdbseqres[pos]):
#X         if seq['res'][pos] != '-' and ((fasta.id[18:22],fasta.id[23]) not in pdbseqres[pos]):
          try:
##X            if aliseqres[pos] != '-' and ((fasta.id[18:22],fasta.id[23]) not in pdbseqres[pos]):
            if alifastaseq[pos] != '-' and ((fasta.id[18:22],fasta.id[23]) not in pdbseqres[pos]):
              pdbseqres[pos].append((fasta.id[18:22],fasta.id[23]))
              readflag = 1  ##X
            else:  ##X
              readflag = 0  # '-' in position or PDBID:CHAIN already read  ##X
          except IndexError:
            print 'INDEX_ERROR:', fasta.id, '\n', fasta.seq, '\n', ''.join(seq['res']), '\n', minindex, '\n', maxindex, '\n'
            break
          else:
#X            if pos >= len(fasta.seq):
              if pos >= len(alifastaseq):
                if PDBss['PDBseqres'][pos] != '.':
                  PDBss['PDBseqres'][pos] = '%s,%s' % (PDBss['PDBseqres'][pos],'-')
                elif PDBss['PDBseq'][pos] != '.':
                  PDBss['PDBseqres'][pos] = '-'
##X              else:
              elif readflag == 1:
                if PDBss['PDBseqres'][pos] != '.':
#X                PDBss['PDBseqres'][pos] = '%s,%s' % (PDBss['PDBseqres'][pos],fasta.seq[pos])
                  PDBss['PDBseqres'][pos] = '%s,%s' % (PDBss['PDBseqres'][pos],alifastaseq[pos])
                elif PDBss['PDBseq'][pos] != '.':
  #          elif fasta.seq[pos] != '-':
#X                PDBss['PDBseqres'][pos] = fasta.seq[pos]
                  PDBss['PDBseqres'][pos] = alifastaseq[pos]
##X            else:  ##X
##X              readflag = 0  # '-' in position or PDBID:CHAIN already read  ##X
#          if PDBss['PDBseqres'][pos] != '.' and ((fasta.id[18:22],fasta.id[23]) not in pdbidchain[pos]):
#            PDBss['PDBseqres'][pos] = '%s,%s' % (PDBss['PDBseqres'][pos],fasta.seq[pos])
#          elif PDBss['PDBseq'][pos] != '.':
##          elif fasta.seq[pos] != '-':
#            PDBss['PDBseqres'][pos] = fasta.seq[pos]
      elif 'SEQATOM' in fasta.id:
#        continue
#chk        for pos in range(0,len(fasta.seq)):
        for pos in range(minindex,maxindex):
##X          if readflag == 1:  # testing readflag, compared with no readflag in SEQRES; if not working, move all below one indentation left
#            if seq['res'][pos] != '-' and ((fasta.id[18:22],fasta.id[23]) not in pdbseqres[pos]):
            try:
##X              if seq['res'][pos] != '-' and ((fasta.id[18:22],fasta.id[23]) not in pdbseqatom[pos]):
              if alifastaseq[pos] != '-' and ((fasta.id[18:22],fasta.id[23]) not in pdbseqatom[pos]):  ##X
                pdbseqatom[pos].append((fasta.id[18:22],fasta.id[23]))
                readflag = 1  ##X
              else:  ##X
                readflag = 0  ##X
            except IndexError:  ##X
              print 'INDEX_ERROR', '\n', fasta.seq, '\n', ''.join(seq['res']).strip("-"), '\n', fastaseq, '\n', seqres, '\n', alifastaseq, '\n', aliseqres, '\n', minindex, '\n', maxindex, '\n', pos  ##X
              break  ##X
            else:  ##X
##X                if pos  >= len(fasta.seq):
                if pos >= len(alifastaseq): ##X
                  if PDBss['PDBseqatom'][pos] != '.':
                    PDBss['PDBseqatom'][pos] = '%s,%s' % (PDBss['PDBseqatom'][pos],'-')
                  elif PDBss['PDBseq'][pos] != '.':
                    PDBss['PDBseqatom'][pos] = '-'
##X                else:
                elif readflag == 1:  ##X
                  if PDBss['PDBseqatom'][pos] != '.':
##X                    PDBss['PDBseqatom'][pos] = '%s,%s' % (PDBss['PDBseqatom'][pos],fasta.seq[pos])
                    PDBss['PDBseqatom'][pos] = '%s,%s' % (PDBss['PDBseqatom'][pos],alifastaseq[pos])  ##X
                  elif PDBss['PDBseq'][pos] != '.':
#              elif fasta.seq[pos] != '-':
##X                    PDBss['PDBseqatom'][pos] = fasta.seq[pos]
                    PDBss['PDBseqatom'][pos] = alifastaseq[pos]  ##X
##X            except IndexError:
##X              print 'INDEX_ERROR', '\n', fasta.seq, '\n', ''.join(seq['res']).strip("-"), '\n', fastaseq, '\n', seqres, '\n', alifastaseq, '\n', aliseqres, '\n', minindex, '\n', maxindex, '\n', pos
##X              break
      elif 'secstr' in fasta.id:
        continue
#chk        for pos in range(0,len(fasta.seq)):
        for pos in range(minindex,maxindex):
          if readflag == 1:
            if PDBss['PDBss'][pos] != '.':
              PDBss['PDBss'][pos] = '%s,%s' % (PDBss['PDBss'][pos],fasta.seq[pos])
            elif PDBss['PDBseq'][pos] != '.':
             PDBss['PDBss'][pos] = fasta.seq[pos]
      elif 'disorder' in fasta.id:
        continue
#        for pos in range(0,len(fasta.seq)):
        for pos in range(minindex,maxindex):
          if readflag == 1:
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

#print(testreadSIFTSparse(inSIFTSparse,ELMpos,seq))
#exit()
PDBss = readSIFTSparse(inSIFTSparse,ELMpos,seq)
#print type(PDBss)

results = seq.copy()  # merge tables
results.update(PDBss)
results.update(predictions)

printTable(results)

