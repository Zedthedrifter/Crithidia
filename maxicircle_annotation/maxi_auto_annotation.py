#!/usr/bin/env python3


#%% MODULES TO IMPORT 

from __future__ import division
import matplotlib.pyplot as plt
import os
import re
import sys
import numpy
import argparse
import subprocess as sbp

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def read_alignments(atxt):
  query,subject=[],[]
  align=open(atxt)
  for l in align:
    if l.startswith('Q'):
      l=[e for e in l.split(' ') if e != ''] 
      query.append(l[2])
    elif l.startswith('S'):
      l=[e for e in l.split(' ') if e != ''] 
      subject.append(l[2])
  query=''.join(query)
  subject=''.join(subject)
  return(query,subject)

def read_lm_annotation():
  lm=[(614,1774,'12S rRNA'),
(1810,2419,'9S rRNA'),
(2493,2792,'ND8'),
(3114,2793,'ND9'),
(3394,3089,'MURF5'),
(3407,4583,'ND7'),
(4587,5467,'COX3'),
(5478,6582,'CYB'),
(6633,7246,'A6'),
(8584,7251,'ND2'),
(8546,8709,'GR3'),
(9614,8673,'ND1'),
(9623,10251,'COX2'),
(10243,11321,'MURF2'),
(12961,11312,'COX1'),
(13224,13013,'GR4'),
(13299,14611,'ND4'),
(14825,14594,'ND3'),
(14797,15002,'RPS12'),
(15063,16834,'ND5')]
  lmdict={item[2]:{'start':item[0],'end':item[1]} for item in lm}
  return(lmdict)

def annotate_cf(cf,lm,lmdict,addcf=2399,addlm=27): #addlm/cf: where the alignments start
  cfdict={}
  for k in lmdict:
    cfcount,lmcount=addcf,addlm
    for c,l in zip(cf,lm):
      if c!='-':
        cfcount+=1
      if l!='-':
        lmcount+=1
      if lmcount==lmdict[k]['start']:
        cfstart=cfcount
      if lmcount==lmdict[k]['end']:
        cfend=cfcount
    cfdict[k]={'start':cfstart,'end':cfend}
  return(cfdict)
    
def output_cfunedited(maxi,cfdict):
  maxi=SeqIO.read(maxi,'fasta')
  records=[SeqRecord(maxi.seq[min([cfdict[k]['start'],cfdict[k]['end']])-1:max([cfdict[k]['start'],cfdict[k]['end']])],id=k,name=k,description=f"{cfdict[k]['start']}-{cfdict[k]['end']}") for k in cfdict]
  records=[]
  for k in cfdict:
    if cfdict[k]['start']<cfdict[k]['end']:
      records.append(SeqRecord(maxi.seq[cfdict[k]['start']-1:cfdict[k]['end']],id=k,name=k,description=f"plus: {cfdict[k]['start']}-{cfdict[k]['end']}"))
    elif cfdict[k]['start']>cfdict[k]['end']:
      records.append(SeqRecord(maxi.seq[cfdict[k]['end']-1:cfdict[k]['start']].reverse_complement(),id=k,name=k,description=f"minus: {cfdict[k]['start']}-{cfdict[k]['end']}"))
  return(records)

def main():
  atxt='aligned_maxicircles.fasta'
  maxi='Cf.maxicircles.fasta'
  cf,lm=read_alignments(atxt)
  lmdict=read_lm_annotation()
  cfdict=annotate_cf(cf,lm,lmdict)
  records=output_cfunedited(maxi,cfdict)
  SeqIO.write(records,'Cf.unedited.mRNAs.fasta','fasta')
  for k in cfdict:
    print(f"drawarrow({cfdict[k]['start']},{cfdict[k]['end']},{k})")
main()