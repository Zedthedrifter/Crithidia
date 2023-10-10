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


def output_lm_align(lmed,lmun,cfun):
  f=open('tmp.ed_un.txt','w')
  for k in lmed:
    if lmed[k].seq.upper().replace('U','')==lmun[k].seq.upper().replace('T',''):
      print(f"{k} validated")
    else:
      print(f"{k} nonT base unmatched")
    f.write(f">{k}_lmed\t{lmed[k].seq}\n>{k}_lmun\t{lmed[k].seq.replace('u',' ')}\n>{k}_cfun\t{cfun[k].seq.replace('T','U')}\n\n")

def main():
  lmed=SeqIO.to_dict(SeqIO.parse('edited_mRNA_small_u.fasta','fasta'))
  lmun=SeqIO.to_dict(SeqIO.parse('Lm_unedited_mRNAs.DNA.fasta','fasta'))
  cfun=SeqIO.to_dict(SeqIO.parse('Cf.unedited.mRNAs.fasta','fasta'))
  output_lm_align(lmed,lmun,cfun)
  
main()