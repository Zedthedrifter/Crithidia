#!/usr/bin/env python3


#%% MODULES TO IMPORT 

from __future__ import division
import matplotlib.pyplot as plt
import os
import re
import sys
import numpy as np
import argparse
import subprocess as sbp

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def reorient(input,output,csb3):
  mdict=SeqIO.to_dict(SeqIO.parse(input,'fasta'))
  for k in mdict:
    if csb3 not in mdict[k].seq:
      mdict[k].seq=mdict[k].seq.reverse_complement() #reorient
      a=mdict[k].seq.find(csb3)
      #print(k,len(mdict[k].seq),a)
  records=[SeqRecord(mdict[k].seq,id=mdict[k].id,name=mdict[k].name,description=mdict[k].description) for k in mdict]
  #SeqIO.write(records, output, "fasta")
  return(mdict)
  
def align(mdict,output,csb3):
  count=0
  for k in mdict:
    #check if two CSB3s are present
    if len(re.findall(csb3,str(mdict[k].seq)))!=2:
      print(f"{k}: only {len(re.findall(csb3,str(mdict[k].seq)))} CSB3 motif is found")
    realign_motif='TTAATATATAGA'  #this motif only appear once
    motif1=mdict[k].seq.find(realign_motif) 
    seq=mdict[k].seq[motif1:]+mdict[k].seq[:motif1]
    mdict[k].seq=seq
    print(len(seq))
  records=[SeqRecord(mdict[k].seq,id=f"Cf_mO{i}",name='',description='') for i,k in enumerate(mdict)]
  SeqIO.write(records,output, "fasta")
  #print(count)
  
def read_coverage(coverage):  
  with open(coverage) as f:
    keys=next(f).strip('\n').split('\t')
    cdict={l.strip('\n').split('\t')[0]:{k:v for k,v in zip(keys,l.strip('\n').split('\t'))} for l in f}
  numreads = np.array([cdict[k]['numreads'] for k in cdict])
  meandepth=np.array([cdict[k]['meandepth'] for k in cdict])
  maxidepth=cdict['Cf_k99_contigk129_1263_len24398']['meandepth']
  mini_permaxi=[round(float(md)/float(maxidepth),3) for md in meandepth]
  percentage=[round(md/sum(mini_permaxi)*100,2) for md in mini_permaxi]
  print(percentage)
  print(mini_permaxi)
  print(cdict)
  plt.pie(mini_permaxi,labels = cdict.keys())
  plt.show()
  plt.pie(numreads,labels = cdict.keys())
  plt.show()
  plt.pie(meandepth,labels = cdict.keys())
  plt.show()
  return(cdict)
read_coverage('coverage.txt')

def main(input,output):
  csb3='GGGGTTGGTGT'
  mdict=reorient(input,output,csb3)
  align(mdict,output,csb3)

#%% MAIN
#if __name__ == '__main__':
#  parser = argparse.ArgumentParser(description='extract circularized contigs.')
#  parser.add_argument('input', help='Name of the input file', metavar='input')
#  parser.add_argument('output', help='Name of the output file', metavar='output')
 # options = parser.parse_args()
  
#  main(options.input,options.output)


def conserved(mdict,output='Cf.noncircularized.fasta'):
  csb3rc='ACACCAACCCC'
  for k in mdict:
    if csb3rc in mdict[k].seq:
      mdict[k].seq=mdict[k].seq.reverse_complement() #reorient
    a=mdict[k].seq.find(csb3)
    motif1=mdict[k].seq.find(realign_motif)
    if motif1!=-1:
      mdict[k].seq=mdict[k].seq[motif1:]+mdict[k].seq[:motif1]
  records=[SeqRecord(mdict[k].seq,id=mdict[k].id,name=mdict[k].name,description=mdict[k].description) for k in mdict]
  SeqIO.write(records, output, "fasta")
  

def process_blast(bfile):
  query,subject='',''
  f=open('top_hits.txt','w')
  for l in open(bfile):
    l2=l.split('\t')
    if l2[0] != query :
      f.write(f"\n{l}")
      query=l2[0]
      subjuct=l2[1]
    elif (l2[0]==query and l2[1]==subject) :
      f.write(f"{l}")
      query=l2[0]
      subjuct=l2[1]
  f.close()
  
def get_unique_uncirc(mdict,bfile,output='Cf.noncircularized.uniq.fasta'):
  match100=sbp.check_output(f"grep '100.000' {bfile} |cut -f 1|uniq",shell=True).decode().strip('\n').split('\n')
  records=[SeqRecord(mdict[k].seq,id=mdict[k].id,name=mdict[k].name,description=mdict[k].description) for k in mdict if (k not in match100)]
  SeqIO.write(records, output, "fasta")
      
def pick_csb3(mdict,output='tmp.Cf.csb3contigs.fasta'):
  csb3rc='ACACCAACCCC'
  records=[]
  for k in mdict:
    if len(mdict[k].seq)>=500:
      if csb3rc in mdict[k].seq:
        mdict[k].seq=mdict[k].seq.reverse_complement() #reorient
        records.append(SeqRecord(mdict[k].seq,id=mdict[k].id,name=mdict[k].name,description=mdict[k].description))
        a=mdict[k].seq.find(csb3)
        print(a,len(mdict[k].seq))
      elif csb3 in mdict[k].seq:
        records.append(SeqRecord(mdict[k].seq,id=mdict[k].id,name=mdict[k].name,description=mdict[k].description))
        a=mdict[k].seq.find(csb3)
        print(a,len(mdict[k].seq))
    motif1=mdict[k].seq.find(realign_motif)
    if motif1!=-1:
      mdict[k].seq=mdict[k].seq[motif1:]+mdict[k].seq[:motif1]
  print(len(records))
  SeqIO.write(records, output, "fasta")
  return(records)
  
def base_content(mini_dict,window,figw=20,figh=20):
  fig,ax = plt.subplots(1,1,figsize=(figw,figh))
  base_color={'A':'green','T':'red','C':'blue','G':'orange'}
  for m in mini_dict:
    seq=mini_dict[m].seq
    for base in base_color.keys():
      bc=[seq[i:i+window].count(base)/window for i in range(0,len(seq),window)]
      ax.plot(range(0,len(seq),window),bc,color=base_color[base])
  plt.show()
    
def base_frequence(aligned_dict,w=1,figw=20,figh=20): #used aligned sequences
  fig,ax = plt.subplots(1,1,figsize=(figw,figh))
  base_color={'A':'green','T':'red','C':'blue','G':'orange'}
  for k in aligned_dict:
    aligned_dict[k].seq=aligned_dict[k].seq.upper()
  l=len(list(aligned_dict.values())[0].seq)
  bdict={i:[aligned_dict[m].seq[i] for m in aligned_dict] for i in range(l)}
  for base in base_color.keys():
    bc=[bdict[i].count(base)/len(bdict[i]) for i in bdict]
    bc=[sum(bc[i:i+w])/w for i in range(0,l,w)]
    ax.plot(range(0,l,w),bc,color=base_color[base])
  plt.show()