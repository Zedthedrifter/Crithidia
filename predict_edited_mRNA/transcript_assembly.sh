#!/bin/bash
#assembly of RNA-seq data using RNASpades

export PATH=$PATH:/home/zed/script_collections:/home/zed/.local/bin/:/usr/bin:/bin
now=$(date '+%y-%m-%d-%H-%M')
log_file=log_${now}.txt
stats_file=stats_${now}.txt

function rnaspades {

acc=$1
r=$2

rnaspades.py -s $r -o tmp
mv tmp/hard_filtered_transcripts.fasta ${acc}_hard_filtered_transcripts.fasta
mv tmp/soft_filtered_transcripts.fasta ${acc}_soft_filtered_transcripts.fasta
mv tmp/transcripts.fasta ${acc}_transcripts.fasta
rm -r tmp

}

function basic_map {

ref_file=$1
ref_name=$2
strain=$3
r=$4

bowtie2-build -f ${ref_file} ${ref_name}
 /usr/bin/bowtie2-align-s --wrapper basic-0 -x ${ref_name} -p 20 --threads 16 --phred33 --very-sensitive-local --quiet --time -S ${ref_name}_${strain}.bam -U ${r} 
  samtools view -Sbh -F 4 ${ref_name}_${strain}.bam|samtools sort -o ${ref_name}_${strain}.sorted.bam
  samtools index ${ref_name}_${strain}.sorted.bam
}


function main {
accs='SRR9289107	SRR9289105 SRR9289103	SRR9289101	SRR9289100 SRR9289099 SRR9289098	SRR9289097 SRR9289096' #Cf transcriptomic data	
  function spades_iter {
    accs=("$@")
    for acc in "${accs[@]}"
    do
    r=/home/zed/disk1/Crithidia/download_data/${acc}_dn.fastq.gz
    rnaspades ${acc} ${r}
    done
    }

#spades_iter $accs
#rnaspades Cf /home/zed/disk1/Crithidia/download_data/Cf.dn.fastq.gz
#basic_map Cf.unedited.mRNAs.fasta Cf_un Cf /home/zed/disk1/Crithidia/download_data/Cf.dn.fastq.gz
basic_map /home/zed/disk1/Crithidia/minicircle_annotation/In_files/Cf.edited_cDNAs.fasta Cf_ed Cf /home/zed/disk1/Crithidia/download_data/Cf.dn.fastq
}

main


