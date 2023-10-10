#!/bin/bash
#download files from SRA with SRR accession numbers

export PATH=$PATH:/home/zed/script_collections:/home/zed/.local/bin/:/usr/bin:/bin
now=$(date '+%y-%m-%d-%H-%M')
log_file=log_${now}.txt
stats_file=stats_${now}.txt

<<<<<<< HEAD
#comment add on master branch
=======
#comment1 added on branch 1
<<<<<<< HEAD
#comment2 on branch1
<<<<<<< HEAD
>>>>>>> branch1
=======

#comment3 branch1
<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> branch1
=======
#comment1 branch2
#
#comment1 branch3
>>>>>>> branch3

=======
#comment1 branch2
#comment2 branch2
>>>>>>> branch2
=======
#comment3 branch1
#comment1 branch4
#comment2 branch4
>>>>>>> branch4
function download {
accs=("$@") #use this argument when passing a list , otherwise only the first of the array will be considered
printf "%s\t%s\t%s\t%s\n" 'read_file' 'total_reads' 'CSB3_containing_reads' 'percentage' >>  $stats_file

for acc in "${accs[@]}"

do

echo $acc
prefetch $acc -O .
fastq-dump --split-e ${acc}.sra #split paired end reads/or single reads
#clean up intermediates
rm ${acc}.sra
#unnecessary if using T-aligner#
gzip ${acc}*.fastq

done

}

#for transcriptomic data, remove things mapped to the nuclear genome
function remove_nuclear_contamination {

reference=$1
r=$2
strain=$3

  /usr/bin/bowtie2-align-s --wrapper basic-0 -x ${reference} -p 20 --threads 16 --phred33 --very-sensitive-local --quiet --time -S ${reference}_${strain}.bam -U ${r} 
  samtools view -Sbh -f 4 -o ${reference}_${strain}.unmapped.bam ${reference}_${strain}.bam
  rm ${reference}_${strain}.bam
  samtools fastq -0 ${strain}_dn.fastq.gz -n ${reference}_${strain}.unmapped.bam
  rm ${reference}_${strain}.unmapped.bam

}


function main {
#Cf  transcriptomic data 
#mosquito: SRR9289107 SRR9289105 SRR9289103	
#adherent: SRR9289101	SRR9289100 SRR9289099 
#swimming: SRR9289098	SRR9289097 SRR9289096
accs=' SRR9289101' 
ref_file=/home/zed/disk1/Crithidia/download_data/Cf_chromosome_only.fasta
ref_name=Cf

download $accs
remove_nuclear_contamination ${ref_name} ${accs}.fastq.gz ${accs}

}


main

#RNAseq remove reads mapped to nuclear genomes
#if the files are downloaded already
function RNAseq {

ref_file=/home/zed/disk1/Crithidia/download_data/Cf_chromosome_only.fasta
ref_name=Cf
#bowtie2-build -f ${ref_file} ${ref_name}
for r in $(ls /home/zed/disk1/Crithidia/download_data/*.all.fastq.gz)
    do
    strain=$(echo ${r}|awk -F'/' '{print $NF}');strain=${strain/.all.fastq.gz/};
    remove_nuclear_contamination ${ref_name} ${r} ${strain}
    done

}

#RNAseq
