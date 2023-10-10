#!/bin/bash
#SNPs calling
export PATH=$PATH:/home/zed/script_collections:/home/zed/.local/bin/:/usr/bin:/bin
now=$(date '+%y-%m-%d-%H-%M')
log_file=log_${now}.txt
stats_file=stats_${now}.txt


function komics_assembly {

log_file=$1
r1=$2
r2=$3
strain=$4
stats_file=$5
kmers=$(echo 89 99 119 129 139 149)

  #echo -e "start KOMICS assembly of ${strain} with kmer size ${kmers}\nr1= ${r1}\nr2= ${r2}" >> ${log_file}
  #assemble minicircle: higher k
  for k in ${kmers}
  do
    komics assemble --threads 16 --kmin ${k} --kmax ${k} ${strain}_k${k} ${r1} ${r2}
    #circularize
    komics circularize ${strain}_k${k} tmp.${strain}_k${k}.csb3contigs.fasta
    #collect all contigs for polishing
    cat  tmp.${strain}_k${k}.circularized.fasta >> tmp.${strain}.circularized.fasta
  done
extract_circ.py tmp.${strain}.circular.fasta tmp.${strain}.circularized.fasta
}

#completeness assessment
function completeness_assessment {

ref=$1
strain=$2
r=$3

  /usr/bin/bowtie2-align-s --wrapper basic-0 -x ${ref} -p 20 --threads 16 --phred33 --very-sensitive --quiet --time --rg-id ${strain} --rg SM:${strain} --rg PL:'ILLUMINA' -S map_mini.bam -U ${r}
  all_reads=$(samtools view map_mini.bam |wc -l)
  mapped=$(samtools view -F 4 map_mini.bam |wc -l)
  mapped_q=$(samtools view -F 4 -q 10 map_mini.bam|wc -l)
  percentmapped=$(echo "scale=2 ; ${mapped} / ${all_reads} * 100" | bc)
  percentmapped_q=$(echo "scale=2 ; ${mapped_q} / ${all_reads} * 100" | bc)
  #printf "Total reads:\n%s\nMapped reads:\n%s %s%%\nMapped reads with quality > 10:\n%s %s%%\n" ${all_reads} ${mapped} ${percentmapped} ${mapped_q} ${percentmapped_q} >> ${log_file}
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s" ${strain} ${minicircles} ${all_reads} ${mapped} ${percentmapped} ${mapped_q} ${percentmapped_q}>>  $stats_file
 #CSB3 reads
  all_reads=$(samtools view map_mini.bam |egrep -c 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC')
  mapped=$(samtools view -F 4 map_mini.bam |egrep -c 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC')
  mapped_q=$(samtools view -F 4 -q 10 map_mini.bam |egrep -c 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC')
  percentmapped=$(echo "scale=2 ; ${mapped} / ${all_reads} * 100" | bc)
  percentmapped_q=$(echo "scale=2 ; ${mapped_q} / ${all_reads} * 100" | bc)
  #printf "Total CSB3 reads:\n%s\nMapped CSB3 reads:\n%s %s%%\nMapped CSB3 reads with quality > 10:\n%s %s%%\n" ${all_reads} ${mapped} ${percentmapped} ${mapped_q} ${percentmapped_q} >> ${log_file}
  printf "\t%s\t%s\t%s\t%s\t%s\n" ${all_reads} ${mapped} ${percentmapped} ${mapped_q} ${percentmapped_q} >>  $stats_file
  samtools view -Sbh -F 4 map_mini.bam|samtools sort -o mini.mapped.sorted.bam
  samtools index mini.mapped.sorted.bam
  samtools view mini.mapped.sorted.bam -Sbh -f 2|bamtools filter -tag XM:'<=1'|samtools sort -o ${strain}.mapped.sorted.paired.XM1.bam
  samtools index ${strain}.mapped.sorted.paired.XM1.bam
  #filter with bamtools (quite helpful in a lot of cases)
  #bamtools filter -tag XM:'<=2' -in $ourdir/${reference}_${strain}.ms.bam -out $ourdir/${reference}_${strain}.msf.bam #<= 2 mismatches
  #bamtools filter -tag XM:0 -in ${reference}_${strain}.ms.bam -out ${reference}_${strain}.msf.bam #0 mismatches
  
}

function bam_stats {

bamfile=$1
  all_reads=$(samtools view $bamfile |wc -l)
  mapped=$(samtools view -F 4 $bamfile |wc -l)
  mapped_q=$(samtools view -F 4 -q 10 $bamfile|wc -l)
  percentmapped=$(echo "scale=2 ; ${mapped} / ${all_reads} * 100" | bc)
  percentmapped_q=$(echo "scale=2 ; ${mapped_q} / ${all_reads} * 100" | bc)
  #printf "Total reads:\n%s\nMapped reads:\n%s %s%%\nMapped reads with quality > 10:\n%s %s%%\n" ${all_reads} ${mapped} ${percentmapped} ${mapped_q} ${percentmapped_q} >> ${log_file}
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s" ${strain} ${minicircles} ${all_reads} ${mapped} ${percentmapped} ${mapped_q} ${percentmapped_q}>>  $stats_file
 #CSB3 reads
  all_reads=$(samtools view $bamfile |egrep -c 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC')
  mapped=$(samtools view -F 4 $bamfile |egrep -c 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC')
  mapped_q=$(samtools view -F 4 -q 10 $bamfile |egrep -c 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC')
  percentmapped=$(echo "scale=2 ; ${mapped} / ${all_reads} * 100" | bc)
  percentmapped_q=$(echo "scale=2 ; ${mapped_q} / ${all_reads} * 100" | bc)
  #printf "Total CSB3 reads:\n%s\nMapped CSB3 reads:\n%s %s%%\nMapped CSB3 reads with quality > 10:\n%s %s%%\n" ${all_reads} ${mapped} ${percentmapped} ${mapped_q} ${percentmapped_q} >> ${log_file}
  printf "\t%s\t%s\t%s\t%s\t%s\n" ${all_reads} ${mapped} ${percentmapped} ${mapped_q} ${percentmapped_q} >>  $stats_file

}

function mafft_multi_alignment {

infile=$1

mafft --genafpair --thread 20 --maxiterate 1000 ${infile}> aligned_${infile}
}

#kmer=31 is the maximum, as the algorithm is different
function velvet_assembly {

kmer=$1
r1=$2
r2=$3

#mkdir paired_${kmer}
velveth paired_${kmer} ${kmer} -fastq -shortPaired -separate ${r1} ${r2} #pair end in seperate files
velvetg paired_${kmer}
}


#use Spades
function spades {

r1=$1
r2=$2

#spades.py --plasmid --careful -t 16 -k 89 99 109 119 127 -o spades -1 $r1 -2 $r2 #assembling only plasmids from WGS data sets
spades.py -k 21,33,55,77 --careful -1 $r1 -2 $r2 -o spades_small_kmer #scaffolds.fasta contains resulting scaffolds (recommended for use as resulting sequences)
#komics circularize ${strain}_k${k} tmp.${strain}_k${k}.csb3contigs.fasta
}


function main {

#list input data
printf "strain\tminicircle_nb\ttotal_reads\tmapped_reads\tmapped_reads_percentage\tmapped_q10\tmapped_q10_percentage\ttotal_csb3\tmapped_csb3\tmapped_csb3_percentage\tmapped_csb3_q10\tmapped_csb3_q10_percentage\n" >> $stats_file
r=/home/zed/disk1/Crithidia/download_data/Cf.dn.fastq.gz
strain=Cf
#run komics assembly and completeness assessment
  completeness_assessment maxi ${strain} ${r}
  #mafft_multi_alignment Cf.minicircles.fasta #Cf.circularized.fasta
  #vsearch --cluster_fast Cf.noncircularized.uniq.fasta --id 0.99 --consout test.fasta
  #bam_stats Cf.circularized.mapped.sorted.mx1.bam
  
}


main