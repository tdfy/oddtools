#BSUB -J process2
#BSUB -o process2.%J.out
#BSUB -e process2.%J.error
#BSUB -n 26
#BSUB -M 204800
#BSUB -R "span[hosts=1] rusage [mem=204800]"

echo -e "Sample\tMapped\tUnMapped\tCentrometric\tPeritelomeric\tMean_Kmer\tSD_Kmer\tMax_Kmer\tMotif kmer2\tMotif kmer2" >> results.tsv

paste TARG | while read var;

do

var2=${var%.sam.bam}

MAPPED=`/home/yoderto/miniconda3/bin/samtools view -c "${var}"`

#/home/yoderto/miniconda3/bin/samtools view -b -f 4 "${var}" > unmapped."${var2}".bam
UNMAP=`/home/yoderto/miniconda3/bin/samtools view -c unmapped."${var2}".bam`


#/home/yoderto/miniconda3/bin/samtools view -b -h -L centro.bed "${var}" > centro."${var2}".bam
CENTRO=`/home/yoderto/miniconda3/bin/samtools view -c centro."${var2}".bam`

#/home/yoderto/miniconda3/bin/samtools view -b -h -L pertelo.38.bed "${var}" > telo."${var2}".bam
TELO=`/home/yoderto/miniconda3/bin/samtools view -c telo."${var2}".bam`


kmer1=`/home/yoderto/miniconda3/bin/samtools view "${var}" | grep -E 'TTAGGGTTAGGG|AATCCCAATCCC' | wc -l`

kmer2=`/home/yoderto/miniconda3/bin/samtools view "${var}" | /home/yoderto/miniconda3/bin/agrep -2 TTAGGGTTAGGGTTAGGG | wc -l`

kmer2c=`/home/yoderto/miniconda3/bin/samtools view "${var}" | /home/yoderto/miniconda3/bin/agrep -2 AATCCCAATCCCAATCCC | wc -l`

(( kmerplus = $kmer2 + $kmer2c ))

grep -n -o  'TTAGGG' /home/yoderto/telo/fastq/"${var2}".fastq | sort -n | uniq -c | cut -d : -f 1 > ${var2}.txt

grep -n -o 'AATCCC' /home/yoderto/telo/fastq/"${var2}".fastq | sort -n | uniq -c | cut -d : -f 1 >> ${var2}.txt

awk '(NR>1) && ($1 > 1) ' ${var2}.txt > ${var2}1.txt

AVG=`awk '{ total += $1 } END { print total/NR }' ${var2}1.txt`

SD=`awk '{sum+=$1; sumsq+=$1*$1} END {print sqrt(sumsq/NR - (sum/NR)^2)}' ${var2}1.txt`

max=`awk 'BEGIN{a=   0}{if ($1>0+a) a=$1} END{print a}' ${var2}1.txt`

echo -e "$var2\t$MAPPED\t$UNMAP\t$CENTRO\t$TELO\t$AVG\t$SD\t$max\t$kmer1\t$kmerplus" >> results.tsv


done;
