#BSUB -J process2
#BSUB -o process2.%J.out
#BSUB -e process2.%J.error
#BSUB -n 26
#BSUB -M 204800
#BSUB -R "span[hosts=1] rusage [mem=204800]"

echo -e "Sample\tMapped\tUnMapped\tCentrometric\tPeritelomeric\tMean_Kmer\tSD_Kmer\tMax_Kmer\tMotif kmer2\tMotif kmer2" >> results.tsv

paste T2T.TARG | while read var;

do

var2=${var%.sam.bam}

var3=${var%_T2T.sam.bam}

MAPPED=`/home/yoderto/miniconda3/bin/samtools view -c "${var}"`

/home/yoderto/miniconda3/bin/samtools view -b -f 4 "${var}" > unmapped."${var2}".bam
UNMAP=`/home/yoderto/miniconda3/bin/samtools view -c unmapped."${var2}".bam`


#/home/yoderto/miniconda3/bin/samtools view -b -h -L centro.bed "${var}" > centro."${var2}".bam
CENTRO= 'null' # `/home/yoderto/miniconda3/bin/samtools view -c centro."${var2}".bam`

/home/yoderto/miniconda3/bin/samtools view -b -h -L T2T.telo.bed "${var}" > telo."${var2}".bam

/home/yoderto/miniconda3/bin/samtools view telo."${var2}".bam | cut -f 10 | awk '{ print length }' >> telo.bases_"$var2".txt

TELO=`/home/yoderto/miniconda3/bin/samtools view -c telo."${var2}".bam`

kmer1=`/home/yoderto/miniconda3/bin/samtools view "${var}" | grep -E 'TTAGGGTTAGGG|AATCCCAATCCC' | wc -l`

/home/yoderto/miniconda3/bin/samtools view "${var}" | grep -E 'TTAGGGTTAGGG|AATCCCAATCCC' | cut -f 10 | awk '{ print length }' >> 2mer_"$var2".txt

kmer2=`/home/yoderto/miniconda3/bin/samtools view "${var}" | /home/yoderto/miniconda3/bin/agrep -2 TTAGGGTTAGGGTTAGGG | wc -l`

samtools view "${var}" | /home/yoderto/miniconda3/bin/agrep -2 TTAGGGTTAGGGTTAGGG | cut -f 10 | awk '{ print length }' >> 3mer_"$var2".txt

kmer2c=`/home/yoderto/miniconda3/bin/samtools view "${var}" | /home/yoderto/miniconda3/bin/agrep -2 AATCCCAATCCCAATCCC | wc -l`

samtools view "${var}" | /home/yoderto/miniconda3/bin/agrep -2 AATCCCAATCCCAATCCC | cut -f 10 | awk '{ print length }' >> 3mer_"$var2".txt

(( kmerplus = $kmer2 + $kmer2c ))

grep -n -o  'TTAGGG' /home/yoderto/telo/fastq/"${var3}".fastq | sort -n | uniq -c | cut -d : -f 1 > ${var2}.txt

grep -n -o 'AATCCC' /home/yoderto/telo/fastq/"${var3}".fastq | sort -n | uniq -c | cut -d : -f 1 >> ${var2}.txt

awk '(NR>1) && ($1 > 1) ' ${var2}.txt > ${var2}1.txt

AVG=`awk '{ total += $1 } END { print total/NR }' ${var2}1.txt`

SD=`awk '{sum+=$1; sumsq+=$1*$1} END {print sqrt(sumsq/NR - (sum/NR)^2)}' ${var2}1.txt`

max=`awk 'BEGIN{a=   0}{if ($1>0+a) a=$1} END{print a}' ${var2}1.txt`

echo -e "$var2\t$MAPPED\t$UNMAP\t$CENTRO\t$TELO\t$AVG\t$SD\t$max\t$kmer1\t$kmerplus" >> T2T.results.tsv


done;
