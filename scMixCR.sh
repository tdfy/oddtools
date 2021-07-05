~/software/bamtofastq_linux --locus=chr7:142299011-142813287  --nthreads=10 possorted_genome_bam.bam  ~/workspace/tcr/

seqtk trimfq bamtofastq_S1_L001_R1_001.fastq.gz -b 44 > trim_bamtofastq_S1_L001_R1_001.fastq

mixcr align -s hs -p rna-seq -OallowPartialAlignments=true -OsaveOriginalReads=true trim_bamtofastq_S1_L001_R1_001.fastq.gz bamtofastq_S1_L001_R2_001.fastq.gz alignments.vdjca

mixcr assemblePartial alignments.vdjca alignmentsRescued_1.vdjca

mixcr assemblePartial alignmentsRescued_1.vdjca alignmentsRescued_2.vdjca

mixcr extend alignmentsRescued_2.vdjca alignmentsRescued_2_extended.vdjca

mixcr assemble -a alignmentsRescued_2_extended.vdjca clones.clna

mixcr exportClones -c TRA -o -t clones.clna clones.txt

mixcr exportReadsForClones --id 2  clones.clna reads_of_my_clones.fastq.gz
