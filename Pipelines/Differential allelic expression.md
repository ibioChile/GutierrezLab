 # Differential allelic expression analysis
 
 This pipeline explains how to screen genes for allele-specific expression among different conditions. The workflow starts with the generation of RNA libraries from the genotypes Arabidopsis Landsberg erecta & Llagostera. Libraries are generated from parental and hybrid genotypes. Overall, the following pipeline follows the next steps:
 
 1. Assembly of the transcriptome of Landsberg.
 2. Calling of discriminant SNPs between Landsberg erecta and Llagostera.
 3. Differential expression among parental genotypes.
 4. Estimate allele-specific expression among hybrids.
 
 
 ## 1. Assembly of the transcriptome of Landsberg.
 
1.1 Trimming adapters and regions of poor quality.

	java -jar trimmomatic-0.39.jar PE -phred33 R1.fastq.gz R2.fastq.gz R1.trim.fastq.gz R1.unpaired.fastq.gz R2.trim.fastq.gz R2.unpaired.fastq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:10:30 MINLEN:50 AVGQUAL:25

1.2 Gather all treatments from the parental genotype of Landsberg and assemble in transcripts using Trinity v.2.8.5 software.

```
#Combine all paired reads
cat LeLe*._1.fq.gz.trim.fil.pair.gz > LeLe_1.all.trim.pair.fq
cat LeLe*._2.fq.gz.trim.fil.pair.gz > LeLe_2.all.trim.pair.fq

#Combine unpaired reads
cat LeLe*._1.fq.gz.trim.fil.unpair.gz > LeLe_1.all.trim.unpair.fq
cat LeLe*._2.fq.gz.trim.fil.unpair.gz > LeLe_2.all.trim.unpair.fq

#Combine paired and unpaired reads
cat LeLe_1.all.trim.pair.fq LeLe_1.all.trim.unpair.fq > LeLe_1.all.trim.fq
cat LeLe_2.all.trim.pair.fq LeLe_2.all.trim.unpair.fq > LeLe_2.all.trim.fq

#Run Trinity
mkdir trinity_out/
Trinity --seqType fq  --left R1_all.trim.fq --right R2_all.trim.fq --CPU 20 --no_normalize_reads --output trinity_out/ --max_memory 200G
```

1.3 Remove redundant transcripts using CD-HIT to cluster transcripts that share 99% similarity.

```
cd-hit-est -i Trinity.fasta -o Trinity.cdhit99.fasta -c 0.99
```

1.4 Evaluate contig statistics with TrinityStats.

```
TrinityStats.pl Trinity.cdhit99.fasta > Trinity.cdhit99.fasta.stats
```

1.4 Evaluate completeness using BUSCO.

```
python3 BUSCO_plants.py -i Trinity.cdhit99.fasta -o Trinity.cdhit99.fasta.BUSCO -l /projects2/software/opt/BUSCO-3.0/datasets/embryophyta_odb9/ -m tran -c 20 -f
```

1.5 Assess the quality of transcriptome using TransRate.

```
transrate --assembly Trinity.cdhit99.fasta --left R2.trim.fastq.gz --right R1.trim.fastq.gz --threads 32
```

##  2. Calling of discriminant SNPs between Landsberg erecta and Llagostera.


2.1 Map Landsberg erecta and Llagostera reads (paired and unpaired) to the assembled transcriptome of Landsberg erecta with bwa.

```
bwa index Trinity.cdhit99.fasta
bwa mem -t 20 -o LeLe.all.trim.sam Trinity.cdhit99.fasta LeLe_1.all.trim.fq LeLe_2.all.trim.fq
bwa mem -t 20 -o LlLl.all.trim.sam Trinity.cdhit99.fasta LlLl_1.all.trim.fq LlLl_2.all.trim.fq
```

2.2 Process sam file, fix mates and remove duplicates.

```
samtools view -S -b LeLe.all.trim.sam > LeLe.all.trim.bam
samtools sort -n  --threads 20 -o LeLe.all.trim.sorted.bam LeLe.all.trim.bam
samtools fixmate -m LeLe.all.trim.sorted.bam LeLe.all.trim.fixmate.bam
samtools sort LeLe.all.trim.fixmate.bam -o LeLe.all.trim.fixmate.sorted.bam
samtools markdup -r LeLe.all.trim.fixmate.sorted.bam LeLe.all.trim.fixmate.sorted.dedup.bam
samtools index LeLe.all.trim.fixmate.sorted.dedup.bam
```
Repeat these steps for reads of Llagostera mapping to the transcriptome of Landsberg.


2.3 Call set of variants using FreeBayes.

```
samtools faidx Trinity.cdhit99.fasta

samtools depth LlLl.all.trim.fixmate.sorted.dedup.bam | coverage_to_regions.py Trinity.cdhit99.fasta.fai 10000 > targets.bed 

picard AddOrReplaceReadGroups I=LeLe.all.trim.fixmate.sorted.dedup.bam O=LeLe.all.trim.fixmate.sorted.dedup.RG.bam RGSM=LeLe.all RGPL=illumina RGLB=LeLe RGPU=LeLe RGID=1

picard AddOrReplaceReadGroups I=LlLl.all.trim.fixmate.sorted.dedup.bam O=LlLl.all.trim.fixmate.sorted.dedup.RG.bam RGSM=LlLl.all RGPL=illumina RGLB=LlLl RGPU=LlLl RGID=2

freebayes-parallel targets.bed 44  -f Trinity.cdhit99.fasta --haplotype-length 0 --standard-filters --min-alternate-fraction 0.05 -p 2 --pooled-discrete --pooled-continuous LeLe.all.trim.fixmate.sorted.dedup.RG.bam LlLl.all.trim.fixmate.sorted.dedup.RG.bam > fb.LeLe.LlLl.vcf
```



 
 
