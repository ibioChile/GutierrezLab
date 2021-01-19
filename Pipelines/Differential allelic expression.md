 # Differential allelic expression analysis
 
 This pipeline explains how to screen genes for allele-specific expression among different conditions. The workflow starts with the generation of RNA libraries from the genotypes Arabidopsis Landsberg erecta & Llagostera. Libraries are generated from parental and hybrid genotypes. Overall, the following pipeline follows the next steps:
 
 1. Assembly of the transcriptome of Landsberg.
 2. Calling of discriminant SNPs between Landsberg erecta and Llagostera.
 3. Differential expression among parental genotypes.
 4. Estimate allele-specific expression among hybrids.
 
 
 ## 1. Assembly of the transcriptome of Landsberg.
 
1.1 Trimming adapters and regions of poor quality.

	java -jar trimmomatic-0.39.jar PE -phred33 R1.fastq.gz R2.fastq.gz R1.trim.fastq.gz R1.unpaired.fastq.gz R2.trim.fastq.gz R2.unpaired.fastq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:10:30 MINLEN:50 AVGQUAL:25

1.2 All treatments from the parental genotype of Landsberg were gathered and assembled in transcripts using Trinity v.2.8.5 software.

```
cat R1.trim.fastq.gz R1.unpaired.fastq.gz > R1_all.trim.fq
cat R2.trim.fastq.gz R2.unpaired.fastq.gz > R2_all.trim.fq
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




 
 
