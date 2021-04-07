 # Differential allelic expression analysis
 
This pipeline explains how to screen genes for allele-specific expression among different conditions. The workflow starts with the generation of RNA libraries from the genotypes Arabidopsis Landsberg erecta & Llagostera. Libraries are generated from parental and hybrid genotypes. In this case, we will use the genome of [Arabidopsis ecotype Col-0](https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas) as reference for reads mapping. Overall, the following pipeline follows the next steps:
 
 2. Calling of discriminant SNPs between Landsberg erecta and Llagostera.
 3. RNA-seq mapping and differential expression among parental genotypes.
 4. Estimate allele-specific expression among hybrids.
 


##  1. Calling of discriminant SNPs between Landsberg erecta and Llagostera.

1.1 Trimming adapters and regions of poor quality.

	java -jar trimmomatic-0.39.jar PE -phred33 R1.fastq.gz R2.fastq.gz R1.trim.fastq.gz R1.unpaired.fastq.gz R2.trim.fastq.gz R2.unpaired.fastq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:10:30 MINLEN:50 AVGQUAL:25
	
1.2 Combine all paired reads from samples of pure Landsberg and Llagostera.

```
cat LeLe*._1.fq.gz.trim.fil.pair.gz > LeLe_1.all.trim.pair.fq
cat LeLe*._2.fq.gz.trim.fil.pair.gz > LeLe_2.all.trim.pair.fq
cat LlLl*._1.fq.gz.trim.fil.pair.gz > LeLe_1.all.trim.pair.fq
cat LlLl*._2.fq.gz.trim.fil.pair.gz > LeLe_2.all.trim.pair.fq
```

1.3 Map Landsberg erecta and Llagostera reads (paired and unpaired) with HiSat2 to the genome of Col-0 (***Tomás***).
 
```
bwa index Trinity.cdhit99.fasta
bwa mem -t 20 -o LeLe.all.trim.sam Trinity.cdhit99.fasta LeLe_1.all.trim.fq LeLe_2.all.trim.fq
bwa mem -t 20 -o LlLl.all.trim.sam Trinity.cdhit99.fasta LlLl_1.all.trim.fq LlLl_2.all.trim.fq
```

1.4 Process sam file, fix mates and remove duplicates.

```
samtools view -S -b LeLe.all.trim.sam > LeLe.all.trim.bam
samtools sort -n  --threads 20 -o LeLe.all.trim.sorted.bam LeLe.all.trim.bam
samtools fixmate -m LeLe.all.trim.sorted.bam LeLe.all.trim.fixmate.bam
samtools sort LeLe.all.trim.fixmate.bam -o LeLe.all.trim.fixmate.sorted.bam
samtools markdup -r LeLe.all.trim.fixmate.sorted.bam LeLe.all.trim.fixmate.sorted.dedup.bam
samtools index LeLe.all.trim.fixmate.sorted.dedup.bam
```

Repeat these steps for reads of Llagostera mapping to Col-0.


1.5 Call set of variants using FreeBayes.

```
#Index genome of Col-0
samtools faidx TAIR10_Chr.all.fa

samtools depth LlLl.all.trim.fixmate.sorted.dedup.bam | coverage_to_regions.py TAIR10_Chr.all.fa.fai 10000 > targets.bed 

picard AddOrReplaceReadGroups I=LeLe.all.trim.fixmate.sorted.dedup.bam O=LeLe.all.trim.fixmate.sorted.dedup.RG.bam RGSM=LeLe.all RGPL=illumina RGLB=LeLe RGPU=LeLe RGID=1

picard AddOrReplaceReadGroups I=LlLl.all.trim.fixmate.sorted.dedup.bam O=LlLl.all.trim.fixmate.sorted.dedup.RG.bam RGSM=LlLl.all RGPL=illumina RGLB=LlLl RGPU=LlLl RGID=2

samtools index LeLe.all.trim.fixmate.sorted.dedup.RG.bam

samtools index LlLl.all.trim.fixmate.sorted.dedup.RG.bam

freebayes-parallel targets.bed 44  -f TAIR10_Chr.all.fa --haplotype-length 0 --standard-filters --min-alternate-fraction 0.05 -p 2 --pooled-discrete --pooled-continuous LeLe.all.trim.fixmate.sorted.dedup.RG.bam LlLl.all.trim.fixmate.sorted.dedup.RG.bam > fb.LeLe.LlLl.vcf
```

1.6 Use bcftools to filter the resulting VCF file.

```
bcftools view -i 'DP > 10 & SAF > 2 & SAR > 2 & RPR > 1 & RPL > 1'  fb.LeLe.LlLl.vcf  > fb.LeLe.LlLl.filt.vcf
```

1.7 Create a dictionary with the transcriptome with GATK.

```
gatk CreateSequenceDictionary R=TAIR10_Chr.all.fa
```

1.8 Use the SelectVariants walker in GATK to keep only monomorphic variants for Llagostera and for Landsberg erecta.

```
gatk3 -T SelectVariants -V fb.LeLe.LlLl.filt.vcf -R TAIR10_Chr.all.fa  -select 'vc.getGenotype("LeLe.all").isHomRef()' -o fb.LeLe.LlLl.filt.homref-1.vcf

gatk3 -T SelectVariants -V fb.LeLe.LlLl.filt.homref-1.vcf -R TAIR10_Chr.all.fa  -select 'vc.getGenotype("LlLl.all").isHomVar()' -o fb.LeLe.LlLl.filt.homvar-2.vcf

gatk3 -T SelectVariants -V fb.LeLe.LlLl.filt.vcf -R TAIR10_Chr.all.fa  -select 'vc.getGenotype("LlLl.all").isHomRef()' -o fb.LeLe.LlLl.filt.homvar-1.vcf

gatk3 -T SelectVariants -V fb.LeLe.LlLl.filt.homvar-1.vcf -R TAIR10_Chr.all.fa  -select 'vc.getGenotype("LeLe.all").isHomVar()' -o fb.LeLe.LlLl.filt.homref-2.vcf
```

1.9 Zip and index vcf files.

```
bgzip fb.LeLe.LlLl.filt.homref-2.vcf
bgzip fb.LeLe.LlLl.filt.homvar-2.vcf
 
tabix fb.LeLe.LlLl.filt.homref-2.vcf.gz
tabix fb.LeLe.LlLl.filt.homvar-2.vcf.gz
```

1.10 Concat vcf files

```
bcftools concat fb.LeLe.LlLl.filt.homvar-2.vcf.gz fb.LeLe.LlLl.filt.homref-2.vcf.gz -o fb.LeLe.LlLl.filt.homref.homvar.vcf -a
```

1.11  To keep only high-confidence sites, filter out sites with AO < 20 and AO / DP < 0.99.

```
grep "#" fb.LeLe.LlLl.filt.homref.homvar.vcf > fb.LeLe.LlLl.filt.final-1.vcf
grep "#" fb.LeLe.LlLl.filt.homref.homvar.vcf > fb.LeLe.LlLl.filt.final-2.vcf
grep -v "#" fb.LeLe.LlLl.filt.homref.homvar.vcf > temp.vcf

while read line; do AO=$(echo $line | awk '{print $11}' | cut -d ":" -f3); DP=$(echo $line | awk '{print $11}' | cut -d ":" -f4); ratio=$(echo $AO $DP | awk '{print $1/$2}'); if [[ $AO > 19 && $ratio > 0.99 ]]; then echo "$line" >> fb.LeLe.LlLl.filt.final-1.vcf ;fi;  done < temp.vcf

while read line; do AO=$(echo $line | awk '{print $10}' | cut -d ":" -f3); DP=$(echo $line | awk '{print $10}' | cut -d ":" -f4); ratio=$(echo $AO $DP | awk '{print $1/$2}'); if [[ $AO > 19 && $ratio > 0.99 ]]; then echo "$line" >> fb.LeLe.LlLl.filt.final-2.vcf ;fi;  done < temp.vcf 
```

1.12 Zip and index vcf files.

```
bgzip fb.LeLe.LlLl.filt.final-1.vcf
bgzip fb.LeLe.LlLl.filt.final-2.vcf
 
tabix fb.LeLe.LlLl.filt.final-1.vcf.gz
tabix fb.LeLe.LlLl.filt.final-2.vcf.gz
```

1.13 Concat vcf files and remove duplicates

```
bcftools concat fb.LeLe.LlLl.filt.final-1.vcf.gz fb.LeLe.LlLl.filt.final-2.vcf.gz -o fb.LeLe.LlLl.filt.final.vcf -a
bgzip fb.LeLe.LlLl.filt.final.vcf
tabix fb.LeLe.LlLl.filt.final.vcf.gz
```

1.14 Remove non-biallelic variants.

```
vcftools --vcf fb.LeLe.LlLl.filt.final.dedup.vcf --min-alleles 2 --max-alleles 2 --recode --out fb.LeLe.LlLl.filt.final.nonbi
```

## 2. RNA-seq mapping and differential expression among parental genotypes.

2.1 In order to avoid mapping biases to the reference when estimating allele-specific expression, built a custom pseudogenome using the FastaAlternateReferenceMaker walker in GATK.

```
gatk IndexFeatureFile -I fb.LeLe.LlLl.filt.final.nonbi.recode.vcf

gatk FastaAlternateReferenceMaker -R TAIR10_Chr.all.fa -V fb.LeLe.LlLl.filt.final.nonbi.recode.vcf --snp-mask fb.LeLe.LlLl.filt.final.nonbi.recode.vcf --snp-mask-priority -O fb.LeLe.LlLl.filt.final.mask.fasta

#Since FastaAlternateReferenceMaker adds a number before each sequence:
cat fb.LeLe.LlLl.filt.final.mask.fasta | sed 's/>.*Chr/>Chr/' | sed 's/:.*//' > fb.LeLe.LlLl.filt.final.mask.fixed.fasta
```

2.2  Map reads from each library against the pseudogenome using HiSat2 (***Tomás***).

```
# build bowtie2 index:
bowtie2-build fb.LeLe.LlLl.filt.final.mask.fixed.fasta fb.LeLe.LlLl.filt.final.mask.fixed.bowtie2index

# Align reads:
bowtie2 -x fb.LeLe.LlLl.filt.final.mask.fixed.bowtie2index -1 $1_1.fq.gz.trim.fil.pair.gz -2 $1_2.fq.gz.trim.fil.pair.gz  -p 20 -S $1.sam >$1.log
```

2.3. Mark duplicates with samtools.

```
samtools view -S -b  $1.sam > $1.bam;
samtools sort -n -o $1.sorted.bam $1.bam;
samtools fixmate -m $1.sorted.bam $1.fixmate.bam;
samtools sort $1.fixmate.bam -o $1.fixmate.sorted.bam;
samtools markdup -r $1.fixmate.sorted.bam $1.fixmate.sorted.dedup.bam;
samtools sort -o $1.sorted.bam $1.fixmate.sorted.dedup.bam;
samtools index $1.fixmate.sorted.dedup.bam;
```

2.4 Quantify expression levels.

```
Prepare a pseudo Rsubread annotation file with the transcript (as chromosomes).

cat  fb.LeLe.LlLl.filt.final.mask.fixed.rsem.transcripts.fa |perl -ane 'chomp($_); print "$_\t"; ' | sed 's/>/\n>/g' |perl -ane '@l=split(/\t/); chomp(@l); $a=length($l[1]); print  "$l[0]\t$l[0]\t1\t$a\t+\n";' |sed 's/>//g' |grep TR>fb.LeLe.LlLl.filt.final.mask.fixed.rsem.transcripts.rsubread

# in R:

library(Rsubread)
ann<-read.table("fb.LeLe.LlLl.filt.final.mask.fixed.rsem.transcripts.rsubread",sep="\t")
colnames(ann)<-c("GeneID","Chr","Start","End", "Strand")
lista.bam<-dir(pattern=".bam$")
c0<-featureCounts(lista.bam,annot.ext=ann,allowMultiOverlap=T,isPairedEnd=T,nthreads=20,strandSpecific=0)
write.table(c0$counts,"conteos.por.transcrito.txt",sep="\t",col.names=NA,quote=F)
```

## 3. Estimate allele-specific expression among hybrids.

3.1 Assign groups to bam files and index.

```
for file in *.dedup.bam; do base=${file##*/}; picard AddOrReplaceReadGroups I=$file O=bam_dedup/${base%.*}.RG.bam RGSM=LeLl RGPL=illumina RGLB=LeLl RGPU=LeLl RGID=1 VALIDATION_STRINGENCY=LENIENT; done

for file in bam_dedup/*.RG.bam; do samtools index $file; done
```

3.2 Use the ASEReadCounter walker in GATK to retrieve Llagostera and Landsberg erecta allele counts at discriminant SNPs 

```
for file in bam_dedup/*.RG.bam; do gatk3 -T ASEReadCounter -I $file -R TAIR10_Chr.all.fa -sites fb.LeLe.LlLl.filt.final.nonbi.recode.vcf -o $file.out -minDepth 30 -mmq 40 -mbq 20 -U ALLOW_N_CIGAR_READS; done
```
 
4.3 Analize ASEReadCounter output files with this [R script](https://github.com/ibioChile/GutierrezLab/blob/master/scripts/VCF_processing2.R).
