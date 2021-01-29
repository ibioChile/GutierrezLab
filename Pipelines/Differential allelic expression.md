 # Differential allelic expression analysis
 
 This pipeline explains how to screen genes for allele-specific expression among different conditions. The workflow starts with the generation of RNA libraries from the genotypes Arabidopsis Landsberg erecta & Llagostera. Libraries are generated from parental and hybrid genotypes. Overall, the following pipeline follows the next steps:
 
 1. Assembly of the transcriptome of Landsberg.
 2. Calling of discriminant SNPs between Landsberg erecta and Llagostera.
 3. RNA-seq mapping and differential expression among parental genotypes.
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

samtools sort LeLe.all.trim.fixmate.sorted.dedup.RG.bam

samtools sort LlLl.all.trim.fixmate.sorted.dedup.RG.bam

freebayes-parallel targets.bed 44  -f Trinity.cdhit99.fasta --haplotype-length 0 --standard-filters --min-alternate-fraction 0.05 -p 2 --pooled-discrete --pooled-continuous LeLe.all.trim.fixmate.sorted.dedup.RG.bam LlLl.all.trim.fixmate.sorted.dedup.RG.bam > fb.LeLe.LlLl.vcf
```

2.4 Use bcftools to filter the resulting VCF file.

```
bcftools view -i 'DP > 10 & SAF > 2 & SAR > 2 & RPR > 1 & RPL > 1'  fb.LeLe.LlLl.vcf  > fb.LeLe.LlLl.filt.vcf
```

2.5 Create a dictionary with the transcriptome with GATK.

```
gatk CreateSequenceDictionary R=Trinity.cdhit99.fasta
```

2.6 Use the SelectVariants walker in GATK to keep only monomorphic variants for an alternative allele in Llagostera and for the reference allele in Landsberg erecta.

```
gatk3 -T SelectVariants -V fb.LeLe.LlLl.filt.vcf -R Trinity.cdhit99.fasta  -select 'vc.getGenotype("LeLe.all").isHomRef()' -o fb.LeLe.LlLl.filt.homref-1.vcf

gatk3 -T SelectVariants -V fb.LeLe.LlLl.filt.homref-1.vcf -R Trinity.cdhit99.fasta  -select 'vc.getGenotype("LlLl.all").isHomVar()' -o fb.LeLe.LlLl.filt.homvar-2.vcf
```

2.7  To keep only high-confidence sites, filter out sites with AO < 20 and AO / DP < 0.99.

```
grep "#" fb.LeLe.LlLl.filt.homvar-2.vcf > fb.LeLe.LlLl.filt.final.vcf
grep -v "#" fb.LeLe.LlLl.filt.homvar-2.vcf > temp.vcf

while read line; do AO=$(echo $line | awk '{print $11}' | cut -d ":" -f3); DP=$(echo $line | awk '{print $11}' | cut -d ":" -f4); ratio=$(echo $AO $DP | awk '{print $1/$2}'); if [[ $AO > 19 && $ratio > 0.99 ]]; then echo "$line" >> fb.LeLe.LlLl.filt.final.vcf ;fi;  done < temp.vcf
```

## 3. RNA-seq mapping and differential expression among parental genotypes.

3.1 In order to avoid mapping biases to the reference when estimating allele-specific expression, built a custom pseudotranscriptome using the FastaAlternateReferenceMaker walker in GATK.

```
gatk IndexFeatureFile -I fb.LeLe.LlLl.filt.final.vcf

gatk FastaAlternateReferenceMaker -R ../../LeLe/trinity_out/Trinity.cdhit99.fasta -V fb.LeLe.LlLl.filt.final.vcf --snp-mask fb.LeLe.LlLl.filt.final.vcf --snp-mask-priority -O fb.LeLe.LlLl.filt.final.mask.fasta

#Since FastaAlternateReferenceMaker adds a number before each sequence:
cat fb.LeLe.LlLl.filt.final.mask.fasta | sed 's/>.*TRI/>TRI/'| sed 's/:.*//' > fb.LeLe.LlLl.filt.final.mask.fixed.fasta
```

3.2  Map reads from each library against the pseudotranscriptome using bowtie2.

```
# build bowtie2 index:
bowtie2-build fb.LeLe.LlLl.filt.final.mask.fixed.fasta fb.LeLe.LlLl.filt.final.mask.fixed.bowtie2index

# Align reads:
bowtie2 -x fb.LeLe.LlLl.filt.final.mask.fixed.bowtie2index -1 $1_1.fq.gz.trim.fil.pair.gz -2 $1_2.fq.gz.trim.fil.pair.gz  -p 20 -S $1.sam >$1.log
```

3.3. Mark duplicates with samtools.

```
samtools view -S -b  $1.sam > $1.bam;
samtools sort -n -o $1.sorted.bam $1.bam;
samtools fixmate -m $1.sorted.bam $1.fixmate.bam;
samtools sort $1.fixmate.bam -o $1.fixmate.sorted.bam;
samtools markdup -r $1.fixmate.sorted.bam $1.fixmate.sorted.dedup.bam;
samtools sort -o $1.sorted.bam $1.fixmate.sorted.dedup.bam;
samtools index $1.fixmate.sorted.dedup.bam;
```

3.4 Quantify expression levels.

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

3.5 Annotate transcripts using blast.

```
# Download and extract cDNA file from arabidosis.:
wget https://www.arabidopsis.org/download_files/Sequences/Araport11_blastsets/Araport11_genes.201606.cdna.new.fasta.gz
gzip -d Araport11_genes.201606.cdna.new.fasta.gz

# Make a blast database for nucleotide blast:
makeblastdb -in Araport11_genes.201606.cdna.new.fasta  -input_type fasta -dbtype nucl -title Araport11_genes.201606.cdna.new.blast

# Run blast:
blastall -p blastn -i ../fb.LeLe.LlLl.filt.final.mask.fixed.rsem.transcripts.fa -d Araport11_genes.201606.cdna.new.fasta -m 8 -a 20 -o salida.m8

# sort file and select the best match:
cat salida.m8 | sort -k1,1  -k12nr,12  >salida.m8.sort

cat salida.m8.sort |perl seleccionaprimeros.pl | perl -ane '@l=split(/\t/); chomp(@l);@g=split(/\./,$l[1]); print join "\t",@l,$g[0]; print "\n";'  >topscore.blast.txt
```

## 4. Estimate allele-specific expression among hybrids.

4.1 Assign groups to bam files and index.

```
for file in *.dedup.bam; do base=${file##*/}; picard AddOrReplaceReadGroups I=$file O=bam_dedup/${base%.*}.RG.bam RGSM=LeLl RGPL=illumina RGLB=LeLl RGPU=LeLl RGID=1 VALIDATION_STRINGENCY=LENIENT; done

for file in bam_dedup/*.RG.bam; do samtools index $file; done
```

4.2 Use the ASEReadCounter walker in GATK to retrieve Llagostera and Landsberg erecta allele counts at discriminant SNPs 

```
for file in bam_dedup/*.RG.bam; do gatk3 -T ASEReadCounter -I $file -R Trinity.cdhit99.fasta -sites fb.LeLe.LlLl.filt.final.vcf -o $file.out -minDepth 30 -mmq 40 -mbq 20 -U ALLOW_N_CIGAR_READS; done
```
 
 
