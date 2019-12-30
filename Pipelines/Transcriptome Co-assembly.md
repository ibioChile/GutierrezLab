## Pipeline overview

![Transcriptome Co-Assembly Pipeline](https://user-images.githubusercontent.com/53570955/71528170-619f9600-28bd-11ea-802b-5597542c3e9e.png)

## QC/QA

Total number of reads per species were counted. Overall quality of raw reads was obtained as follows:
```bash
fastqc -o reads *fastq.gz
multiqc reads/
```

To remove low quality bases, reads were trimmed with Trimommatic with the following parameters:
```bash
java -jar /projects2/software/opt/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 R1.fastq.gz R2.fastq.gz R1.trim.fastq.gz R1.unpaired.fastq.gz R2.trim.fastq.gz R2.unpaired.fastq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:10:30 MINLEN:50 AVGQUAL:25
```
> Beware of Illumina adapters: if they are present in the sample, they can be removed with ILLUMINACLIP option. In our case, adapters were removed by the Illumina machine.

Quality of trimmed reads was re-evaluated as follows:
```bash
fastqc -o trimmed *trim.fastq.gz
multiqc trimmed/
```

> Multiqc generates a report in html format that you can visually inspect to evaluate the average quality of your libraries. 

## De novo transcriptome assembly

Three different tools will be used for de novo transcriptome assembly: Bridger, Trinity and RNASpades. 

Running Bridger:
```bash
export LD_LIBRARY_PATH=/usr/local/boost1.57/lib:$LD_LIBRARY_PATH
export CPPFLAGS="-I/usr/local/boost1.57/"
export LIBS="-L/usr/local/boost1.57/lib"

perlbrew switch perl5.24.0_with
/projects2/software/opt/Bridger_r2014-12-01/Bridger.pl --seqType fq --left R1.trim.fastq.gz --right R2.trim.fastq.gz --CPU 20 -o BridgerOut --min_kmer_coverage 5 --SS_lib_type RF
```

Running Trinity:
```
cat R1.trim.fastq.gz R1.unpaired.fastq.gz > R1_all.trim.fq
cat R2.trim.fastq.gz R2.unpaired.fastq.gz > R2_all.trim.fq
mkdir trinity_out/
Trinity --seqType fq  --left R1_all.trim.fq --right R2_all.trim.fq --CPU 20 --no_normalize_reads --output trinity_out/ --max_memory 200G
```

Running RNASpades

```
mmkdir out_spades/
rnaspades.py --ss-rf --pe1-1 R1.trim.fastq.gz --pe1-2 R2.trim.fastq.gz --pe1-s R1.unpaired.fastq.gz --pe2-s R2.unpaired.fastq.gz -o out_spades/ -t 32 -m 500
```

> Notice that these commands consider a paired-end (PE) library.

Importantly, Bridger is not isoform-aware, it outputs the consensus transcript (longest possible isoform). This is important to consider in downstream analysis.

## De novo transcriptome assembly evaluation

De novo transcriptomes for non-model species can be evaluated using different approaches. 

First, we evaluated contig statistics with TrinityStats:
```bash
awk 'BEGIN{count=0}{if($1~/>/){{count++}{print ">comp"count"_c0__""'"$1"'"}}else{print $0}}' asm.bridger.fasta > asm.bridger.fasta_rename # renaming Bridger output to Trinity-like format
/home2/programs/trinityrnaseq_r20131110/util/TrinityStats.pl asm.bridger.fasta_rename > asm.bridger.fasta_rename.stats
```

To evaluate completeness, we used BUSCO:
```bash
# Running BUSCO using plant database
python3 /usr/local/BUSCO_v1.1b1/plant_early_release/BUSCO_plants.py -i asm.bridger.fasta -o asm.bridger.fasta.BUSCO -l /usr/local/BUSCO_v1.1b1/plant_early_release/plantae -m trans -c 20 -f
# Running BUSCO using complete database
python3 /projects2/software/opt/BUSCO-3.0/scripts/run_BUSCO.py -i FILENAME.fasta -o DIRECTORY -l /projects2/software/opt/BUSCO-3.0/datasets/embryophyta_odb9/ -m tran -c 20 -f
```

## Co-assembly

The script used to combined multiple assemblies is part of this [publication](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1406-x#MOESM3) and can be found [here]( https://static-content.springer.com/esm/art%3A10.1186%2Fs12859-016-1406-x/MediaObjects/12859_2016_1406_MOESM3_ESM.pl). We named this script ```Suplementary_file.pl```.

```
conda install -c bioconda cd-hit 
conda install -c bioconda transdecoder=3.0.0 

mkdir raw_assemblies/
cp ../salida.inicial/transcripts.fasta ../BridgerOut_paired/Bridger.fasta ../trinity_out/Trinity.fasta raw_assemblies/

mv Bridger.fasta alerce_assembly_bridger.fasta
mv Trinity.fasta alerce_assembly_trinity.fasta
mv transcripts.fasta alerce_assembly_spades.fasta

perl Suplementary_file.pl -i raw_assemblies/ -w Concatenated_assembly -Cd /home/pcamejo/anaconda2/envs/trinity/bin/ -Tr /home/pcamejo/anaconda2/envs/trinity/bin/
```
## Annotation

