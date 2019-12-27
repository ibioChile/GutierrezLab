## Pipeline overview

<img src="https://user-images.githubusercontent.com/7208981/57542790-23f7f800-7342-11e9-9a01-88da4aeda87c.png" width="400">

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
/projects2/software/opt/Bridger_r2014-12-01/Bridger.pl --seqType fq --left R1.trimmed.fastq --right R2.trimmed.fastq --CPU 20 -o BridgerOut --min_kmer_coverage 5 --SS_lib_type RF
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
> CEGMA is not mantained anymore by the Korf Lab.

Reads used for the assembly can be mapped back to the assembly. The number of reads that can be correctly mapped to a transcript indicates the amount of information incorporated in the final de novo assembly. These statistics can be obtained with transrate:
```bash
transrate --assembly asm.bridger.fasta --left R1.trimmed.fastq --right R2.trimmed.fastq
```

## Assemblies merging

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

## Protein prediction and clustering

From assembled transcripts, protein sequence was predicted using Transdecoder. This software finds all 6 possible open reading frames within the transcript and translates the coding sequence into aminoacids. 

Proteins were obtained with the following command:
```bash
perl /projects2/software/opt/trinityrnaseq_r20131110/trinity-plugins/transdecoder/TransDecoder -t asm.bridger.fasta --workdir transdecoder --cd_hit_est /usr/local/bin/cd-hit-est
```

Then, CD-HIT was used to group together similar proteins and remove from the dataset small proteins contained within a larger protein version. CD-HIT was run as follows:
```bash
cd-hit -i asm.bridger.fasta.transdecoder.pep -o asm.bridger.fasta.transdecoder.cdhit90.pep -T 0 -M 0 -c 0.9
```

## Removing known non-plant species

Proteins matching non-plants species were removed from the proteomes, leaving plants proteins as well as unclassified/unknown proteins in the final dataset. 

- Proteomes were compared to Uniprot SwissProt using BLASTP
- Alignments with e-value < 1e-10 and identity > 50% were selected
- Taxa corresponding to the subject was added
- Protein IDs matching a "Viridiplantae" sequences were saved in a separate file
- Protein IDs matching something different than "Viridiplantae" were saved in a separate file
- All protein IDs were saved in a separate file
- Protein IDs not matching any species were obtained by substracting total protein IDs and protein IDs matching and organism different than "Viridiplantae"

These steps can be achieved as follows:
```bash
blastall -p blastp -a 30 -m 8 -i asm.bridger.fasta.transdecoder.cdhit90.pep -o asm.bridger.fasta.transdecoder.cdhit90.pep.uniprot -d /projects2/databases/UniprotKB/uniprot_sprot.fasta
awk -F"\t" '{if($11 < 1e-10) print $0}' asm.bridger.fasta.transdecoder.cdhit90.pep.uniprot | awk -F"\t" '{if($3 > 50.00) print $0}' | sort -k1,1 -u | cut -f1,2 > tempfile1
join -1 2 -2 1 -a 1 <(sort -k2,2 tempfile1) <(sort -k1,1 /projects1/dcsoto/databases/swissprot_taxa.tab) | awk '{print $2"\t"$1"\t"$3"\t"$5}' > tempfile2
grep -w "Viridiplantae" tempfile2 | cut -f1 > tempfile3 # plant sequences
grep -w -v "Viridiplantae" tempfile2 | cut -f1 > tempfile4 # no plant sequences
grep ">" asm.bridger.fasta.transdecoder.cdhit90.pep | cut -d" " -f1 > tempfile5
grep -F -v -f tempfile4 tempfile5 > tempfile6
perl /projects1/dcsoto/code/extract-seq-id.pl tempfile6 asm.bridger.fasta.transdecoder.cdhit90.pep > asm.bridger.fasta.transdecoder.cdhit90.filtered.pep
```
> Note: Sequence ID names may differ and renaming steps might be needed. 

Proteome size was calculated before and after this steps, and total amount of proteins from non-plants species was reported. 

## Big Plant Pipeline dataset

The final set of predicted proteins (plant+unclassified) was sent to our collaborators for the Big Plant pipeline. This dataset included 32 Atacama plants and 32 "sister" species from the California Desert, selected according to their phylogenetic distance from Atacama species as well as their availability. These dataset comes from 3 different sources SRA, Dryad and 1KP. Some of these species underwent the same pipeline as Atacama species, while others were retrieved as assemblies. We also included 6 model organisms and crops. Proteomes for these species were downloaded from their respective websites or JGI.

## Annotation

Transcripts were annotated using InterProScan (functional annotation based on protein domains prediction), and E2P2 (pathway annotation), and Gene Ontology (GO). 

InterPro domains were found as follow:
```bash
/projects2/software/opt/interproscan-5.18-57.0/interproscan.sh -t n -i asm.bridger.cdhit90.fasta -iprlookup -goterms -pa -f tsv
/projects2/software/opt/interproscan-5.18-57.0/interproscan.sh -i asm.bridger.cdhit90.fasta.transdecoder.pep -iprlookup -goterms -pa -f tsv
```

Pathway annotation was performed as follow:
```bash
cd /projects2/software/opt/e2p2-3.0/
./runE2P2.v3.0.py -i asm.bridger.cdhit90.fasta.transdecoder.pep -o asm.bridger.cdhit90.fasta.transdecoder.pep.e2p2v3
```

GO terms were transferred to transcripts using 3 sources of information: _A. thaliana_ hits, Uniprot SwissProt hits and InterProScan. 

1) _A. thaliana_ GO terms:
```bash
blastall -p blastp -a 20 -m 8 -i asm.bridger.fasta.transdecoder.pep.cdhit -o asm.bridger.fasta.transdecoder.pep.cdhit.athal -d /projects1/dcsoto/databases/TAIR10_pep_20110103_representative_gene_model_updated
paste <(awk  -F"[|]m." '{print $1}' asm.bridger.fasta.transdecoder.pep.cdhit.athal) <(cut -f2- asm.bridger.fasta.transdecoder.pep.cdhit.athal) | sed 's/cds.//g' > tempfile1
awk '{if($3>50 && $11<1e-10){print $0}}' tempfile1| sort -k11,11g | sort -u -k1,1 > tempfile2
join -t $'\t' -1 2 -2 1 <(cut -d"." -f1 tempfile2 |sort -k2,2) <(sort -k1,1 /projects1/dcsoto/databases/ATH_GO.tab2) | cut -f2- > asm.bridger.fasta.annot.go.arath
```

2) Uniprot GO terms:
```bash
awk '{if($3>50 && $11<1e-10){print $0}}' asm.bridger.fasta.transdecoder.pep.cdhit.uniprot |cut -f1,14 |sed 's/cds.//g' |sed 's/; /\t/g' |grep "GO:" > tempfile
paste <(cut -d"|" -f1 tempfile) <(cut -f2- tempfile) > asm.bridger.fasta.annot.go.uniprot
```

3) InterPro GO terms
```bash
cut -f1,14 asm.bridger.cdhit90.fasta.transdecoder.pep.tsv | awk '{if($2!=""){print $0}}' > tempfile1 # InterProScan output
paste <(awk  -F"[|]m." '{print $1}' tempfile1) <(cut -f2 tempfile1 |sed 's/|/\t/g') |sed 's/cds.//g' |sort -u > tempfile2
python /projects1/dcsoto/code/columnizador.py -i tempfile2 -o tempfile3
sort -u tempfile3 > tempfile4
python /projects1/dcsoto/code/filalizador.py -i tempfile4 -o asm.bridger.fasta.annot.go.interproscan
```

All three GO annotation approaches were merged and formatted for GoStats:
```bash
sort asm.bridger.fasta.annot.go.arath asm.bridger.fasta.annot.go.interproscan asm.bridger.fasta.annot.go.swissprot > tempfile1
python /projects1/dcsoto/code/columnizador.py -i tempfile1 -o tempfile2
sort -u tempfile2 | awk '{print $2"\tISS\t"$1}' > asm.bridger.fasta.annot.go.all.gostats
```

## Transcript abundance analysis

This analysis was performed using manually calculated transcripts per million (TPM). However, I strongly recommend using a software like Salmon for transcripts quantification. This software uses a probabilistic model to obtain adjusted TPMs that consider the complexities of multiple isoforms and multimapping reads, as well as RNA-seq biases.

Steps for this analysis:
- Calculate TPMs using a RNA-seq quantification software
- Compare ranking of relative abundance in different species
- Explore annotations for transcripts with high relative abundance

