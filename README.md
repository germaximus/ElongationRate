# ElongationRate

**Prerequisites:**  
[cutadapt 2.5](https://cutadapt.readthedocs.io/en/stable/index.html)  
[STAR-2.7.2b](https://github.com/alexdobin/STAR)  
[bowtie 1.2.3](http://bowtie-bio.sourceforge.net/index.shtml)
[gffread utility](http://ccb.jhu.edu/software/stringtie/gff.shtml)  
[blast+ 2.9.0](https://blast.ncbi.nlm.nih.gov/)

Transcriptome (polyA captured mRNA-seq) samples were sequenced in PE100 mode on Illumina sequencer. Libraries were prepared with Nextera kit.
Raw sequencing files are available from [GEO]().

### Preparing genome annotation and index files
Mouse genomic sequences and annotation files (GRCm38.p6) were downloaded from the [NCBI repository](http://ftp.ncbi.nih.gov/genomes/M_musculus/). 
To obtain genome assembly, download fasta files of individual chromosomes from ```Assembled_chromosomes/seq/``` folder and concatenate them in the ascending order (omit mitochondrial chromosome and sex chromosomes). Edit chromosome names to match annotation names in gff3 (for example convert ```>ref|NC_000067.6|``` to ```>NC_000067.6```) and say thanks to the NIH staff for letting you do it.  

| files               | MD5 check sum (unzipped)         | Description                                               |
| ------------------- |:--------------------------------:| ----------------------------------------------------------|
| GRCm38.p6.rna.fa    | e90edf06df1dade6f60126c694d58ec6 | RNA in fasta format, coding + noncoding                   |
| GRCm38.p6.genome.fa | 632b214701fb60537eb5f9bf1bf37983 | Genome sequence (nuclear genome only, no sex chromosomes) |
| GRCm38.p6.gbk       | c18a112a3b16049f3091f862c2b32024 | RNA in gene bank format, coding + noncoding               |
| GRCm38.p6.gff3      | e39f2eaf16e41e62b4e6564f36a5c437 | Genome annotation                                         | 


<details><summary><b>Customizing genome annotation</b></summary>  

**Customize genome annotation**  
Annotation of extrachromosomal contigs and sex chromosomes was omitted. 'Gnomon' (Predicted) records from gff file were also omitted and only 'RefSeq' and 'BestRefSeq' (manually curated) kept. Perl and R scripts are included in the GitHub repository.   
```bash
Discard_extrachromosomal_annotation.pl GRCm38.p6.gff3 >GRCm38.p6.custom.gff
Discard_gnomon_annotation.pl >GRCm38.p6.Refseq.gff	# automatically takes GRCm38.p6.custom.gff as an input
```
**Remove non-coding RNA genes**, leave only coding genes with their mRNA, transcript, exon, and CDS children. Fix the gff annotation from previous script by matching gene coordinates with the childern coordinates (occured due to removal of Gnomon features).
```bash
Discard_noncoding_annotation.R
```

**Convert annotation from GFF3 to GTF format**  
```bash
gffread GRCm38.p6.Refseq.coding.gff -T -o GRCm38.p6.Refseq.coding.gtf
# -T          - convert gff/gtf
```
</details>


<details><summary><b>Extract ORF sequences for translation rate estimation</b></summary>  

**Fetch all mRNA records**  
Run ```mRNA_extractor.pl```. First, it takes ```GRCm38.p6.gbk``` and extracts all RefSeq records for every gene including CDS, 5UTR, 3UTR lengths and a sequence. Then, it selects the single RefSeq record as the longest isoform. Sometimes, the ORF lengths of two isoforms are equal, in that case the longest isoform is selected based on the UTR length with 5UTR taking precedence over 3UTR. The script also trims mRNAs by 100 nucleotides flanking CDS.  If 5UTR and/or 3UTR are shorter than 100 nt, it raises a "flag".  

```bash
 perl mRNA_extractor.pl /path/GRCm38.p6.gbk
 # creates an output file named temp3
```
Fill missing 5UTR and 3UTRs with genomic sequences in cases when they are shorter than 100 nt.  
```bash
perl mRNA_genome_filler.pl 
# requires requires temp3 from the previous step in the same folder
# outputs mRNA_100.fasta file
```

mRNA_100.fasta file contains transcripts that can share high degree of homology. It is beneficial to eliminate highly similar transcripts prior to engaging to the main ribo-seq analysis. Run nucleotide blast in all vs. all mode

```bash
# build a database with local sequences
makeblastdb -in mRNA_100.fasta -title "mRNA_100" -dbtype nucl
# blast all sequences against each other
blastn -task blastn -num_threads 4 -outfmt 6 -evalue 0.001 -db mRNA_100.fasta -query mRNA_100.fasta -out blast_result.txt   
```

Extract non-redundant genes from ```blast_result.txt```. Selected blast parameters are not very strict and often assign a good score to a pair of genes that are not too similar.   
```bash
BLASTNprocessor.pl blast_result.txt
# outputs mRNA_100uniq.fasta file
```
</details>

<details><summary><b>Pre-process ribosomal profiling reads</b></summary>  

| Table of index sequences used to multiplex libraries |||
| Index             | 8-nt barcode sequence   | Mice where it was used             |
| ------------------|:-----------------------:| -----------------------------------|
| Ribo-seq Index 1  | TCGCCTTA                |  19-month old mice                 |
| Ribo-seq Index 2  | CTAGTACG                |  19-month old mice                 |
| Ribo-seq Index 3  | TTCTGCCT                |  19-month old mice                 |
| Ribo-seq Index 4  | GCTCAGGA                |  19-month old mice                 | 
| Ribo-seq Index 5  | AGGAGTCC                |  19-month old mice                 |
| Ribo-seq Index 6  | CATGCCTA                |  19-month old mice                 |
| Ribo-seq Index 7  | GTAGAGAG                |  19-month old mice                 |
| Ribo-seq Index 8  | CCTCTCTG                |  19-month old mice                 |
| Ribo-seq Index 9  | AGCGTAGC                |  19-month old mice                 |
| Ribo-seq Index 10 | TCCTCTAC                |  19-month old mice                 |
| Ribo-seq Index 11 | CCTGAGAT                |  19-month old mice                 |
| Ribo-seq Index 11 | TAGCGAGT                |  19-month old mice                 |  



Libraries of 19-month old mice were prepared with custom 8-nt barcodes and sequenced at Novogene in 150 PE mode. Ribosomal footprints are short, therefore only the forward read file (R1) is needed and the R2 file can be discarded.  

```bash
cutadapt -u 1 -m 23 -a AGATCGGAAGAGCACACGTCT --discard-untrimmed
```


<\details>





<details><summary><b>Mapping ribosomal footprints to unique ORFs</b></summary>  
 
Build a Bowtie index out of ```mRNA_100uniq.fasta```  
```bash
bowtie-build  ./bowtie/genomes/mRNA_100uniq.fasta ./bowtie/Mouse_indices/mRNA_100uniq

```
 
<\details>
