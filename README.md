# ElongationRate

Data analysis for the study [Gerashchenko MV, Peterfi Z, Yim SH, Gladyshev VN. Translation elongation rate varies among organs and decreases with age. Nucleic Acids Res. 2021 Jan 25;49(2):e9. PMID: 33264395; PMCID: PMC7826258](https://doi.org/10.1093/nar/gkaa1103)      

**Prerequisites:**  
[cutadapt 2.5](https://cutadapt.readthedocs.io/en/stable/index.html)  
[STAR-2.7.2b](https://github.com/alexdobin/STAR)  
[bowtie 1.2.3](http://bowtie-bio.sourceforge.net/index.shtml)
[gffread utility](http://ccb.jhu.edu/software/stringtie/gff.shtml)  
[blast+ 2.9.0](https://blast.ncbi.nlm.nih.gov/)

Transcriptome (polyA captured mRNA-seq) samples were sequenced in PE100 mode on Illumina sequencer. Libraries were prepared with Nextera kit.
Raw sequencing files are available from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112223).

### Preparing genome annotation and index files
Mouse genomic sequences and annotation files (GRCm38.p6) were downloaded from the [NCBI repository](https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/all_assembly_versions/GCF_000001635.26_GRCm38.p6/). 
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

Two sets of indexed primers were used for library multiplexing. One set has 6-nt barcodes and the other 8-nt barcodes.  
<details><summary><b>Table of custom 6-nt index sequences used to multiplex libraries</b></summary>  
 

</details>

<details><summary><b>Table of custom 8-nt index sequences used to multiplex libraries</b></summary>  

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
| Ribo-seq Index 12 | TAGCGAGT                |  19-month old mice                 |  
</details>


<details><summary><b>Ribo-seq of livers from 19-month old mice</b></summary>  
 
Liver Ribo-seq libraries of 19-month old mice were prepared with custom 8-nt barcodes and sequenced at Novogene in 150 PE mode. Total 12 libraries were pooled together and sequenced on a single lane. Ribosomal footprints are short, therefore only the forward read file (R1) is needed and the R2 file can be discarded. It is also more convenient to trim and remove rRNA reads from the pooled library before demultiplexing. If you are reproducing this analysis by downloading raw data from NCBI GEO repository, then demultiplexing was already done, but you still have to trimm and remove ribosomal contaminants from individual files.   

```bash
cutadapt -j 25 -u 1 -m 23 -M 40 -a AGATCGGAAGAGCACACGTCT --discard-untrimmed -o trimmed.fastq input.fastq
# j - number of processors
# u - delete first nucleotide of each read
# m - minimum length after adapter trimming
# M - maximum length after adapter trimming
```

Remove ribosomal contaminants  
```bash
 bowtie -p 20 -v 2 --un genomic.fastq ../bowtie-1.2.3/Mouse_indices/rmtRNA trimmed.fastq >/dev/null
```
Demultiplex pooled library into 12 individual samples  
```bash
perl BarcodeSplitter_8nt.pl genomic.fastq
```
</details>

</details>

<details><summary><b>Mapping ribosomal footprints to unique ORFs</b></summary>  
 
Build a Bowtie index out of ```mRNA_100uniq.fasta``` and ribosomal RNA    
```bash
bowtie-build  ./bowtie/genomes/mRNA_100uniq.fasta ./bowtie/Mouse_indices/mRNA_100uniq
bowtie-build  ./bowtie/genomes/Mouse_rmtRNA.fasta ./bowtie/Mouse_indices/rmtRNA
```
Align ribosomal fotprints against ```mRNA_100uniq.fasta```  
```bash
 bowtie -p 20 -v 2 -m 1 --norc --max /sample/redundant.fastq /bowtie-1.2.3/Mouse_indices/mRNA_100uniq /sample/genomic.fastq >uniq.bwt
```
</details>

<details><summary><b>Quick and dirty way of making Ribosome coverage plots</b></summary>  
 
 Run custom perl scipts to calculate ORFs coverage profiles and a metaprofile for every sample.   
 ```bash
 perl Coverage.pl uniq.bwt
 # requires mRNA_100uniq.fastq in the same folder with uniq.bwt
 perl Coverage_processor.pl 2000 start *.coverage
 ```
</details>

<details><summary><b>R-friendly ribosome coverage plots</b></summary>  
 
 Run custom perl scipts to calculate ORFs coverage profiles for every sample.   
 ```bash
 perl Coverage.pl uniq.bwt
 # requires mRNA_100uniq.fastq in the same folder with uniq.bwt
 ```
Transfer coverage files to a separate folder, give them appropriate names, for instance "MI26Li.coverage" and organize a txt table with sample names and description factors for subsequent analysis. Below is the example of a table I got:  

<details><summary>Table</summary> 

| Sample coverage file        | timepoint (sec)   |    organ     |   age  |
| --------------------------- |:-----------------:| :----------: | :----: |
|   MI26K.coverage            |   0               |   kidney     |   3    |
|   MI26Li.coverage           |   0               |   liver      |   3    |
|   MI26SKM.coverage          |   0               |   skeletal   |   3    |
|   MI27K.coverage            |   0               |   kidney     |   3    |
|   MI27Li.coverage           |   0               |   liver      |   3    |
|   MI27SKM.coverage          |   0               |   skeletal   |   3    |
|   MI28Li.coverage           |   0               |   liver      |   3    |
|   MI29Li.coverage           |   0               |   liver      |   3    |
|   MI43K.coverage            |   30              |   kidney     |   3    |
|   MI43Li.coverage           |   30              |   liver      |   3    |
|   MI43SKM.coverage          |   30              |   skeletal   |   3    |
|   MI44K.coverage            |   45              |   kidney     |   3    |
|   MI44Li.coverage           |   45              |   liver      |   3    |
|   MI45Li.coverage           |   60              |   liver      |   3    |
|   MI45K.coverage            |   60              |   kidney     |   3    |
|   MI50K.coverage            |   15              |   kidney     |   3    |
|   MI50Li.coverage           |   15              |   liver      |   3    |
|   MI51K.coverage            |   15              |   kidney     |   3    | 
|   MI51Li.coverage           |   15              |   liver      |   3    |
|   MI51SKM.coverage          |   15              |   skeletal   |   3    |
|   MI52K.coverage            |   30              |   kidney     |   3    |
|   MI52Li.coverage           |   30              |   liver      |   3    |
|   MI52SKM.coverage          |   30              |   skeletal   |   3    |
|   MI53K.coverage            |   45              |   kidney     |   3    |
|   MI53Li.coverage           |   45              |   liver      |   3    |
|   MI53SKM.coverage          |   45              |   skeletal   |   3    |
|   MI54Li.coverage           |   60              |   liver      |   3    |
|   MI54K.coverage            |   60              |   kidney     |   3    |
|   MI70SKM.coverage          |   15              |   skeletal   |   3    |
|   MI71K.coverage            |   15              |   kidney     |   3    |
|   MI71Li.coverage           |   15              |   liver      |   3    |
|   MI71SKM.coverage          |   15              |   skeletal   |   3    |
|   MI72K.coverage            |   30              |   kidney     |   3    |
|   MI72Li.coverage           |   30              |   liver      |   3    |
|   MI72SKM.coverage          |   30              |   skeletal   |   3    |
|   MI73K.coverage            |   45              |   kidney     |   3    |
|   MI74Li.coverage           |   60              |   liver      |   3    |
|   MI74K.coverage            |   60              |   kidney     |   3    |
|   MI75K.coverage            |   15              |   kidney     |   3    |
|   MI76K.coverage            |   30              |   kidney     |   3    |
|   MI77K.coverage            |   45              |   kidney     |   3    |
|   MI77Li.coverage           |   45              |   liver      |   3    |
|   MI77SKM.coverage          |   45              |   skeletal   |   3    |
|   MI78Li.coverage           |   60              |   liver      |   3    |
|   MI78K.coverage            |   60              |   kidney     |   3    |
|   MI13-16poolLi.coverage    |   300             |   liver      |   3    |
|   MI13-16poolK.coverage     |   300             |   kidney     |   3    |
|   MI13-16poolSKM.coverage   |   300             |   skeletal   |   3    |
|   MI106Li.coverage          |   30              |   liver      |   18   |
|   MI107Li.coverage          |   45              |   liver      |   18   |
|   MI108Li.coverage          |   15              |   liver      |   18   |
|   MI109Li.coverage          |   30              |   liver      |   18   |
|   MI111Li.coverage          |   15              |   liver      |   18   |
|   MI113Li.coverage          |   45              |   liver      |   18   |
|   MI118Li.coverage          |   30              |   liver      |   18   |
|   MI119Li.coverage          |   15              |   liver      |   18   |
|   MI120Li.coverage          |   45              |   liver      |   18   |
|   MI122Li.coverage          |   30              |   liver      |   18   |
|   MI125Li.coverage          |   45              |   liver      |   18   |
|   MI127Li.coverage          |   15              |   liver      |   18   |  
 

</details>
 
</details>

<details><summary><b>Main Analysis</b></summary>  
 
Main analysis including statistical analysis, plots, and interactive visualization tools can be accesses by opening ```elongationRate.R``` file in Rstudio and proceeding from there.    
</details>

<details><summary><b>Extras</b></summary>  

Extra bits of analysis that are not part of the manuscript  

### Epas1 upstream reading frames  
Epas-1 gene has one annotated uORF and several putative ORFs base on the 5' UTR sequence.  
Harringtonine treatment highlights UTR regions capable of translation initiation.  



</details>
