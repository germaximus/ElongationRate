# ElongationRate

**Prerequisites:**  
[cutadapt 2.5](https://cutadapt.readthedocs.io/en/stable/index.html)  
[STAR-2.7.2b](https://github.com/alexdobin/STAR)  
[gffread utility](http://ccb.jhu.edu/software/stringtie/gff.shtml)  

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


<details><summary><b>Extract ORF sequences</b></summary>  

**Fetch all mRNA records**  
Run ```mRNA_extractor.pl```. First, it takes ```GRCm38.p6.gbk``` and extracts all RefSeq records for every gene including CDS, 5UTR, 3UTR lengths and a sequence. Then, it selects the single RefSeq record as the longest isoform. Sometimes, the ORF lengths of two isoforms are equal, in that case the longest isoform is selected based on the UTR length with 5UTR taking precedence over 3UTR. The script also trims mRNAs by 100 nucleotides flanking CDS.  If 5UTR and/or 3UTR are shorter than 100 nt, it raises a "flag".  

```bash
 perl mRNA_extractor.pl /path/GRCm38.p6.gbk
 # creates an output file named temp3
```
Fill missing 5UTR and 3UTRs with genomic sequences in cases when they are shorter than 100 nt.  
```bash
perl mRNA_genome_filler.pl 
# requires requires temp3 and individual chromosomal *.fa files in the same folder.
```

</details>



