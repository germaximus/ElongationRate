#!/usr/bin/perl
# calculates coverage for every gene in the bowtie output.
# Usage: perl Coverage.pl [index] [input] [output]
# Where:
# [index] - mRNA_100uniq.fasta index file used to create Bowtie index. Better put it in the same folder with the script.
# [input] - uniq.bwt file containing Bowtie output in native format
# [output] - "coverage" filename with or w/o path

my %genes;
my @gene_order;

#initialize index hash. Due to memory requirements I have to run it on a server
open (INDEX, "<$ARGV[0]") || die "no index file found";
while ($line = <INDEX>)
  {  chomp($line);
     $name = substr($line,1);
     $line = <INDEX>; chomp($line);
     $length = length($line);
     push(@gene_order, $name);
     for(my $i=0; $i < $length; $i++)    {   $genes{$name}[$i] = 0;    }
  }
close (INDEX);


open (INPUT, "<$ARGV[1]"); open (OUT, ">$ARGV[2]"); 
while ($line = <INPUT>)
  {
     @line = split /\t/,$line;
     
     for(my $i = $line[3]; $i < ($line[3] + length($line[4])); $i++)
       {
         $genes{$line[2]}[$i]++;
       }
  }
foreach $name (@gene_order)
   {
      print OUT "$name\t";
      print OUT join("\t",@{$genes{$name}})."\n";
   }
close(INPUT); close(OUT);
