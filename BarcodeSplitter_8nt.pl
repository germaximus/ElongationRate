#!/usr/bin/perl
open (FILE, "<$ARGV[0]");
open (FILE1, ">Index1.fq");
open (FILE2, ">Index2.fq");
open (FILE3, ">Index3.fq");
open (FILE4, ">Index4.fq");
open (FILE5, ">Index5.fq");
open (FILE6, ">Index6.fq");
open (FILE7, ">Index7.fq");
open (FILE8, ">Index8.fq");
open (FILE9, ">Index9.fq");
open (FILE10, ">Index10.fq");
open (FILE11, ">Index11.fq");
open (FILE12, ">Index12.fq");
open (FILE13, ">noMatch.fq");

my $index1 = "TAAGGCGA";
my $index2 = "CGTACTAG";
my $index3 = "AGGCAGAA";
my $index4 = "TCCTGAGC";
my $index5 = "GGACTCCT";
my $index6 = "TAGGCATG";
my $index7 = "CTCTCTAC";
my $index8 = "CAGAGAGG";
my $index9 = "GCTACGCT";
my $index10 = "GTAGAGGA";
my $index11 = "ATCTCAGG";
my $index12 = "ACTCGCTA";

my %index1;
my %index2;
my %index3;
my %index4;
my %index5;
my %index6;
my %index7;
my %index8;
my %index9;
my %index10;
my %index11;
my %index12;

for (my $i=0; $i < 8; $i++) {
   
   my $mutated_index1 = $index1;
   my $mutated_index2 = $index2;
   my $mutated_index3 = $index3;
   my $mutated_index4 = $index4;
   my $mutated_index5 = $index5;
   my $mutated_index6 = $index6;
   my $mutated_index7 = $index7;
   my $mutated_index8 = $index8;
   my $mutated_index9 = $index9;
   my $mutated_index10 = $index10;
   my $mutated_index11 = $index11;
   my $mutated_index12 = $index12;
   
   substr($mutated_index1,$i,1,'A'); $index1{$mutated_index1} = 0;
   substr($mutated_index1,$i,1,'T'); $index1{$mutated_index1} = 0;
   substr($mutated_index1,$i,1,'G'); $index1{$mutated_index1} = 0;
   substr($mutated_index1,$i,1,'C'); $index1{$mutated_index1} = 0;
   substr($mutated_index1,$i,1,'N'); $index1{$mutated_index1} = 0;
   
   substr($mutated_index2,$i,1,'A'); $index2{$mutated_index2} = 0;
   substr($mutated_index2,$i,1,'T'); $index2{$mutated_index2} = 0;
   substr($mutated_index2,$i,1,'G'); $index2{$mutated_index2} = 0;
   substr($mutated_index2,$i,1,'C'); $index2{$mutated_index2} = 0;
   substr($mutated_index2,$i,1,'N'); $index2{$mutated_index2} = 0;
   
   substr($mutated_index3,$i,1,'A'); $index3{$mutated_index3} = 0;
   substr($mutated_index3,$i,1,'T'); $index3{$mutated_index3} = 0;
   substr($mutated_index3,$i,1,'G'); $index3{$mutated_index3} = 0;
   substr($mutated_index3,$i,1,'C'); $index3{$mutated_index3} = 0;
   substr($mutated_index3,$i,1,'N'); $index3{$mutated_index3} = 0;
   
   substr($mutated_index4,$i,1,'A'); $index4{$mutated_index4} = 0;
   substr($mutated_index4,$i,1,'T'); $index4{$mutated_index4} = 0;
   substr($mutated_index4,$i,1,'G'); $index4{$mutated_index4} = 0;
   substr($mutated_index4,$i,1,'C'); $index4{$mutated_index4} = 0;
   substr($mutated_index4,$i,1,'N'); $index4{$mutated_index4} = 0;
   
   substr($mutated_index5,$i,1,'A'); $index5{$mutated_index5} = 0;
   substr($mutated_index5,$i,1,'T'); $index5{$mutated_index5} = 0;
   substr($mutated_index5,$i,1,'G'); $index5{$mutated_index5} = 0;
   substr($mutated_index5,$i,1,'C'); $index5{$mutated_index5} = 0;
   substr($mutated_index5,$i,1,'N'); $index5{$mutated_index5} = 0;
   
   substr($mutated_index6,$i,1,'A'); $index6{$mutated_index6} = 0;
   substr($mutated_index6,$i,1,'T'); $index6{$mutated_index6} = 0;
   substr($mutated_index6,$i,1,'G'); $index6{$mutated_index6} = 0;
   substr($mutated_index6,$i,1,'C'); $index6{$mutated_index6} = 0;
   substr($mutated_index6,$i,1,'N'); $index6{$mutated_index6} = 0;
   
   substr($mutated_index7,$i,1,'A'); $index7{$mutated_index7} = 0;
   substr($mutated_index7,$i,1,'T'); $index7{$mutated_index7} = 0;
   substr($mutated_index7,$i,1,'G'); $index7{$mutated_index7} = 0;
   substr($mutated_index7,$i,1,'C'); $index7{$mutated_index7} = 0;
   substr($mutated_index7,$i,1,'N'); $index7{$mutated_index7} = 0;
   
   substr($mutated_index8,$i,1,'A'); $index8{$mutated_index8} = 0;
   substr($mutated_index8,$i,1,'T'); $index8{$mutated_index8} = 0;
   substr($mutated_index8,$i,1,'G'); $index8{$mutated_index8} = 0;
   substr($mutated_index8,$i,1,'C'); $index8{$mutated_index8} = 0;
   substr($mutated_index8,$i,1,'N'); $index8{$mutated_index8} = 0;
   
   substr($mutated_index9,$i,1,'A'); $index9{$mutated_index9} = 0;
   substr($mutated_index9,$i,1,'T'); $index9{$mutated_index9} = 0;
   substr($mutated_index9,$i,1,'G'); $index9{$mutated_index9} = 0;
   substr($mutated_index9,$i,1,'C'); $index9{$mutated_index9} = 0;
   substr($mutated_index9,$i,1,'N'); $index9{$mutated_index9} = 0;

   substr($mutated_index10,$i,1,'A'); $index10{$mutated_index10} = 0;
   substr($mutated_index10,$i,1,'T'); $index10{$mutated_index10} = 0;
   substr($mutated_index10,$i,1,'G'); $index10{$mutated_index10} = 0;
   substr($mutated_index10,$i,1,'C'); $index10{$mutated_index10} = 0;
   substr($mutated_index10,$i,1,'N'); $index10{$mutated_index10} = 0;

   substr($mutated_index11,$i,1,'A'); $index11{$mutated_index11} = 0;
   substr($mutated_index11,$i,1,'T'); $index11{$mutated_index11} = 0;
   substr($mutated_index11,$i,1,'G'); $index11{$mutated_index11} = 0;
   substr($mutated_index11,$i,1,'C'); $index11{$mutated_index11} = 0;
   substr($mutated_index11,$i,1,'N'); $index11{$mutated_index11} = 0;

   substr($mutated_index12,$i,1,'A'); $index12{$mutated_index12} = 0;
   substr($mutated_index12,$i,1,'T'); $index12{$mutated_index12} = 0;
   substr($mutated_index12,$i,1,'G'); $index12{$mutated_index12} = 0;
   substr($mutated_index12,$i,1,'C'); $index12{$mutated_index12} = 0;
   substr($mutated_index12,$i,1,'N'); $index12{$mutated_index12} = 0;
   
}


while($line1 = <FILE>){
   chomp $line1;
   $seq = <FILE>;
   $plus = <FILE>;
   $score  = <FILE>;
       
   @line1 = split /\s/, $line1;
   $index = $line1[1];
   $index = substr($index, -8);
 
   if(exists($index1{$index})){
        print FILE1 $line1."\n";
        print FILE1 $seq;
        print FILE1 $plus;
        print FILE1 $score;
   }

   elsif(exists($index2{$index})){
      print FILE2 $line1."\n";
      print FILE2 $seq;
      print FILE2 $plus;
      print FILE2 $score;
   }

   elsif(exists($index3{$index})){
      print FILE3 $line1."\n";
      print FILE3 $seq;
      print FILE3 $plus;
      print FILE3 $score;
   }
   
   elsif(exists($index4{$index})){
      print FILE4 $line1."\n";
      print FILE4 $seq;
      print FILE4 $plus;
      print FILE4 $score;
   }
   elsif(exists($index5{$index})){
      print FILE5 $line1."\n";
      print FILE5 $seq;
      print FILE5 $plus;
      print FILE5 $score;
   }
   elsif(exists($index6{$index})){
      print FILE6 $line1."\n";
      print FILE6 $seq;
      print FILE6 $plus;
      print FILE6 $score;
   }
   elsif(exists($index7{$index})){
      print FILE7 $line1."\n";
      print FILE7 $seq;
      print FILE7 $plus;
      print FILE7 $score;
   }
   elsif(exists($index8{$index})){
      print FILE8 $line1."\n";
      print FILE8 $seq;
      print FILE8 $plus;
      print FILE8 $score;
   }
   elsif(exists($index9{$index})){
      print FILE9 $line1."\n";
      print FILE9 $seq;
      print FILE9 $plus;
      print FILE9 $score;
   }
   elsif(exists($index10{$index})){
      print FILE10 $line1."\n";
      print FILE10 $seq;
      print FILE10 $plus;
      print FILE10 $score;
   }
   elsif(exists($index11{$index})){
      print FILE11 $line1."\n";
      print FILE11 $seq;
      print FILE11 $plus;
      print FILE11 $score;
   }
   elsif(exists($index12{$index})){
      print FILE12 $line1."\n";
      print FILE12 $seq;
      print FILE12 $plus;
      print FILE12 $score;
   }
   else{
      print FILE13 $line1."\n";
      print FILE13 $seq;
      print FILE13 $plus;
      print FILE13 $score;
   }
}
close(FILE);
close(FILE1);
close(FILE2);
close(FILE3);
close(FILE4);
close(FILE5);
close(FILE6);
close(FILE7);
close(FILE8);
close(FILE9);
close(FILE10);
close(FILE11);
close(FILE12);
close(FILE13);