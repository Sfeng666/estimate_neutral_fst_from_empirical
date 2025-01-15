#!/usr/bin/perl -w
#use strict;

# my @chrs = ('chrX','chr2L','chr2R','chr3L','chr3R');
my @chrs = ('chrX');

my $PopA = ('FR');
my $PopB = ('ZI');

my $PopATargetN = 18;
my $PopBTargetN = 18;

my $OutputFile = 'SNP_FSTRey_Dist_FR_ZI.txt';

#declarations
my $c = 0;
my $i = 0;
my $j = 0;
my $site = 0;
my $CountFileA = '';
my $CountFileB = '';
my $PopANumAlleles = 0;
my $PopBNumAlleles = 0;
my $TotalNumAlleles = 0;
my $PopACurrentN = 0;
my $PopBCurrentN = 0;
my $random = 0;
my $MajorAllele = 0;
my $MinorAllele = 0;
my $MajorFreq1 = 0;
my $MajorFreq2 = 0;
my $MinorFreq1 = 0;
my $MinorFreq2 = 0;
my $SharedNum = 0;
my $NumA = 0;
my $FracNum = 0;
my $FracDen = 0;
my $frac = 0;
my $WholeNum = 0;
my $DenFracNum = 0;
my $DenFrac = 0;
my $WholeDen = 0;
my $SampleSize1 = $PopATargetN;
my $SampleSize2 = $PopBTargetN;
my $FST = 0;

my @FSTDistAoA = ();
my @line = ();
my @AlleleCounts = ();
my @SampledCounts = ();

#Make FST output array of arrays
for ($i = 0; $i < 1001; $i++){
  push @line, 0;
}
for ($c = 0; $c < @chrs; $c++){
  push @FSTDistAoA, [ @line ];
}

# Open output files for each chromosome arm to store FST values
my %out_files;
foreach my $chr (@chrs) {
    open($out_files{$chr}, '>', "fst_$chr.txt") or die "Cannot open output file for $chr: $!";
}

#chromosome arm loop
for ($c = 0; $c < @chrs; $c++){
  $CountFileA = 'AllCounts_' . $PopA . '_NoInv_' . $chrs[$c] . '_rmlowrecomb_syn_seg.txt';
  open A, "<../../data/output/$CountFileA" or die "can not open $CountFileA\n";
  $CountFileB = 'AllCounts_' . $PopB . '_NoInv_' . $chrs[$c] . '_rmlowrecomb_syn_seg.txt';
  open B, "<../../data/output/$CountFileB" or die "can not open $CountFileB\n";
  $site = 0;
  
#Read in data for one row/site from both input count files
  while (<A>){
    @AlleleCounts = ();
    chomp;
    @line = split;
    push @AlleleCounts, @line;
    $_ = (<B>);
    chomp;
    @line = split;
    push @AlleleCounts, @line;
    last if (@AlleleCounts != 8);
    $site++;

#Determine whether sufficient sample size and initially variable (otherwise skip to next site)
    $PopACurrentN = $AlleleCounts[0] + $AlleleCounts[1] + $AlleleCounts[2] + $AlleleCounts[3];
    $PopBCurrentN = $AlleleCounts[4] + $AlleleCounts[5] + $AlleleCounts[6] + $AlleleCounts[7];
    next if ($PopACurrentN < $PopATargetN);
    next if ($PopBCurrentN < $PopBTargetN);
    for ($i = 0; $i < 3.5; $i++){
      if (($AlleleCounts[$i] + $AlleleCounts[$i+4]) > 0){
    $TotalNumAlleles++;
      }
    }
    next if ($TotalNumAlleles < 1.5);

# skip down-sampling and let SampledCounts equal AlleleCounts
    @SampledCounts = @AlleleCounts;

# #Down-sample and confirm if we sampled two alleles (otherwise skip to next site)
#     @SampledCounts = (0,0,0,0,0,0,0,0);
#     for ($i = 0; $i < $PopATargetN; $i++){
#       $random = int(rand($PopACurrentN));
#       if ($random < $AlleleCounts[0]){
#     $SampledCounts[0]++;
#     $AlleleCounts[0]--;
#     $PopACurrentN--;
#       }
#       elsif ($random < ($AlleleCounts[0] + $AlleleCounts[1])){
#     $SampledCounts[1]++;
#     $AlleleCounts[1]--;
#     $PopACurrentN--;
#       }
#       elsif ($random < ($AlleleCounts[0] + $AlleleCounts[1] + $AlleleCounts[2])){
#     $SampledCounts[2]++;
#     $AlleleCounts[2]--;
#     $PopACurrentN--;
#       }
#       else{
#     $SampledCounts[3]++;
#     $AlleleCounts[3]--;
#     $PopACurrentN--;
#       }
#     }
#     for ($i = 0; $i < $PopBTargetN; $i++){
#       $random = int(rand($PopBCurrentN));
#       if ($random < $AlleleCounts[4]){
#     $SampledCounts[4]++;
#     $AlleleCounts[4]--;
#     $PopBCurrentN--;
#       }
#       elsif ($random < ($AlleleCounts[4] + $AlleleCounts[5])){
#     $SampledCounts[5]++;
#     $AlleleCounts[5]--;
#     $PopBCurrentN--;
#       }
#       elsif ($random < ($AlleleCounts[4] + $AlleleCounts[5] + $AlleleCounts[6])){
#     $SampledCounts[6]++;
#     $AlleleCounts[6]--;
#     $PopBCurrentN--;
#       }
#       else{
#     $SampledCounts[7]++;
#     $AlleleCounts[7]--;
#     $PopBCurrentN--;
#       }
#     }
#     $TotalNumAlleles = 0;
#     for ($i = 0; $i < 3.5; $i++){
#       if (($SampledCounts[$i] + $SampledCounts[$i+4]) > 0){
#     $TotalNumAlleles++;
#       }
#     }
#     next if ($TotalNumAlleles != 2);
      
#Define major and minor allele frequencies
    $MajorAllele = -1;
    $MinorAllele = -1;
    for ($i = 0; $i < 3.5; $i++){
      if ((($SampledCounts[$i]+$SampledCounts[$i+4]) / ($PopATargetN + $PopBTargetN)) >= 0.5){
    $MajorAllele = $i;
    last;
      }
    }
    die if ($MajorAllele == -1);
    for ($i = 0; $i < 3.5; $i++){
      next if ($i == $MajorAllele);
      if (($SampledCounts[$i]+$SampledCounts[$i+4]) > 0){
    $MinorAllele = $i;
    last;
      }
    }
    die if ($MinorAllele == -1);
    # next if (($SampledCounts[$MinorAllele] + $SampledCounts[$MinorAllele+4]) < 1.5);
    $MajorFreq1 = $SampledCounts[$MajorAllele] / $PopATargetN;
    $MajorFreq2 = $SampledCounts[$MajorAllele+4] / $PopBTargetN;
    $MinorFreq1 = $SampledCounts[$MinorAllele] / $PopATargetN;
    $MinorFreq2 = $SampledCounts[$MinorAllele+4] / $PopBTargetN;

    $MajorFreq1 = 1;
    $MajorFreq2 = 17/18;
    $MinorFreq1 = 0;
    $MinorFreq2 = 1/18;
    $SampleSize1 = 18;
    $SampleSize2 = 18;

#Get Reynolds FST
    $SharedNum = $SampleSize1 * ($MinorFreq1 -  ($MinorFreq1 ** 2)) + $SampleSize2 * ($MinorFreq2 - ($MinorFreq2 ** 2));
    # $NumA = ($MinorFreq1 - $MinorFreq2) ** 2;
    $NumA = (($MinorFreq1 - $MinorFreq2) ** 2 + ($MajorFreq1 - $MajorFreq2) ** 2)/2;
    $FracNum = (($SampleSize1 + $SampleSize2)/2) * $SharedNum;
    $FracDen = $SampleSize1 * $SampleSize2 * ((($SampleSize1 + $SampleSize2)/2) - 1);
    $frac = $FracNum / $FracDen;
    $WholeNum = $NumA - $frac;
    $DenFracNum = ($SampleSize1 * $SampleSize2 - ($SampleSize1 + $SampleSize2)/2) * $SharedNum;
    $DenFrac = $DenFracNum / $FracDen;
    $WholeDen = $NumA + $DenFrac;
    if ($WholeDen != 0){
      $FST = $WholeNum / $WholeDen;
    }
    else{
      $FST = 0;
    }

#Add to output array
    # if ($FST < 0){
    #   $FST = 0;
    # }
    # Write FST value to the corresponding output file
    my $out_fh = $out_files{$chrs[$c]};
    print $out_fh "$FST\n";
    print "On $chrs[$c] at site $site, FST is $FST, MajorFreq1 is $MajorFreq1,  MinorFreq1 is $MinorFreq1, MajorFreq2 is $MajorFreq2, MinorFreq2 is $MinorFreq2\n";

    $FST = int($FST * 1000);
    $FSTDistAoA[$c][$FST]++;
    # print "On $chrs[$c] at site $site, FST is 0.$FST\n";


  }	
  close A;
  close B;
}

# Close all output files
foreach my $chr (@chrs) {
    close($out_files{$chr});
}

#Output
open O, ">$OutputFile" or die;
print O "FSTBinMax";
for ($c = 0; $c < @chrs; $c++){
  print O "\t$chrs[$c]";
}
for ($i = 0; $i < @{$FSTDistAoA[0]}; $i++){
  $j = $i * 0.001;
  print O "\n$j";
  for ($c = 0; $c < @chrs; $c++){
    print O "\t$FSTDistAoA[$c][$i]";
  }
}
print O "\n";
close O;