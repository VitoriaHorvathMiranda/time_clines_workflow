#!/usr/bin/perl -w
use strict;
use Getopt::Std;

#Get command line arguments
our($opt_f, $opt_w, $opt_s, $opt_a, $opt_r, $opt_b);
getopt('fwsarb');

unless (defined($opt_f) && defined($opt_w) && defined($opt_s) && defined($opt_a) && defined($opt_r) && defined($opt_b)){
    print <<EOF;

GO enrichment analysis

It requires the following arguments:
-f [region file name]
-w [window statistic file prefix]
-s [statistic column to include in window stat file]
-a [exon file prefix]
-r [number of replicates to simulate]
-b [batch number, optional to allow parallel processing]

EOF
exit 1;
}

#Get command line arguments
my $OutliersFile = $opt_f;
my $WinStatFilePrefix = $opt_w;
my $StatColumn = $opt_s;
my $ExonFilePrefix = $opt_a; #assumes full name of annotation file is something like annot_chr2L.txt
my $BatchNum = $opt_b;
my $PermutationReps = $opt_r;

#GO enrichment analysis using all genes overlapping each outlier region center

my $GenesForAllWindowsFile = 'windows_overlap_genes.txt'; #will use if it exists, or else create it to save time in the future

my $TotalWinColumn = 3; #which column from outlier region file gives the total number of windows in each outlier region
my $FbgnColumn = 4;  #which column from outlier region file are the Fbgns in, starting from column 0

my $InputHeaderRows = 1;  #how many header rows are present in the input file before the data
my $MissingDataCode = 'NA'; # originally -999

my $GeneGOCatFile = '../resources/GO_gene_cats_parents19052025.txt';

my $GOCatDescFile = '../resources/GO_desc12_condensed_parents.txt';

my $OutputFile = $OutliersFile;
$OutputFile =~ s/\.tsv//;
$OutputFile = $OutputFile . '_GO' . $BatchNum . '.txt';

my @chrs = ('2L','2R','3L','3R','X');  #keep in alphatbetically sorted order


my $c = 0;
my $f = 0;
my $i = 0;
my $file = '';
my $center = 0;
my @line = ();
my @ChrList = ();
#my @PosList = ();
my @CDSAoA = ();
my @AllOverlapGenes = ();
my @AllOverlapFbgns = ();
my @AllWinStarts = ();
my @AllWinStops = ();
my @AllOverlapFbgnAoA = ();
my @OutlierRegionWindowLengths = ();

my @GeneOrientationList = ();
my @GenePositionList = ();
my $j = 0;
my $k = 0;
my $g = 0;
my $h = 0;
#my $NearestExon = -1;
my $NearestDistance = 1000000000;
my $distance = 0;
#my @GenesWithin = ();
my @UpstreamGenes = ();
my @DownstreamGenes = ();
my @OverlapGenes = ();
my $AlreadyListed = 0;
my $orientation = "";
my $GenePosition = 0;
my $UpstreamLimit = 0;
my $DownstreamLimit = 0;
my $direction = 0;
my $GeneList = "";
my @OverlapFbgns = ();
my $FbgnList = "";
my $ClosestGene = "";
my $ClosestDist = 100000000;
my $start = 0;

#Test if GenesForAllWindowsFile exists.  If so, read it.  If not, evaluate genes for each window.
if (-e $GenesForAllWindowsFile){
  open C, "<$GenesForAllWindowsFile" or die "can not open $GenesForAllWindowsFile\n";
  while (<C>){
    chomp;
    last if m/^$/;
    @line = split;
    push @ChrList, $line[0];
    push @AllWinStarts, $line[1];
    push @AllWinStops, $line[2];
    splice @line, 0, 3;
    push @AllOverlapFbgnAoA, [ @line ];
#    push @AllOverlapGeneAoA, $line[3];
  }
  close C;
print "Obtained genes associated with each window from file\n";
}
else{

#For each chromosome arm, get locations of all windows starts and stops.  Skip windows with missing data
for ($c = 0; $c < @chrs; $c++){
  $file = $WinStatFilePrefix . $chrs[$c] . '_.tsv';
  open W, "<$file" or die "can not open windows stat file $file\n";
  for ($i = 0; $i < $InputHeaderRows; $i++){
    scalar (<W>);
  }
  while (<W>){
    chomp;
    last if m/^$/;
    @line = split;
    next if ($line[$StatColumn] eq $MissingDataCode);
    push @ChrList, $line[0];
    push @AllWinStarts, $line[1];
    push @AllWinStops, $line[2];
  }
  close W;

#Read in all coding sequence blocks
  $file = $ExonFilePrefix . $chrs[$c] . '.txt';
  open E, "<$file" or die "can not open exon file $file\n";
  while (<E>){
    chomp;
    last if m/^$/;
    @line = split;
    $line[2] =~ s/c/C/;
    push @CDSAoA, [ @line ];
  }
  close E;
}

#For each window's start and stop positions, find all overlapping Fbgns.  
#If >1, write in the form: "FBgn0025378/FBgn0025391".
#Obtain all genes that overlap the full outlier region, plus the next exon to each side
for ($i = 0; $i < @AllWinStarts; $i++){
    @OverlapFbgns = ();
    @OverlapGenes = ();
    print "Evaluating closest gene for window $i.  $ChrList[$i]:$AllWinStarts[$i]\n";
   for ($j = $start; $j < @CDSAoA; $j++){
#Make sure comparing same arm
      next if ($CDSAoA[$j][2] lt $ChrList[$i]);
      last if ($CDSAoA[$j][2] gt $ChrList[$i]);
      next if ($AllWinStarts[$i] > $CDSAoA[$j][4]);
#if we hit our first exon that starts after the region start, go back and get genes with exons that overlap this boundary, or else the closest exon before it
      if (@OverlapFbgns == 0){
	$start = $j;
	$ClosestGene = -1;
	$ClosestDist = 1000000000;
	for ($k = $j - 1; $k >= 0; $k--){
	  if (($AllWinStarts[$i] < $CDSAoA[$k][5]) && ($CDSAoA[$k][2] eq $ChrList[$i])){
	    for ($g = 0; $g < @OverlapFbgns; $g++){
	      if ($CDSAoA[$k][0] eq $OverlapFbgns[$g]){
		splice @OverlapFbgns, $g, 1;
		splice @OverlapGenes, $g, 1;
		last;
	      }
	    }
	    push @OverlapFbgns, $CDSAoA[$k][0];
	    push @OverlapGenes, $CDSAoA[$k][1];
	    next;
	  }
	  if ( (($AllWinStarts[$i] - $CDSAoA[$k][5]) < $ClosestDist) && ($CDSAoA[$k][2] eq $ChrList[$i])){
	    $ClosestGene = $k;
	    $ClosestDist = $AllWinStarts[$i] - $CDSAoA[$k][5];
	  }
	  if ( ($k <= ($j-500)) || ($k == 0) || ($CDSAoA[$k][2] lt $ChrList[$i]) ){
	    if ( (@OverlapFbgns == 0) && ($ClosestGene > -0.5)){
	       push @OverlapFbgns, $CDSAoA[$ClosestGene][0];
	       push @OverlapGenes, $CDSAoA[$ClosestGene][1];
	     }
	    last;
	  }
	}
	$ClosestGene = -1;
	$ClosestDist = 1000000000;
      }
#Gather all genes with exons that fall completely within the outlier region
      if ($CDSAoA[$j][5] < $AllWinStops[$i]){
	for ($g = 0; $g < @OverlapFbgns; $g++){
	  if ($CDSAoA[$j][0] eq $OverlapFbgns[$g]){
	    splice @OverlapFbgns, $g, 1;
	    splice @OverlapGenes, $g, 1;
	    last;
	  }
	}
	push @OverlapFbgns, $CDSAoA[$j][0];
	push @OverlapGenes, $CDSAoA[$j][1];
	next;
      }
#when we reach the first exon that starts after the region stop, include this gene if there were none overlapping the outlier region stop, then terminate the loop
      elsif ($CDSAoA[$j][4] > $AllWinStops[$i]){
	if ($ClosestDist > 0){
	  for ($g = 0; $g < @OverlapFbgns; $g++){
	    if ($CDSAoA[$j][0] eq $OverlapFbgns[$g]){
	      splice @OverlapFbgns, $g, 1;
	      splice @OverlapGenes, $g, 1;
	      last;
	    }
	  }
	  push @OverlapFbgns, $CDSAoA[$j][0];
	  push @OverlapGenes, $CDSAoA[$j][1];
	}
	last;
      }
#include genes with exons that overlap the region stop
      else{
	for ($g = 0; $g < @OverlapFbgns; $g++){
	  if ($CDSAoA[$j][0] eq $OverlapFbgns[$g]){
	    splice @OverlapFbgns, $g, 1;
	    splice @OverlapGenes, $g, 1;
	    last;
	  }
	}
	push @OverlapFbgns, $CDSAoA[$j][0];
	push @OverlapGenes, $CDSAoA[$j][1];
	$ClosestDist = 0;
      }
    }
###
#print "Genes found: @OverlapFbgns\n";
###
    push @AllOverlapFbgnAoA, [ @OverlapFbgns ];
#    push @AllOverlapGeneAoA, [ @OverlapGenes ];
  }
#write this information to file, so we can skip the above analysis if we run GO enrichment for these windows again
open C, ">$GenesForAllWindowsFile";
for ($i = 0; $i < @ChrList; $i++){
  print C "$ChrList[$i]\t$AllWinStarts[$i]\t$AllWinStops[$i]";
  for ($j = 0; $j < @{$AllOverlapFbgnAoA[$i]}; $j++){
    print C "\t$AllOverlapFbgnAoA[$i][$j]";
  }
  print C "\n";
}
close C;
print "Recorded genes associated with each window to file (for future reference).\n";
}

#Transfer to new arrays for GO enrichment analysis
my $l = 0;
my $match = 0;
my @AllFbgnsRepeats = ();
my @AllFbgnsUnique = ();
my @AllGenesUnique = ();
my @GOpreAoA = ();
my @GOAoA = ();
my @AllGORepeats = ();
my @AllGOUnique = ();
my @AllGOCounts = ();
my @OutlierFbgnAoA = ();
my @RegionGOCats = ();
my @OutlierGOList = ();
my @OutlierGOCounts = ();
for ($i = 0; $i < @AllOverlapFbgnAoA; $i++){
  for ($j = 0; $j < @{$AllOverlapFbgnAoA[$i]}; $j++){
    push @AllFbgnsRepeats, $AllOverlapFbgnAoA[$i][$j];
  }
}
@AllFbgnsRepeats = sort{ lc($a) cmp lc($b) } @AllFbgnsRepeats;
$i = @AllFbgnsRepeats;
print "AllFbgnsRepeats has $i\n";
###
#for ($i = 0; $i < @AllFbgnsRepeats; $i++){
#  print "$AllFbgnsRepeats[$i]\n";
#}
###
push @AllFbgnsUnique, $AllFbgnsRepeats[0];
for ($i = 1; $i < @AllFbgnsRepeats; $i++){
  if ($AllFbgnsRepeats[$i] ne $AllFbgnsRepeats[$i-1]){
    push @AllFbgnsUnique, $AllFbgnsRepeats[$i];
  }
}
#@AllFbgnsUnique = @AllFbgnsRepeats;
#for ($i = 1; $i < @AllFbgnsUnique; $i++){
#  if ($AllFbgnsUnique[$i] eq $AllFbgnsUnique[$i-1]){
#    splice @AllFbgnsUnique, $i, 1;
#    $i--;
#  }
#}
$i = @AllFbgnsUnique;
print "AllFbgnsUnique has $i\n";

#Read GO category list into pre-AoA
open G, "<$GeneGOCatFile" or die "Can't open $GeneGOCatFile\n";
while (<G>){
  chomp;
  last if m/^$/;
  @line = split;
  for ($i = 0; $i < @AllFbgnsUnique; $i++){
    if ($line[0] eq $AllFbgnsUnique[$i]){
      push @GOpreAoA, [ @line ];
      last;
    }
  }
}
close G;
$i = @GOpreAoA;
print "Found GO categories for $i genes\n";

#Build a GO AoA that corresponds to indexing of AllFbgnsUnique
for ($i = 0; $i < @AllFbgnsUnique; $i++){
  $match = 0;
  for ($j = 0; $j < @GOpreAoA; $j++){
    if ($AllFbgnsUnique[$i] eq $GOpreAoA[$j][0]){
      push @AllGenesUnique, $GOpreAoA[$j][1];
      @line = ();
      for ($k = 2; $k < @{$GOpreAoA[$j]}; $k++){
	push @line, $GOpreAoA[$j][$k];
#	$AlreadyListed = 0;
#	for ($l = 1; $l < @line; $l++){
#	  if ($GOpreAoA[$j][$k] eq $line[$l]){
#	    $AlreadyListed = 1;
#	    last;
#	  }
#	}
#	if ($AlreadyListed == 0){
#	  push @line, $GOpreAoA[$j][$k];
#	}
      }
      push @GOAoA, [ @line ];
      $match = 1;
      last;
    }
  }
  if ($match == 0){
    @line = ();
    push @line, 'GO:9999999';
    push @AllGenesUnique, '-';
    push @GOAoA, [ @line ];
  }
}
#
if (@GOAoA != @AllFbgnsUnique){
  $i = @AllFbgnsUnique;
  print "Error:  AllFbgnsUnique has $i entries,";
  $i = @GOAoA;
  print "but GOAoA has $i entries\n";
  die;
}
if (@AllGenesUnique != @AllFbgnsUnique){
  $i = @AllFbgnsUnique;
  print "Error:  AllFbgnsUnique has $i entries,";
  $i = @AllGenesUnique;
  print "but AllGenesUnique has $i entries\n";
  die;
}

#Make arrays with all GO terms (first with repeats, then two more arrays with unique terms and their corresponding counts)
for ($i = 0; $i < @GOAoA; $i++){
  for ($j = 0; $j <  @{$GOAoA[$i]}; $j++){
    push @AllGORepeats, $GOAoA[$i][$j];
  }
}
@AllGORepeats = sort{ lc($a) cmp lc($b) } @AllGORepeats;
$i = @AllGORepeats;
print "Set of all genes has a total of $i GO listings,\n";

push @AllGOUnique, $AllGORepeats[0];
push @AllGOCounts, 1;
for ($i = 1; $i < @AllGORepeats; $i++){
  if ($AllGORepeats[$i] eq $AllGORepeats[$i-1]){
    $AllGOCounts[-1]++;
  }
  else{
    push @AllGOUnique, $AllGORepeats[$i];
    push @AllGOCounts, 1;
  }
}
$i = @AllGOUnique;
print "representing $i distinct GO categories.\n";

#Read in the list of outlier genes
#make OutlierFbgnAoA (one row of genes per outlier region) so we can count each GO category max 1 time per outlier region
$k = 0;
open F, "<$OutliersFile" or die "can't open $OutliersFile\n";
scalar (<F>);
while (<F>){
  chomp;
  last if m/^$/;
  @line = split;
  push @OutlierRegionWindowLengths, $line[$TotalWinColumn];
  $_ = $line[$FbgnColumn];
  @line = split(/\//);
  push @OutlierFbgnAoA, [ @line ];
  $k += @line;
###
#  print "@line\n";
###
}
close F;
$i = @OutlierFbgnAoA;
print "Found $k outlier genes for $i regions,\n";

#Eliminate duplicate genes in OutlierFbgnAoA (remove from the outlier region containing more genes)
for ($i = 0; $i < @OutlierFbgnAoA; $i++){
  for ($j = 0; $j < @{$OutlierFbgnAoA[$i]}; $j++){
    $match = 0;
    for ($k = $i + 1; $k < @OutlierFbgnAoA; $k++){
      for ($l = 0; $l < @{$OutlierFbgnAoA[$k]}; $l++){
	if ($OutlierFbgnAoA[$i][$j] eq $OutlierFbgnAoA[$k][$l]){
	  if (@{$OutlierFbgnAoA[$i]} < @{$OutlierFbgnAoA[$k]}){
	    @line = @{$OutlierFbgnAoA[$k]};
	    splice @line, $l, 1;
	    @{$OutlierFbgnAoA[$k]} = @line;
	    last;
	  }
	  else{
	    @line = @{$OutlierFbgnAoA[$i]};
	    splice @line, $j, 1;
	    @{$OutlierFbgnAoA[$i]} = @line;
	    $match = 1;
	    last;
	  }
	}
      }
      if ($match == 1){
	$j--;
	last;
      }
    }
  }
}
$k = 0;
for ($i = 0; $i < @OutlierFbgnAoA; $i++){
  $k += @{$OutlierFbgnAoA[$i]};
}
#print "After removing duplicates of genes between outlier regions, $k genes remain\n";

#For each region, obtain the list of GO categories associated with overlapping genes.  
#Eliminate duplicate GO categories in same outlier region (avoid false positive GO enrichment results due to clusters of functionally related genes)
for ($i = 0; $i < @OutlierFbgnAoA; $i++){
  @RegionGOCats = ();
  for ($j = 0; $j < @{$OutlierFbgnAoA[$i]}; $j++){
    for ($k = 0; $k < @GOAoA; $k++){
      if ($AllFbgnsUnique[$k] eq $OutlierFbgnAoA[$i][$j]){
	for ($l = 0; $l < @{$GOAoA[$k]}; $l++){
	  push @RegionGOCats, $GOAoA[$k][$l];
	}
      }
    }
  }
  @RegionGOCats  = sort{ lc($a) cmp lc($b) } @RegionGOCats;
  for ($j = 1; $j < @RegionGOCats; $j++){
    if ($RegionGOCats[$j] eq $RegionGOCats[$j-1]){
      splice @RegionGOCats, $j, 1;
      $j--;
    }
  }
  push @OutlierGOList, @RegionGOCats;
#  push @OutlierGOMatrix, [ @RegionGOCats ];
}
$i = @OutlierGOList;
print "with a total of $i GO listings\n";
for ($i = 0; $i < @AllGOUnique; $i++){
  push @OutlierGOCounts, 0;
}
@OutlierGOList = sort{ lc($a) cmp lc($b) } @OutlierGOList;    
$j = 0;
for ($i = 0; $i < @OutlierGOList; $i++){
  while ($OutlierGOList[$i] ne $AllGOUnique[$j]){
    $j++;
  }
  $OutlierGOCounts[$j]++;
}
print "Empirical GO category counts established.\n";

#PERMUTATION - Sample random outlier regions for p values
my $r = 0;
my $random = 0;
my $FirstWin = 0;
my $LastWin = 0;
my $WinLength = 0;
my @RegionFbgns = ();
my @SampledFbgnAoA = ();
my @SampledGOList = ();
my @SampledGOCounts = ();
my @PValues = ();
for ($i = 0; $i < @AllGOUnique; $i++){
  push @PValues, 0;
}

for ($r = 0; $r < $PermutationReps; $r++){
  @SampledFbgnAoA = ();
  @SampledGOList = ();
  @SampledGOCounts = ();  
#get sampled genes from a randomly located outlier region of random window length (based on empirical lengths)
  while (@SampledFbgnAoA < @OutlierFbgnAoA){
    $FirstWin = int(rand @AllOverlapFbgnAoA);
    $random = int(rand @OutlierRegionWindowLengths);
    $WinLength = $OutlierRegionWindowLengths[$random];
    @RegionFbgns = ();
    $LastWin = ($FirstWin + $WinLength) - 1;
    next if (( !(defined($ChrList[$LastWin]))) || ($ChrList[$FirstWin] ne $ChrList[$LastWin]));
    for ($i = $FirstWin; $i <= $LastWin; $i++){
      @line = @{$AllOverlapFbgnAoA[$i]};
      push @RegionFbgns, @line;
    }
    @RegionFbgns =  sort{ lc($a) cmp lc($b) } @RegionFbgns;
    for ($i = 1; $i < @RegionFbgns; $i++){
      if ($RegionFbgns[$i] eq $RegionFbgns[$i-1]){
	splice @RegionFbgns, $i, 1;
	$i--;
      }
    }
    push @SampledFbgnAoA, [ @RegionFbgns ];
  }
#Eliminate duplicate genes in SampledFbgnAoA (remove from the outlier region containing more genes)
  for ($i = 0; $i < @SampledFbgnAoA; $i++){
    for ($j = 0; $j < @{$SampledFbgnAoA[$i]}; $j++){
      $match = 0;
      for ($k = $i + 1; $k < @SampledFbgnAoA; $k++){
	for ($l = 0; $l < @{$SampledFbgnAoA[$k]}; $l++){
	  if ($SampledFbgnAoA[$i][$j] eq $SampledFbgnAoA[$k][$l]){
	    if (@{$SampledFbgnAoA[$i]} < @{$SampledFbgnAoA[$k]}){
	      @line = @{$SampledFbgnAoA[$k]};
	      splice @line, $l, 1;
	      @{$SampledFbgnAoA[$k]} = @line;
	      last;
	    }
	    else{
	      @line = @{$SampledFbgnAoA[$i]};
	      splice @line, $j, 1;
	      @{$SampledFbgnAoA[$i]} = @line;
	      $match = 1;
	      last;
	    }
	  }
	}
	if ($match == 1){
	  $j--;
	  last;
	}
      }
    }
  }

#get GO categories from sampled regions
  for ($i = 0; $i < @SampledFbgnAoA; $i++){
    @RegionGOCats = ();
    for ($j = 0; $j < @{$SampledFbgnAoA[$i]}; $j++){
      for ($k = 0; $k < @GOAoA; $k++){
	if ($AllFbgnsUnique[$k] eq $SampledFbgnAoA[$i][$j]){
	  for ($l = 0; $l < @{$GOAoA[$k]}; $l++){
	    push @RegionGOCats, $GOAoA[$k][$l];
	  }
	}
      }
    }
    @RegionGOCats  = sort{ lc($a) cmp lc($b) } @RegionGOCats;
    for ($j = 1; $j < @RegionGOCats; $j++){
      if ($RegionGOCats[$j] eq $RegionGOCats[$j-1]){
	splice @RegionGOCats, $j, 1;
	$j--;
      }
    }
    push @SampledGOList, @RegionGOCats;
  }

  for ($i = 0; $i < @AllGOUnique; $i++){
    push @SampledGOCounts, 0;
  }
  @SampledGOList = sort{ lc($a) cmp lc($b) } @SampledGOList;    
  $j = 0;
  for ($i = 0; $i < @SampledGOList; $i++){
    while ($SampledGOList[$i] ne $AllGOUnique[$j]){
      $j++;
    }
    $SampledGOCounts[$j]++;

  }

  
#Use a P-value array to add 1/reps each time the resampled data has more counts for a GO category than the empirical data did.
  for ($i = 0; $i < @OutlierGOCounts; $i++){
    if ($SampledGOCounts[$i] >= $OutlierGOCounts[$i]){
      $PValues[$i] += (1/$PermutationReps);
    }
  }
  print "Done with resampled data set $r\n";
}

#For each GO category with one or more outliers, look up the common gene names of the outliers
my $GOCatOutliers = '';
my @GOCatOutlierList = ();
my @OutlierFbgns = ();
my @OutlierGenes = ();
my @OutlierGOAoA = ();
for ($i = 0; $i < @OutlierFbgnAoA; $i++){
  for ($j = 0; $j < @{$OutlierFbgnAoA[$i]}; $j++){
    push @OutlierFbgns, $OutlierFbgnAoA[$i][$j];
  }
}
for ($i = 0; $i < @OutlierFbgns; $i++){
  for ($j = 0; $j < @AllFbgnsUnique; $j++){
    if ($OutlierFbgns[$i] eq $AllFbgnsUnique[$j]){
      push @OutlierGenes, $AllGenesUnique[$j];
      @line = @{$GOAoA[$j]};
      push @OutlierGOAoA, [ @line ];
      last;
    }
    if ($j == (@AllFbgnsUnique - 1)){
      print "could not find $OutlierFbgns[$i] in AllFbgnsUnique\n";
    }
  }
}
$i = @OutlierFbgns;
$j = @OutlierGOAoA;
if ($i != $j){
  die "wrong number of entries in OutlierGOAoA ($j instead of $i)\n";
}
for ($i = 0; $i < @OutlierGOCounts; $i++){
  $GOCatOutliers = '';
  if ( ($OutlierGOCounts[$i] > 0) && ($i < (@OutlierGOCounts - 1))){
    for ($j = 0; $j < @OutlierGOAoA; $j++){
      for ($k = 0; $k < @{$OutlierGOAoA[$j]}; $k++){
	if ($AllGOUnique[$i] eq $OutlierGOAoA[$j][$k]){
	  next if ( defined($GOCatOutliers) && ($GOCatOutliers =~ m/$OutlierGenes[$j]/));
	  if ( defined($GOCatOutliers) && (length($GOCatOutliers)) > 0 ){
	    $GOCatOutliers = $GOCatOutliers . ',' . $OutlierGenes[$j];
	  }
	  else{
	    $GOCatOutliers = $OutlierGenes[$j];
	  }
	}
      }
    } 
    push @GOCatOutlierList, $GOCatOutliers;
  }
  else{
    $_ = '';
    push @GOCatOutlierList, $_;
  }
}

#Look up descriptions for each GO category
my @ontologies = ();
my @descriptions = ();
my @GODescAoA = ();
my $desc = '';
open D, "<$GOCatDescFile" or die "can not open $GOCatDescFile\n";
while (<D>){
  chomp;
  last if m/^$/;
  @line = split;
  push @GODescAoA, [ @line ];
}
close D;
for ($i = 0; $i < @AllGOUnique; $i++){
  for ($j = 0; $j < @GODescAoA; $j++){
    if ($AllGOUnique[$i] eq $GODescAoA[$j][0]){
      push @ontologies, $GODescAoA[$j][1];
      $desc = $GODescAoA[$j][2];
      for ($k = 3; $k < @{$GODescAoA[$j]}; $k++){
	$desc = $desc . ' ' . $GODescAoA[$j][$k];
      }
      push @descriptions, $desc;
      last;
    }
    if ($j == (@GODescAoA - 1)){
      $_ = '';
      push @ontologies, $_;
      push @descriptions, $_;
    }
  }
}

#Put all output in one AoA, sort by P value, and send to output file
my @OutputAoA = ();
for ($i = 0; $i < @AllGOUnique; $i++){
  @line = ();
  push @line, $AllGOUnique[$i];
  push @line, $ontologies[$i];
  push @line, $descriptions[$i];
  push @line, $OutlierGOCounts[$i];
  push @line, $AllGOCounts[$i];
  push @line, $PValues[$i];
  push @line, $GOCatOutlierList[$i];
  push @OutputAoA, [ @line ];
}
#@OutputAoA = sort {$a->[5] cmp $b->[5]} @OutputAoA;
open O, ">$OutputFile";
for ($i = 0; $i < @OutputAoA; $i++){
  for ($j = 0; $j < @{$OutputAoA[$i]}; $j++){
    print O $OutputAoA[$i][$j];
    if ($j < (@{$OutputAoA[$i]} - 1)){
      print O "\t";
    }
    else{
      print O "\n";
    }
  }
}
close O;
