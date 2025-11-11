#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use File::Basename;
use File::Spec;
#merges output files from GO_overlap_parallel into one output files

#Get command input (GO_overlap_parallel.pl output prefix)
our($opt_i);
getopt('i');


#my $InputFileKey = 'KF_new_PBE_regions_01_GO';
my $InputFileKey = $opt_i;
#my $OutputFile = $InputFileKey . '_merged.txt';


my ($basename, $dir) = fileparse($InputFileKey);
$dir = File::Spec->rel2abs($dir);

my $OutputFile = File::Spec->catfile($dir, $basename . '_merged.txt');

print "base name $basename \n";
print "dir name $dir \n";
print "output name $OutputFile \n";

my $i = 0;
my $j = 0;
my $PValue = 0;
my @line = ();


#get sim peaks from input file(s)
my @OutputAoA = ();
my @AllFiles = ();
my @InputFiles = ();


opendir DIR, $dir or die "couldn't open directory $dir\n";
@AllFiles = readdir(DIR);
closedir DIR;
for ($i = 0; $i < @AllFiles; $i++){
  if (($AllFiles[$i] =~ m/$basename/) && ($AllFiles[$i] =~ m/txt/)){
    unless ($AllFiles[$i] =~ m/merged\.txt$/) {
      push @InputFiles, "$dir/$AllFiles[$i]"; 
    }
  }
}

$j = @AllFiles;
print "Found $j files in $dir\n";


$i = @InputFiles;
print "Found $i input files\n";

open I, "<$InputFiles[0]" or die "can't open $InputFiles[0]: $!\n";
while (<I>){
  chomp;
  last if m/^$/;
  @line = split("\t",$_);
  push @OutputAoA, [ @line ];
}
close I;

$i = @OutputAoA;
print "Found $i GO categories in first permutation file\n";


for ($i = 1; $i < @InputFiles; $i++){
  open I, "<$InputFiles[$i]" or die;
  $j = 0;
  while (<I>){
    chomp;
    last if m/^$/;
    @line = split("\t",$_);
    $PValue = $line[5];
    $OutputAoA[$j][5] = $OutputAoA[$j][5] + $PValue;
    $j++;
  }
  close I;
}

for ($i = 0; $i < @OutputAoA; $i++){
  $OutputAoA[$i][5] = $OutputAoA[$i][5] / @InputFiles;
}

open O, ">$OutputFile";
for ($i = 0; $i < @OutputAoA; $i++){
  for ($j = 0; $j < @{$OutputAoA[$i]}; $j++){
    print O $OutputAoA[$i][$j];
    if ($j == (@{$OutputAoA[$i]} - 1)){
      print O "\n";
    }
    else{
      print O "\t";
    }
  }
}
close O;


