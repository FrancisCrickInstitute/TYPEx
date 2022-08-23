#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use File::Find::Rule;
use File::Path qw(make_path);

my ($inDir, $outDir, $run) = @ARGV;
die "USAGE: basename($0) <inDir> <outDir> <run>"
  unless defined defined $inDir && defined $outDir;

my $delim = ",";
# my $analysisDir = join("/", $ENV{"HOME"}, "labwd/analyses/imc/nextflow/");
my @features = ("MeanIntensity");
 print($outDir, "\n");
if(! -d "$outDir" && ! -e "$outDir" && ! -l $outDir )  {
  make_path("$outDir");
}

my (%stats);
print "$inDir\n";
my @inDirs=File::Find::Rule->maxdepth(1)
  ->directory
  ->in("$inDir");
print $#inDirs + 1, "\n";

for my $runDir (@inDirs) {
  print "$runDir\n";
  my @inFiles = File::Find::Rule->file()
                            ->name('Cells.csv')
                            ->in("$runDir/results/segmentation");
  print $#inFiles + 1, ' files found in ', $runDir, "\n";

  foreach my $inFile (@inFiles) {

    print 'File ', $inFile, "\n";
    my $sample=$inFile;
    $sample =~ s/.*segmentation.(.*).Cells.csv/$1/;
    print "Sample ", $sample, "\n";

    # next if( -f join("/", $outDir, "MeanIntensity", $run, "$sample.txt");
    $sample =~ s/\//-/;
    open IN, $inFile or die "Could not open $inFile:$!\n";
    my @header = split /$delim/, <IN>;
    $header[$#header] =~ s/\n//;
    @header = map {s/Intensity_//; $_} @header;
    @header = map {s/Location_Center/LocationCenter/; $_} @header;
    #print join(",", @header), "\n";
    my @features = @header;
    map {s/_.*//; $_} @features;
    my @markers = @header;
    map { $_=~s/^[_]+//; $_}  @markers;
    map {s/^.*_([^_]+)_c1$/$1/; $_} @markers;
    my ($object_col) = grep { $header[$_] eq "ObjectNumber" } (0 .. $#header);
    print $object_col, "\n";

    my (%in_hash, %markers);
    while(<IN>) {
      my @columns = split /$delim/;
      my $object=$columns[$object_col];
      for (my $i = 0; $i < $#columns; $i++) {
        next if($header[$i] =~ "Location_Max|Location_CenterMass");
        die "Duplicate IDs $sample $features[$i] $markers[$i] $object $header[$i] \n"
          if exists $in_hash{$features[$i]}{$object}{$markers[$i]};
        $in_hash{$features[$i]}{$object}{$markers[$i]} = $columns[$i];
        $markers{$features[$i]}{$markers[$i]} = 1;
      }
    }
    close(IN);
    print "Read $inFile \n";

    # SPLIT BY FEATURE, PANEL,
    # summarize number of cells per file,
    # number of cells with NA values, or zero values?
    foreach my $feature (keys %in_hash)	{
      if(! -d join("/", $outDir, $feature)) {
        mkdir(join("/", $outDir, $feature));
      }
      if(! -d join("/", $outDir, $feature, $run)) {
        mkdir(join("/", $outDir, $feature, $run));
      }
      my @columns = keys %{ $markers{$feature}};
      my $fileOut = join("/", $outDir, $feature, $run, "$sample.txt");
      open OUT, ">". $fileOut or die "Could not open $fileOut:$!\n";
      print OUT join($delim, "ObjectNumber", @columns), "\n";
      my @cells = sort keys %{$in_hash{$feature}};
      foreach my $cell (@cells) {
        print OUT join($delim, $cell,
          map {
            if(exists($in_hash{$feature}{$cell}{$_})) {
              $in_hash{$feature}{$cell}{$_}
            } else {
              ""
            }
          } @columns), "\n";
        }
        close(OUT);
    }
  }
}

sub unique {
  my @array = @_;
  my %hash;
  foreach my $element (@array)  {
    $hash{$element} = 1;
  }
  return(keys %hash);
}
