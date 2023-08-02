#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use File::Find::Rule;
use File::Path qw(make_path);

my ($inFile, $core, $feature, $imcyto) = @ARGV;
die "USAGE: basename($0) <inFile> <core> <feature> [<imcyto>]"
	unless defined $inFile &&  defined $core && defined $feature;

if(not defined $imcyto) {
	$imcyto = 1;
}

my $delim = ",";
if($imcyto eq 'false') {
	$delim = "\t";
}
print $delim, "\n";
print 'File ', $inFile, "\n";
  
$core =~ s/.*segmentation.(.*).Cells.csv/$1/;
print "core ", $core, "\n";

$core =~ s/\//-/;
open IN, $inFile or die "Could not open $inFile:$!\n";
my @header = split /$delim/, <IN>;
$header[$#header] =~ s/\n//;

if($imcyto eq 'true')	{
	@header = map {s/Intensity_//; $_} @header;
	@header = map {s/Location_Center/LocationCenter/; $_} @header;
	# simple segmentaiton
	@header = map {s/mean_intensity/MeanIntensity/i; $_} @header;
	@header = map {s/centroid-?0/LocationCenter_X/; $_} @header;
	@header = map {s/centroid-?1/LocationCenter_Y/; $_} @header;
	@header = map {s/label/ObjectNumber/; $_} @header;
	@header = map {s/area/AreaShape_Area/; $_} @header;
} else {
	print('deep-imcyto independent run');
	# Input format: "ObjectNumber"  "imagename"     "Center_X"       "Center_Y" 
		# "area" [ marker 1 \t marker 2 ... marker n]
	@header = map {
		print join("\t", $_, $header[$_], "\n");
			if($_ == 0) {
				"ObjectNumber"
			} elsif ($_ == 1) {
				"imagename"
			} elsif (($_ == 2 || $_ == 3) && not $header[$_] =~ /LocationCenter/) {
				"LocationCenter_$header[$_]"
			} elsif( $_ == 4) {
				"AreaShape_$header[$_]"
			} else {
				"MeanIntensity_$header[$_]"
			}
		} (0 .. $#header);
	
}

# print join(",", @header), "\n";
my @features = @header;
map {s/_.*//; $_} @features;
my @markers = @header;
map { $_=~s/^[_]+//; $_}  @markers;
map {s/^.*_([^_]+)_c1$/$1/; $_} @markers;
my ($object_col) = grep { $header[$_] eq "ObjectNumber" } (0 .. $#header);
my ($imageID_col) = grep { $header[$_] eq "imagename" } (0 .. $#header);
# print join("\t", $object_col, $imageID_col, "\n");
print join(",", @header, "\n");

my (%in_hash, %markers);
while(<IN>) {
	$_ =~ s/\n//;
      my @columns = split /$delim/;
      my $object=$columns[$object_col];
			my $imageID;
			if(! defined $imageID_col) {
				$imageID=$core;
			} else {
				$imageID=$columns[$imageID_col];
			}
			
      for (my $i = 0; $i <= $#columns; $i++) {
        next if($header[$i] =~ "Location_Max|Location_CenterMass");
        die "Duplicate IDs $core $features[$i] $markers[$i] $object $header[$i] \n"
          if exists $in_hash{$imageID}{$features[$i]}{$object}{$markers[$i]};
        $in_hash{$imageID}{$features[$i]}{$object}{$markers[$i]} = $columns[$i];
        $markers{$features[$i]}{$markers[$i]} = 1;
      }
 }
close(IN);

print "Read $inFile \n";
foreach my $imageID (keys %in_hash) {

    # SPLIT BY FEATURE, PANEL,
    # summarize number of cells per file,
    # number of cells with NA values, or zero values?
    # %in_hash

		if(! -d $feature) {
        mkdir($feature);
    }
		my @columns = keys %{ $markers{$feature}};
		my $fileOut = join("/", $feature, "$imageID.txt");
		open OUT, ">". $fileOut or die "Could not open $fileOut:$!\n";
		print OUT join("\t", "ObjectNumber", @columns), "\n";
		my @cells = sort keys %{$in_hash{$imageID}{$feature}};
		
		foreach my $cell (@cells) {
			print OUT join("\t", $cell,
				map {
					if(exists($in_hash{$imageID}{$feature}{$cell}{$_})) {
						$in_hash{$imageID}{$feature}{$cell}{$_}
					} else {
						""
					}
				} @columns), "\n";
		}
		close(OUT);
}

sub unique {
	
  my @array = @_;
  my %hash;
  foreach my $element (@array)  {
    $hash{$element} = 1;
  }
  return(keys %hash);
}
