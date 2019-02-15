#!/usr/bin/perl -w
use strict;

# Si 1 : alors on change le signe sinon non.
my $changeSign = $ARGV[0];
my $input = $ARGV[1];

open(IN,$input);
while(<IN>){
    chomp;
    my ($chr,
	$start,$end,
	$name,
	$score,
	$strand,
	$thickStart,$thickEnd,
	$color,
	$nBlocks,
	$blockSizes,
	$blockStarts) = split(/\t/);

    my @blockSizes = split(",",$blockSizes);
    my @blockStarts = split(",",$blockStarts);

    if($changeSign eq "y"){
	$strand = ($strand eq "+") ? "-":"+";
    }
    
    for(my $i=1;$i<$nBlocks;$i++){
	my $fivePrime = $blockStarts[$i-1]+$blockSizes[$i-1]+$start;
	my $threePrime = $blockStarts[$i]+$start;
	print $chr."\t".$fivePrime."\t".$threePrime."\tINTRON\t".$score."\t".$strand."\n";
    }
}
close(IN);
