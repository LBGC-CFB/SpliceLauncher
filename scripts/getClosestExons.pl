#!/usr/bin/perl

use strict;
use warnings;

# Argument 1: Les exons des NM au format BED
# Argument 2: Le fichier de jonctions  annotées et comptées

my ($captureExonsNMBed,$annotFile) = @ARGV;


my %exonsStarts;
my %exonsEnds;

open(CAPTURE,$captureExonsNMBed);
while(<CAPTURE>){
    #chr2    215590369       215593732       NM_001282543_exon_0_0_chr2_215590370_r  0       -
    chomp;
    my ($chr,$start,$end,$annot) = split(/\t/);
    my @exon = split(/_/,$annot);
    
    $exonsStarts{"NM_".$exon[1]}{$start} = "NM_".$exon[1]."_".($exon[3]+1)."_".$start."_".$end;
    $exonsEnds{"NM_".$exon[1]}{$end} = "NM_".$exon[1]."_".($exon[3]+1)."_".$start."_".$end;
}

close(CAPTURE);


# Les introns détectés et qui sont sur un NM, avec les comptages
open(ANNOT,$annotFile);
while(<ANNOT>){
    #chr13   32889804        32890558        +       NM_000059       8       1       1047    2106
    chomp;
    my ($chr,$start,$end,$strand,$nm,$c1,$c2,$c3,$c4) = split(/\t/);
    
    my $dist1 = 1000000000;
    my $closestPrevExon = "";
    my $closestNextExon = "";
    
    for my $exonEnd (keys(%{$exonsEnds{$nm}})){
	if(abs($start-$exonEnd)<$dist1){
	    $dist1 = abs($exonEnd-$start);
	    $closestPrevExon = $exonsEnds{$nm}{$exonEnd};
	}
    }
    my $dist2 = 1000000000;
    for my $exonStart (keys (%{$exonsStarts{$nm}})){
	if(abs($exonStart-$end)<$dist2){
	    $dist2 = abs($exonStart-$end);
	    $closestNextExon = $exonsStarts{$nm}{$exonStart};
	}
    }

    #if($dist1<200 || $dist2<200){
	print $chr."\t".$start."\t".$end."\t".$strand."\t".$nm."\t".$closestPrevExon."\t".$closestNextExon."\t".$c1."\t".$c2."\t".$c3."\t".$c4."\n";
    #}
}
close(ANNOT);
