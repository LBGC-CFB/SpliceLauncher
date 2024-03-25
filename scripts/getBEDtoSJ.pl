#!/usr/bin/perl

use strict;
use warnings;

# Argument 1: Les exons des NM au format BED

my ($SJ_file) = @ARGV;

my %exonsStarts;
my %exonsEnds;

open(SJ,$SJ_file);
while(<SJ>){
    #chr1    185351    185490    2    2    1    3
    chomp;
    my ($chr,$start,$end,$strand,$motif,$known,$count,$count2,$OverHang) = split(/\t/);
    $strand =~s/2/-/;
    $strand =~s/1/+/;
    $start-=1;
    if($count ne "0"){
        print $chr."\t".$start."\t".$end."\t".$count."\t255\t".$strand."\n";
    }
}

close(SJ);
