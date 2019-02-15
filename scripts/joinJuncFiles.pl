#!/usr/bin/perl

use strict;
use warnings;

my @files = @ARGV;

my %allJuncs;

for my $f (@files){
    open(IN,$f);
    while(<IN>){
	my ($chr,$start,$end,$strand,$NM,$intron,$count) = split(/\t/);

	$allJuncs{$start."_".$end."_".$strand}{data} = 
	{
	    chr => $chr,
	    start=>$start,
	    end=> $end,
	    strand=>$strand,
	    nm=>$NM,
	    intron=>$intron,
	};
	$allJuncs{$start."_".$end."_".$strand}{count}{$f}=$count;
    }
    close(IN);
}


print "chr\tstart\tend\tstrand\tNM\tintron";
for my $f (@files){
    print "\t".$f;
}
print "\n";


for my $key (keys(%allJuncs)){
    my $data = $allJuncs{$key}{data};
    print   $data->{chr}."\t".$data->{start}."\t".$data->{end}."\t".$data->{strand}."\t".$data->{nm}."\t".$data->{intron};
    
    for my $f (@files){
	if(defined $allJuncs{$key}{count}{$f}){
	    print  "\t".$allJuncs{$key}{count}{$f};
	}else{
	    print  "\t0";
	}
    }
    print "\n";
}
