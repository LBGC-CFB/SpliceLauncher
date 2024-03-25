#!/usr/bin/perl

use strict;
use warnings;

my @files = @ARGV;

my %allJuncs;

for my $f (@files){
    open(IN,$f);
    while(<IN>){
	my ($chr,$start,$end,$strand,$gene,$count) = split(/\t/);

	$allJuncs{$start."_".$end."_".$strand."_".$gene}{data} = 
	{
	    chr => $chr,
	    start=>$start,
	    end=> $end,
	    strand=>$strand,
	    gene=>$gene
	};
	$allJuncs{$start."_".$end."_".$strand."_".$gene}{count}{$f}=$count;
    }
    close(IN);
}


print "chr\tstart\tend\tstrand\tgene";
for my $f (@files){
    print "\t".$f;
}
print "\n";


for my $key (keys(%allJuncs)){
    my $data = $allJuncs{$key}{data};
    print   $data->{chr}."\t".$data->{start}."\t".$data->{end}."\t".$data->{strand}."\t".$data->{gene};
    
    for my $f (@files){
	if(defined $allJuncs{$key}{count}{$f}){
	    my $junc_count = $allJuncs{$key}{count}{$f};
		$junc_count =~s/\n//;
	    print  "\t".$junc_count;
	}else{
	    print  "\t0";
	}
    }
    print "\n";
}
