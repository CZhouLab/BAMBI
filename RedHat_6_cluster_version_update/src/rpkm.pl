#! /usr/bin/perl 
use strict;

my $usage="$0 reads_count(e.g.HTSeq output) gene_len_file Unique/or_total_mapped_reads > of\nNote: reads_count file should have two columns: the first column is the id, and the 2nd column is the counts of reads mapped to each id\n";
if(@ARGV<3){die $usage;}
my $total_reads_num=$ARGV[2];

#10^9*No.of-reads-mapped-to-exons/(exon-length*total-unique-mapped-reads)

#TGFBR3  205
open GLEN, $ARGV[1];
chomp(my @glens=<GLEN>);
close GLEN;
my %g2len=();
for (my $i=0;$i<@glens; $i++){
	my @fields=split /\t/, $glens[$i];
	my $id=$fields[0];
	my $len=$fields[1];
	if(not exists $g2len{$id}){
		$g2len{$id}=$len;
	}else {warn "Duplicated $id\n";}
}


open IF, $ARGV[0];
my $line;

while(chomp($line=<IF>)){
	my @fields=split /\t/, $line;
	my $gene=$fields[0];
	my $score=$fields[1];
	if($g2len{$gene}>0){
		my $rpkm=1000*1000*1000*$score / ($g2len{$gene} * $total_reads_num);
		print "$gene\t$rpkm\n";
	}else{
		warn "$gene does not exist in FILE $ARGV[1]\n";

	}
}

close IF;
