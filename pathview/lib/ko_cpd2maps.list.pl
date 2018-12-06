#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
my $usage=<<USAGE;

Usage:    perl $0 --gene_data ko.table.csv --cpd_data cpd.table.csv --ko2map ko2maps.xls --cpd2map cpd2maps.xls --outdir ./ --maps_list select_maps.list --prefix pathview

Function: select maps that contain the genes and cpds

Command: --gene_data Gene Data accepts data matrices in tab- or comma-delimited format (xls, csv)
		 --cpd_data Compound Data accepts data matrices in tab- or comma-delimited format (xls, csv)
		 --ko2map ko2maps.xls 
		 --cpd2map cpd2maps.xls 
		 --outdir ./ 
		 --maps_list output file [default "select_maps.list"] 
		 
Author:   yangfenglong

USAGE


my %opt = (
	ko2maps=>"$Bin/ko_cpd2maps/ko2maps.xls",
	cpd2maps=>"$Bin/ko_cpd2maps/cpd2maps.xls",
	outdir=> './', maps_list=>'select_maps.list',
);
GetOptions(
    \%opt,"gene_data:s","cpd_data:s","ko2map:s","cpd2map:s","outdir:s","maps_list:s",
);

die $usage unless (exists $opt{gene_data} || exists $opt{cpd_data});
(-d $opt{outdir}) || `mkdir -p $opt{outdir}`;

my @maps;
if($opt{gene_data}){
	my %h = split/\s+/,`less $opt{ko2maps}`;
	for(`less $opt{gene_data}`){
		if($_=~/^(K\d{5})/){ 
		    my $id=$1; 
		    $h{$id} &&  push @maps, (split/,/,$h{$id});
		}
	}
}

if($opt{cpd_data}){
	my %h = split/\s+/,`less $opt{cpd2maps}`;
	for(`less $opt{cpd_data}`){
		if($_=~/^(C\d{5})/){
		    my $id=$1;
		    $h{$id} &&  push @maps, (split/,/,$h{$id});
		}
	}
}

my %count;
my @uniq = grep {++$count{$_}<2} @maps;
open(MAP,">$opt{outdir}/$opt{maps_list}");
my $mapids = join("\n",@uniq);
print MAP $mapids,"\n";
close MAP;
print join("\n",@uniq);
