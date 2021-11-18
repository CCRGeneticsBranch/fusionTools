#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;

my $source_gtf;
my $output_gtf;
my $canonical_file;
my $canonical_gtf;
my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

Options:

  -s  <string>  source GTF
  -o  <string>  output GTF
  -c  <string>  canonical output file
  -g  <string>  canonical output GTF
  
__EOUSAGE__



GetOptions (
  's=s' => \$source_gtf,
  'o=s' => \$output_gtf,
  'c=s' => \$canonical_file,
  'g=s' => \$canonical_gtf  
);

open(IN_FILE, "$source_gtf") or die "Cannot open file $source_gtf";
open(OUT_GTF, ">$output_gtf") or die "Cannot open file $output_gtf";
open(OUT_C, ">$canonical_file") or die "Cannot open file $canonical_file";
open(OUT_C_GTF, ">$canonical_gtf") or die "Cannot open file $canonical_gtf";
print OUT_C "Symbol\tEnsembl\tCanonical\tMANE\n";
while(<IN_FILE>) {
	chomp;
	my @s=split(/\t/);
	my ($gn)=$_=~/gene_name "(.*?)".*/;
	my ($tid)=$_=~/transcript_id "(.*?)\..*/;
	if ($tid) {
		$_ =~ s/gene_name "(.*?)"/gene_name "$1:\Q${tid}\E"/;
		print OUT_GTF "$_\n";
		my $is_canonical = (/Ensembl_canonical/)? 'Y': 'N';
		my $is_mane = (/MANE_Select/)? 'Y': 'N';;
		if ($is_canonical eq "Y") {
			print OUT_C_GTF "$_\n";
		}
		if ($is_canonical eq "Y" || $is_mane eq "Y") {
			if ($s[2] eq "transcript") {
				print OUT_C"$gn\t$tid\t$is_canonical\t$is_mane\n";
			}
		}
		if (/MANE_Select/) {
		}
	}
}
close(OUT_GTF);
close(OUT_C_GTF);
close(OUT_C);
close(IN_FILE);