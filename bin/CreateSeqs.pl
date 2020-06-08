#!/usr/bin/perl

#####Documentation######
$num_args = $#ARGV + 1;
if ($num_args != 2) {
  print "\nUsage: CreateSeqs.pl input_fasta input_table [STDOUT]\n\n";
  print "\tinput_fasta: Fasta file with sequence to convert\n";
  print "\tinput_table: AA conversions for hashing\n";
  print "\n\nNote: This program relies is meant for only internal usage\n";
  exit;
}
############
use warnings;
#use strict; 
##Reading and verifying input commands
my $Fasta=$ARGV[0];
my $Ref=$ARGV[1];
my %AAtab;
my %hash;

##Reading and hashing genomic sequence
my $string="";

open(op,$Ref) or die "cannot open Ref file\n";

while($line=<op>){

chomp($line);

($codon, $AA) = split /\t/, $line;
    push @{$AAtab{$AA}}, $codon;
	$hash{$codon}=$AA;

}

close(op);

open(op, $Fasta) or die "cannot open Fasta file\n";

while($line=<op>){

chomp($line);

if(($line =~ m/^>/)){next;}

@chars = split //, $line;

}

close(op);

@start = ("AAAA");


for($i=0; $i <= (scalar(@chars) - 1) ; $i=$i+3){

@start = append_AA_end($chars[$i].$chars[$i+1].$chars[$i+2],@start);
}

foreach $seqi (@start){
print $seqi."\n"
}

exit;


sub append_AA_end {

@res=();
($codon, @seq) = @_;
@array = @{ $AAtab{ $hash{$codon} } };

foreach $amimi (@array){
foreach $part (@seq){
push (@res, $part.$amimi);
}
}

return(@res);
}





