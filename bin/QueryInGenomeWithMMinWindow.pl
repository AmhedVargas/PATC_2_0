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

##Reading and verifying input commands
my $Fasta=$ARGV[0];
my $Ref=$ARGV[1];
my @dict;


##Reading and hashing genomic sequence
my $string="";
open(op,$Ref) or die "cannot open Ref file\n";
while($line=<op>){
chomp($line);
if(($line =~ m/^>/)){next;}
$string=$string.$line
}
close(op);

$Genomic = $string;
$Revcomp = reverse $Genomic;
$Revcomp =~ tr/ACGTacgt/TGCAtgca/;

#Hash both strands
#print "\nIndexing genome\n";
my %hash;
for($i=0;$i<=(length($Genomic)-$k);$i++){
$hash{substr($Genomic,$i,$k)}++;
$hash{substr($Revcomp,$i,$k)}++;
}

##Read each query and perform actions
my $Qseq;
my @array;
my $sarray;
my $val;
open(op,$Query) or die "cannot open query file\n";
#print "\nReading sequences\n";
##my $con=0;
while($line=<op>){
##$con++;
chomp($line);
if($line =~ m/^>/){$Qseq=$line}else{
if(length($line) != $k){die $Qseq." query sequence has a different size (".length($Qseq).") than k-mer-length (".$k.").\n";}
#if($hash{$line}>1){next;}
#$val=$hash{$line};
#my $flag=1;
#print $Qseq."\t".$line."\t".$val."\n"}}

##If mm larger than 0 go for array approach otherwise print hashed value
$val = 0;
if($mm == 0){
  $val=$hash{$line};
  }else{
    #print "Line to test".$line."\n";
    @array = produceSeqs($line,$pwin1,$pwin2,$mm);
    foreach $finito(@array){
      #print "\nHola:".$finito."\n";
      $val = $val + $hash{$finito}; 
      if(length($finito) != $k){die "Amhed you did something weird with the sequence ".$line.' producing '.$finito}}
  }

print $Qseq."\t".$line."\t".$val."\n";
#if(flag ==1){print $Qseq."\t".$line."\t".$val."\n";}
} #Stop if not header
##if(!($con % 5)){print "\n".$con." sequences analyzed\n";}
} ##Stop reading file
close(op);

exit;

###Sub function: Produce $line $start $end $number of mismatches
sub produceSeqs {
  #print "\nA\n";
  @seq=split('',$_[0]);
  #print @seq;
  #print "\n";
  $left="";
  $right="";
  if($_[1] > 0){$left= $left.join('',@seq[0..($_[1]-1)])}
  if(($_[2]+1) < scalar(@seq)){$right= $right.join('',@seq[($_[2]+1)..(scalar(@seq)-1)])}
  #print "Left:".$left."\n";
  #print "Right:".$right."\n";
  @temp = repseq($_[3],@seq[$_[1]..$_[2]]);
  @seqs=();
  foreach $tem (@temp){push(@seqs,$left.$tem.$right)}
  #print "\nSeqsA\n";
  #print @seqs;
  return(@seqs);
}

sub repseq{
  #print "\nB Args\n";
  @args=@_;
  #print @args;
  $rep=$args[0];
  #print "\nRepetitions:".$rep;
  @letters=@args[1..(scalar(@args)-1)];
  #print "\n Args again\n";
  #print @args;
  #print "\nSize:".scalar(@args);
  #print "\nMaxIndex:".$#args;
  #print "\nletters:\n";
  #print @letters;
  #@dict=("A","T","C","G");
  @res=();
  push(@res,join('',@letters));
  my %rest;
  #$size = keys %rest;
  #print "\nAvantResB\n";
  #print $size;
  while($rep != 0){
    @res=perseq(@res);
    $rep--;
  }
  foreach $s (@res){$rest{$s}++}
  #print "\nResB\n";
  #$size = keys %rest;
  #print $size;
  return(keys(%rest));
}

sub perseq{
  #print "\nC\n";
  @fin=();
  foreach $sequence (@_){
    #print "\nSequence inside C:\n";
    #print "-".$sequence."-";
    @pat=split('',$sequence);
    for($j=0; $j< scalar(@pat);$j++){
      @tmep=@pat;
      $tmep[$j]="A";
      push(@fin,join("",@tmep));
      $tmep[$j]="T";
      push(@fin,join("",@tmep));
      $tmep[$j]="C";
      push(@fin,join("",@tmep));
      $tmep[$j]="G";
      push(@fin,join("",@tmep));
    }
  }
  return(@fin);
}
