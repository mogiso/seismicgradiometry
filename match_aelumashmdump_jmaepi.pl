#!/usr/bin/env perl

$arrayresultfile = $ARGV[0];
$jmahypofile = $ARGV[1];
$matchresultfile = $ARGV[2];
$nonmatchresultfile = $ARGV[3];


open IN, "<", $arrayresultfile;
while(<IN>){
  chomp $_;
  $_ =~ s/^\s+/(.*?)\s+$/$1/;
  @tmp = split /\s+/, $_;
  push @arrayresult, $_;
  push @arraydatetime, $tmp[0];
  push @arrayepilon, $tmp[3];
  push @arrayepilat, $tmp[6];
}
close IN;


