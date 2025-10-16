#!/usr/bin/env perl

$out = $ARGV[0];

$lon_diff_max = 2.0;
$lat_diff_max = 2.0;
$appvel_max = 20.0;

open OUT, ">", $out;
for($i = 1; $i <= $#ARGV; $i++){
  open IN, "<", $ARGV[$i];
  while(<IN>){
    chomp $_;
    $_ =~ s/^\s+(.*?)\s+$/$1/;
    @tmp = split /\s+/, $_;
    $lon_diff = $tmp[7] + $tmp[8];
    $lat_diff = $tmp[9] + $tmp[10];
    $ot_diff = $tmp[11] + $tmp[12];
    if($lon_diff <= $lon_diff_max && $lat_diff <= $lat_diff_max && $tmp[16] <= $appvel_max){
      print OUT "$_\n";
    }
  }
  close IN;
}
close OUT;
    

