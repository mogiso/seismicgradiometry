#!/usr/bin/env perl 

$in = $ARGV[0];

$filesize = -s $in;
$file_read = 0;
open IN, "<", $in;
while($file_read <= $filesize){
  read IN, $buf, 4;
  $lon = unpack "f", $buf;
  read IN, $buf, 4;
  $lat = unpack "f", $buf;
  read IN, $buf, 4;
  $deg = unpack "f", $buf;
  read IN, $buf, 4;
  $appvel = unpack "f", $buf;
  read IN, $buf, 4;
  $min_correlation = unpack "f", $buf;
  read IN, $buf, 4;
  $array_maxamp = unpack "f", $buf;
  read IN, $buf, 4;
  $array_lta = unpack "f", $buf;


  $file_read = $file_read + 4 * 7;
  print stdout "$lon $lat $deg $appvel $min_correlation $array_maxamp $array_lta\n";
}
close IN;
 
