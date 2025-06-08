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
  $file_read = $file_read + 4 * 4;
  print stdout "$lon $lat $deg $appvel\n";
}
close IN;
 
