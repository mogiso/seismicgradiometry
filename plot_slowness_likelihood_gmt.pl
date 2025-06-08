#!/usr/bin/env perl

$out = $ARGV[0];
$in_vector = $ARGV[1];
$in_likelihood = $ARGV[2];

$size_x = 15;
$lon_w = 122.0;
$lon_e = 149.0;
$lat_s = 23;
$lat_n = 47;
$annot_lon = "a5";
$annot_lat = "a5";
$appvelcpt = "appvel.cpt";
$vector_len = "0.5c";

system "gmt set PS_LINE_JOIN round";

system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -K -P -X3c -Y4c > $out";

system "gmt pscoast -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -W0.8p,black -Dh \\
                    -Bpx$annot_lon -Bpy$annot_lat -BWSen -O -K -P >> $out";

open OUT, " | gmt psxy -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -SVb0.15/0.3/0.2 -C$appvelcpt -W0.8p,black -O -K -P >> $out";
$filesize = -s $in_vector;
$fileread = 0;
open IN, "<", $in_vector;
while($fileread < $filesize){
  read IN, $lon, 4;
  read IN, $lat, 4;
  read IN, $az,  4;
  read IN, $apv, 4;
  $lon = unpack "f", $lon;
  $lat = unpack "f", $lat;
  $az  = unpack "f", $az;
  $apv = unpack "f", $apv;
  print OUT "$lon $lat $apv $az $vector_len\n";
  $fileread = $fileread + 4 * 4;
}
close IN;
close OUT;


system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -P >> $out";

unlink "gmt.conf";
unlink "gmt.history";

