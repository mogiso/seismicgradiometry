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
$annot_appvelcpt = "+l\"Apparent velocity (km/s)\"";
$appvelcpt_y = -1.0;
$vector_len = "0.5c";
$vector_shape = "0.12/0.23/0.2";

$likelihoodcpt = "likelihood.cpt";
$annot_likelihoodcpt = "+l\"Normalized likelihood\"";
$likelihoodcpt_y = -3.0;
$particle_size = 0.2;
$symbolpen = "0.4p,black";

$cpt_x = $size_x / 2.0;
$cpt_len = $size_x;
$cpt_width = "0.3ch";


system "gmt set PS_LINE_JOIN round";
system "gmt set LABEL_FONT 12p,Helvetica,black";

system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -K -P -X3c -Y7c > $out";

system "gmt psxy $in_likelihood -bif -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
                                -Sc$particle_size -C$likelihoodcpt -W$symbolpen -O -K -P >> $out";

open OUT, " | gmt psxy -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -SVb$vector_shape -C$appvelcpt -W$symbolpen -O -K -P >> $out";
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

system "gmt pscoast -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -W0.5p,black -Di \\
                    -Bpx$annot_lon -Bpy$annot_lat -BWSen -O -K -P >> $out";

system "gmt psscale -Dx$cpt_x/$appvelcpt_y/$cpt_len/$cpt_width -C$appvelcpt -B$annot_appvelcpt -O -K -P >> $out";
system "gmt psscale -Dx$cpt_x/$likelihoodcpt_y/$cpt_len/$cpt_width -C$likelihoodcpt -B$annot_likelihoodcpt -Q -O -K -P >> $out";

system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -P >> $out";

unlink "gmt.conf";
unlink "gmt.history";

