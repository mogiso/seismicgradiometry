#!/usr/bin/env perl
use Math::Trig 'rad2deg';

$out = $ARGV[0];
$in_dir = $ARGV[1];
$yr = $ARGV[2];
$mo = $ARGV[3];
$dy = $ARGV[4];
$hr = $ARGV[5];
$mi = $ARGV[6];
$sc = $ARGV[7];

$datetime = "${yr}-${mo}-${dy}T${hr}:${mi}:${sc}";
$epicenterlist = "${in_dir}/${yr}${mo}${dy}_epicenter.txt";
$in_vector = "${in_dir}/${yr}${mo}${dy}${hr}${mi}${sc}_array_obs.dat";
$in_likelihood = "${in_dir}/${yr}${mo}${dy}${hr}${mi}${sc}_likelihood_particle.dat";

print stderr "$datetime\n";
print stderr "$epicenterlist\n";

open IN, "<", $epicenterlist;
while(<IN>){
  if($_ =~ /^$datetime/){
    print stdout "$_";
    @tmp = split /\s+/, $_;
    $eplon = $tmp[3];
    $eplat = $tmp[6];
    $eplon_err = $tmp[7];
    $eplat_err = $tmp[8];
  }
}
close IN;


$size_x = 15;
#$lon_w = 122.0;
#$lon_e = 149.0;
#$lat_s = 23;
#$lat_n = 47;
$lon_w = 126.5;
$lon_e = 134.5;
$lat_s = 27.5;
$lat_n = 34.5;
$size_y = `echo $lon_w $lat_n | gmt mapproject -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n | awk '{print $2}'`;
@tmp = split /\s+/, $size_y;
$size_y = $tmp[1];
$txt_x = 0.1;
$txt_y = $size_y - 0.2;
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
$particle_size = 0.15;
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
  read IN, $buf, 4; #min_correlation
  read IN, $buf, 4; #array_maxamp
  read IN, $buf, 4; #array_lta
  $lon = unpack "f", $lon;
  $lat = unpack "f", $lat;
  $az  = unpack "f", $az;
  $apv = unpack "f", $apv;
  print OUT "$lon $lat $apv $az $vector_len\n";
  $fileread = $fileread + 4 * 7;
}
close IN;
close OUT;

open OUT, " | gmt psxy -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -W1p,black -Sa0.4c -Gwhite -Exy+p1p,gray -O -K -P >> $out";
print OUT "$eplon $eplat $eplon_err $eplat_err\n";
close OUT;

system "gmt pscoast -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -W0.5p,black -Di \\
                    -Bpx$annot_lon -Bpy$annot_lat -BWSen -O -K -P >> $out";

open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -F+f+a+j -Gwhite -O -K -P >> $out";
print OUT "$txt_x $txt_y 14p,Helvetica,black 0 LT $datetime\n";
close OUT;


system "gmt psscale -Dx$cpt_x/$appvelcpt_y/$cpt_len/$cpt_width -C$appvelcpt -B$annot_appvelcpt -O -K -P >> $out";
system "gmt psscale -Dx$cpt_x/$likelihoodcpt_y/$cpt_len/$cpt_width -C$likelihoodcpt -B$annot_likelihoodcpt -Q -O -K -P >> $out";

system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -P >> $out";

system "gmt psconvert $out -Tg -A";

unlink "gmt.conf";
unlink "gmt.history";

