#!/usr/bin/env perl

$out = $ARGV[0];
$out_time = $ARGV[1];
$epi = $ARGV[2];
$header = $ARGV[3];

$size_x = 13.0;
#$lon_w = 120.0;
#$lon_e = 149.0;
#$lat_s = 22.5;
#$lat_n = 48.0;
##Kyushu
#$lon_w = 126.5;
#$lon_e = 134.5;
#$lat_s = 27.5;
#$lat_n = 34.5;
#Tokara
$lon_w = 128.5;
$lon_e = 130.5;
$lat_s = 29.0;
$lat_n = 31.5;


$size_y = `echo $lon_w $lat_n | gmt mapproject -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n`;
chomp $size_y;
@tmp = split /\s+/, $size_y;
$size_y = $tmp[1];
$txt_x = 0.1;
$txt_y = $size_y - 0.2;

$date_cpt = "../date.cpt";
$cpt_x = $size_x / 2.0;
$cpt_y = -1.2;
$cpt_len = $size_x;
$cpt_width = "0.3ch";
$annot_cpt = "a3hg1h+l\"Time (hour)\"";

$annot_lon = "a1.0";
$annot_lat = "a1.0";
$clon = 137.7122;
$clat = 37.5144;
$horizon = 90.0;
$symbolsize = "0.5";
$symbolsize_circle = "0.2";

if($header > 0){
  $header = "-H${header}";
}else{
  $header = "";
}

system "gmt set PS_LINE_JOIN round";
system "gmt set GMT_HISTORY false";

system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -K -P -X3c -Y4c > $out";

system "gmt pscoast -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Df -Bpx${annot_lon} -Bpy${annot_lat} -BWSen \\
                    -W0.5p,black -O -K -P >> $out";

open OUT, " | gmt psxy -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
                       -Sa$symbolsize -W0.6p,black -C$date_cpt ${header} -O -K -P >> $out";
$count = 0;
$datemin = substr $ARGV[4], 0, 8;
$datemax = substr $ARGV[$#ARGV], 0, 8;
for($i = 3; $i <= $#ARGV; $i++){
  open IN, "<", $ARGV[$i];
  while(<IN>){
    chomp $_;
    @tmp = split /\s+/, $_;
    if(($tmp[3] >= $lon_w && $tmp[3] <= $lon_e && $tmp[6] >= $lat_s && $tmp[6] <= $lat_n)){
      $count++;
      print OUT "$tmp[3] $tmp[6] $tmp[0]\n";
    }
  }
  close IN;
}
close OUT;

open OUT, " | gmt psxy -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Sc$symbolsize_circle -W0.4p,black -C$date_cpt -O -K -P >> $out";
open IN, "<", $epi;
while(<IN>){
  chomp $_;
  $yr = substr $_, 1, 4;
  $mo = substr $_, 5, 2;
  $dy = substr $_, 7, 2;
  $hr = substr $_, 9, 2;
  $mi = substr $_, 11, 2;
  $sc = substr $_, 13, 4;
  $sc = $sc / 100.0;
  $lat = substr $_, 22, 2;
  $lat_m = substr $_, 24, 4;
  $lat = $lat + $lat_m / 100.0 / 60.0;
  $lon = substr $_, 33, 3;
  $lon_m = substr $_, 36, 4;
  $lon = $lon + $lon_m / 100.0 / 60.0;
  if(($lon >= $lon_w && $lon <= $lon_e && $lat >= $lat_s && $lat <= $lat_n)){
    print OUT "$lon $lat ${yr}-${mo}-${dy}T${hr}:${mi}:${sc}\n";
  }
}
close IN;
close OUT;


open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -F+f+a+j -Gwhite -O -K -P >> $out";
print OUT "$txt_x $txt_y 14p,Helvetica,black 0 LT ${datemin}-${datemax} N=${count}\n";
close OUT;

system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -C$date_cpt -B${annot_cpt} -O -K -P >> $out";
system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -P >> $out";

system "gmt psconvert $out -Tg -A";


##M-T plot
$time_min = "2025-07-05T00:00";
$time_max = "2025-07-06T00:00";
$annot_x = "a6hf1h+l\"Time (hour)\"";
$mag_min = 2.5;
$mag_max = 5.8;
$annot_y = "a1+l\"Magnitude\"";
$size_x = 10;
$size_y = $size_x * 9.0 / 16.0;

system "gmt psbasemap -JX${size_x}T/$size_y -R$time_min/$time_max/$mag_min/$mag_max \\
                      -Bpx${annot_x} -Bpy${annot_y} -BWSen -K -P -X3c -Y4c > $out_time";
open IN, "<", $epi;
while(<IN>){
  chomp $_;
  $yr = substr $_, 1, 4;
  $mo = substr $_, 5, 2;
  $dy = substr $_, 7, 2;
  $hr = substr $_, 9, 2;
  $mi = substr $_, 11, 2;
  $sc = substr $_, 13, 4;
  $sc = $sc / 100.0;
  $lat = substr $_, 22, 2;
  $lat_m = substr $_, 24, 4;
  $lat = $lat + $lat_m / 100.0 / 60.0;
  $lon = substr $_, 33, 3;
  $lon_m = substr $_, 36, 4;
  $lon = $lon + $lon_m / 100.0 / 60.0;
  $mag = substr $_, 52, 2;
  $mag = $mag / 10.0;
  if(($lon >= $lon_w && $lon <= $lon_e && $lat >= $lat_s && $lat <= $lat_n)){
    open OUT, " | gmt psxy -JX${size_x}T/$size_y -R$time_min/$time_max/$mag_min/$mag_max -W0.8p,black -O -K -P >> $out_time";
    print OUT "${yr}-${mo}-${dy}T${hr}:${mi}:${sc} $mag_min\n";
    print OUT "${yr}-${mo}-${dy}T${hr}:${mi}:${sc} $mag\n";
    close OUT;
    open OUT, " | gmt psxy -JX${size_x}T/$size_y -R$time_min/$time_max/$mag_min/$mag_max -W0.4p,red -Sc0.25 -O -K -P >> $out_time";
    print OUT "${yr}-${mo}-${dy}T${hr}:${mi}:${sc} $mag\n";
    close OUT;
  }
}
close IN;

open OUT, " | gmt psxy -JX${size_x}T/$size_y -R$time_min/$time_max/$mag_min/$mag_max \\
                       -Sa$symbolsize -W0.5p,black -Gwhite -N -O -K -P >> $out_time";
$count = 0;
$datemin = substr $ARGV[4], 0, 8;
$datemax = substr $ARGV[$#ARGV], 0, 8;
for($i = 3; $i <= $#ARGV; $i++){
  open IN, "<", $ARGV[$i];
  while(<IN>){
    chomp $_;
    @tmp = split /\s+/, $_;
    if(($tmp[3] >= $lon_w && $tmp[3] <= $lon_e && $tmp[6] >= $lat_s && $tmp[6] <= $lat_n)){
      $count++;
      print OUT "$tmp[0] $mag_max\n";
    }
  }
  close IN;
}
close OUT;

system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -P >> $out_time";
system "gmt psconvert $out_time -Tg -A";

unlink "gmt.history";
unlink "gmt.conf";

