#!/usr/bin/env perl

$out = $ARGV[0];
$gscsv = $ARGV[1];
$header = $ARGV[2];

$size_x = 13.0;
$lon_w = 120.0;
$lon_e = 149.0;
$lat_s = 22.5;
$lat_n = 48.0;
$annot_lon = "a5.0";
$annot_lat = "a5.0";
$clon = 137.7122;
$clat = 37.5144;
$horizon = 180.0;
$symbolsize = "0.4";

$txt_x = 0.1;
$txt_y = $size_x + 0.1;

if($header > 0){
  $header = "-H${header}";
}else{
  $header = "";
}


$date_cpt = "../date.cpt";
$cpt_x = $size_x / 2.0;
$cpt_y = -0.6;
$cpt_len = $size_x;
$cpt_width = "0.3ch";
$annot_cpt = "a3hg1h+l\"Time (hour)\"";

system "gmt set PS_LINE_JOIN round";
system "gmt set GMT_HISTORY false";

system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -K -P -X3c -Y4c > $out";

system "gmt pscoast -JE$clon/$clat/$horizon/$size_x -Rg -B90 -Dc -A10000 -W0.5p,black -O -K -P >> $out";

open OUT, " | gmt psxy -JE$clon/$clat/$horizon/$size_x -Rg -Sa$symbolsize -W0.5p,black -C$date_cpt ${header} -O -K -P >> $out";
$count = 0;
$datemin = substr $ARGV[3], 0, 8;
$datemax = substr $ARGV[$#ARGV], 0, 8;
for($i = 3; $i <= $#ARGV; $i++){
  open IN, "<", $ARGV[$i];
  while(<IN>){
    chomp $_;
    @tmp = split /\s+/, $_;
    if(!($tmp[3] >= $lon_w && $tmp[3] <= $lon_e && $tmp[6] >= $lat_s && $tmp[6] <= $lat_n)){
      $count++;
      print OUT "$tmp[3] $tmp[6] $tmp[0]\n";
    }
  }
  close IN;
}
close OUT;

open OUT, " | gmt psxy -JE$clon/$clat/$horizon/$size_x -Rg -Sc0.2c -W0.3p,black -C$date_cpt -H1 -O -K -P >> $out";
open IN, "<", $gscsv;
while(<IN>){
  chomp $_;
  @tmp = split /,/, $_;
  if(!($tmp[2] >= $lon_w && $tmp[2] <= $lon_e && $tmp[1] >= $lat_s && $tmp[1] <= $lat_n)){
    print OUT "$tmp[2] $tmp[1] $tmp[0]\n";
  }
}
close IN;
close OUT;

open OUT, " | gmt pstext -JX$size_x/$size_x -R0/$size_x/0/$size_x -F+f+a+j -Gwhite -N -O -K -P >> $out";
print OUT "$txt_x $txt_y 14p,Helvetica,black 0 LB ${datemin}-${datemax} N=${count}\n";
close OUT;


system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -C$date_cpt -B${annot_cpt} -O -K -P >> $out";
system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -P >> $out";

system "gmt psconvert $out -Tg -A";


unlink "gmt.history";
unlink "gmt.conf";

