#!/usr/bin/env perl
use Math::Trig qw(pi rad2deg);

$in = $ARGV[0];
$out = $ARGV[1];
$in2 = $ARGV[2];


$min_x = -700;
$max_x = 500;
$min_y = -400;
$max_y = 700;

$size_x = 12;
$size_y = 11;

$cpt = "amplitude_gradiometry.cpt";

$txt_x = 0.0;
$txt_y = $size_y + 0.5;
$cpt_x = $size_x / 2.0;
$cpt_y = -1.0;
$cpt_len = $size_x;
$cpt_width = "0.3ch";
$slowness_length = 0.3;
$decimate = 30;

system "gmt set PS_LINE_JOIN round";

system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -K -P -X3c -Y4c > $out";

system "gmt grdimage $in -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -C$cpt -O -K -P >> $out";
system "gmt psbasemap -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Bpxa100f50 -Bpya100f50 -BWSen -O -K -P >> $out";
open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -O -K -P >> $out";
  print OUT "$txt_x $txt_y 13p,Helvetica,black 0 LB $in\n";
close OUT;

if(-f $in2){
  $filesize = -s $in2;
  $read_filesize = 0;
  open OUT, " | gmt psxy -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y \\
                -Svb0.08c/0.12c/0.125c -W0.5p,black -O -K -P >> $out";
  open IN, "<", $in2;
  while($read_filesize < $filesize){
    read IN, $buf, 4;
    $x_east = unpack "f", $buf;
    read IN, $buf, 4;
    $y_north = unpack "f", $buf;
    read IN, $buf, 4;
    $slowness_x = unpack "f", $buf;
    read IN, $buf, 4;
    $slowness_y = unpack "f", $buf;
    $read_filesize = $read_filesize + 4 * 4;
    if($slowness_x != 0.0 && $slowness_y != 0.0 && $x_east % $decimate == 0 && $y_north % $decimate == 0){
      $slowness = sqrt($slowness_x * $slowness_x + $slowness_y * $slowness_y);
      $direction = rad2deg(atan2($slowness_y, $slowness_x));
      #print stdout "$x_east $y_north $direction $slowness_length $slowness\n";
      print OUT "$x_east $y_north $direction $slowness_length $slowness\n";
    }
  }
  close IN;
  close OUT;
}
  
  
  

system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -C$cpt -O -K -P >> $out";


system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -P >> $out";

system "gmt psconvert $out -Tg -A";

unlink "gmt.history";
unlink "gmt.conf";

