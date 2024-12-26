#!/usr/bin/env perl
# Copyright 2022 Masashi Ogiso (masashi.ogiso@gmail.com)
# Released under the MIT license.
# see https://opensource.org/licenses/MIT

use Math::Trig qw(pi rad2deg);
use Parallel::ForkManager;
$MAX_PROCESS = 16;

$in_dir = $ARGV[0];
$index_begin = $ARGV[1];
$index_end = $ARGV[2];
$coastline = $ARGV[3];

$ref_yr = "2022";
$ref_mo = "1";
$ref_dy = 15;
$ref_hh = 9;
$ref_mm = 30;
$ref_ss = 0;
$dt = 60;

$ref_sec = $ref_hh * 3600 + $ref_mm * 60 + $ref_ss;

#S-net ref: 142.5E, 38.25N
$dgrid_x = 20;
$dgrid_y = 20;
$min_x = -350;
$max_x = 350;
$min_y = -600;
$max_y = 600;
$size_x = 4.8;
$size_y = 8.22;
$annot = "a200f100";

#DONET ref: 135.75E, 33.2N
#$dgrid_x = 10;
#$dgrid_y = 10;
#$min_x = -150;
#$max_x = 150;
#$min_y = -100;
#$max_y = 100;
#$size_x = 4.8;
#$size_y = 3.2;
#$annot = "a100f50";

$dx = $size_x + 0.8;


$ngrid_x = int(($max_x - $min_x) / $dgrid_x) + 1;
$ngrid_y = int(($max_y - $min_y) / $dgrid_y) + 1;

$cpt_vel = "app_vel_err.cpt";
$cpt_az = "az_err.cpt";

$txt_x = 0.15;
$txt_y = $size_y - 0.15;
$cpt_x = $size_x / 2.0;
$cpt_y = -1.4;
$cpt_len = $size_x;
$cpt_width = "0.3";
$slowness_length = 0.35;
$decimate = 1;

$txt_x2 = -1.9;
$txt_x3 = -0.2;
$txt_y2 = $size_y + 0.3;

$symbolsize = $size_x / $ngrid_x;

system "gmt set PS_LINE_JOIN round";
system "gmt set FONT_LABEL 11p,Helvetica";
system "gmt set FONT_ANNOT_PRIMARY 9p,Helvetica";
system "gmt set MAP_LABEL_OFFSET 5p";

for($index = $index_begin; $index <= $index_end; $index++){
  push @index_array, $index;
}

$pm = new Parallel::ForkManager ($MAX_PROCESS);

foreach $index (@index_array){

  $pm->start and next;

  $time_index = sprintf "%04d", $index;
  $current_sec = $ref_sec + ($index - 1) * $dt;
  $current_yr = $ref_yr;
  $current_mo = $ref_mo;
  $current_dy = $ref_dy;
  if($current_sec >= 24 * 60 * 60){
    $current_sec = $current_sec - 24 * 60 * 60;
    $current_dy = $current_dy + 1;
  }
  $current_hh = int($current_sec / 3600);
  $current_mm = int(($current_sec - $current_hh * 3600) / 60);
  $current_ss = $current_sec - 3600 * $current_hh - 60 * $current_mm;
  $current_hh = sprintf "%02d", $current_hh;
  $current_mm = sprintf "%02d", $current_mm;
  $current_ss = sprintf "%02d", $current_ss;

  $in2 = "$in_dir/slowness_gradiometry_${time_index}.dat";
  $out = "$in_dir/velocity_gradiometry_err_${time_index}.ps";
  print stderr "$out\n";

  if(-f $in2){
    $filesize = -s $in2;
    print stderr "$filesize\n";
    open IN, "<", $in2;
    $read_filesize = 0;
    while($read_filesize < $filesize){
      read IN, $buf, 4;
      push @x_east, (unpack "f", $buf);
      read IN, $buf, 4;
      push @y_north, (unpack "f", $buf);
      read IN, $buf, 4;
      push @slowness_x, (unpack "f", $buf);
      read IN, $buf, 4;
      push @slowness_y, (unpack "f", $buf);
      read IN, $buf, 4;
      push @sigma_slowness_x, (unpack "f", $buf);
      read IN, $buf, 4;
      push @sigma_slowness_y, (unpack "f", $buf);
      read IN, $buf, 4;
      push @ampterm_x, (unpack "f", $buf);
      read IN, $buf, 4;
      push @ampterm_y, (unpack "f", $buf);
      read IN, $buf, 4;
      push @sigma_ampterm_x, (unpack "f", $buf);
      read IN, $buf, 4;
      push @sigma_ampterm_y, (unpack "f", $buf);
      $read_filesize = $read_filesize + 4 * 10;
      }
      close IN;
    }
  
  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -K -P -X3c -Y6c > $out";

  system "gmt psbasemap -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Bpx${annot}+l\"Easting (km)\" \\
                        -Bpy${annot}+l\"Northing (km)\" -BWSen -O -K -P >> $out";
  if (-f $coastline){
    system "gmt psxy $coastline -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -W0.8p,black -O -K -P >> $out";
  }

  if (-f $in2){
    open OUT, " | gmt psxy -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Ss$symbolsize -W0.5p,black \\
                           -C$cpt_vel -O -K -P >> $out";
    for($i = 0; $i <= $#x_east; $i++){
      $slowness_sq = ($slowness_x[$i] ** 2 + $slowness_y[$i] ** 2);
      $velocity = sqrt(1.0 / $slowness_sq) * 1000.0;
        $velocity_err = 1000.0 / $slowness_sq * 
                      sqrt(($sigma_slowness_x[$i] ** 2 * $slowness_x[$i] ** 2 + $sigma_slowness_y[$i] ** 2 * $slowness_y[$i] ** 2)
                           / $slowness_sq);
        $velocity_err = $velocity_err / $velocity * 100.0;
        print OUT "$x_east[$i] $y_north[$i] $velocity_err\n";
      #print stdout "$x_east[$i] $y_north[$i] $velocity_err $velocity\n";
    }
    close OUT;
  }

  open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -Gwhite -O -K -P >> $out";
    print OUT "$txt_x $txt_y 12p,Helvetica,black 0 LT $current_yr/$current_mo/$current_dy $current_hh:$current_mm\n";
    print OUT "$txt_x2 $txt_y2 14p,Helvetica,black 0 LB (a)\n";
  close OUT;


  system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/${cpt_width}ch -Ef \\
                      -Ba10+l\"Relative error of velocity (\%)\" -C$cpt_vel -O -K -P >> $out";

  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -K -P -X$dx >> $out";

  system "gmt psbasemap -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Bpx${annot}+l\"Easting (km)\" -Bpy${annot} \\
                        -BSwen -O -K -P >> $out";
  if (-f $coastline){
    system "gmt psxy $coastline -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -W0.8p,black -O -K -P >> $out";
  }

  if (-f $in2){
    open OUT, " | gmt psxy -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Ss$symbolsize -W0.5p,black \\
                           -C$cpt_az -O -K -P >> $out";
    for($i = 0; $i <= $#x_east; $i++){
      $slowness_sq = $slowness_x[$i] ** 2 + $slowness_y[$i] ** 2;
      $az_err = 1.0 / $slowness_sq * 
                      sqrt(($sigma_slowness_x[$i] * $slowness_y[$i]) ** 2 + ($sigma_slowness_y[$i] * $slowness_x[$i]) ** 2);
      $az_err = rad2deg($az_err);
      print OUT "$x_east[$i] $y_north[$i] $az_err\n";
      #print stdout "$x_east[$i] $y_north[$i] $az_err\n";
    }
    close OUT;
  }


  open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -O -K -P >> $out";
    print OUT "$txt_x3 $txt_y2 14p,Helvetica,black 0 LB (b)\n";
  close OUT;

  system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/${cpt_width}ch -Ef -B+l\"Azimuth error (\\232)\" -C$cpt_az -O -K -P >> $out";

  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -P >> $out";
  system "gmt psconvert $out -Tg -A";

  $pm->finish;
}
$pm->wait_all_children;

unlink "gmt.history";
unlink "gmt.conf";

