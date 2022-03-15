#!/usr/bin/env perl
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
$ref_hh = 18;
$ref_mm = 30;
$ref_ss = 0;
$dt = 60;

$ref_sec = $ref_hh * 3600 + $ref_mm * 60 + $ref_ss;

$min_x = -150;
$max_x = 150;
$min_y = -100;
$max_y = 100;

$size_x = 8;
$size_y = 5;
$dx = $size_x + 0.8;


$cpt = "amplitude_gradiometry_OBP.cpt";
$app_vel_cpt = "app_vel_gradiometry_OBP.cpt";

$txt_x = 0.0;
$txt_y = $size_y + 0.3;
$cpt_x = $size_x / 2.0;
$cpt_y = -1.7;
$cpt_len = $size_x;
$cpt_width = "0.3ch";
$slowness_length = 0.3;
$decimate = 1;
$dgrid_x = 10;
$dgrid_y = 10;

system "gmt set PS_LINE_JOIN round";
system "gmt set FONT_LABEL 13p,Helvetica";
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

  $in = "$in_dir/amplitude_gradiometry_${time_index}.grd";
  $in2 = "$in_dir/slowness_gradiometry_${time_index}.dat";
  $out = $in;
  $out =~ s/\.grd$/\.ps/;
  print stderr "$out\n";

  
  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -K -P -X3c -Y6c > $out";

  system "gmt grdimage $in -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -C$cpt -O -K -P >> $out";
  system "gmt psbasemap -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Bpxa100f50+l\"Easting (km)\" \\
                        -Bpya100f50+l\"Northing (km)\" -BWSen -O -K -P >> $out";
  if (-f $coastline){
    system "gmt psxy $coastline -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -W0.8p,black -O -K -P >> $out";
  }
  open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -O -K -P >> $out";
    print OUT "$txt_x $txt_y 14p,Helvetica,black 0 LB $current_yr/$current_mo/$current_dy $current_hh:$current_mm:$current_ss\n";
  close OUT;

  if (-f "station_location.txt"){
    system "gmt psxy station_location.txt -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Sc0.1 -W0.5p,black -O -K -P >> $out";
  }

  system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -B+l\"Amplitude (hPa)\" -C$cpt -O -K -P >> $out";

  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -K -P -X$dx >> $out";

  system "gmt psbasemap -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Bpxa100f50+l\"Easting (km)\" -Bpya100f50 \\
                        -BSwen -O -K -P >> $out";
  if (-f $coastline){
    system "gmt psxy $coastline -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -W0.8p,black -O -K -P >> $out";
  }

  if(-f $in2){
    $filesize = -s $in2;
    $read_filesize = 0;
    open OUT, " | gmt psxy -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y \\
                  -Svb0.08c/0.12c/0.125c -W0.5p,black -C$app_vel_cpt -O -K -P >> $out";
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
      read IN, $buf, 4;
      $sigma_slowness_x = unpack "f", $buf;
      read IN, $buf, 4;
      $sigma_slowness_y = unpack "f", $buf;
      read IN, $buf, 4;
      $ampterm_x = unpack "f", $buf;
      read IN, $buf, 4;
      $ampterm_y = unpack "f", $buf;
      read IN, $buf, 4;
      $sigma_ampterm_x = unpack "f", $buf;
      read IN, $buf, 4;
      $sigma_ampterm_y = unpack "f", $buf;
      $read_filesize = $read_filesize + 4 * 10;
      if($slowness_x != 0.0 && $slowness_y != 0.0 && 
         $x_east % ($decimate * $dgrid_x) == 0 && $y_north % ($decimate * $dgrid_y) == 0){
        $app_vel = 1.0 / sqrt($slowness_x * $slowness_x + $slowness_y * $slowness_y) * 1000.0;
        $direction = rad2deg(atan2($slowness_y, $slowness_x));
        #print stdout "$x_east $y_north $direction $slowness_length $slowness\n";
        print OUT "$x_east $y_north $app_vel $direction $slowness_length\n";
      }
    }
    close IN;
    close OUT;

  }

  if (-f "station_location.txt"){
    open OUT, " | gmt psxy -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Sc0.1 -C$app_vel_cpt -W0.5p,black -O -K -P >> $out";
    open IN, "<", "station_location.txt";
    while(<IN>){
      chomp $_;
      $_ =~ s/^\s*(.*?)\s*$/$1/;
      @tmp = split /\s+/, $_;
      $app_vel = sqrt(9.8 * $tmp[4] * 1000.0);
      print OUT "$tmp[0] $tmp[1] $app_vel\n";
    }
    close IN;
    close OUT;
  }

  system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -Ba50+l\"Velocity (m/s)\" -C$app_vel_cpt -O -K -P >> $out";
  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -P >> $out";
  system "gmt psconvert $out -Tg -A";

  $pm->finish;
}
$pm->wait_all_children;

unlink "gmt.history";
unlink "gmt.conf";

