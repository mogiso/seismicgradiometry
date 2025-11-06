#!/usr/bin/env perl
use Math::Trig qw(pi rad2deg deg2rad);
use Parallel::ForkManager;
$MAX_PROCESS = 10;

$in_dir = $ARGV[0];
$index_begin = $ARGV[1];
$index_end = $ARGV[2];
$coastline = $ARGV[3];
$simulation = $ARGV[4];
$tidestation = $ARGV[5];

$ref_yr = "2025";
$ref_mo = "7";
$ref_dy = 29;
$ref_hh = 0;
$ref_mm = 0;
$ref_ss = 0;
$dt = 60;
$ref_sec = $ref_hh * 3600 + $ref_mm * 60 + $ref_ss;

if($simulation == 1){
  $ref_sec = 0;
}
if($simulation == 0){
  $ref_sec = -(86400 + 8 * 3600 + 24 * 60 + 52.0); 
}

#S-net ref: 142.5E, 38.25N
$dgrid_x = 20;
$dgrid_y = 20;
$min_x = -350;
$max_x = 350;
$min_y = -600;
$max_y = 600;
$size_x = 3.8;
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

$dx = $size_x + 0.7;


$ngrid_x = int(($max_x - $min_x) / $dgrid_x) + 1;
$ngrid_y = int(($max_y - $min_y) / $dgrid_y) + 1;


$cpt = "amplitude_gradiometry_OBP.cpt";
if($simulation == 0){
  $cpt = "amplitude_gradiometry_OBP.cpt";
}
if($simulation == 1){
  $cpt = "amplitude_gradiometry_OBP_simulation.cpt";
}
#$cpt = "amplitude_gradiometry_OBP_simulation_240116.cpt";
#$cpt = "amplitude_gradiometry_OBP_simulation.cpt";
#$cpt = "amplitude_gradiometry_soratena.cpt";
#$cpt = "amplitude_gradiometry_OBP_231202.cpt";
$app_vel_cpt = "app_vel_gradiometry_OBP.cpt";
#$ampterm_cpt = "ampterm_OBP.cpt";
$ampterm_cpt = "ampterm_OBP_simulation2.cpt";
$tsunami_vel_grd = "tsunami_velocity_etopo.grd";

if (-f $tidestation){
  open IN, "<", $tidestation;
  while(<IN>){
    chomp $_;
    $_ =~ s/^\s*(.*?)\s*$/$1/;
    @tmp = split /\s+/, $_;
    push @tide_lon, $tmp[0];
    push @tide_lat, $tmp[1];
    push @tide_xeast, $tmp[2];
    push @tide_ynorth, $tmp[3];
  }
  close IN;
}

$txt_x = 0.1;
$txt_y = $size_y - 0.1;
$cpt_x = $size_x / 2.0;
$cpt_y = -1.4;
$cpt_len = $size_x;
$cpt_width = "0.3ch";
$slowness_length = 0.35;
$decimate = 2;

$txt_x2 = -1.4;
$txt_x3 = -0.2;
$txt_y2 = $size_y + 0.3;

system "gmt set PS_LINE_JOIN round";
system "gmt set FONT_LABEL 9p,Helvetica";
system "gmt set FONT_ANNOT_PRIMARY 9p,Helvetica";
system "gmt set MAP_LABEL_OFFSET 5p";
system "gmt set GMT_AUTO_DOWNLOAD off";
system "gmt set GMT_HISTORY false";

for($index = $index_begin; $index <= $index_end; $index++){
  push @index_array, $index;
}

$pm = new Parallel::ForkManager ($MAX_PROCESS);

foreach $index (@index_array){

  $pid = $pm->start and next;


  $time_index = sprintf "%04d", $index;
  $current_sec = $ref_sec + ($index - 1) * $dt;

  if($simulation != 1 && $simulation != 0){
    $current_yr = $ref_yr;
    $current_mo = $ref_mo;
    $current_dy = $ref_dy;
    while($current_sec >= 24 * 60 * 60){
      $current_sec = $current_sec - 24 * 60 * 60;
      $current_dy = $current_dy + 1;
    }
    if($current_dy > 31){
      $current_mo = $current_mo + 1;
      $current_dy = $current_dy - 31;
    }
  }

  $current_hh = int($current_sec / 3600);
  $current_mm = int(($current_sec - $current_hh * 3600) / 60);
  $current_ss = $current_sec - 3600 * $current_hh - 60 * $current_mm;
  $current_hh = sprintf "%02d", $current_hh;
  $current_mm = sprintf "%02d", $current_mm;
  $current_ss = sprintf "%02d", $current_ss;

  $in = "$in_dir/amplitude_gradiometry_${time_index}.grd";
  $in2 = "$in_dir/slowness_gradiometry_${time_index}.dat";
  $in3 = "$in_dir/velocity_ratio_${time_index}.grd";
  $out = $in;
  $out =~ s/\.grd$/\.ps/;
  print stderr "$out\n";

  
  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -K -X2c -Y6c -P > $out";

  system "gmt grdimage $in -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -C$cpt -O -K >> $out";
  system "gmt psbasemap -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Bpx${annot}+l\"Easting (km)\" \\
                        -Bpy${annot}+l\"Northing (km)\" -BWSen -O -K >> $out";
  if (-f $coastline){
    system "gmt psxy $coastline -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -W0.4p,black -O -K >> $out";
  }
  if (-f $tidestation){
    open OUT, " | gmt psxy -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Si0.25 -W0.5p,whitesmoke -Gblack -O -K >> $out";
    for($i = 0; $i <= $#tide_xeast; $i++){
      print OUT "$tide_xeast[$i] $tide_ynorth[$i]\n";
    }
    close OUT;
  }
    if($simulation == 1){
      open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -M -N -F+f+a+j -Gwhite -O -K >> $out";
        print OUT "> $txt_x $txt_y 9p,Helvetica,black 0 LT 10p 2.23c l\nTime: ${current_hh}hr ${current_mm}m\nSynthetic\n";
      close OUT;
    }elsif($simulation == 0){
      open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -M -N -F+f+a+j -Gwhite -O -K >> $out";
        print OUT "> $txt_x $txt_y 9p,Helvetica,black 0 LT 10p 2.46c l\nTime: ${current_hh}hr ${current_mm}m\nObservation\n";
      close OUT;
    }else{
      open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -Gwhite -O -K >> $out";
        print OUT "$txt_x $txt_y 9p,Helvetica,black 0 LT $current_yr/$current_mo/$current_dy $current_hh:$current_mm\n";
      close OUT;
    } 
  open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -Gwhite -O -K >> $out";
    print OUT "$txt_x2 $txt_y2 12p,Helvetica,black 0 LB (a)\n";
  close OUT;

  if (-f "$in_dir/station_location.txt"){
    system "gmt psxy $in_dir/station_location.txt -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Sc0.1 -W0.5p,black -O -K >> $out";
  }

  if($simulation == 0){
    system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -B+l\"Amplitude (hPa)\" -C$cpt -O -K >> $out";
  }elsif($simulation == 1){
    system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -B+l\"Amplitude (m)\" -C$cpt -O -K >> $out";
  }else{
    system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -B+l\"Amplitude (hPa)\" -C$cpt -O -K >> $out";
  }

  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -K -X$dx >> $out";

  system "gmt psbasemap -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Bpx${annot}+l\"Easting (km)\" -Bpy${annot} \\
                        -BSwen -O -K >> $out";
  if (-f $coastline){
    system "gmt psxy $coastline -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -W0.4p,black -O -K >> $out";
  }

  #if (-f "$in_dir/station_location.txt"){
  #  open OUT, " | gmt psxy -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Sc0.12 -W0.5p,black -O -K >> $out";
  #  open IN, "<", "$in_dir/station_location.txt";
  #  while(<IN>){
  #    chomp $_;
  #    $_ =~ s/^\s*(.*?)\s*$/$1/;
  #    @tmp = split /\s+/, $_;
  #    $app_vel = sqrt(9.8 * $tmp[4] * 1000.0);
  #    print OUT "$tmp[0] $tmp[1] $app_vel\n";
  #  }
  #  close IN;
  #  close OUT;
  #}

  if(-f $in2){
    $filesize = -s $in2;
    $read_filesize = 0;
    @x_east_tmp = ();
    @y_north_tmp = ();
    @amp_radiation = ();
    @theta_radiation = ();
    open OUT, " | gmt psxy -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y \\
                  -SVb0.12c/0.16c/0.16c -W0.5p,black -C$app_vel_cpt -O -K >> $out";
    open IN, "<", $in2;
    $count_x = 0;
    $count_y = 0;
    while($read_filesize < $filesize){
      $count_x++;
      if($count_x == $ngrid_x){
        $count_y++;
        $count_x = 0;
      }
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
      if($slowness_x != 0.0 && $slowness_y != 0.0){
        $app_vel = 1.0 / sqrt($slowness_x * $slowness_x + $slowness_y * $slowness_y) * 1000.0;
        $direction = atan2($slowness_x, $slowness_y);
        $amp_radiation_tmp = $ampterm_x * sin($direction) + $ampterm_y * cos($direction);
        $amp_radiation_tmp = $amp_radiation_tmp * 100.0;
        $theta_radiation_tmp = $ampterm_x * cos($direction) - $ampterm_y * sin($direction);
        $theta_radiation_tmp = $theta_radiation_tmp * 100.0;
        if($amp_radiation_tmp != 0.0){
          push @x_east_tmp, $x_east;
          push @y_north_tmp, $y_north;
          push @amp_radiation, $amp_radiation_tmp;
          push @theta_radiation, $theta_radiation_tmp;
        }
        next if (int(($x_east - $min_x) / $dgrid_x + 0.5) % $decimate != 0);
        next if (int(($y_north - $min_y) / $dgrid_y + 0.5) % $decimate != 0);
        $direction = rad2deg($direction);
        print OUT "$x_east $y_north $app_vel $direction $slowness_length\n";
      }
    }
    close IN;
    close OUT;

  }
  if (-f $tidestation){
    open OUT, " | gmt psxy -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Si0.25 -W0.5p,whitesmoke -Gblack -O -K >> $out";
    for($i = 0; $i <= $#tide_xeast; $i++){
      print OUT "$tide_xeast[$i] $tide_ynorth[$i]\n";
    }
    close OUT;
  }

  open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -O -K >> $out";
    print OUT "$txt_x3 $txt_y2 14p,Helvetica,black 0 LB (b)\n";
  close OUT;

  system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -B+l\"Apparent velocity (m/s)\" -C$app_vel_cpt -O -K >> $out";



  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -K -X$dx >> $out";
  #if(-f $in3){
  if(@x_east_tmp){
    $grdfile = "${time_index}_tmp.grd";
    open OUT, " | gmt xyz2grd -G$grdfile -I$dgrid_x/$dgrid_y -R$min_x/$max_x/$min_y/$max_y -di";
    for($i = 0; $i <= $#x_east_tmp; $i++){
      print OUT "$x_east_tmp[$i] $y_north_tmp[$i] $amp_radiation[$i]\n";
    }
    close OUT;
    #system "gmt grdimage $in3 -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -C$vel_ratio_cpt -O -K >> $out";
    system "gmt grdimage $grdfile -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -C$ampterm_cpt -O -K >> $out";
    unlink "$grdfile";
   }
  if (-f $tsunami_vel_grd){
    system "gmt grdcontour $tsunami_vel_grd -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y \\
                                            -Ccontour.txt -O -K >> $out";
  }

  system "gmt psbasemap -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Bpx${annot}+l\"Easting (km)\" \\
                        -Bpy${annot}+l\"Northing (km)\" -BwSen -O -K >> $out";
  #if (-f "$in_dir/station_location.txt"){
  #  system "gmt psxy $in_dir/station_location.txt -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Sc0.1 -W0.5p,dimgray -O -K >> $out";
  #}
  if (-f $coastline){
    system "gmt psxy $coastline -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -W0.4p,black -O -K >> $out";
  }
  if (-f $tidestation){
    open OUT, " | gmt psxy -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Si0.25 -W0.5p,whitesmoke -Gblack -O -K >> $out";
    for($i = 0; $i <= $#tide_xeast; $i++){
      print OUT "$tide_xeast[$i] $tide_ynorth[$i]\n";
    }
    close OUT;
  }
  open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -Gwhite -O -K >> $out";
    print OUT "$txt_x3 $txt_y2 14p,Helvetica,black 0 LB (c)\n";
  close OUT;

  system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -Ba1f0.5g0.5+l\"Geometrical spreading (10\@+-2\@+ km\@+-1\@+)\" \\
                       -C$ampterm_cpt -O -K >> $out";


  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -K -X$dx >> $out";
  if(@x_east_tmp){
    $grdfile = "${time_index}_tmp.grd";
    open OUT, " | gmt xyz2grd -G$grdfile -I$dgrid_x/$dgrid_y -R$min_x/$max_x/$min_y/$max_y -di";
    for($i = 0; $i <= $#x_east_tmp; $i++){
      print OUT "$x_east_tmp[$i] $y_north_tmp[$i] $theta_radiation[$i]\n";
    }
    close OUT;
    #system "gmt grdimage $in3 -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -C$vel_ratio_cpt -O -K >> $out";
    system "gmt grdimage $grdfile -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -C$ampterm_cpt -O -K >> $out";
    unlink "$grdfile";
   }
  if (-f $tsunami_vel_grd){
    system "gmt grdcontour $tsunami_vel_grd -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y \\
                                            -Ccontour.txt -A100 -O -K >> $out";
  }

  system "gmt psbasemap -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Bpx${annot}+l\"Easting (km)\" \\
                        -Bpy${annot}+l\"Northing (km)\" -BwSen -O -K >> $out";
  #if (-f "$in_dir/station_location.txt"){
  #  system "gmt psxy $in_dir/station_location.txt -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Sc0.1 -W0.5p,dimgray -O -K >> $out";
  #}
  if (-f $coastline){
    system "gmt psxy $coastline -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -W0.4p,black -O -K >> $out";
  }
  if (-f $tidestation){
    open OUT, " | gmt psxy -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Si0.25 -W0.5p,whitesmoke -Gblack -O -K >> $out";
    for($i = 0; $i <= $#tide_xeast; $i++){
      print OUT "$tide_xeast[$i] $tide_ynorth[$i]\n";
    }
    close OUT;
  }
  open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -Gwhite -O -K >> $out";
    print OUT "$txt_x3 $txt_y2 14p,Helvetica,black 0 LB (d)\n";
  close OUT;

  system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -Ba1f0.5g0.5+l\"Radiation pattern (10\@+-2\@+ km\@+-1\@+)\" \\
                       -C$ampterm_cpt -O -K >> $out";

  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O >> $out";
  system "gmt psconvert $out -Tg -A";

  $pm->finish;
}
$pm->wait_all_children;

unlink "gmt.history";
unlink "gmt.conf";

