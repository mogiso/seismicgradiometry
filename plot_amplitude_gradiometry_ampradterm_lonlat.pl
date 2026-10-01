#!/usr/bin/env perl
use Math::Trig qw(pi rad2deg deg2rad);
use Parallel::ForkManager;
$MAX_PROCESS = 10;

$in_dir = $ARGV[0];
$index_begin = $ARGV[1];
$index_end = $ARGV[2];
$simulation = $ARGV[3];
$tidestation = $ARGV[4];
$period = $ARGV[5];

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
$dgrid_lon = 0.2;
$dgrid_lat = 0.2;
$grdlon_w = 131.0;
$grdlon_e = 147.0;
$grdlat_s = 30.0;
$grdlat_n = 44.0;
$lon_w = 131.0;
$lon_e = 147.0;
$lat_s = 30.0;
$lat_n = 44.0;
$size_x = 8.0;
$size_y = `echo $lon_e $lat_n | gmt mapproject -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n`;
@tmp = split /\s+/, $size_y;
$size_y = $tmp[1];
$annot = "a5f1";

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

$dx = $size_x + 1.0;
$dy = $size_y + 3.5;


$ngrid_lon = int(($grdlon_e - $grdlon_w) / $dgrid_lon + 0.5) + 1;
$ngrid_lat = int(($grdlat_n - $grdlat_s) / $dgrid_lat + 0.5) + 1;


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
$ampterm_cpt = "ampterm_OBP_simulation.cpt";
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
$cpt_y = -1.0;
$cpt_len = $size_x;
$cpt_width = "0.25ch";
$slowness_length = 0.35;
$decimate = 2;

$txt_x2 = -0.7;
$txt_x3 = -0.2;
$txt_y2 = $size_y + 0.5;
$txt_x_period = $size_x - 0.1;
$txt_y_period = 0.125;

system "gmt set PS_LINE_JOIN round";
system "gmt set FONT_LABEL 9p,Helvetica";
system "gmt set FONT_ANNOT 9p,Helvetica";
system "gmt set MAP_LABEL_OFFSET 5p";
system "gmt set GMT_AUTO_DOWNLOAD off";
system "gmt set GMT_HISTORY false";
system "gmt set MAP_FRAME_PEN thick,black";

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

  
  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -K -X2c -Y19c -P > $out";

  system "gmt grdimage $in -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -C$cpt -O -K >> $out";
  system "gmt psbasemap -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Bpx${annot}+l\"Longitude\" \\
                        -Bpy${annot}+l\"Latitude\" -BWSen -O -K >> $out";
  system "gmt pscoast -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Dh -W0.4p,black -O -K -P >> $out";
  if (-f $tidestation){
    open OUT, " | gmt psxy -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Si0.25 -W0.5p,whitesmoke -Gblack -O -K >> $out";
    for($i = 0; $i <= $#tide_lon; $i++){
      print OUT "$tide_lon[$i] $tide_lat[$i]\n";
    }
    close OUT;
  }
  if($simulation == 1){
    open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -M -N -F+f+a+j -Gwhite -O -K >> $out";
      print OUT "> $txt_x $txt_y 9p,Helvetica,black 0 LT 9p 2.23c l\nTime: ${current_hh}hr ${current_mm}m\nSynthetic\n";
    close OUT;
  }elsif($simulation == 0){
    open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -M -N -F+f+a+j -Gwhite -O -K >> $out";
      print OUT "> $txt_x $txt_y 9p,Helvetica,black 0 LT 9p 2.44c l\nTime: ${current_hh}hr ${current_mm}m\nObservation\n";
    close OUT;
  }else{
    open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -Gwhite -O -K >> $out";
      print OUT "$txt_x $txt_y 9p,Helvetica,black 0 LT $current_yr/$current_mo/$current_dy $current_hh:$current_mm\n";
    close OUT;
  } 
  if($period){
    open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -Gwhite -O -K >> $out";
      print OUT "$txt_x_period $txt_y_period 9p,Helvetica,black 0 RB Period: $period\n";
    close OUT;
  }
  open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -Gwhite -O -K >> $out";
    print OUT "$txt_x2 $txt_y2 12p,Helvetica,black 0 LB (a)\n";
  close OUT;

  if (-f "$in_dir/station_location.txt"){
    system "gmt psxy $in_dir/station_location.txt -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
                                                  -Sc0.1 -W0.5p,black -O -K >> $out";
  }

  if($simulation == 0){
    system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -B+l\"Amplitude (hPa)\" -C$cpt -O -K >> $out";
  }elsif($simulation == 1){
    system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -B+l\"Amplitude (m)\" -C$cpt -O -K >> $out";
  }else{
    system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -B+l\"Amplitude (hPa)\" -C$cpt -O -K >> $out";
  }

  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -K -X$dx >> $out";

  system "gmt psbasemap -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Bpx${annot}+l\"Lontigude\" -Bpy${annot} \\
                        -BSwen -O -K >> $out";
  system "gmt pscoast -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Dh -W0.4p,black -O -K -P >> $out";

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
    @lon_tmp = ();
    @lat_tmp = ();
    @amp_radiation = ();
    @theta_radiation = ();
    open OUT, " | gmt psxy -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
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
      $lon = unpack "f", $buf;
      read IN, $buf, 4;
      $lat = unpack "f", $buf;
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
        #if($amp_radiation_tmp != 0.0){
        #  push @lon_tmp, $lon;
        #  push @lat_tmp, $lat;
        #  push @amp_radiation, $amp_radiation_tmp;
        #  push @theta_radiation, $theta_radiation_tmp;
        #}
        next if (int(($lon - $lon_w) / $dgrid_lon + 0.5) % $decimate != 0);
        next if (int(($lat - $lat_s) / $dgrid_lat + 0.5) % $decimate != 0);
        $direction = rad2deg($direction);
        print OUT "$lon $lat $app_vel $direction $slowness_length\n";
      }
    }
    close IN;
    close OUT;

  }
  if (-f $tidestation){
    open OUT, " | gmt psxy -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Si0.25 -W0.5p,whitesmoke -Gblack -O -K >> $out";
    for($i = 0; $i <= $#tide_xeast; $i++){
      print OUT "$tide_lon[$i] $tide_lat[$i]\n";
    }
    close OUT;
  }

  open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -O -K >> $out";
    print OUT "$txt_x3 $txt_y2 14p,Helvetica,black 0 LB (b)\n";
  close OUT;

  system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -B+l\"Apparent velocity (m/s)\" -C$app_vel_cpt -O -K >> $out";



  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -K -X-$dx -Y-$dy >> $out";
  #if(@lon_tmp){
    $grdfile = "geospread_gradiometry_${time_index}.grd";
    #$grdfile = "${time_index}_tmp.grd";
    #open OUT, " | gmt xyz2grd -G$grdfile -I$dgrid_lon/$dgrid_lat -R$lon_w/$lon_e/$lat_s/$lat_n -di";
    #for($i = 0; $i <= $#lon_tmp; $i++){
    #  print OUT "$lon_tmp[$i] $lat_tmp[$i] $amp_radiation[$i]\n";
    #}
    close OUT;
    #system "gmt grdimage $in3 -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -C$vel_ratio_cpt -O -K >> $out";
    system "gmt grdimage $grdfile -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -C$ampterm_cpt -O -K >> $out";
    #unlink "$grdfile";
  #}
  #if (-f $tsunami_vel_grd){
  #  system "gmt grdcontour $tsunami_vel_grd -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y \\
  #                                          -Ccontour.txt -O -K >> $out";
  #}

  system "gmt psbasemap -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Bpx${annot}+l\"Longitude\" \\
                        -Bpy${annot}+l\"Latitude\" -BWSen -O -K >> $out";
  system "gmt pscoast -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Dh -W0.4p,black -O -K -P >> $out";
  #if (-f "$in_dir/station_location.txt"){
  #  system "gmt psxy $in_dir/station_location.txt -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Sc0.1 -W0.5p,dimgray -O -K >> $out";
  #}
  #if (-f $coastline){
  #  system "gmt psxy $coastline -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -W0.4p,black -O -K >> $out";
  #}
  if (-f $tidestation){
    open OUT, " | gmt psxy -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Si0.25 -W0.5p,whitesmoke -Gblack -O -K >> $out";
    for($i = 0; $i <= $#tide_lon; $i++){
      print OUT "$tide_lon[$i] $tide_lat[$i]\n";
    }
    close OUT;
  }
  open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -Gwhite -O -K >> $out";
    print OUT "$txt_x2 $txt_y2 14p,Helvetica,black 0 LB (c)\n";
  close OUT;

  system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -Ba1f0.5g0.5+l\"Geometrical spreading (10\@+-2\@+ km\@+-1\@+)\" \\
                       -C$ampterm_cpt -O -K >> $out";


  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -K -X$dx >> $out";
  #if(@lon_tmp){
    $grdfile = "radpattern_gradiometry_${time_index}.grd";
    #$grdfile = "${time_index}_tmp.grd";
    #open OUT, " | gmt xyz2grd -G$grdfile -I$dgrid_lon/$dgrid_lat -R$lon_w/$lon_e/$lat_s/$lat_n -di";
    #for($i = 0; $i <= $#lon_tmp; $i++){
    #  print OUT "$lon_tmp[$i] $lat_tmp[$i] $theta_radiation[$i]\n";
    #}
    #close OUT;
    #system "gmt grdimage $in3 -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -C$vel_ratio_cpt -O -K >> $out";
    system "gmt grdimage $grdfile -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -C$ampterm_cpt -O -K >> $out";
    #unlink "$grdfile";
  #}
  #if (-f $tsunami_vel_grd){
  #  system "gmt grdcontour $tsunami_vel_grd -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y \\
  #                                          -Ccontour.txt -A100 -O -K >> $out";
  #}

  system "gmt psbasemap -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Bpx${annot}+l\"Easting (km)\" \\
                        -Bpy${annot}+l\"Northing (km)\" -BwSen -O -K >> $out";
  system "gmt pscoast -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Dh -W0.4p,black -O -K -P >> $out";
  #if (-f "$in_dir/station_location.txt"){
  #  system "gmt psxy $in_dir/station_location.txt -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Sc0.1 -W0.5p,dimgray -O -K >> $out";
  #}
  #if (-f $coastline){
  #  system "gmt psxy $coastline -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -W0.4p,black -O -K >> $out";
  #}
  if (-f $tidestation){
    open OUT, " | gmt psxy -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Si0.25 -W0.5p,whitesmoke -Gblack -O -K >> $out";
    for($i = 0; $i <= $#tide_xeast; $i++){
      print OUT "$tide_lon[$i] $tide_lat[$i]\n";
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

