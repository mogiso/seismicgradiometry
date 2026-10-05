#!/usr/bin/env perl
use Math::Trig qw(pi rad2deg deg2rad);
use Parallel::ForkManager;
$MAX_PROCESS = 10;

$in_dir = $ARGV[0];
$index_begin = $ARGV[1];
$index_end = $ARGV[2];
$period = $ARGV[3];
$simulation = $ARGV[4];

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

$dx = $size_x + 0.7;

$dx = $size_x + 1.0;
$dy = $size_y + 3.5;


$ngrid_lon = int(($grdlon_e - $grdlon_w) / $dgrid_lon + 0.5) + 1;
$ngrid_lat = int(($grdlat_n - $grdlat_s) / $dgrid_lat + 0.5) + 1;


$cpt_a = "ampterm_OBP_simulation.cpt";
$cpt_a_err = "gradiometry_relerr.cpt";
$cpt_p = "slowness_gradiometry.cpt";
$cpt_p_err = "gradiometry_relerr.cpt";


$txt_x = 0.1;
$txt_y = $size_y - 0.1;
$cpt_x = $size_x / 2.0;
$cpt_y = -1.0;
$cpt_len = $size_x;
$cpt_width = "0.25ch";
$slowness_length = 0.35;

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

  $in = "$in_dir/slowness_gradiometry_${time_index}.dat";
  $out = "$in_dir/gradiometry_axay_err_${time_index}.ps";
  $out2 = "$in_dir/gradiometry_pxpy_err_${time_index}.ps";
  print stderr "$out $out2\n";


  $read_filesize = 0;
  @lon = ();
  @lat = ();
  @slowness_x = ();
  @slowness_y = ();
  @sigma_slowness_x = ();
  @sigma_slowness_y = ();
  @ampterm_x = ();
  @ampterm_y = ();
  @sigma_ampterm_x = ();
  @sigma_ampterm_y = ();
  if (-f $in){
    open IN, "<", $in;
    $filesize = -s $in;
    while($read_filesize < $filesize){
      read IN, $buf, 4;
      push @lon, (unpack "f", $buf);
      read IN, $buf, 4;
      push @lat, (unpack "f", $buf);
      read IN, $buf, 4;
      push @slowness_x, (unpack "f", $buf);
      read IN, $buf, 4;
      push @slowness_y, (unpack "f", $buf);
      read IN, $buf, 4;
      push @sigma_slowness_x, sqrt(unpack "f", $buf);
      read IN, $buf, 4;
      push @sigma_slowness_y, sqrt(unpack "f", $buf);
      read IN, $buf, 4;
      push @ampterm_x, ((unpack "f", $buf) * 100.0);
      read IN, $buf, 4;
      push @ampterm_y, ((unpack "f", $buf) * 100.0);
      read IN, $buf, 4;
      push @sigma_ampterm_x, (sqrt(unpack "f", $buf) * 100.0);
      read IN, $buf, 4;
      push @sigma_ampterm_y, (sqrt(unpack "f", $buf) * 100.0);
      $read_filesize = $read_filesize + 4 * 10;
    }
    close IN;
  }
  $grdfile = "_tmp_${time_index}.grd";

  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -K -X2c -Y18c -P > $out";
  ##(a): ax
  open OUT, " | gmt xyz2grd -G$grdfile -R$grdlon_w/$grdlon_e/$grdlat_s/$grdlat_n -I$dgrid_lon/$dgrid_lat -di";
  for($i = 0; $i <= $#lon; $i++){
    print OUT "$lon[$i] $lat[$i] $ampterm_x[$i]\n";
  }
  close OUT;
  
  system "gmt grdimage $grdfile -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -C$cpt_a -O -K >> $out";
  unlink "$grdfile";
  system "gmt pscoast -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Df -W1p,black \\
                      -Bpx${annot} -Bpy${annot} -BWSen -O -K >> $out";
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
    system "gmt psxy $in_dir/station_location.txt -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Sc0.1 -W0.5p,black -O -K >> $out";
  }
  system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -B+l\"A\@-x\@- (10\@+-2\@+ km\@+-1\@+)\" -C$cpt_a -O -K >> $out";

  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -K -X$dx >> $out";

  ##(b): ay
  open OUT, " | gmt xyz2grd -G$grdfile -R$grdlon_w/$grdlon_e/$grdlat_s/$grdlat_n -I$dgrid_lon/$dgrid_lat -di";
  for($i = 0; $i <= $#lon; $i++){
    print OUT "$lon[$i] $lat[$i] $ampterm_y[$i]\n";
  }
  close OUT;
  
  system "gmt grdimage $grdfile -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -C$cpt_a -O -K >> $out";
  unlink "$grdfile";
  system "gmt pscoast -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Df -W1p,black \\
                      -Bpx${annot} -Bpy${annot} -BwSen -O -K >> $out";
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
    print OUT "$txt_x2 $txt_y2 12p,Helvetica,black 0 LB (b)\n";
  close OUT;

  if (-f "$in_dir/station_location.txt"){
    system "gmt psxy $in_dir/station_location.txt -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Sc0.1 -W0.5p,black -O -K >> $out";
  }
  system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -B+l\"A\@-y\@- (10\@+-2\@+ km\@+-1\@+)\" -C$cpt_a -O -K >> $out";


  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -K -X-$dx -Y-$dy >> $out";

  ##(c): ax
  open OUT, " | gmt xyz2grd -G$grdfile -R$grdlon_w/$grdlon_e/$grdlat_s/$grdlat_n -I$dgrid_lon/$dgrid_lat -di";
  for($i = 0; $i <= $#lon; $i++){
    $relval = abs($sigma_ampterm_x[$i] / $ampterm_x[$i]) * 100.0;
    print OUT "$lon[$i] $lat[$i] $relval\n";
  }
  close OUT;
  
  system "gmt grdimage $grdfile -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -C$cpt_a_err -O -K >> $out";
  unlink "$grdfile";
  system "gmt pscoast -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Df -W1p,black \\
                      -Bpx${annot} -Bpy${annot} -BWSen -O -K >> $out";
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
    print OUT "$txt_x2 $txt_y2 12p,Helvetica,black 0 LB (c)\n";
  close OUT;

  if (-f "$in_dir/station_location.txt"){
    system "gmt psxy $in_dir/station_location.txt -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Sc0.1 -W0.5p,black -O -K >> $out";
  }
  system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -B+l\"A\@-x\@- Relative error (%)\" -C$cpt_a_err -O -K >> $out";


  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -K -X$dx >> $out";

  ##(d): ay error
  open OUT, " | gmt xyz2grd -G$grdfile -R$grdlon_w/$grdlon_e/$grdlat_s/$grdlat_n -I$dgrid_lon/$dgrid_lat -di";
  for($i = 0; $i <= $#lon; $i++){
    $relval = abs($sigma_ampterm_y[$i] / $ampterm_y[$i]) * 100.0;
    print OUT "$lon[$i] $lat[$i] $relval\n";
  }
  close OUT;
  
  system "gmt grdimage $grdfile -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -C$cpt_a_err -O -K >> $out";
  unlink "$grdfile";
  system "gmt pscoast -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Df -W1p,black \\
                      -Bpx${annot} -Bpy${annot} -BwSen -O -K >> $out";
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
    print OUT "$txt_x2 $txt_y2 12p,Helvetica,black 0 LB (d)\n";
  close OUT;

  if (-f "$in_dir/station_location.txt"){
    system "gmt psxy $in_dir/station_location.txt -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Sc0.1 -W0.5p,black -O -K >> $out";
  }
  system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -B+l\"A\@-y\@- Relative error (%)\" -C$cpt_a_err -O -K >> $out";

  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -P >> $out";
  system "gmt psconvert $out -Tg -A";

  ##slowness
  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -K -X2c -Y18c -P > $out2";
  ##(a): px
  open OUT, " | gmt xyz2grd -G$grdfile -R$grdlon_w/$grdlon_e/$grdlat_s/$grdlat_n -I$dgrid_lon/$dgrid_lat -di";
  for($i = 0; $i <= $#lon; $i++){
    print OUT "$lon[$i] $lat[$i] $slowness_x[$i]\n";
  }
  close OUT;
  
  system "gmt grdimage $grdfile -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -C$cpt_p -O -K >> $out2";
  unlink "$grdfile";
  system "gmt pscoast -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Df -W1p,black \\
                      -Bpx${annot} -Bpy${annot} -BWSen -O -K >> $out2";
  if($simulation == 1){
    open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -M -N -F+f+a+j -Gwhite -O -K >> $out2";
      print OUT "> $txt_x $txt_y 9p,Helvetica,black 0 LT 9p 2.23c l\nTime: ${current_hh}hr ${current_mm}m\nSynthetic\n";
    close OUT;
  }elsif($simulation == 0){
    open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -M -N -F+f+a+j -Gwhite -O -K >> $out2";
      print OUT "> $txt_x $txt_y 9p,Helvetica,black 0 LT 9p 2.44c l\nTime: ${current_hh}hr ${current_mm}m\nObservation\n";
    close OUT;
  }else{
    open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -Gwhite -O -K >> $out2";
      print OUT "$txt_x $txt_y 9p,Helvetica,black 0 LT $current_yr/$current_mo/$current_dy $current_hh:$current_mm\n";
    close OUT;
  } 
  if($period){
    open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -Gwhite -O -K >> $out2";
      print OUT "$txt_x_period $txt_y_period 9p,Helvetica,black 0 RB Period: $period\n";
    close OUT;
  }
  open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -Gwhite -O -K >> $out2";
    print OUT "$txt_x2 $txt_y2 12p,Helvetica,black 0 LB (a)\n";
  close OUT;

  if (-f "$in_dir/station_location.txt"){
    system "gmt psxy $in_dir/station_location.txt -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Sc0.1 -W0.5p,black -O -K >> $out2";
  }
  system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -B+l\"p\@-x\@- (km\@+-1\@+)\" -C$cpt_p -O -K >> $out2";

  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -K -X$dx >> $out2";

  ##(b): py
  open OUT, " | gmt xyz2grd -G$grdfile -R$grdlon_w/$grdlon_e/$grdlat_s/$grdlat_n -I$dgrid_lon/$dgrid_lat -di";
  for($i = 0; $i <= $#lon; $i++){
    print OUT "$lon[$i] $lat[$i] $slowness_y[$i]\n";
  }
  close OUT;
  
  system "gmt grdimage $grdfile -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -C$cpt_p -O -K >> $out2";
  unlink "$grdfile";
  system "gmt pscoast -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Df -W1p,black \\
                      -Bpx${annot} -Bpy${annot} -BwSen -O -K >> $out2";
  if($simulation == 1){
    open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -M -N -F+f+a+j -Gwhite -O -K >> $out2";
      print OUT "> $txt_x $txt_y 9p,Helvetica,black 0 LT 9p 2.23c l\nTime: ${current_hh}hr ${current_mm}m\nSynthetic\n";
    close OUT;
  }elsif($simulation == 0){
    open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -M -N -F+f+a+j -Gwhite -O -K >> $out2";
      print OUT "> $txt_x $txt_y 9p,Helvetica,black 0 LT 9p 2.44c l\nTime: ${current_hh}hr ${current_mm}m\nObservation\n";
    close OUT;
  }else{
    open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -Gwhite -O -K >> $out2";
      print OUT "$txt_x $txt_y 9p,Helvetica,black 0 LT $current_yr/$current_mo/$current_dy $current_hh:$current_mm\n";
    close OUT;
  } 
  if($period){
    open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -Gwhite -O -K >> $out2";
      print OUT "$txt_x_period $txt_y_period 9p,Helvetica,black 0 RB Period: $period\n";
    close OUT;
  }
  open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -Gwhite -O -K >> $out2";
    print OUT "$txt_x2 $txt_y2 12p,Helvetica,black 0 LB (b)\n";
  close OUT;

  if (-f "$in_dir/station_location.txt"){
    system "gmt psxy $in_dir/station_location.txt -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Sc0.1 -W0.5p,black -O -K >> $out2";
  }
  system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -B+l\"p\@-y\@- (km\@+-1\@+)\" -C$cpt_p -O -K >> $out2";


  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -K -X-$dx -Y-$dy >> $out2";

  ##(c): ax
  open OUT, " | gmt xyz2grd -G$grdfile -R$grdlon_w/$grdlon_e/$grdlat_s/$grdlat_n -I$dgrid_lon/$dgrid_lat -di";
  for($i = 0; $i <= $#lon; $i++){
    $relval = abs($sigma_slowness_x[$i] / $slowness_x[$i]) * 100.0;
    print OUT "$lon[$i] $lat[$i] $relval\n";
  }
  close OUT;
  
  system "gmt grdimage $grdfile -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -C$cpt_p_err -O -K >> $out2";
  unlink "$grdfile";
  system "gmt pscoast -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Df -W1p,black \\
                      -Bpx${annot} -Bpy${annot} -BWSen -O -K >> $out2";
  if($simulation == 1){
    open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -M -N -F+f+a+j -Gwhite -O -K >> $out2";
      print OUT "> $txt_x $txt_y 9p,Helvetica,black 0 LT 9p 2.23c l\nTime: ${current_hh}hr ${current_mm}m\nSynthetic\n";
    close OUT;
  }elsif($simulation == 0){
    open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -M -N -F+f+a+j -Gwhite -O -K >> $out2";
      print OUT "> $txt_x $txt_y 9p,Helvetica,black 0 LT 9p 2.44c l\nTime: ${current_hh}hr ${current_mm}m\nObservation\n";
    close OUT;
  }else{
    open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -Gwhite -O -K >> $out2";
      print OUT "$txt_x $txt_y 9p,Helvetica,black 0 LT $current_yr/$current_mo/$current_dy $current_hh:$current_mm\n";
    close OUT;
  } 
  if($period){
    open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -Gwhite -O -K >> $out2";
      print OUT "$txt_x_period $txt_y_period 9p,Helvetica,black 0 RB Period: $period\n";
    close OUT;
  }
  open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -Gwhite -O -K >> $out2";
    print OUT "$txt_x2 $txt_y2 12p,Helvetica,black 0 LB (c)\n";
  close OUT;

  if (-f "$in_dir/station_location.txt"){
    system "gmt psxy $in_dir/station_location.txt -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Sc0.1 -W0.5p,black -O -K >> $out2";
  }
  system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -B+l\"p\@-x\@- Relative error (%)\" -C$cpt_a_err -O -K >> $out2";


  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -K -X$dx >> $out2";

  ##(d): ay error
  open OUT, " | gmt xyz2grd -G$grdfile -R$grdlon_w/$grdlon_e/$grdlat_s/$grdlat_n -I$dgrid_lon/$dgrid_lat -di";
  for($i = 0; $i <= $#lon; $i++){
    $relval = abs($sigma_slowness_y[$i] / $slowness_y[$i]) * 100.0;
    print OUT "$lon[$i] $lat[$i] $relval\n";
  }
  close OUT;
  
  system "gmt grdimage $grdfile -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -C$cpt_p_err -O -K >> $out2";
  unlink "$grdfile";
  system "gmt pscoast -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Df -W1p,black \\
                      -Bpx${annot} -Bpy${annot} -BwSen -O -K >> $out2";
  if($simulation == 1){
    open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -M -N -F+f+a+j -Gwhite -O -K >> $out2";
      print OUT "> $txt_x $txt_y 9p,Helvetica,black 0 LT 9p 2.23c l\nTime: ${current_hh}hr ${current_mm}m\nSynthetic\n";
    close OUT;
  }elsif($simulation == 0){
    open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -M -N -F+f+a+j -Gwhite -O -K >> $out2";
      print OUT "> $txt_x $txt_y 9p,Helvetica,black 0 LT 9p 2.44c l\nTime: ${current_hh}hr ${current_mm}m\nObservation\n";
    close OUT;
  }else{
    open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -Gwhite -O -K >> $out2";
      print OUT "$txt_x $txt_y 9p,Helvetica,black 0 LT $current_yr/$current_mo/$current_dy $current_hh:$current_mm\n";
    close OUT;
  } 
  if($period){
    open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -Gwhite -O -K >> $out2";
      print OUT "$txt_x_period $txt_y_period 9p,Helvetica,black 0 RB Period: $period\n";
    close OUT;
  }
  open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -N -F+f+a+j -Gwhite -O -K >> $out2";
    print OUT "$txt_x2 $txt_y2 12p,Helvetica,black 0 LB (d)\n";
  close OUT;

  if (-f "$in_dir/station_location.txt"){
    system "gmt psxy $in_dir/station_location.txt -JM$size_x -R$lon_w/$lon_e/$lat_s/$lat_n -Sc0.1 -W0.5p,black -O -K >> $out2";
  }
  system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -B+l\"p\@-y\@- Relative error (%)\" -C$cpt_p_err -O -K >> $out2";

  system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -P >> $out2";
  system "gmt psconvert $out2 -Tg -A";
  $pm->finish;
}
$pm->wait_all_children;

unlink "gmt.history";
unlink "gmt.conf";

