#!/usr/bin/env perl
use Math::Trig qw(pi rad2deg deg2rad);

$in_geospreaddiffgrd = $ARGV[0];
$in_radtermdiffgrd = $ARGV[1];
$in_geospreadradtermdiff_txt = $ARGV[2];
$coastline = $ARGV[3];
$obpstation = $ARGV[4];
$tidestation = $ARGV[5];
$out = $ARGV[6];


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

$dx = $size_x + 0.8;

$ngrid_x = int(($max_x - $min_x) / $dgrid_x) + 1;
$ngrid_y = int(($max_y - $min_y) / $dgrid_y) + 1;

$tsunami_vel_grd = "./tsunami_velocity_etopo.grd";
$cpt = "amplitude_diff.cpt";
if (-f $obpstation){
  open IN, "<", $obpstation;
  while(<IN>){
    chomp $_;
    $_ =~ s/^\s*(.*?)\s*$/$1/;
    @tmp = split /\s+/, $_;
    push @obpstation_xeast, $tmp[0];
    push @obpstation_ynorth, $tmp[1];
    push @obpstation_lon, $tmp[2];
    push @obpstation_lat, $tmp[3];
  }
  close IN;
}

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

$txt_x2 = -1.4;
$txt_x3 = -0.2;
$txt_y2 = $size_y + 0.3;

system "gmt set PS_LINE_JOIN round";
system "gmt set FONT_LABEL 9p,Helvetica";
system "gmt set FONT_ANNOT_PRIMARY 9p,Helvetica";
system "gmt set MAP_LABEL_OFFSET 5p";
system "gmt set GMT_AUTO_DOWNLOAD off";
system "gmt set GMT_HISTORY false";


open OUT, " | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -K -X2c -Y15c -P > $out";
close OUT;

system "gmt grdimage $in_geospreaddiffgrd -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -C$cpt -O -K >> $out";
if (-f $tsunami_vel_grd){
  system "gmt grdcontour $tsunami_vel_grd -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y \\
                                          -Ccontour.txt -O -K >> $out";
}

system "gmt psbasemap -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -Bpx${annot}+l\"Easting (km)\" \\
                      -Bpy${annot}+l\"Northing (km)\" -BWSen -O -K >> $out";
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
  print OUT "$txt_x2 $txt_y2 14p,Helvetica,black 0 LB (a)\n";
close OUT;

system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -B+l\"Difference (10\@+-2\@+ km\@+-1\@+)\" \\
                     -C$cpt -O -K >> $out";


system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -K -X$dx >> $out";
system "gmt grdimage $in_radtermdiffgrd -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -C$cpt -O -K >> $out";
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
  print OUT "$txt_x3 $txt_y2 14p,Helvetica,black 0 LB (b)\n";
close OUT;

system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len/$cpt_width -B+l\"Difference (10\@+-2\@+ km\@+-1\@+)\" \\
                     -C$cpt -O -K >> $out";


$size_x = 3.8;
$size_y = 3.8;
$txt_x = 0.1;
$txt_y = $size_y - 0.1;

$txt_x2 = -1.4;
$txt_x3 = -0.2;
$txt_y2 = $size_y + 0.3;
$bin = 0.2;
$min_x = -1.6;
$max_x = 1.6;
$min_y = 0;
$max_y = 160;
$annot_x = "a0.8+l\"Difference\"";
$annot_y = "a40+l\"Counts\"";
system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -K -X-$dx -Y-7.5c >> $out";

open OUT, " | gmt pshistogram -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -T$bin -Ggray \\
                              -Bx$annot_x -By$annot_y -BWSen -F -W0.5p,black -O -K -P >> $out";
open IN, "<", $in_geospreadradtermdiff_txt;
while(<IN>){
  chomp $_;
  $_ =~ s/^\s*(.*?)\s*$/$1/;
  @tmp = split /\s+/, $_;
  print OUT "$tmp[0]\n";
}
close IN;
close OUT;
open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -F+f+a+j -N -O -K -P >> $out";
print OUT "$txt_x2 $txt_y2 14p,Helvetica,black 0 LB (c)\n";
close OUT;

system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -K -X$dx >> $out";
open OUT, " | gmt pshistogram -JX$size_x/$size_y -R$min_x/$max_x/$min_y/$max_y -T$bin -Ggray \\
                              -Bx$annot_x -By$annot_y -BwSen -F -W0.5p,black -O -K -P >> $out";
open IN, "<", $in_geospreadradtermdiff_txt;
while(<IN>){
  chomp $_;
  $_ =~ s/^\s*(.*?)\s*$/$1/;
  @tmp = split /\s+/, $_;
  print OUT "$tmp[1]\n";
}
close IN;
close OUT;
open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -F+f+a+j -N -O -K -P >> $out";
print OUT "$txt_x3 $txt_y2 14p,Helvetica,black 0 LB (d)\n";
close OUT;





system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O >> $out";
system "gmt psconvert $out -Tg -A";


unlink "gmt.history";
unlink "gmt.conf";

