#!/usr/bin/env perl
# Copyright 2022 Masashi Ogiso (masashi.ogiso@gmail.com)
# Released under the MIT license.
# see https://opensource.org/licenses/MIT

use Math::Trig qw(pi rad2deg);

$in = $ARGV[0];

#S-net ref: 142.5E, 38.25N
#$dgrid_x = 20;
#$dgrid_y = 20;
#$min_x = -350;
#$max_x = 350;
#$min_y = -600;
#$max_y = 600;

#DONET ref: 135.75E, 33.2N
$dgrid_x = 10;
$dgrid_y = 10;
$min_x = -150;
$max_x = 150;
$min_y = -100;
$max_y = 100;


$ngrid_x = int(($max_x - $min_x) / $dgrid_x) + 1;
$ngrid_y = int(($max_y - $min_y) / $dgrid_y) + 1;

$decimate = 1;

if(-f $in){
  $filesize = -s $in;
  print stderr "$filesize\n";
  open IN, "<", $in;
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
  
$count_x = 0;
$count_y = 0;
for($i = 0; $i <= $#x_east; $i++){
  $count_x++;
  if($count_x == $ngrid_x){
    $count_y++;
    $count_x = 0;
  }
  $slowness_sq = ($slowness_x[$i] ** 2 + $slowness_y[$i] ** 2);
  $velocity = sqrt(1.0 / $slowness_sq) * 1000.0;
  $velocity = sprintf "%.4f", $velocity;
  $velocity_err = 1000.0 / $slowness_sq * 
                  sqrt(($sigma_slowness_x[$i] ** 2 * $slowness_x[$i] ** 2 + $sigma_slowness_y[$i] ** 2 * $slowness_y[$i] ** 2)
                        / $slowness_sq);
  $velocity_err = $velocity_err / $velocity * 100.0;
  $velocity_err = sprintf "%.2f", $velocity_err;
  $az = 90.0 - rad2deg(atan2($slowness_y[$i], $slowness_x[$i]));
  $az = sprintf "%.2f", $az;
  $az_err = 1.0 / $slowness_sq * 
                  sqrt(($sigma_slowness_x[$i] * $slowness_y[$i]) ** 2 + ($sigma_slowness_y[$i] * $slowness_x[$i]) ** 2);
  $az_err = rad2deg($az_err);
  $az_err = sprintf "%.2f", $az_err;
  if($count_x % $decimate == 0 && $count_y % $decimate == 0){
    print stdout "$x_east[$i] $y_north[$i] $velocity $velocity_err $az $az_err\n";
  }
}


