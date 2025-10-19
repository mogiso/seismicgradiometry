#!/usr/bin/env perl
use Math::Trig 'deg2rad sin cos'; 

##calculate cosine similarity

$in = $ARGV[0];
$in_dir = $ARGV[1];

open IN, "<", $in;
while(<IN>){
  chomp $_;
  push @arrayepicenter, $_;
}
close IN;

for($i = 0; $i <= $#arrayepicenter; $i++){
  $yr = substr $arrayepicenter[$i], 0, 4;
  $mo = substr $arrayepicenter[$i], 5, 2;
  $dy = substr $arrayepicenter[$i], 8, 2;
  $hr = substr $arrayepicenter[$i], 11, 2;
  $mi = substr $arrayepicenter[$i], 14, 2;
  $sc = substr $arrayepicenter[$i], 17, 2;
  $arrrayobsfile = "$in_dir/${yr}${mo}${dy}${hr}${mi}${sc}_array_obs.dat";
  @tmp = split /\s+/, $arrayepicenter[$i];
  $arraynum = $tmp[14];

  die "No $arrayobsfile\n" if (!-f $arrayobsfile);
  $filesize = -s $arrayobsfile;
  die "Filesize match error\n" if $filesize != 4 * 7 * $arraynum;

  $read_fiiesize = 0;
  open IN, "<", $arrayobsfile;
  while($read_filesize < $filesize){
    read IN, $buf, 4;
    push @arraylon, (unpack "f", $buf);  
    read IN, $buf, 4;
    push @arraylat, (unpack "f", $buf);  
    read IN, $buf, 4;
    push @arrayaz, (unpack "f", $buf) * deg2rad;  
    read IN, $buf, 4;
    push @arrayapv, (unpack "f", $buf);  
    read IN, $buf, 4;
    read IN, $buf, 4;
    read IN, $buf, 4;
    $read_filesize = $read_filesize + 4 * 7;
  }
  close IN;

  $min_cos_similarity = 1.0;
  for($k = 0; $k <= $arraynum - 1; $k++){
    for($j = $k + 1; $j <= $arraynum; $k++){
      $cos_similarity_tmp = cos($arrayaz[$k]) * cos($arrayaz[$j]) + sin($arrayaz[$k]) * sin($arrayaz[$j]);
      $min_cos_similarity = $cos_similarity if($cos_similarity_tmp <= $min_cos_similarity);
    }
  }

  print stdout "$arrayepicenter[$i] $min_cos_similarity\n";
}


    

  
  
