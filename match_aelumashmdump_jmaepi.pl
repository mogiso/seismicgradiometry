#!/usr/bin/env perl
use Math::Trig qw(tan atan acos pi great_circle_distance great_circle_direction deg2rad);

$arrayresultfile = $ARGV[0];
$jmahypofile = $ARGV[1];
$matchresultfile = $ARGV[2];
$nonmatchresultfile = $ARGV[3];

$sec_diff_threshold = 150.0;

open IN, "<", $arrayresultfile;
while(<IN>){
  chomp $_;
  $_ =~ s/^\s+(.*?)\s+$/$1/;
  @tmp = split /\s+/, $_;
  push @arrayresult, $_;
  push @arraydatetime, $tmp[0];
  push @arrayhypolon, $tmp[3];
  push @arrayhypolat, $tmp[6];
}
close IN;

open IN, "<", $jmahypofile;
while(<IN>){
  chomp $_;
  push @jmahypodatetime, (substr $_, 1, 16);
  $jmalat = substr $_, 21, 3;
  $jmalat_m = substr $_, 24, 4;
  push @jmahypolat, ($jmalat + $jmalat_m / 100.0 / 60.0);
  $jmalon = substr $_, 32, 4;
  $jmalon_m = substr $_, 36, 4;
  push @jmahypolon, ($jmalon + $jmalon_m / 100.0 / 60.0);
  $jmamag = substr $_, 52, 2;
  push @jmahypomag, $jmamag * 0.1;
} 
close IN;


open OUT, ">", $matchresultfile;
open OUT2, ">", $nonmatchresultfile;
@match_flag = ();
for($i = 0; $i <= $#arraydatetime; $i++){
  ($arrayyr, $arraymo, $arraydy, $arrayhr, $arraymi, $arraysc) = split /[-T:]/, $arraydatetime[$i];
  &ymd2jday($arrayyr, $arraymo, $arraydy, $arrayjday);
  $arraysec_from_day = $arrayhr * 3600.0 + $arraymi * 60.0 + $arraysc;


  $match_flag[$i] = -1;
  @distance = ();
  for($j = 0; $j <= $#jmahypodatetime; $j++){
    @point1 = ();
    $jmayr = substr $jmahypodatetime[$j], 0, 4;
    next if ($jmayr != $arrayyr);
    $jmamo = substr $jmahypodatetime[$j], 4, 2;
    $jmady = substr $jmahypodatetime[$j], 6, 2;
    $jmahr = substr $jmahypodatetime[$j], 8, 2;
    $jmami = substr $jmahypodatetime[$j], 10, 2;
    $jmasc = substr $jmahypodatetime[$j], 12, 4;
    &ymd2jday($jmayr, $jmamo, $jmady, $jmajday);
    $jday_diff = $arrayjday - $jmajday;
    next if (abs($jday_diff) > 1);

    $jmasec_from_day = $jmahr * 3600.0 + $jmami * 60.0 + $jmasc / 100.0;
    $sec_diff = $arraysec_from_day - $jmasec_from_day;
    $sec_diff = $sec_diff + $jday_diff * 86400.0;
    next if (abs($sec_diff) > $sec_diff_threshold);

    $lat0 = geograph2geocen($arrayhypolat[$i]);
    @point0 = NESW($arrayhypolon[$i], $lat0);
    $lat1 = geograph2geocen($jmahypolat[$j]);
    @point1 = NESW($jmahypolon[$j], $lat1);
    $distance[$j] = sprintf "%.3f", (great_circle_distance(@point0, @point1, 6370.291));
    next if ($distance[$j] > $distance_threshold);

    $match_flag[$i] = $j;
    print stderr "$sec_diff $arraydatetime[$i] $jmahypodatetime[$j]\n";
  }

  if($match_flag[$i] != -1){
    print stdout "$arrayresult[$i] $jmahypodatetime[$match_flag[$i]] $jmahypolon[$match_flag[$i]] $jmahypolat[$match_flag[$i]] $jmahypomag[$match_flag[$i]] $distance[$match_flag[$i]]\n";
  }else{
    print OUT2 "$arrayresult[$i]\n";
  }
}
close OUT;
close OUT2;






##距離計算のためのサブルーチン
sub NESW{
  deg2rad($_[0]), deg2rad(90 - $_[1])
}
##地理緯度を地心緯度に変換
sub geograph2geocen{

  ##楕円体の定数(GRS80)
  my $r1 = 6378137.0;
  my $f = 1.0 / 298.257222101;
  my $r2 = $r1 * (1.0 - $f);

  return atan(($r2 / $r1) ** 2 * tan($_[0] * pi / 180.0)) * 180.0 / pi;
}
sub leapyear{
  if($_[0] % 400 == 0){
    $_[1] = 1;
  }elsif($_[0] % 100 == 0){
    $_[1] = 0;
  }elsif($_[0] % 4 == 0){
    $_[1] = 1;
  }else{
    $_[1] = 0;
  }
}
sub jday2ymd{
  my $leap;
  my $jday = $_[0];
  my $year = $_[1];
  my $mm   = $_[2];
  my $dd   = $_[3];

  &leapyear($year, $leap);
  if($leap){
    $dom[1] = 29;
  }else{
    $dom[1] = 28;
  }

  @dom = (31, $dom[1], 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);
  my $jday_tmp = $dom[0];
  for(my $i = 1; $i <= 12; $i++){
    if($jday <= $jday_tmp){
      $mm = $i - 1;
      last;
    }
    $jday_tmp = $jday_tmp + $dom[$i];
  }
  $dd = $jday - ($jday_tmp - $dom[$mm]);

  $_[2] = $mm + 1;
  $_[3] = $dd;
} 
  
sub ymd2jday{
  my $leap;
  my $year = $_[0];
  my $mm   = $_[1];
  my $dd   = $_[2];
  my $jday;

  &leapyear($year, $leap);
  if($leap){
    $dom[1] = 29;
  }else{
    $dom[1] = 28;
  }

  @dom = (31, $dom[1], 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);
  $jday = $dd;
  if($mm - 1 > 0){
    for(my $i = 0; $i < $mm - 1; $i++){
      $jday = $jday + $dom[$i];
    }
  }
  $_[3] = $jday;
} 
