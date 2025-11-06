#!/usr/bin/env perl

$out = $ARGV[0];

##plot likelihood function

##(a): g_0
$c0 = 0.7;
$c1 = 20.0;

$size_x = 6;
$size_y = 6;
$txt_x = -1.0;
$txt_y = $size_y + 0.5;
$dx = $size_x + 2.5;
$dy = $size_y + 3.0;
$x_min = 0;
$x_max = 60;
$y_min = 0.0;
$y_max = 1.0;
$annot_x = "a20f10";
$label_x = "N\@%12%\@-\161\@-\@%%\@+obs\@+";
$annot_y = "a0.2";
$label_y = "g\@-0\@-(N\@%12%\@-\161\@-\@%%\@+obs\@+)";
$dxval = 0.1;

@xvalarray = (1.0, 10.0, 30.0, 50.0);
@penarray = ("1p,black", "1p,dimgray", "1p,darkgray", "1p,lightgray");

system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -K -P -X3c -Y14c > $out";
open OUT, " | gmt psxy -JX$size_x/$size_y -R$x_min/$x_max/$y_min/$y_max -W1.2p,black -O -K -P >> $out";
$xval = $x_min;
while($xval <= $x_max){
  $yval = 1.0 - $c0 * exp(-$xval * $xval * 0.5 / $c1 / $c1);
  print OUT "$xval $yval\n";
  #print stdout "$xval $yval\n";
  $xval = $xval + $dxval;
}
close OUT;

#c0=0.5
#$c0 = 0.5;
#open OUT, " | gmt psxy -JX$size_x/$size_y -R$x_min/$x_max/$y_min/$y_max -W0.6p,black,- -O -K -P >> $out";
#$xval = $x_min;
#while($xval <= $x_max){
#  $yval = 1.0 - $c0 * exp(-$xval * $xval * 0.5 / $c1 / $c1);
#  print OUT "$xval $yval\n";
#  #print stdout "$xval $yval\n";
#  $xval = $xval + $dxval;
#}
#close OUT;

system "gmt psbasemap -JX$size_x/$size_y -R$x_min/$x_max/$y_min/$y_max \\
                      -Bpx${annot_x} -Bpy$annot_y -BWSen \\
                      -O -K -P >> $out";

open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -F+f+a+j -N -O -K -P >> $out";
  print OUT "$txt_x $txt_y 16p,Helvetica,black 0 LB (a)\n";
  print OUT "3.0 -1.0     16p,Helvetica,black 0 CT N\@-\@%12%\161\@%%\@-\@+obs\@+\n";
  print OUT "-1.8 3.0     16p,Helvetica,black 90 CT g\@-0\@-(N\@-\@%12%\161\@%%\@-\@+obs\@+)\n";
close OUT;

system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -K -P -X$dx >> $out";


$x_min = -3.14;
$x_max = 3.14;
$annot_x = "a1.57f0.785";
$label_x = "\@%12%\161\@%%\@-obs\@-\-@%12%\161\@%%\@-cal\@-";
$annot_y = "a0.2";
$label_y = "L\@%12%\@-\161\@-\@%%";
$c2 = 15.0 * 3.14159265 / 180.0;

$dxval = 0.005;
for($i = 0; $i <= $#xvalarray; $i++){
  $pen = $penarray[$i];
  print stdout "$i $pen\n";
  open OUT, " | gmt psxy -JX$size_x/$size_y -R$x_min/$x_max/$y_min/$y_max -W${pen} -O -K -P >> $out";
  $xval = $x_min;
  while($xval <= $x_max){
    $yval_tmp = 1.0 - $c0 * exp(-$xvalarray[$i] * $xvalarray[$i] * 0.5 / $c1 / $c1);
    $yval = (1.0 - $yval_tmp) * exp(-0.5 * $xval * $xval / $c2 / $c2) + $yval_tmp;
    print OUT "$xval $yval\n";
    #print stdout "$xval $yval\n";
    $xval = $xval + $dxval;
    }
  close OUT;
}

system "gmt psbasemap -JX$size_x/$size_y -R$x_min/$x_max/$y_min/$y_max \\
                      -Bpx${annot_x}+l\"$label_x\" -Bpy$annot_y+l\"$label_y\" -BWsen \\
                      -O -K -P >> $out";

open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -F+f+a+j -N -O -K -P >> $out";
  print OUT "0   -0.35 16p,Helvetica,black 0 CT \055\@%12%\160\@%%\n";
  print OUT "1.5 -0.35 16p,Helvetica,black 0 CT \055\@%12%\160\@%%\@:12p:\/2\@::\n";
  print OUT "3.0 -0.35 12p,Helvetica,black 0 CT 0\n";
  print OUT "4.5 -0.35 16p,Helvetica,black 0 CT \@%12%\160\@%%\@:12p:\/2\@::\n";
  print OUT "6.0 -0.35 16p,Helvetica,black 0 CT \@%12%\160\@%%\n";
  print OUT "3.0 -1.0  16p,Helvetica,black 0 CT \@%12%\161\@%%\@-obs\@-\055\@%12%\161\@%%\@-cal\@- (radian)\n";
  print OUT "$txt_x $txt_y 16p,Helvetica,black 0 LB (b)\n";
close OUT;
open OUT, " | gmt pstext -JX$size_x/$size_y -R$x_min/$x_max/$y_min/$y_max -F+f+a+j -N -O -K -P >> $out";
  print OUT "-1.57  0.29  12p,Helvetica,black 0 CT N\@-\@%12%\161\@%%\@-\@+obs\@+=$xvalarray[0]\n";
  print OUT "-1.57  0.41  12p,Helvetica,dimgray 0 CB N\@-\@%12%\161\@%%\@-\@+obs\@+=$xvalarray[1]\n";
  print OUT "-1.57  0.76  12p,Helvetica,darkgray 0 CT N\@-\@%12%\161\@%%\@-\@+obs\@+=$xvalarray[2]\n";
  print OUT "-1.57  0.95  12p,Helvetica,lightgray 0 CT N\@-\@%12%\161\@%%\@-\@+obs\@+=$xvalarray[3]\n";
close OUT;


system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -K -P -X-$dx -Y-$dy >> $out";

##(c)
$x_min = 1;
$x_max = 1000;
$y_min = 0.0;
$y_max = 1.0;
$annot_x = "a1p";
$label_x = "\@%12%\104\@%%\@-obs\@- (km)";
$annot_y = "a0.2";
$label_y = "g\@-0\@-(\@%12%\104\@%%\@-obs\@-)";
$dxval = 1.0;
@xvalarray = (10.0, 50.0, 200.0, 600.0);

$c0 = 0.7;
$c1 = log(100.0);
open OUT, " | gmt psxy -JX${size_x}l/$size_y -R$x_min/$x_max/$y_min/$y_max -W1.2p,black -O -K -P >> $out";
$xval = $x_min;
while($xval <= $x_max){
  $xtmp = log($xval);
  $yval = 1.0 - $c0 * exp(-$xtmp * $xtmp * 0.5 / $c1 / $c1);
  print OUT "$xval $yval\n";
  #print stdout "$xval $yval\n";
  $xval = $xval + $dxval;
}
close OUT;

#c0=0.5
#$c0 = 0.5;
#open OUT, " | gmt psxy -JX$size_x/$size_y -R$x_min/$x_max/$y_min/$y_max -W0.6p,black,- -O -K -P >> $out";
#$xval = $x_min;
#while($xval <= $x_max){
#  $yval = 1.0 - $c0 * exp(-$xval * $xval * 0.5 / $c1 / $c1);
#  print OUT "$xval $yval\n";
#  #print stdout "$xval $yval\n";
#  $xval = $xval + $dxval;
#}
#close OUT;

system "gmt psbasemap -JX${size_x}l/$size_y -R$x_min/$x_max/$y_min/$y_max \\
                      -Bpx${annot_x}+l\"$label_x\" -Bpy$annot_y+l\"$label_y\" -BWSen \\
                      -O -K -P >> $out";

open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -F+f+a+j -N -O -K -P >> $out";
  print OUT "$txt_x $txt_y 16p,Helvetica,black 0 LB (c)\n";
close OUT;

system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -K -P -X$dx >> $out";


$x_min = -300;
$x_max = 300;
$annot_x = "a100";
$annot_y = "a0.2";
$label_y = "L\@-T0\@-";
$c2 = 60.0;

for($i = 0; $i <= $#xvalarray; $i++){
  $pen = $penarray[$i];
  print stdout "$i $pen\n";
  open OUT, " | gmt psxy -JX$size_x/$size_y -R$x_min/$x_max/$y_min/$y_max -W${pen} -O -K -P >> $out";
  $xval = $x_min;
  while($xval <= $x_max){
    $yval_tmp = 1.0 - $c0 * exp(-log($xvalarray[$i]) * log($xvalarray[$i]) * 0.5 / $c1 / $c1);
    $yval = (1.0 - $yval_tmp) * exp(-0.5 * $xval * $xval / $c2 / $c2) + $yval_tmp;
    print OUT "$xval $yval\n";
    #print stdout "$xval $yval\n";
    $xval = $xval + $dxval;
    }
  close OUT;
}

system "gmt psbasemap -JX$size_x/$size_y -R$x_min/$x_max/$y_min/$y_max \\
                      -Bpx${annot_x} -Bpy$annot_y+l\"$label_y\" -BWSen \\
                      -O -K -P >> $out";

open OUT, " | gmt pstext -JX$size_x/$size_y -R0/$size_x/0/$size_y -F+f+a+j -N -O -K -P >> $out";
  print OUT "3.0 -0.8  16p,Helvetica,black 0 CT T\@-0\@-\@+obs\@+\055T\@-0\@-\@+cal\@+ (s)\n";
  print OUT "$txt_x $txt_y 16p,Helvetica,black 0 LB (d)\n";
close OUT;
open OUT, " | gmt pstext -JX$size_x/$size_y -R$x_min/$x_max/$y_min/$y_max -F+f+a+j -N -O -K -P >> $out";
  print OUT "-299  0.38  12p,Helvetica,black 0 LT \@%12%\104\@%%\@-obs\@-=$xvalarray[0]\n";
  print OUT "-299  0.5  12p,Helvetica,dimgray 0 LT \@%12%\104\@%%\@-obs\@-=$xvalarray[1]\n";
  print OUT "-299  0.63  12p,Helvetica,darkgray 0 LT \@%12%\104\@%%\@-obs\@-=$xvalarray[2]\n";
  print OUT "-299  0.72  12p,Helvetica,lightgray 0 LT \@%12%\104\@%%\@-obs\@-=$xvalarray[3]\n";
close OUT;


system "cat /dev/null | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O -P >> $out";
system "gmt psconvert -A -Tg $out";


unlink "gmt.conf";
unlink "gmt.history";

