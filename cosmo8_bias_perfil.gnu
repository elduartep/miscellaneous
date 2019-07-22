
reset
set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,20'

load 'params_best-fit_bd_0_6.dat'

set xr[0:10]
set yr[-0.9:0.6]
set mxtics 2
set mytics 4
set key b r
set pointsize 0.5
set style fill transparent solid 0.2 noborder
set xlabel 'r/r_v'
set ylabel '{/Symbol x}_{vm}'

set label '{/Symbol L}CDM' at 8,0.4

#  plot 'densidad_n_esferas_512x256_03.txt' u 1:2 lt 1 w p t ''
#replot 'densidad_n_esferas_512x256_03.txt' u 1:($2+$3):($2-$3) lt 1 with filledcurves t 'r_{0.2}=1.7 Mpc / h'
#re
  plot 'densidad_n_esferas_512x256_03.txt' u 1:4 lt 2 w p t ''
replot 'densidad_n_esferas_512x256_03.txt' u 1:($4+$5):($4-$5) lt 2 with filledcurves t 'r_{0.2}=2.3 Mpc / h'
replot 'densidad_n_esferas_512x256_03.txt' u 1:6 lt 3 w p t ''
replot 'densidad_n_esferas_512x256_03.txt' u 1:($6+$7):($6-$7) lt 3 with filledcurves t 'r_{0.2}=3.0 Mpc / h'
replot 'densidad_n_esferas_512x256_03.txt' u 1:8 lt 4 w p t ''
replot 'densidad_n_esferas_512x256_03.txt' u 1:($8+$9):($8-$9) lt 4 with filledcurves t 'r_{0.2}=4.0 Mpc / h'
replot 'densidad_n_esferas_512x256_03.txt' u 1:10 lt 5 w p t ''
replot 'densidad_n_esferas_512x256_03.txt' u 1:($10+$11):($10-$11) lt 5 with filledcurves t 'r_{0.2}=5.2 Mpc / h'
replot 'densidad_n_esferas_512x256_03.txt' u 1:12 lt 6 w p t ''
replot 'densidad_n_esferas_512x256_03.txt' u 1:($12+$13):($12-$13) lt 6 with filledcurves t 'r_{0.2}=6.9 Mpc / h'
replot 'densidad_n_esferas_512x256_03.txt' u 1:14 lt 7 w p t ''
replot 'densidad_n_esferas_512x256_03.txt' u 1:($14+$15):($14-$15) lt 7 with filledcurves t 'r_{0.2}=9.0 Mpc / h'
replot 'densidad_n_esferas_512x256_03.txt' u 1:16 lt 8 w p t ''
replot 'densidad_n_esferas_512x256_03.txt' u 1:($16+$17):($16-$17) lt 8 with filledcurves t 'r_{0.2}=12. Mpc / h'

#dc=0.867554
#f(x)=-dc*(1.-(x/u1)**u2)/(1.0+(x/u3)**u4)
#g(x,A,B)=exp(-A*x**-B)
#a=1.3
#b=0.98
##replot 'lin_matter_correlation_fuction.txt' u ($1/r0):(f($1/r0*b)+$2*b0*g($1/r0,A0,B0)) w l lc 1 lt 3 t ''
#replot 'lin_matter_correlation_fuction.txt' u ($1/r1):(f($1/r1*b)+$2*b1*g($1/r1,A1*a,B1)) w l lc 2 lt 3 t ''
#replot 'lin_matter_correlation_fuction.txt' u ($1/r2):(f($1/r2*b)+$2*b2*g($1/r2,A2*a,B2)) w l lc 3 lt 3 t ''
#replot 'lin_matter_correlation_fuction.txt' u ($1/r3):(f($1/r3*b)+$2*b3*g($1/r3,A3*a,B3)) w l lc 4 lt 3 t ''
#replot 'lin_matter_correlation_fuction.txt' u ($1/r4):(f($1/r4*b)+$2*b4*g($1/r4,A4*a,B4)) w l lc 5 lt 3 t ''
#replot 'lin_matter_correlation_fuction.txt' u ($1/r5):(f($1/r5*b)+$2*b5*g($1/r5,A5*a,B5)) w l lc 6 lt 3 t ''
#replot 'lin_matter_correlation_fuction.txt' u ($1/r6):(f($1/r6*b)+$2*b6*g($1/r6,A6*a,B6)) w l lc 7 lt 3 t ''
#replot 'lin_matter_correlation_fuction.txt' u ($1/r7):(f($1/r7*b)+$2*b7*g($1/r7,A7*a,B7)) w l lc 8 lt 3 t ''

dc=0.867554
f(x,a,b)=exp(-(a*x)**-b)
g(x)=(x/u1)**u2/(1.+(x/u3)**u4)

#replot 'lin_matter_correlation_fuction.txt' u ($1/r0):(f($1/r0,A0,B0)*(dc+$2*b0)-dc+g($1/r0)) w l lc 1 lt 3 t ''
replot 'lin_matter_correlation_fuction.txt' u ($1/r1):(f($1/r1,A1,B1)*(dc+$2*b1)-dc+g($1/r1)) w l lc 2 lt 3 t ''
replot 'lin_matter_correlation_fuction.txt' u ($1/r2):(f($1/r2,A2,B2)*(dc+$2*b2)-dc+g($1/r2)) w l lc 3 lt 3 t ''
replot 'lin_matter_correlation_fuction.txt' u ($1/r3):(f($1/r3,A3,B3)*(dc+$2*b3)-dc+g($1/r3)) w l lc 4 lt 3 t ''
replot 'lin_matter_correlation_fuction.txt' u ($1/r4):(f($1/r4,A4,B4)*(dc+$2*b4)-dc+g($1/r4)) w l lc 5 lt 3 t ''
replot 'lin_matter_correlation_fuction.txt' u ($1/r5):(f($1/r5,A5,B5)*(dc+$2*b5)-dc+g($1/r5)) w l lc 6 lt 3 t ''
replot 'lin_matter_correlation_fuction.txt' u ($1/r6):(f($1/r6,A6,B6)*(dc+$2*b6)-dc+g($1/r6)) w l lc 7 lt 3 t ''
replot 'lin_matter_correlation_fuction.txt' u ($1/r7):(f($1/r7,A7,B7)*(dc+$2*b7)-dc+g($1/r7)) w l lc 8 lt 3 t ''

replot 0 lc -1 t ''
set output 'densidad_vs_bias_lupa.png'
replot
reset

