
reset
set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,22'

load 'best-fit_bias_0_6.dat'
load '../l/params_best-fit_bd_0_6.dat'

set xr[0:7.5]
set yr[-0.9:0.6]
set mxtics 2
set mytics 4
set key b r
set pointsize 0.5
set style fill transparent solid 0.2 noborder
set xlabel 'r/r_v'
set ylabel '{/Symbol x}_{vm}'

set label '{/Symbol L}CDM' at 6,0.4

#  plot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:4 lt 1 w p t ''
#replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($4+$5):($4-$5) lt 1 with filledcurves t 'r_{0.2}=1.7 Mpc / h'
#re
  plot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:6 lt 1 w p t ''
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($6+$7):($6-$7) lt 1 with filledcurves t 'r_{0.2}=2.3 Mpc / h'
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:8 lt 2 w p t ''
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($8+$9):($8-$9) lt 2 with filledcurves t 'r_{0.2}=3.0 Mpc / h'
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:10 lt 3 w p t ''
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($10+$11):($10-$11) lt 3 with filledcurves t 'r_{0.2}=4.0 Mpc / h'
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:12 lt 4 w p t ''
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($12+$13):($12-$13) lt 4 with filledcurves t 'r_{0.2}=5.2 Mpc / h'
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:14 lt 6 w p t ''
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($14+$15):($14-$15) lt 6 with filledcurves t 'r_{0.2}=6.9 Mpc / h'
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:16 lt 7 w p t ''
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($16+$17):($16-$17) lt 7 with filledcurves t 'r_{0.2}=9.0 Mpc / h'
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:18 lt 8 w p t ''
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($18+$19):($18-$19) lt 8 with filledcurves t 'r_{0.2}=12. Mpc / h'


f(x)=-0.867554*(1.-(x/Pbd2)**Pbd0)/(1.+(x/Pbd3)**Pbd1)
#tail(x,r)=exp(-Pbd4*r**(-Pbd6)*x**-(Pbd5))
tail(x,r)=exp(-Pbd4*r**(-1)*x**-(Pbd5))



#replot 'lin_matter_correlation_fuction.txt' u ($1/r1):(f($1/r1)+$2*b1*tail($1/r1,r1)) w l lc 1 lt 3 t ''
replot 'lin_matter_correlation_fuction.txt' u ($1/r2):(f($1/r2)+$2*b2*tail($1/r2,r2)) w l lc 1 lt 3 t '' 
replot 'lin_matter_correlation_fuction.txt' u ($1/r3):(f($1/r3)+$2*b3*tail($1/r3,r3)) w l lc 2 lt 3 t '' 
replot 'lin_matter_correlation_fuction.txt' u ($1/r4):(f($1/r4)+$2*b4*tail($1/r4,r4)) w l lc 3 lt 3 t '' 
replot 'lin_matter_correlation_fuction.txt' u ($1/r5):(f($1/r5)+$2*b5*tail($1/r5,r5)) w l lc 4 lt 3 t '' 
replot 'lin_matter_correlation_fuction.txt' u ($1/r6):(f($1/r6)+$2*b6*tail($1/r6,r6)) w l lc 6 lt 3 t '' 
replot 'lin_matter_correlation_fuction.txt' u ($1/r7):(f($1/r7)+$2*b7*tail($1/r7,r7)) w l lc 7 lt 3 t '' 
replot 'lin_matter_correlation_fuction.txt' u ($1/r8):(f($1/r8)+$2*b8*tail($1/r8,r8)) w l lc 8 lt 3 t '' 

replot 0 lc -1 t ''
set output 'densidad_vs_bias_ll_.png'
replot
reset




#########################################################

#########################################################


reset
set terminal pngcairo size 1000,750 enhanced font 'TimesNewRoman,18'

load 'best-fit_bias_0_6.dat'
load '../l/params_best-fit_bd_0_6.dat'

set xr[0:7.5]
set yr[-0.3:0.6]
set mxtics 2
set mytics 4
set key t r
set pointsize 0.5
set style fill transparent solid 0.2 noborder
set xlabel 'r/r_v'
set ylabel '{/Symbol x}_{vm}   -   {/Symbol r}_{v} / {~{/Symbol r}{.7-}}_{m}'


ff(x)=-0.867554*(1.-(x/Pbd2)**Pbd0)/(1.+(x/Pbd3)**Pbd1)

#  plot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($4-ff($1)) lt 1 w p t ''
#replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($4-ff($1)+($5**2+$57**2)**0.5):($4-ff($1)-($5**2+$57**2)**0.5) lt 1 with filledcurves t 'r_{0.2}=1.7 Mpc / h'
#re
  plot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($6-ff($1)) lt 1 w p t ''
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($6-ff($1)+($7**2+$57**2)**0.5):($6-ff($1)-($7**2+$57**2)**0.5) lt 1 with filledcurves t 'r_{0.2}=2.3 Mpc / h'
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($8-ff($1)) lt 2 w p t ''
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($8-ff($1)+($9**2+$57**2)**0.5):($8-ff($1)-($9**2+$57**2)**0.5) lt 2 with filledcurves t 'r_{0.2}=3.0 Mpc / h'
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($10-ff($1)) lt 3 w p t ''
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($10-ff($1)+($11**2+$57**2)**0.5):($10-ff($1)-($11**2+$57**2)**0.5) lt 3 with filledcurves t 'r_{0.2}=4.0 Mpc / h'
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($12-ff($1)) lt 4 lw 2 w p t ''
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($12-ff($1)+($13**2+$57**2)**0.5):($12-ff($1)-($13**2+$57**2)**0.5) lt 4 with filledcurves t 'r_{0.2}=5.2 Mpc / h'
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($14-ff($1)) lt 5 w p t ''
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($14-ff($1)+($15**2+$57**2)**0.5):($14-ff($1)-($15**2+$57**2)**0.5) lt 5 with filledcurves t 'r_{0.2}=6.9 Mpc / h'
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($16-ff($1)) lt 7 w p t ''
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($16-ff($1)+($17**2+$57**2)**0.5):($16-ff($1)-($17**2+$57**2)**0.5) lt 7 with filledcurves t 'r_{0.2}=9.0 Mpc / h'
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($18-ff($1)) lt 8 w p t ''
replot 'densidad_n_esferas_512x256_03_plus_5.txt' u 1:($18-ff($1)+($19**2+$57**2)**0.5):($18-ff($1)-($19**2+$57**2)**0.5) lt 8 with filledcurves t 'r_{0.2}=12. Mpc / h'


#f(x)=-0.867554*(1.-(x/Pbd2)**Pbd0)/(1.+(x/Pbd3)**Pbd1)
f(x)=0

tail(x,r)=exp(-(Pbd4*(r/6.9)**Pbd5*x)**-(Pbd6*(r/6.9)**Pbd7))


#replot 'lin_matter_correlation_fuction.txt' u ($1/r1):(f($1/r1)+$2*b1*tail($1/r1,r1)) w l lc 1 lt 3 t ''
replot 'lin_matter_correlation_fuction.txt' u ($1/r2):(f($1/r2)+$2*b2*tail($1/r2,r2)) w l lc 1 lt 3 t '' 
replot 'lin_matter_correlation_fuction.txt' u ($1/r3):(f($1/r3)+$2*b3*tail($1/r3,r3)) w l lc 2 lt 3 t '' 
replot 'lin_matter_correlation_fuction.txt' u ($1/r4):(f($1/r4)+$2*b4*tail($1/r4,r4)) w l lc 3 lt 3 t '' 
replot 'lin_matter_correlation_fuction.txt' u ($1/r5):(f($1/r5)+$2*b5*tail($1/r5,r5)) w l lc 4 lt 3 t '' 
replot 'lin_matter_correlation_fuction.txt' u ($1/r6):(f($1/r6)+$2*b6*tail($1/r6,r6)) w l lc 5 lw 3 t '' 
replot 'lin_matter_correlation_fuction.txt' u ($1/r7):(f($1/r7)+$2*b7*tail($1/r7,r7)) w l lc 7 lt 3 t '' 
replot 'lin_matter_correlation_fuction.txt' u ($1/r8):(f($1/r8)+$2*b8*tail($1/r8,r8)) w l lc 8 lt 3 t '' 

replot 0 lc -1 t ''
set output 'densidad_vs_bias_l.png'
replot
reset






reset
set terminal pngcairo size 1000,750 enhanced font 'TimesNewRoman,18'

load '../l/params_best-fit_bd_0_6.dat'
set xzeroaxis
set samples 400

set xr[0:7.5]
set yr[-0.9:0.6]
set mxtics 2
set mytics 4
set key t r
set pointsize 0.5
set style fill transparent solid 0.2 noborder
set xlabel 'r/r_v'
set ylabel '{/Symbol r}_{v} / {~{/Symbol r}{.7-}}_{m}'


f(x)=-0.867554*(1.-(x/Pbd2)**Pbd0)/(1.+(x/Pbd3)**Pbd1)


plot f(x) lc -1 t 'Universal Void Profile'


set output 'universal_profile.png'
replot
reset


