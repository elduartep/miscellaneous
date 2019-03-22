reset
set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'log_{10}|f_{R0}|'
set label 'z_{SSB}=3' at -4.7,0.9
unset key
set xtics 0.2
set mxtics 2
plot 'lista2_0_1_0_1_0.txt' u 2:(exp(-0.5*($1-1.892142e+02))) w l dt 2 lc 1 lw 2 t 'abundance'
replot 'lista2_0_1_1_0_0.txt' u 2:(exp(-0.5*($1-1.471204e+03))) w l dt 3 lc 2 lw 3 t 'density profile'
replot 'lista2_0_1_0_0_1.txt' u 2:(exp(-0.5*($1-1.124549e+01))) w l dt 5 lc 4 lw 2 t 'bias'
replot 'lista2_0_1_1_1_1.txt' u 2:(exp(-0.5*($1-1.726016e+03))) w l dt 1 lc 7 lw 2 t 'all'
set xr[-4.800000e+00:-2.800000e+00]
set output 'histo_fD.png'
replot
reset


reset
set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'log_{10}|f_{R0}|'
set label 'z_{SSB}=2' at -5.9,0.9
unset key
set xtics 0.2
set mxtics 2
plot 'lista2_1_1_0_1_0.txt' u 2:(exp(-0.5*($1-8.312710e+01))) w l dt 2 lc 1 lw 2 t 'abundance'
replot 'lista2_1_1_1_0_0.txt' u 2:(exp(-0.5*($1-9.553256e+02))) w l dt 3 lc 2 lw 3 t 'density profile'
replot 'lista2_1_1_0_0_1.txt' u 2:(exp(-0.5*($1-4.306938e+00))) w l dt 5 lc 4 lw 2 t 'bias'
replot 'lista2_1_1_1_1_1.txt' u 2:(exp(-0.5*($1-1.096380e+03))) w l dt 1 lc 7 lw 2 t 'all'
set xr[-6.000000e+00:-4.000000e+00]
set output 'histo_fB.png'
replot
reset


reset
set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel '|f_{R0}| {/Symbol \264} 10^6'
set label 'z_{SSB}=1' at 2,0.9
set key t r
set xtics 1
set mxtics 5
plot 'lista2_2_1_0_1_0.txt' u (10.**$2*10**6):(exp(-0.5*($1-1.589104e+01))) w l dt 2 lc 1 lw 2 t 'abundance'
replot 'lista2_2_1_1_0_0.txt' u (10.**$2*10**6):(exp(-0.5*($1-7.274093e+02))) w l dt 3 lc 2 lw 3 t 'density profile'
replot 'lista2_2_1_0_0_1.txt' u (10.**$2*10**6):(exp(-0.5*($1-2.824014e+00))) w l dt 5 lc 4 lw 2 t 'bias'
replot 'lista2_2_1_1_1_1.txt' u (10.**$2*10**6):(exp(-0.5*($1-7.556652e+02))) w l dt 1 lc 7 lw 2 t 'all'
set xr[:5]
set output 'histo_fA.png'
replot
reset


reset
set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'log_{10}|f_{R0}|'
set label 'z_{SSB}=0' at 0.0,0.9
unset key
set xtics 0.5
set mxtics 5
set xlabel '|f_{R0}| {/Symbol ´} 10^8'
set label 'z_{SSB}=0' at 0.0,0.9
set format x '%1.0t'
set xtics 1e-8
 set mxtics 2
plot 'lista2_3_1_0_1_0.txt' u (10.**$2):(exp(-0.5*($1-1.069846e+01))) w l dt 2 lc 1 lw 2 t 'abundance'
replot 'lista2_3_1_1_0_0.txt' u (10.**$2):(exp(-0.5*($1-5.837314e+02))) w l dt 3 lc 2 lw 3 t 'density profile'
replot 'lista2_3_1_0_0_1.txt' u (10.**$2):(exp(-0.5*($1-4.558067e+00))) w l dt 5 lc 4 lw 2 t 'bias'
replot 'lista2_3_1_1_1_1.txt' u (10.**$2):(exp(-0.5*($1-5.991741e+02))) w l dt 1 lc 7 lw 2 t 'all'
set xr[0.000000e+00:7.000000e-08]
set output 'histo_fl.png'
replot
reset


reset
set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'z_{SSB}'
set label '|f_{R0}|=10^{0}' at 0.1,0.9
unset key
set xtics 0.2
set mxtics 2
plot 'lista2_3_2_0_1_0.txt' u 2:(exp(-0.5*($1-2.212768e+01))) w l dt 2 lc 1 lw 2 t 'abundance'
replot 'lista2_3_2_1_0_0.txt' u 2:(exp(-0.5*($1-5.869694e+02))) w l dt 3 lc 2 lw 3 t 'density profile'
replot 'lista2_3_2_0_0_1.txt' u 2:(exp(-0.5*($1-5.043077e+00))) w l dt 5 lc 4 lw 2 t 'bias'
replot 'lista2_3_2_1_1_1.txt' u 2:(exp(-0.5*($1-6.141470e+02))) w l dt 1 lc 7 lw 2 t 'all'
set xr[0.000000e+00:1.500000e+00]
set output 'histo_sl.png'
replot
reset


reset
set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'z_{SSB}'
set label '|f_{R0}|=10^{-6}' at 0.1,0.9
set key t r
set xtics 0.2
set mxtics 2
plot 'lista2_4_2_0_1_0.txt' u 2:(exp(-0.5*($1-1.407469e+01))) w l dt 2 lc 1 lw 2 t 'abundance'
replot 'lista2_4_2_1_0_0.txt' u 2:(exp(-0.5*($1-6.659247e+02))) w l dt 3 lc 2 lw 3 t 'density profile'
replot 'lista2_4_2_0_0_1.txt' u 2:(exp(-0.5*($1-4.292453e+00))) w l dt 5 lc 4 lw 2 t 'bias'
replot 'lista2_4_2_1_1_1.txt' u 2:(exp(-0.5*($1-7.006072e+02))) w l dt 1 lc 7 lw 2 t 'all'
set xr[0.000000e+00:2.000000e+00]
set output 'histo_s6.png'
replot
reset


reset
set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'z_{SSB}'
set label '|f_{R0}|=10^{-5}' at 1.1,0.9
unset key
set xtics 0.2
set mxtics 2
plot 'lista2_5_2_0_1_0.txt' u 2:(exp(-0.5*($1-3.816026e+01))) w l dt 2 lc 1 lw 2 t 'abundance'
replot 'lista2_5_2_1_0_0.txt' u 2:(exp(-0.5*($1-7.356159e+02))) w l dt 3 lc 2 lw 3 t 'density profile'
replot 'lista2_5_2_0_0_1.txt' u 2:(exp(-0.5*($1-6.404226e-01))) w l dt 5 lc 4 lw 2 t 'bias'
replot 'lista2_5_2_1_1_1.txt' u 2:(exp(-0.5*($1-8.523015e+02))) w l dt 1 lc 7 lw 2 t 'all'
set xr[1.000000e+00:3.000000e+00]
set output 'histo_s5.png'
replot
reset


reset
set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'z_{SSB}'
set label '|f_{R0}|=10^{-4}' at 1.7,0.9
unset key
set xtics 0.2
set mxtics 2
plot 'lista2_6_2_0_1_0.txt' u 2:(exp(-0.5*($1-1.405290e+02))) w l dt 2 lc 1 lw 2 t 'abundance'
replot 'lista2_6_2_1_0_0.txt' u 2:(exp(-0.5*($1-1.250966e+03))) w l dt 3 lc 2 lw 3 t 'density profile'
replot 'lista2_6_2_0_0_1.txt' u 2:(exp(-0.5*($1-1.144633e+01))) w l dt 5 lc 4 lw 2 t 'bias'
replot 'lista2_6_2_1_1_1.txt' u 2:(exp(-0.5*($1-1.568164e+03))) w l dt 1 lc 7 lw 2 t 'all'
set xr[1.600000e+00:3.600000e+00]
set output 'histo_s4.png'
replot
reset


