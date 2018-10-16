reset
set terminal pngcairo size 1000,750 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'log_{10}|f_{R0}|'
unset key
set xtics 0.1
set mxtics 2
plot 'lista_0_1_0_1.txt' u 2:(exp(-0.5*($1-2.539107e+01))) w l lw 2 dt 2 lc 1 t 'abundance'
replot 'lista_0_1_1_0.txt' u 2:(exp(-0.5*($1-7.826989e+02))) w l lw 2 dt 5 lc 4 t 'density profile'
replot 'lista_0_1_1_1.txt' u 2:(exp(-0.5*($1-8.084227e+02))) w l lw 2 dt 1 lc 7 t 'both'
set xr[-4.200000e+00:-3.800000e+00]
set output 'histo_f4.png'
replot
reset


reset
set terminal pngcairo size 1000,750 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'log_{10}|f_{R0}|'
unset key
set xtics 0.1
set mxtics 2
plot 'lista_1_1_0_1.txt' u 2:(exp(-0.5*($1-1.750140e+01))) w l lw 2 dt 2 lc 1 t 'abundance'
replot 'lista_1_1_1_0.txt' u 2:(exp(-0.5*($1-6.442387e+02))) w l lw 2 dt 5 lc 4 t 'density profile'
replot 'lista_1_1_1_1.txt' u 2:(exp(-0.5*($1-6.630666e+02))) w l lw 2 dt 1 lc 7 t 'both'
set xr[-5.200000e+00:-4.800000e+00]
set output 'histo_f5.png'
replot
reset


reset
set terminal pngcairo size 1000,750 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'log_{10}|f_{R0}|'
unset key
set xtics 0.1
set mxtics 2
plot 'lista_2_1_0_1.txt' u 2:(exp(-0.5*($1-1.293482e+01))) w l lw 2 dt 2 lc 1 t 'abundance'
replot 'lista_2_1_1_0.txt' u 2:(exp(-0.5*($1-6.267569e+02))) w l lw 2 dt 5 lc 4 t 'density profile'
replot 'lista_2_1_1_1.txt' u 2:(exp(-0.5*($1-6.407556e+02))) w l lw 2 dt 1 lc 7 t 'both'
set xr[-6.200000e+00:-5.800000e+00]
set output 'histo_f6.png'
replot
reset


reset
set terminal pngcairo size 1000,750 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'log_{10}|f_{R0}|'
set key t r
set xtics 0.2
set mxtics 4
set xlabel '|f_{R0}|'
set format x '%2.0t{/Symbol ´}10^{%L}'
set xtics 2e-8
 set mxtics 4
plot 'lista_3_1_0_1.txt' u (10.**$2):(exp(-0.5*($1-1.069846e+01))) w l lw 2 dt 2 lc 1 t 'abundance'
replot 'lista_3_1_1_0.txt' u (10.**$2):(exp(-0.5*($1-5.837314e+02))) w l lw 2 dt 5 lc 4 t 'density profile'
replot 'lista_3_1_1_1.txt' u (10.**$2):(exp(-0.5*($1-5.946130e+02))) w l lw 2 dt 1 lc 7 t 'both'
set xr[0.000000e+00:7.000000e-08]
set output 'histo_fl.png'
replot
reset


reset
set terminal pngcairo size 1000,750 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'z_{SSB}'
set key t r
set xtics 0.1
set mxtics 2
plot 'lista_3_2_0_1.txt' u 2:(exp(-0.5*($1-2.212768e+01))) w l lw 2 dt 2 lc 1 t 'abundance'
replot 'lista_3_2_1_0.txt' u 2:(exp(-0.5*($1-5.869694e+02))) w l lw 2 dt 5 lc 4 t 'density profile'
replot 'lista_3_2_1_1.txt' u 2:(exp(-0.5*($1-6.091039e+02))) w l lw 2 dt 1 lc 7 t 'both'
set xr[0.000000e+00:4.000000e-01]
set output 'histo_sl.png'
replot
reset


reset
set terminal pngcairo size 1000,750 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'z_{SSB}'
unset key
set xtics 0.1
set mxtics 2
plot 'lista_4_2_0_1.txt' u 2:(exp(-0.5*($1-9.405609e+00))) w l lw 2 dt 2 lc 1 t 'abundance'
replot 'lista_4_2_1_0.txt' u 2:(exp(-0.5*($1-6.963311e+02))) w l lw 2 dt 5 lc 4 t 'density profile'
replot 'lista_4_2_1_1.txt' u 2:(exp(-0.5*($1-7.058757e+02))) w l lw 2 dt 1 lc 7 t 'both'
set xr[8.000000e-01:1.200000e+00]
set output 'histo_sA.png'
replot
reset


reset
set terminal pngcairo size 1000,750 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'z_{SSB}'
unset key
set xtics 0.1
set mxtics 2
plot 'lista_5_2_0_1.txt' u 2:(exp(-0.5*($1-2.455454e+01))) w l lw 2 dt 2 lc 1 t 'abundance'
replot 'lista_5_2_1_0.txt' u 2:(exp(-0.5*($1-6.551378e+02))) w l lw 2 dt 5 lc 4 t 'density profile'
replot 'lista_5_2_1_1.txt' u 2:(exp(-0.5*($1-6.804290e+02))) w l lw 2 dt 1 lc 7 t 'both'
set xr[1.800000e+00:2.200000e+00]
set output 'histo_sB.png'
replot
reset


reset
set terminal pngcairo size 1000,750 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'z_{SSB}'
unset key
set xtics 0.1
set mxtics 2
plot 'lista_6_2_0_1.txt' u 2:(exp(-0.5*($1-1.797614e+01))) w l lw 2 dt 2 lc 1 t 'abundance'
replot 'lista_6_2_1_0.txt' u 2:(exp(-0.5*($1-8.241263e+02))) w l lw 2 dt 5 lc 4 t 'density profile'
replot 'lista_6_2_1_1.txt' u 2:(exp(-0.5*($1-8.430150e+02))) w l lw 2 dt 1 lc 7 t 'both'
set xr[2.800000e+00:3.200000e+00]
set output 'histo_sD.png'
replot
reset


