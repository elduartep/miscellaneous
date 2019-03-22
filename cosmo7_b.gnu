reset
set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'log_{10}|f_{R0}|'
set label '|f_{R0}|=10^{-4}' at -3.3,0.9
unset key
set xtics 0.2
set mxtics 2
plot 'lista_0_1_0_1_0.txt' u 2:(exp(-0.5*($1-2.539107e+01))) w l dt 2 lc 1 lw 2 t 'abundance'
replot 'lista_0_1_1_0_0.txt' u 2:(exp(-0.5*($1-7.826989e+02))) w l dt 3 lc 2 lw 3 t 'density profile'
replot 'lista_0_1_0_0_1.txt' u 2:(exp(-0.5*($1-7.242166e+00))) w l dt 5 lc 4 lw 2 t 'bias'
replot 'lista_0_1_1_1_1.txt' u 2:(exp(-0.5*($1-8.156661e+02))) w l dt 1 lc 7 lw 2 t 'all'
set xr[-5.00000e+00:-3.00000e+00]
set output 'histo_f4_b.png'
replot
reset


reset
set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'log_{10}|f_{R0}|'
set label '|f_{R0}|=10^{-5}' at -4.3,0.9
unset key
set xtics 0.2
set mxtics 2
plot 'lista_1_1_0_1_0.txt' u 2:(exp(-0.5*($1-1.750140e+01))) w l dt 2 lc 1 lw 2 t 'abundance'
replot 'lista_1_1_1_0_0.txt' u 2:(exp(-0.5*($1-6.442387e+02))) w l dt 3 lc 2 lw 3 t 'density profile'
replot 'lista_1_1_0_0_1.txt' u 2:(exp(-0.5*($1-5.392222e-01))) w l dt 5 lc 4 lw 2 t 'bias'
replot 'lista_1_1_1_1_1.txt' u 2:(exp(-0.5*($1-6.636076e+02))) w l dt 1 lc 7 lw 2 t 'all'
set xr[-6.00000e+00:-4.00000e+00]
set output 'histo_f5_b.png'
replot
reset


reset
set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'log_{10}|f_{R0}|'
set label '|f_{R0}|=10^{-6}' at -5.3,0.9
unset key
set xtics 0.2
set mxtics 2
plot 'lista_2_1_0_1_0.txt' u 2:(exp(-0.5*($1-1.293482e+01))) w l dt 2 lc 1 lw 2 t 'abundance'
replot 'lista_2_1_1_0_0.txt' u 2:(exp(-0.5*($1-6.267569e+02))) w l dt 3 lc 2 lw 3 t 'density profile'
replot 'lista_2_1_0_0_1.txt' u 2:(exp(-0.5*($1-4.385109e+00))) w l dt 5 lc 4 lw 2 t 'bias'
replot 'lista_2_1_1_1_1.txt' u 2:(exp(-0.5*($1-6.451407e+02))) w l dt 1 lc 7 lw 2 t 'all'
set xr[-7.00000e+00:-5.00000e+00]
set output 'histo_f6_b.png'
replot
reset


reset
set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'log_{10}|f_{R0}|'
set label '{/Symbol L}CDM' at -6,0.4
set key t r
set xtics 0.5
set mxtics 5
set xlabel 'log_{10}|f_{R0}|'
#set format x '%2.t'
set xtics 1
 set mxtics 5
plot 'lista_3_1_0_1_0.txt' u 2:(exp(-0.5*($1-1.069846e+01))) w l dt 2 lc 1 lw 2 t 'abundance'
replot 'lista_3_1_1_0_0.txt' u 2:(exp(-0.5*($1-5.837314e+02))) w l dt 3 lc 2 lw 3 t 'density profile'
replot 'lista_3_1_0_0_1.txt' u 2:(exp(-0.5*($1-4.558067e+00))) w l dt 5 lc 4 lw 2 t 'bias'
replot 'lista_3_1_1_1_1.txt' u 2:(exp(-0.5*($1-5.991741e+02))) w l dt 1 lc 7 lw 2 t 'all'
set xr[-10:-5.5]
set output 'histo_fl_b.png'
replot
reset


reset
set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'z_{SSB}'
set label '{/Symbol L}CDM' at 0.85,0.4
set key t r
set xtics 0.2
set mxtics 2
plot 'lista_3_2_0_1_0.txt' u 2:(exp(-0.5*($1-2.212768e+01))) w l dt 2 lc 1 lw 2 t 'abundance'
replot 'lista_3_2_1_0_0.txt' u 2:(exp(-0.5*($1-5.869694e+02))) w l dt 3 lc 2 lw 3 t 'density profile'
replot 'lista_3_2_0_0_1.txt' u 2:(exp(-0.5*($1-5.043077e+00))) w l dt 5 lc 4 lw 2 t 'bias'
replot 'lista_3_2_1_1_1.txt' u 2:(exp(-0.5*($1-6.141470e+02))) w l dt 1 lc 7 lw 2 t 'all'
set xr[0.000000e+00:10.000000e-01]
set output 'histo_sl_b.png'
replot
reset


reset
set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'z_{SSB}'
set label 'z_{SSB}=1' at 1.7,0.9
unset key
set xtics 0.2
set mxtics 2
plot 'lista_4_2_0_1_0.txt' u 2:(exp(-0.5*($1-9.405609e+00))) w l dt 2 lc 1 lw 2 t 'abundance'
replot 'lista_4_2_1_0_0.txt' u 2:(exp(-0.5*($1-6.963311e+02))) w l dt 3 lc 2 lw 3 t 'density profile'
replot 'lista_4_2_0_0_1.txt' u 2:(exp(-0.5*($1-2.790169e+00))) w l dt 5 lc 4 lw 2 t 'bias'
replot 'lista_4_2_1_1_1.txt' u 2:(exp(-0.5*($1-7.086669e+02))) w l dt 1 lc 7 lw 2 t 'all'
set xr[0.000000e-01:2.00000e+00]
set output 'histo_sA_b.png'
replot
reset


reset
set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'z_{SSB}'
set label 'z_{SSB}=2' at 2.7,0.9
unset key
set xtics 0.2
set mxtics 2
plot 'lista_5_2_0_1_0.txt' u 2:(exp(-0.5*($1-2.455454e+01))) w l dt 2 lc 1 lw 2 t 'abundance'
replot 'lista_5_2_1_0_0.txt' u 2:(exp(-0.5*($1-6.551378e+02))) w l dt 3 lc 2 lw 3 t 'density profile'
replot 'lista_5_2_0_0_1.txt' u 2:(exp(-0.5*($1-2.214233e+00))) w l dt 5 lc 4 lw 2 t 'bias'
replot 'lista_5_2_1_1_1.txt' u 2:(exp(-0.5*($1-6.826433e+02))) w l dt 1 lc 7 lw 2 t 'all'
set xr[1.00000e+00:3.00000e+00]
set output 'histo_sB_b.png'
replot
reset


reset
set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,22' dashlength 1.7
unset ytics
set yr[0:1.05]
set xlabel 'z_{SSB}'
set label 'z_{SSB}=3' at 3.7,0.9
unset key
set xtics 0.2
set mxtics 2
plot 'lista_6_2_0_1_0.txt' u 2:(exp(-0.5*($1-1.797614e+01))) w l dt 2 lc 1 lw 2 t 'abundance'
replot 'lista_6_2_1_0_0.txt' u 2:(exp(-0.5*($1-8.241263e+02))) w l dt 3 lc 2 lw 3 t 'density profile'
replot 'lista_6_2_0_0_1.txt' u 2:(exp(-0.5*($1-6.102131e+00))) w l dt 5 lc 4 lw 2 t 'bias'
replot 'lista_6_2_1_1_1.txt' u 2:(exp(-0.5*($1-8.491182e+02))) w l dt 1 lc 7 lw 2 t 'all'
set xr[2.00000e+00:4.00000e+00]
set output 'histo_sD_b.png'
replot
reset


