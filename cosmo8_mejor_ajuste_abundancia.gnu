reset
set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,22' 
set logscale
set xr [1.7:19]
set yr [1.e-6:5.e-3]
set format y "10^{%L}"
set xlabel 'r_{0.2} [Mpc/h]'
set ylabel 'dn/dln r_{0.2} [Mpc/h]^{-3}'

  plot      'dnv_esferas_512x256_03.txt' u 1:2:($2-$3):($2+$3) with yerrorbars lc 1 t '{/Symbol L}CDM'
replot '../6/dnv_esferas_512x256_03.txt' u 1:2:($2-$3):($2+$3) with yerrorbars lc 2 t '|f_{R0}|=10^{-6}'
replot '../5/dnv_esferas_512x256_03.txt' u 1:2:($2-$3):($2+$3) with yerrorbars lc 3 t '|f_{R0}|=10^{-5}'
replot '../4/dnv_esferas_512x256_03.txt' u 1:2:($2-$3):($2+$3) with yerrorbars lc 4 t '|f_{R0}|=10^{-4}'
replot 'mejor_ajuste_abundancia_l_1.dat' w l lw 1 lc 1 t ''
replot 'mejor_ajuste_abundancia_6_1.dat' w l lw 1 lc 2 t ''
replot 'mejor_ajuste_abundancia_5_1.dat' w l lw 1 lc 3 t ''
replot 'mejor_ajuste_abundancia_4_1.dat' w l lw 1 lc 4 t ''
set output 'best_abundance_f.png' 
replot
set xr [1.7:6]
set yr [3.e-4:5.e-3]
set output 'best_abundance_f_zoom.png'
replot
reset








reset
set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,22'
set logscale
set xr [1.7:19]
set yr [1.e-6:5.e-3]
set format y "10^{%L}"
set xlabel 'r_{0.2} [Mpc/h]'
set ylabel 'dn/dln r_{0.2} [Mpc/h]^{-3}'

  plot      'dnv_esferas_512x256_03.txt' u 1:2:($2-$3):($2+$3) with yerrorbars lc 1 t '{/Symbol L}CDM'
replot '../A/dnv_esferas_512x256_03.txt' u 1:2:($2-$3):($2+$3) with yerrorbars lc 5 t 'z_{SSB}=1'
replot '../B/dnv_esferas_512x256_03.txt' u 1:2:($2-$3):($2+$3) with yerrorbars lc 6 t 'z_{SSB}=2'
replot '../D/dnv_esferas_512x256_03.txt' u 1:2:($2-$3):($2+$3) with yerrorbars lc 7 t 'z_{SSB}=3'
replot 'mejor_ajuste_abundancia_l_2.dat' w l lw 1 lc 1 t ''
replot 'mejor_ajuste_abundancia_A_2.dat' w l lw 1 lc 5 t ''
replot 'mejor_ajuste_abundancia_B_2.dat' w l lw 1 lc 6 t ''
replot 'mejor_ajuste_abundancia_D_2.dat' w l lw 1 lc 7 t ''
set output 'best_abundance_s.png'
replot
set xr [1.7:6]
set yr [3.e-4:5.e-3]
set output 'best_abundance_s_zoom.png'
replot
reset


