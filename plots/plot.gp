set terminal pngcairo enhanced font "Arial,12" size 800,600
set output "plot.png"

# Labels and title
set title  "Eigenvectors"
set xlabel "Atomic site"
set ylabel "Orbital coefficient"

# Grid and axes
set grid
set key top left  # Move legend to the top right
set border 3       # Draw axes on bottom and left only

# Set range (optional, adjust based on data)
#set xrange [0:3]  
set yrange [-0.15:0.2]

# Set styles
set style line 1 lc rgb "#0060ad" lw 1 pt 7 ps 1.5 # dark blue
set style line 2 lc rgb "#56B4E9" lw 1 pt 7 ps 1.5 # light blue
set style line 3 lc rgb "#009E73" lw 2 pt 7 ps 1.5 # green

plot "eigen_vectors_100at_eq_1.out" using 1:2 with lp ls 1 title "{/Symbol a_1} = 0; {/Symbol a_2 = 0}",  "eigen_vectors.out" u 1:2 w lp ls 3 title "{/Symbol a_1} = 0  ; {/Symbol a_2 = -1}  "#, "eigen_vectors_100at_eq_05.out" u 1:2 w lp ls 2 title "{/Symbol b_1} = -1; {/Symbol b_2 = -0.5}", "eigen_vectors_100at_eq_1.out" u 1:2 w p ls 3 title "{/Symbol b_1} = -1  ; {/Symbol b_2 = -1}  "

