set terminal pngcairo enhanced font "Arial,12" size 800,600
set output "plot_rings.png"

# Labels and title
set title  "Eigenvalues"
set xlabel "State"
set ylabel "Energies (A.U.)"

# Grid and axes
set grid
set key top left  # Move legend to the top right
set border 3       # Draw axes on bottom and left only

# Set range (optional, adjust based on data)
#set xrange [0:100]  # Change according to your data
set yrange [-2:2]

# Set styles
set style line 1 lc rgb "#0060ad" lw 1 pt 7 ps 1.5
set style line 2 lc rgb "#00ad00" lw 1 pt 7 ps 1.5

plot "energies_ring_100.out" using 1:2 with p ls 1 title "100 atoms" , \
     "energies_ring_4.out" using 1:2 with p ls 2 title "4 atoms"

