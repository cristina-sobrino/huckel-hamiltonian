set terminal pngcairo enhanced font "Arial,12" size 800,600
set output "plot_vectors.png"

# Labels and title
set title  "Eigenvectors"
set xlabel "Atomic site"
set ylabel "Orbital coefficient"

# Grid and axes
set grid
set key top left  # Move legend to the top right
set border 3       # Draw axes on bottom and left only

# Set range (optional, adjust based on data)
#set xrange [0:3]  # Change according to your data
#set yrange [-0.2:0.2]

# Set styles
set style line 1 lc rgb "#0060ad" lw 2 pt 7 ps 1.5  # Blue
set style line 2 lc rgb "#dd181f" lw 2 pt 7 ps 1.5  # Red
set style line 3 lc rgb "#00ad00" lw 2 pt 7 ps 1.5  # Green


plot "eigen_vectors.out" using 1:2 with lp ls 1 title "First eigenvector ", \
     "eigen_vectors.out" using 1:3 with lp ls 2 lc 2 title "Second eigenvector ", \
     "eigen_vectors.out" using 1:4 with lp ls 1 lc 3 title "Third eigenvector "

