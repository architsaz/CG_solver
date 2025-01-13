#include <stdio.h>
#include <stdlib.h>

// Function to initialize Gnuplot
FILE* initialize_gnuplot(void) {
    FILE *gnuplot = popen("gnuplot -persistent", "w");
    if (!gnuplot) {
        fprintf(stderr, "Error: Could not open pipe to Gnuplot.\n");
        exit(EXIT_FAILURE);
    }

    // Gnuplot configuration
    fprintf(gnuplot, "set terminal pngcairo\n");
    fprintf(gnuplot, "set output 'residual_plot.png'\n");
    fprintf(gnuplot, "set title 'Residual vs. Iteration'\n");
    fprintf(gnuplot, "set xlabel 'Iteration'\n");
    fprintf(gnuplot, "set ylabel 'Residual'\n");
    fprintf(gnuplot, "set grid\n");
    fprintf(gnuplot, "plot '-' with lines title 'Residual'\n");
    fflush(gnuplot);

    return gnuplot;
}

// Function to send data to Gnuplot
void update_gnuplot(FILE *gnuplot, int iteration, double residual) {
    fprintf(gnuplot, "%d %f\n", iteration, residual);
    fflush(gnuplot);
}

// Function to finalize Gnuplot
void finalize_gnuplot(FILE *gnuplot) {
    fprintf(gnuplot, "e\n"); // End of data stream
    fflush(gnuplot);
    pclose(gnuplot);         // Close the pipe
}

// "No operation" functions to disable Gnuplot actions
void noop_update(FILE *gnuplot, int iteration, double residual) {
    (void)gnuplot; (void)iteration; (void)residual; // Do nothing
}

void noop_finalize(FILE *gnuplot) {
    (void)gnuplot; // Do nothing
}