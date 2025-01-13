#ifndef RTGNUPLOT_H
#define RTGNUPLOT_H
    #include <stdio.h>
    FILE* initialize_gnuplot();
    void update_gnuplot(FILE *gnuplot, int iteration, double residual);
    void finalize_gnuplot(FILE *gnuplot);
    void noop_update(FILE *gnuplot, int iteration, double residual);
    void noop_finalize(FILE *gnuplot);
#endif
