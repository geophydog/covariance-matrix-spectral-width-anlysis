#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>
#include <gsl/all.h>
#include "sacio.h"

int pow2_next(int n) {
    float f;
    int m;
    f = log((float)n)/log(2.) + 1;
    m = (int)pow(2., f);
    return m;
}

void no_spa( char *ps) {
    char *pt = ps;
    while ( *ps != '\0' ) {
        if ( *ps != ' ' && *ps != '\n' ) {
            *pt++ = *ps;
        }
        ++ps;
    }
    *pt = '\0';
}

int main (int argc, char *argv[]) {
    int i, j, k, m, n, seg_num, sub_seg_num, sub_seg_npts, seg_npts, sac_npts, fftw_npts, count = 0, fre_index, sta_index = 0;
    double real, imag, fs, fre, t1, t2, fre1, fre2, fre_step, delta, seg_t, sub_seg_t, time, **dat,**fft, ***tmp, tmp1, tmp2, width, time_start, time_end;
    float *data;
    fftw_complex *in, *out;
    fftw_plan p;
    FILE *fin, *fout, *fp;
    char ss[256];
    SACHEAD hd;

    if ( argc != 10 ) {
        fprintf(stderr,"+----------------------------------------------------------------------------------------------------+\n");
        fprintf(stderr,"+   Usage: cov <listfile> <t1> <t2> <fre1> <fre2> <fre_step> <seg_t> <sub_seg_t> <output_file>       +\n");
        fprintf(stderr,"+          <listfile>     [1]   file containing SAC format files                                     +\n");
        fprintf(stderr,"+          <t1>           [2]   begining time of executing program                                   +\n");
        fprintf(stderr,"+          <t2>           [3]   ending time of executing program                                     +\n");
        fprintf(stderr,"+          <fre1>         [4]   low corner frequency limitation                                      +\n");
        fprintf(stderr,"+          <fre2>         [5]   high corner frequency limitation                                     +\n");
        fprintf(stderr,"+          <fre_step>     [6]   frequency step length of executing program                           +\n");
        fprintf(stderr,"+          <seg_t>        [7]   segmentation length of executing program                             +\n");
        fprintf(stderr,"+          <sub_seg_t>    [8]   sub-segmentation length of executing program                         +\n");
        fprintf(stderr,"+          <output_file>  [9]   outputting results file, containing 3 columns:                       +\n");
        fprintf(stderr,"+                               col1: time  col2: frequency  col3: covariance matrix spectral width  +\n");
        fprintf(stderr,"+----------------------------------------------------------------------------------------------------+\n");
        exit(1);
    }

    time_start = clock();
    fin = fopen(argv[1], "r");
    t1 = atof(argv[2]); t2 = atof(argv[3]); fre1 = atof(argv[4]); fre2 = atof(argv[5]); fre_step = atof(argv[6]);
    seg_t = atof(argv[7]); sub_seg_t = atof(argv[8]); fout = fopen(argv[9],"w");
    while ( fgets(ss, 256, fin) ) {
        no_spa(ss);
        data = read_sac(ss, &hd);
        delta = hd.delta;
        count += 1;
    }
    fclose(fin);
/*-------------------------------------------------------------------------------------------------------------------------------------------*/
    gsl_complex z;
    gsl_matrix_complex *cov = gsl_matrix_complex_alloc(count, count);
    gsl_vector *eval = gsl_vector_alloc(count);
    gsl_eigen_herm_workspace *w = gsl_eigen_herm_alloc(count);

    fs = 1./delta; seg_npts = (int)(seg_t/delta); sub_seg_npts = (int)(sub_seg_t/delta);
    fftw_npts = pow2_next(sub_seg_npts); seg_num = (int)((t2-t1)/seg_t); sub_seg_num = (int)(seg_t/sub_seg_t); sac_npts = (int)((t2-t1)/delta);
    printf("seg_num: %d  sub_seg_num: %d  seg_npts: %d  sub_seg_npts: %d  t1: %g  t2:%g  fre1: %g  fre2: %g  fre_step: %g  seg_t: %g  sub_seg_t: %g sac_npts: %d\n", seg_num, sub_seg_num, seg_npts, sub_seg_npts, t1, t2,fre1, fre2, fre_step, seg_t, sub_seg_t, sac_npts);

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftw_npts);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftw_npts);
    dat = (double**) malloc(sizeof(double) * count);
    fft = (double**) malloc(sizeof(double*) * count);
    for ( i = 0; i < count; i ++ ) {
        fft[i] = (double*) malloc(sizeof(double) * 2);
        dat[i] = (double*) malloc(sizeof(double) * sac_npts);
    }

    fin = fopen(argv[1], "r");
    while( fgets(ss, 256, fin) ) {
        no_spa(ss);
        data = read_sac(ss,&hd);
        for ( i = 0; i < sac_npts; i ++ )
            dat[sta_index][i] = (double)data[i];
        sta_index += 1;
    }
    fclose(fin); free(data);

    tmp = (double***) malloc(sizeof(double**) * count);
    for ( i = 0; i < count; i ++ ) {
        tmp[i] = (double**) malloc(sizeof(double*) * count);
        for ( j = 0; j < count; j ++ )
            tmp[i][j] = (double*) malloc(sizeof(double) * 2);
    }

    fre = fre1;
    while ( fre <= fre2 ) {
        printf("Now computing covariance matrix spectral width of frequency %g Hz\n", fre);
        fre_index = (int) (fre*fftw_npts/fs);
        for ( i = 0; i < seg_num; i ++ ) {
            time = t1 + i*seg_t + seg_t/2.;
            real = 0.; imag = 0.;
            for ( m = 0; m < count; m++ )
                for ( n = 0; n < count; n ++ ) {
                    tmp[m][n][0] = 0.; tmp[m][n][1] = 0.;
                }
            for ( j = 0; j < sub_seg_num; j ++ ){
                    sta_index = 0;
                    while ( sta_index < count ) {
                        for ( k = 0; k < fftw_npts; k++ ) {
                            if ( k < sub_seg_npts ) in[k][0] = dat[sta_index][i*seg_npts+j*sub_seg_npts+k];
                            else in[k][0] = 0.;
                            in[k][1] = 0.;
                        }
                            p = fftw_plan_dft_1d(fftw_npts, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
                            fftw_execute(p);
                            fft[sta_index][0] = out[fre_index][0];
                            fft[sta_index][1] = out[fre_index][1];
                            fftw_destroy_plan(p);
                        sta_index += 1;
                    }
                    for ( m = 0; m < count; m ++ )
                        for ( n = 0; n < count; n ++ ) {
                            tmp[m][n][0] += (fft[m][0]*fft[n][0] + fft[m][1]*fft[n][1]);
                            tmp[m][n][1] += (fft[n][0]*fft[m][1] - fft[m][0]*fft[n][1]);
                        }
            }

            for ( m = 0; m < count; m ++ )
                for ( n = 0; n < count; n ++ ) {
                    //printf("%g + %gi    ||   %g + %gi\n", tmp[m][n][0], tmp[m][n][1], tmp[n][m][0], tmp[n][m][1]);
                    //if( n == (count-1) ) printf("\n");
                    GSL_SET_COMPLEX(&z, tmp[m][n][0], tmp[m][n][1]);
                    gsl_matrix_complex_set(cov, m, n, z);
                }
            gsl_eigen_herm(cov, eval, w);
            tmp1 = 0.; tmp2 = 0.;
            for ( m = 0; m < count; m ++ ) {
                double eval_m = gsl_vector_get(eval, m);
                printf("eigenvalue %d : %g\n", m+1, eval_m);
                tmp1 += (m*eval_m); tmp2 += eval_m;
            }
            //printf("tmp1: %g  tmp2: %g\n", tmp1, tmp2);
            width = tmp1/tmp2;
            fprintf(fout,"%g %g %g\n", time, fre, width);
        }
        fre += fre_step;
    }

    printf("EIGENVALUES FOUND FINISHING!!!\n");

    fclose(fout); fftw_free(in); fftw_free(out);
    gsl_vector_free(eval); gsl_eigen_herm_free(w); gsl_matrix_complex_free(cov);
    for ( i = 0; i < count; i ++ ) {
        free(fft[i]); free(dat[i]);
    }
    free(fft); free(dat);

    for ( m = 0; m < count; m ++ ) {
        for ( n = 0; n < count; n ++ ) {
            free(tmp[m][n]);
        }
        free(tmp[m]);
    }
    free(tmp);
    time_end = clock();
    printf("TIme used: %g seconds!!!\n", (time_end-time_start)/CLOCKS_PER_SEC);
    return 0;
}
