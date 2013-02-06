
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <math.h>
#include <lapacke.h>

#include "hdf5.h"
#include "hdf5_hl.h"


void test_dstegr(void);

int main(int argc, char *argv[]) {
    
//    test_dstegr();

    if (argc != 4 + 1) {
        printf("Usage: %s N W K output.h5\n", argv[0]);
        return -1;
    }
    
    lapack_int  N = atoi(argv[1]);
    double      W = atof(argv[2]);
    lapack_int  K = atoi(argv[3]);
    char *outname = argv[4];
    
//    printf("N = %d\nW = %g\nK = %d\noutname = %s\n", N, W, K, outname);

    if (!(0 < K && K <= N && 0 < W && W < 0.5)) {
        printf("The arguments must satisfy 0 < K <= N and 0 < W < 0.5\n");
        return -1;
    }

    double *d = calloc(N, sizeof(double));
    double *e = calloc(N, sizeof(double));
    double *w = calloc(N, sizeof(double));
    double *z = calloc(K * N, sizeof(double));
    lapack_int m   = 0;
    lapack_int ldz = N;
    lapack_int *isuppz = calloc(2 * K, sizeof(lapack_int));

    assert(d && e && w && z && isuppz);

    double cos_two_pi_W = cos(2 * 3.14159265358979323846 * W);
    double x;
    for (int i = 0; i < N; i++) {
        x = (0.5 * (N - 1) - i);
        d[i] = x * x * cos_two_pi_W;
    }

    for (int i = 0; i < N - 1; i++)
        e[i] = 0.5 * (i + 1) * (N - i - 1);

    
/* lapack_int LAPACKE_dstegr( int matrix_order, char jobz, char range, */
/*                            lapack_int n, double* d, double* e, double vl, */
/*                            double vu, lapack_int il, lapack_int iu, */
/*                            double abstol, lapack_int* m, double* w, double* z, */
/*                            lapack_int ldz, lapack_int* isuppz ); */

    printf("Before DSTEGR\n");
    lapack_int info = LAPACKE_dstegr(LAPACK_COL_MAJOR, 'V', 'I', 
                                     N, d, e, 0, 
                                     0, N - K + 1, N,
                                     0, &m, w, z,
                                     ldz, isuppz);
    printf("After DSTEGR\n");
    printf("K = %d\nm = %d\n", K, m);
    
    if (info) {
        printf("Some error occurred in DSTEGR, info = %d\n", info);
        return -1;
    }

    hid_t   file_id;
    hsize_t dims_w[1] = {K};
    hsize_t dims_z[2] = {K, N};
    hsize_t dims_scalar[1] = {1};
    herr_t  status;
    
    file_id = H5Fcreate(outname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    assert(file_id >= 0);
    status  = H5LTmake_dataset_double(file_id, "/eigenvalues",  1, dims_w, w);
    assert(status >= 0);
    status  = H5LTmake_dataset_double(file_id, "/eigenvectors", 2, dims_z, z);
    assert(status >= 0);
    status  = H5LTmake_dataset_double(file_id, "/W", 1, dims_scalar, &W);
    assert(status >= 0);
    status  = H5LTmake_dataset_int(file_id, "/K", 1, dims_scalar, &K);
    assert(status >= 0);
    status  = H5LTmake_dataset_int(file_id, "/N", 1, dims_scalar, &N);
    assert(status >= 0);
    status  = H5Fclose (file_id);
    assert(status >= 0);

    /* free(d);  */
    /* free(e); */
    /* free(w); */
    /* free(z); */
    /* free(isuppz); */

    return 0;
    
}


void test_dstegr(void) {

    int matrix_order = LAPACK_COL_MAJOR;
    char jobz = '\0', range = '\0';
    lapack_int n = 0;
    double *d = NULL, 
           *e = NULL;
    double vl = 0, 
           vu = 0;
    lapack_int il = 4,
               iu = 5;
    double abstol = 0; /* not used */
    lapack_int *m = NULL;
    double *w = NULL;
    double *z = NULL;
    lapack_int ldz = 0;
    lapack_int *isuppz = NULL;
    
    lapack_int info;

    jobz = 'V';
    range = 'I';
    n = 5;
    d = calloc(n, sizeof(double));
    e = calloc(n, sizeof(double));
    w = calloc(n, sizeof(double));
    z = calloc(n * n, sizeof(double));
    m = calloc(1, sizeof(lapack_int));
    ldz = n;
    isuppz = calloc(2 * n, sizeof(lapack_int));

    assert(d && e && w && z && m);

    for (int i = 0; i < n; i++) 
        d[i] = i;
    for (int i = 0; i < n - 1; i++) 
        e[i] = i;
    
    info = LAPACKE_dstegr(matrix_order, jobz, range,
                          n, d, e, 
                          vl, vu, il, iu, abstol, 
                          m, w, z, ldz, isuppz);
    
    printf("info = %d\n", info);
    printf("m    = %d\n", *m);
    printf("w    =\n");
    for (int i = 0; i < *m; i++) 
        printf("\t% f", w[i]);
    printf("\n");
    printf("z    =\n");
    for (int j = 0; j < *m; j++) {
        for (int i = 0; i < n; i++)
            printf("\t% f", z[j * ldz + i]);
        printf("\n");
    }

}
