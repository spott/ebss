#pragma once

//call DSBEVD( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO )
extern "C" { 
    void dsbevd_( const char* jobz, const char* uplo, int* n, int* kd, double* ab, int* ldab, double* w, double* z, int* ldz, double* work, int* lwork, int* iwork, int* liwork, int* info); 
    void dsbev_( const char* jobz, const char* uplo, int* n, int* kd, double* ab, int* ldab, double* w, double* z, int* ldz, double* work, int* info); 
}

int dsbevd(int n, int kd, double *input, double *evalues)
{
    int ldab = kd + 1;
    int lwork = 2 * n;
    int liwork = 1;
    double *work = new double[lwork];
    int *iwork = new int[liwork];
    int info;
    double *null = new double;

    dsbevd_("N","U", &n, &kd, input, &ldab, evalues, null, &n, work, &lwork, iwork, &liwork,  &info);

	delete[] null;
    delete[] work;
    delete[] iwork;

    return info;
}

int dsbevd(int n, int kd, double *input, double *evalues, double *evectors)
{

    int ldab = kd+1;
    int lwork = 1 + 5*n + 2*n*n;
    int liwork = 3 + 5*n;
    double *work = new double[lwork];
    int *iwork = new int[liwork];
    int info;

    dsbevd_("V","U", &n, &kd, input, &ldab, evalues, evectors, &n, work, &lwork, iwork, &liwork,  &info);

    delete[] work;
    delete[] iwork;
    return info;
}

int dsbev(int n, int kd, double *input, double *evalues)
{

    int ldab = kd+1;
    double *work = new double[3*n-2];
    int info;
    double *null = new double;

    dsbev_("V","U", &n, &kd, input, &ldab, evalues, null, &n, work, &info);

	delete[] null;
    delete[] work;
    return info;
}
int dsbev(int n, int kd, double *input, double *evalues, double *evectors)
{

    int ldab = kd+1;
    double *work = new double[3*n-2];
    int info;

    dsbev_("V","U", &n, &kd, input, &ldab, evalues, evectors, &n, work, &info);

    delete[] work;
    return info;
}
