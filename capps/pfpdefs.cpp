//
// Created by nicolet on 31/10/17.
//

#include "pfpdefs.hpp"

void dclean(double *a, msize_t n){
    for (msize_t i = 0; i < n; ++i) {
        a[i] = 0.0;
    }
}

void propagate_error(double* a, double err, msize_t k, msize_t n) {
    for(msize_t i = k; i < n; i++){
        a[i] += err;
    }
}

void println_darray(double *a, msize_t from, msize_t to){
    if(from <= to - 1) {
        printf("[");
        for (msize_t i = from; i < to - 1; i++) {
            printf("%f, ", a[i]);
        }
        printf("%f]\n", a[to - 1]);
    } else {
        printf("[]\n");
    }
}


double mymax(double r,double n) {
    if (r > n)
        return r;
    return n;
}

double myparallelmax(double* data, int ndata) {
    double m;
    m = data[0];
#pragma omp declare reduction \
      (reduce_max:double:omp_out=mymax(omp_out,omp_in)) \
      initializer(omp_priv=0)

#pragma omp parallel for reduction(reduce_max:m)
    for (int idata=0; idata<ndata; idata++)
        m = mymax(m,data[idata]);

    return m;
}


_mts_ mtsjoin(_mts_ r,_mts_ n) {
    // r is the already reduced value
    // n is the new value
    double new_mts;
    new_mts = n.sum + r.mts;
    if (new_mts < n.mts){
        new_mts = n.mts;
    }
    r.sum = r.sum + n.sum;
    r.mts = new_mts;
    return r;
}

_mts_ myparallelmts(double* data, int ndata) {
    _mts_ m;
    m = {data[0], 0};
#pragma omp declare reduction \
      (mts_reduction:_mts_:omp_out=mtsjoin(omp_out,omp_in)) \
      initializer(omp_priv={0,0})

#pragma omp parallel for reduction(mts_reduction:m) num_threads(4)
    for (int idata=0; idata<ndata; idata++){
        m.sum += data[idata];
        m.mts = (m.mts + data[idata] > 0) ? m.mts + data[idata] : 0;
    }

    return m;
}

_mts_ custom_reduce_mts(double* data, int ndata) {
    _mts_ m;
    m = {data[0], 0};

#pragma omp parallel
    {
        _mts_ local_mts = {0, 0};

#pragma omp for nowait
        for (int idata = 0; idata < ndata; idata++) {
            local_mts.sum += data[idata];
            local_mts.mts = (local_mts.mts + data[idata] > 0) ? local_mts.mts + data[idata] : 0;
        }
    }

    return m;
}
