#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <lapacke.h>

void * Malloc(size_t nbytes)
{
    void * ptr = malloc(nbytes);
    if (ptr == 0x0)
    {
        fprintf(stderr,"malloc failed\n");
        exit(-1);
    }
    return ptr;
}

struct rbf_t 
{
    double * wts;
    const double * rpts;
    const double * vals;
    double (*phi)(double,double);
    double r0;
    size_t npts,ndims;
    int norm;
};

typedef struct rbf_t rbf_t;

rbf_t * rbf_init(const double *r, const double *vals,
    double r0, double (*phi)(double,double),
    size_t npts, size_t ndims, int normalized) 
{

    size_t i,j,k;
    const double *ri;
    const double *rj;
    double dx,r2,wsum;
    double phi_;
    double *A;
    int * pivs;
    rbf_t * rbf = (rbf_t*)Malloc(sizeof(rbf_t));
    int e;
        
    rbf->norm = normalized;
    rbf->wts = (double*)Malloc(sizeof(double)*npts);
    rbf->phi = phi;
    // shallow copy that perhaps should be a deep copy
    rbf->rpts = r;
    rbf->vals = vals;
    rbf->npts = npts;
    rbf->ndims = ndims;
    rbf->r0 = r0;
    A = (double*)Malloc(sizeof(double)*npts*npts);
    pivs = (int*)Malloc(sizeof(int)*npts);    
    for (i=0;i<npts;++i)
    {
        wsum = 0.0;
        ri = r+i*ndims;
        for (j=0;j<npts;++j)
        {
            rj = r + j*ndims;
            r2 = 0.0;
            for (k=0;k<ndims;++k)
            {
                dx = ri[k] - rj[k];
                r2 += dx * dx;
            }
            phi_ = (*phi)(sqrt(r2),r0);
            A[i*npts+j] = phi_;  
            wsum += phi_;
        }
        if (normalized) {         
            (rbf->wts)[i] = vals[i] * wsum;
        }else{
            (rbf->wts)[i] = vals[i];        
        }
        pivs[i] = i;
    }
    // solve A * wts = y
    // e = LAPACKE_dgesv(LAPACK_ROW_MAJOR,npts,npts,A,npts,
    //    pivs,rbf->wts,1);
    e = LAPACKE_dgels(LAPACK_ROW_MAJOR,'N',npts,npts,1,A,npts,rbf->wts,1);
    if (e != 0)
    {
        fprintf(stderr,"linear system is not solved!\n");
        exit(-1);
    }
    free(pivs);
    free(A);
    fprintf(stderr,"initialized rbf\n");
    return rbf;
};

void rbf_free(rbf_t *rbf)
{
    free(rbf->wts);
    rbf->phi = 0x0;
    free(rbf);
}

double rbf_interp(rbf_t *rbf,double *rp)
{
    size_t i,j;
    double r2,dx,v;
    double wsum = 0.0;
    double psum = 0.0;
    const double * r;
    const double * wt = rbf->wts;
    const double * ri = rbf->rpts;
    const double * yi = rbf->vals;
    const size_t str = rbf->ndims;
    double (*phi)(double,double) = rbf->phi;
    for (i=0;i<rbf->npts;++i)
    {
        r2 = 0.0;
        r = ri + i * str;
        for (j=0;j<rbf->ndims;++j)
        {
            dx = rp[j] - r[j];
            r2 += dx * dx;
        }
        v = (*phi)(sqrt(r2),rbf->r0);
//        v = 1./sqrt(r2+rbf->r0);
        wsum += wt[i];
        psum += wt[i] * v;
    }
    return rbf->norm ? (psum/wsum):psum;
}

double rbf_inverse_quadratic(const double r,const double r0)
{
    double t = r * r + r0;
    return sqrt(1./t);
}

double rbf_phi_gaussian(const double r,const double r0)
{
    double t = r/r0;
    return exp(-0.5*t*t);
}

void interp(
    double * r,
    const double * y,
    size_t ndim,
    size_t npts,
    double * dx,
    double * x0,
    size_t * nn)
{
    rbf_t * rbf;
    FILE * fp;
    size_t k,j,ninterp;
    double vp;
    double r0 = 1.0;
    size_t * indx = (size_t *)malloc(sizeof(size_t)*ndim);
    double * rp = (double *)malloc(sizeof(size_t)*ndim);
    ninterp = 1;
    for (k=0;k<ndim;++k) ninterp *= nn[k];
    for (k=0;k<ndim;++k) indx[k] = 0;
    fprintf(stderr,"# of interpolated points = %lu\n",ninterp);
    
    rbf = rbf_init((const double*)r,(const double*)y,r0,&rbf_inverse_quadratic,
        npts,ndim,0);

    fp = fopen("interp.dat","w");
    for (k=0;k<ninterp;++k)
    {
        for (j=0;j<ndim;++j)
        {
            rp[j] = x0[j] + indx[j] * dx[j];
            fprintf(stderr," %lg",rp[j]);
        }
        vp = rbf_interp(rbf,rp);
        for (j=0;j<ndim;++j)
        {
            fprintf(fp," %15.10lf,",rp[j]);
        } 
        fprintf(fp," %15.10lf\n",vp);
        fprintf(stderr," %lg\n",vp);
        for (j=ndim;j;)
        {
            --j;
            indx[j] += 1;
            if ( indx[j] == nn[j]) {
                fprintf(fp,"\n");
                indx[j] = 0;
            }else{
                break;
            }
        }
    }
    fclose(fp);
    rbf_free(rbf);
    free(rp);
    free(indx);
}

void interp2d(
    double * r,
    const double * y,
    size_t npts,
    double * dx,
    double * x0,
    size_t * nn)
{
    rbf_t * rbf;
    FILE * fp;
    FILE * fph;
    size_t k,j,ninterp;
    double vp;
    size_t ndim = 2;
    double r0 = 1.0;
    size_t * indx = (size_t *)malloc(sizeof(size_t)*ndim);
    double * rp = (double *)malloc(sizeof(size_t)*ndim);
    ninterp = 1;
    for (k=0;k<ndim;++k) ninterp *= nn[k];
    for (k=0;k<ndim;++k) indx[k] = 0;
    fprintf(stderr,"# of interpolated points = %lu\n",ninterp);
    
    rbf = rbf_init((const double*)r,(const double*)y,r0,&rbf_inverse_quadratic,
        npts,ndim,0);

    fp = fopen("interp.dat","w");
    fph = fopen("interp_heat.dat","w");
    for (k=0;k<ninterp;++k)
    {
        for (j=0;j<ndim;++j)
        {
            rp[j] = x0[j] + indx[j] * dx[j];
            fprintf(stderr," %lg",rp[j]);
        }
        vp = rbf_interp(rbf,rp);
        for (j=0;j<ndim;++j)
        {
            fprintf(fp," %15.10lf,",rp[j]);
            fprintf(fph," %15.10lf",rp[j]);
        } 
        fprintf(fp," %15.10lf\n",vp);
        fprintf(stderr," %lg\n",vp);
        fprintf(fph," %15.10lf\n",vp);
        for (j=ndim;j;)
        {
            --j;
            indx[j] += 1;
            if ( indx[j] == nn[j]) {
                fprintf(fph,"\n");
                indx[j] = 0;
            }else{
                break;
            }
        }
    }
    rbf_free(rbf);
    fclose(fp);
    fclose(fph);
    free(rp);
    free(indx);
}

int main(int argc, char **argv)
{
    size_t i,j;
    double *r;
    double *y;
    size_t npts,ndim;
    size_t * nn;
    double * x0;
    double * dx;
    double * xend;
    size_t str;
    size_t phi_pow = 3;
    int result;
    
    FILE * fp = fopen("data.txt","r");
    fscanf(fp,"%lu %lu",&npts,&ndim);
    fprintf(stderr,"# pts = %lu\n",npts);
    fprintf(stderr,"# dims = %lu\n",ndim);
    if (npts == 0)
    {
        fprintf(stderr,"# points = 0\n");
        exit(-1);
    }
    if (ndim == 0)
    {
        fprintf(stderr,"# dims = 0\n");
        exit(-1);
    }
    y = (double *)malloc(sizeof(double)*npts);
    r = (double *)malloc(sizeof(void*)*npts*ndim);
    nn = (size_t*)malloc(sizeof(size_t)*ndim);
    x0 = (double *)malloc(sizeof(double)*ndim);
    dx = (double *)malloc(sizeof(double)*ndim);    
    xend = (double *)malloc(sizeof(double)*ndim);
    str = ndim;
    for (i=0;i<npts;++i)
    {
        for (j=0;j<ndim;++j) 
        {
            fscanf(fp,"%lg",r+str*i+j);
        }
        fscanf(fp,"%lg",y+i);
    }
    for (j=0;j<ndim;++j)
    {
        fscanf(fp,"%lu",&(nn[j]));
    }
    for (j=0;j<ndim;++j)
    {
        fscanf(fp,"%lg",&(x0[j]));
    }
    for (j=0;j<ndim;++j)
    {
        fscanf(fp,"%lg",&(xend[j]));
    }
    for (j=0;j<ndim;++j)
    {
        dx[j] = (xend[j] - x0[j])/(nn[j]-1);
    }    
    fclose(fp);
    for (j=0;j<ndim;++j) {
        fprintf(stderr,"x0[%lu] x1[%lu] dx[%lu] = %lg %lg %lg\n",j,j,j,x0[j],xend[j],dx[j]);
    }
    free(xend);
    
    if (ndim == 2) {
        interp2d(r,y,npts,dx,x0,nn);
    }else{
        interp(r,y,ndim,npts,dx,x0,nn);    
    }
    free(dx);
    free(x0);
    free(nn); 
    free(r);
    free(y);
    return EXIT_SUCCESS;
}
