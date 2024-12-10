#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>

FILE * Fopen(const char *name, const char *mode)
{
    FILE *fp = fopen(name,mode);
    if (!fp) {
        fprintf(stderr,"could not open %s in mode %s\n",
            name,mode);
        exit(EXIT_FAILURE);    
    }
    return fp;
}

void * Malloc(size_t nbytes)
{
    void * ptr = malloc(nbytes);
    if (!ptr) {
        fprintf(stderr,"malloc failed\n");
        exit(-1);
    }
    return ptr;
}

double shep_phi(const double r,const size_t p) { 
    return pow((1.0/r),p);
}

double shep_interp(
    double * rp, /** point to interpolated */
    double ** ri, /** sample points **/
    const double *vi, /** sample point value */
    const size_t npts,
    const size_t ndim, /** # of sample points */
    const size_t phi_pow,
    int norm)
{
    size_t i,j;
    double dx,dr2,w;
    double wsum = 0.0;
    double psum = 0.0;
    
    for (i=0;i<npts;++i)
    {
        dr2 = 0.0;
        for (j=0;j<ndim;++j) {
            dx = rp[j] - ri[i][j];
            dr2 += dx * dx;
        }
        if ( dr2 < 1.e-15) return vi[i];
        // Shephard's phi function 1/pow(r,phi_pow)
        dr2 = 1./sqrt(dr2);
        w = pow(dr2,phi_pow);
        wsum += w;
        psum += vi[i] * w;
    }
    if ( !norm ) return psum;
    return psum/wsum;
}

void interp2d(
    double ** r,
    const double * y,
    size_t npts,
    double * dx,
    double * x0,
    size_t * nn,
    size_t phi_pow,
    int norm)
{
    size_t ndim = 2;
    FILE * fp = 0x0;
    FILE * fph = 0x0;
    size_t k,j,ninterp;
    double vp;
    size_t * indx = (size_t *)Malloc(sizeof(size_t)*ndim);
    double * rp = (double *)Malloc(sizeof(size_t)*ndim);
    ninterp = 1;
    for (k=0;k<ndim;++k) ninterp *= nn[k];
    for (k=0;k<ndim;++k) indx[k] = 0;
    fprintf(stderr,"# of interpolated points = %lu\n",ninterp);
    fp = Fopen("interp.dat","w");
    fph = Fopen("interp_heat.dat","w");
    for (k=0;k<ninterp;++k)
    {
        for (j=0;j<ndim;++j)
        {
            rp[j] = x0[j] + indx[j] * dx[j];
        }
        vp = shep_interp(rp,r,y,npts,ndim,phi_pow,norm);
        for (j=0;j<ndim;++j)
        {
            fprintf(fp," %15.10lf",rp[j]);
            fprintf(fph," %15.10lf",rp[j]);
        } 
        fprintf(fp," %15.10lf\n",vp);
        fprintf(fph," %15.10lf\n",vp);
        for (j=ndim;j;)
        {
            --j;
            indx[j] += 1;
            if ( indx[j] == nn[j]) {
                indx[j] = 0;
                fprintf(fph,"\n");
            }else{
                break;
            }
        }
    }
    fclose(fph);
    fclose(fp);
    free(rp);
    free(indx);
}


void interp(
    double ** r,
    const double * y,
    size_t ndim,
    size_t npts,
    double * dx,
    double * x0,
    size_t * nn,
    size_t phi_pow,
    int norm)
{
    FILE * fp = 0x0;
    size_t k,j,ninterp;
    double vp;
    size_t * indx = (size_t *)Malloc(sizeof(size_t)*ndim);
    double * rp = (double *)Malloc(sizeof(size_t)*ndim);
    ninterp = 1;
    for (k=0;k<ndim;++k) ninterp *= nn[k];
    for (k=0;k<ndim;++k) indx[k] = 0;
    fprintf(stderr,"# of interpolated points = %lu\n",ninterp);
    fp = Fopen("interp.dat","w");
    for (k=0;k<ninterp;++k)
    {
        for (j=0;j<ndim;++j)
        {
            rp[j] = x0[j] + indx[j] * dx[j];
            fprintf(stderr," %lg",rp[j]);
        }
        vp = shep_interp(rp,r,y,npts,ndim,phi_pow,norm);
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
                indx[j] = 0;
            }else{
                break;
            }
        }
    }
    fclose(fp);
    free(rp);
    free(indx);
}


int main(int argc, char **argv)
{
    size_t i,j;
    double **r;
    double *y;
    size_t npts,ndim;
    size_t * nn;
    double * x0;
    double * dx;
    double * xend;
    size_t phi_pow = 3;
    char *end;
    int result;
    int norm = 1;
    
    if (argc > 1) {
        phi_pow = strtoul(argv[1],&end,10);
        if (argc > 2) norm = atoi(argv[2]);
    }
    FILE * fp = Fopen("data.txt","r");
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
    y = (double *)Malloc(sizeof(double)*npts);
    r = (double**)Malloc(sizeof(void*)*npts);
    for (i=0;i<npts;++i) r[i] = (double*)Malloc(sizeof(double)*ndim);
    nn = (size_t*)Malloc(sizeof(size_t)*ndim);
    x0 = (double *)Malloc(sizeof(double)*ndim);
    dx = (double *)Malloc(sizeof(double)*ndim);    
    xend = (double *)Malloc(sizeof(double)*ndim);
    for (i=0;i<npts;++i)
    {
        for (j=0;j<ndim;++j) 
        {
            fscanf(fp,"%lg",r[i]+j);
        }
        fscanf(fp,"%lg",y+i);
    }
    for (j=0;j<ndim;++j)
    {
        fscanf(fp,"%lu",&(nn[j]));
    }
    fprintf(stderr,"\n");

    for (i=0;i<npts;++i)
    {
        for (j=0;j<ndim;++j) 
        {
            fprintf(stderr," %lg",r[i][j]);
        }
        fprintf(stderr," %lg\n",y[i]);
    }
    for (j=0;j<ndim;++j)
    {
        fscanf(fp,"%lg",&(x0[j]));
    }
    for (j=0;j<ndim;++j)
    {
        fscanf(fp,"%lg",&(xend[j]));
    }
    fprintf(stderr,"\n");
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
        interp2d(r,y,npts,dx,x0,nn,phi_pow,norm);
    }else{
        interp(r,y,ndim,npts,dx,x0,nn,phi_pow,norm);
    }
    free(dx);
    free(x0);
    free(nn); 
    for (i=npts;i;)
    {
        --i;
        free(r[i]);
    }
    free(r);
    free(y);
    return EXIT_SUCCESS;
}
