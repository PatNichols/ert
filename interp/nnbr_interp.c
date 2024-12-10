#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>

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

double nn_interp(
    double * rp, /** point to interpolated */
    double ** ri, /** sample points **/
    const double *vi, /** sample point value */
    double * dist,
    size_t * index,
    const size_t ndim,
    const size_t npts, /** # of sample points */
    size_t num_nbrs)
{
    size_t i,j,ip;
    double dx,dr2,w;
    double p;
    double psum = 0.0;

    for (i=0;i<npts;++i)
    {
        dr2 = 0.0;
        for (j=0;j<ndim;++j) {
            dx = rp[j] - ri[i][j];
            dr2 += dx * dx;
        }
        if ( dr2 < 1.e-15) return vi[i];
        dist[i] = dr2;
        index[i] = i;
    }
//    fprintf(stderr,"\n  ^^^ # pts =  %lu\n",npts);
//    fprintf(stderr,"\n  ^^^ # nbrs =  %lu\n",num_nbrs);
//    for (i=0;i<npts;++i)
//    {
//        fprintf(stderr,"\n *** %lu %lu %lg\n",i,index[i],dist[i]);    
//    }
    for (j=1; j<npts; ++j)
    {
        p=dist[j];
        ip =  index[j];
        i=j-1;
        while (i>=0 && dist[i]>p)
        {
            index[i+1] = index[i];
            dist[i+1]=dist[i];
            --i;
        }
        dist[i+1]=p;
        index[i+1] =  ip;
    }
//    for (i=0;i<npts;++i)
//    {
//        fprintf(stderr,"\n +++ %lu %lu %lg\n",i,index[i],dist[i]);    
//    }
    /** average over n points */
    psum = 0.0;
    for (i=0;i<num_nbrs;++i)
    {
        j = index[i];
        psum += vi[j];
//        fprintf(stderr,"\n  --- %lu %lu %lg\n",i,j,vi[j]);
    }
//    fprintf(stderr," sum = %lg\n",psum);
    return psum/num_nbrs;
}

void interp(
    double ** r,
    const double * y,
    size_t ndim,
    size_t npts,
    double * dx,
    double * x0,
    size_t * nn,
    size_t num_nbrs)
{
    FILE * fp;
    size_t k,j,ninterp;
    double vp;
    size_t * index;
    double * dist;
    size_t * aindx = (size_t *)Malloc(sizeof(size_t)*ndim);
    double * rp = (double *)Malloc(sizeof(size_t)*ndim);
    
    
    ninterp = 1;
    for (k=0;k<ndim;++k) ninterp *= nn[k];
    for (k=0;k<ndim;++k) aindx[k] = 0;
    fprintf(stderr,"# of interpolated points = %lu\n",ninterp);
    if ( npts < num_nbrs) num_nbrs = npts;
    fprintf(stderr,"# of neighbors = %lu\n",num_nbrs);
    dist = (double*)Malloc(npts*sizeof(double));
    index = (size_t*)Malloc(npts * sizeof(size_t));
    fp = fopen("interp.dat","w");
    for (k=0;k<ninterp;++k)
    {
        for (j=0;j<ndim;++j)
        {
            rp[j] = x0[j] + aindx[j] * dx[j];
//            fprintf(stderr," %lg",rp[j]);
        }
        vp = nn_interp(rp,r,y,dist,index,ndim,npts,num_nbrs);
        for (j=0;j<ndim;++j)
        {
            fprintf(fp," %15.10lf,",rp[j]);
        } 
        fprintf(fp," %15.10lf\n",vp);
//        fprintf(stderr," %lg\n",vp);
        for (j=ndim;j;)
        {
            --j;
            aindx[j] += 1;
            if ( aindx[j] == nn[j]) {
                aindx[j] = 0;
            }else{
                break;
            }
        }
    }
    fclose(fp);
    free(index);
    free(dist);
    free(rp);
    free(aindx);
}

void interp2d(
    double ** r,
    const double * y,
    size_t npts,
    double * dx,
    double * x0,
    size_t * nn,
    size_t num_nbrs)
{
    size_t ndim = 2;
    FILE * fph;
    FILE * fp;
    size_t k,j,ninterp;
    double vp;
    size_t * index;
    double * dist;
    size_t * aindx = (size_t *)Malloc(sizeof(size_t)*ndim);
    double * rp = (double *)Malloc(sizeof(size_t)*ndim);
    
    
    ninterp = 1;
    for (k=0;k<ndim;++k) ninterp *= nn[k];
    for (k=0;k<ndim;++k) aindx[k] = 0;
    fprintf(stderr,"# of interpolated points = %lu\n",ninterp);
    if ( npts < num_nbrs) num_nbrs = npts;
    fprintf(stderr,"# of neighbors = %lu\n",num_nbrs);
    dist = (double*)Malloc(npts*sizeof(double));
    index = (size_t*)Malloc(npts * sizeof(size_t));
    fp = fopen("interp.dat","w");
    fph = fopen("interp_heat.dat","w");
    for (k=0;k<ninterp;++k)
    {
        for (j=0;j<ndim;++j)
        {
            rp[j] = x0[j] + aindx[j] * dx[j];
//            fprintf(stderr," %lg",rp[j]);
        }
        vp = nn_interp(rp,r,y,dist,index,ndim,npts,num_nbrs);
        for (j=0;j<ndim;++j)
        {
            fprintf(fp," %15.10lf,",rp[j]);
            fprintf(fph," %15.10lf,",rp[j]);
        } 
        fprintf(fp," %15.10lf\n",vp);
        fprintf(fph," %15.10lf\n",vp);
        for (j=ndim;j;)
        {
            --j;
            aindx[j] += 1;
            if ( aindx[j] == nn[j]) {
                aindx[j] = 0;
                fprintf(fph,"\n");
            }else{
                break;
            }
        }
    }
    fclose(fp);
    fclose(fph);
    free(index);
    free(dist);
    free(rp);
    free(aindx);
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
    size_t num_nbrs = 1;
    char *end;
    int result;
    int norm = 1;
    
    if (argc > 1) {
        num_nbrs = strtoul(argv[1],&end,10);
    }
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
//        fprintf(fp," %lu",nn[j]);
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
    fprintf(stderr,"\n");
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
    if ( ndim == 2) {
        fprintf(stderr,"2 dim interp\n");
        interp2d(r,y,npts,dx,x0,nn,num_nbrs);
    }else{
        interp(r,y,ndim,npts,dx,x0,nn,num_nbrs);
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
