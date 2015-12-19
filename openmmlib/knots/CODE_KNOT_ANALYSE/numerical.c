#include <math.h>
#include <stdio.h>

#define TINY 1.0e-200;

/* help fcts */
                                                                                                                             
void nrerror(char *error_text)
  {
          void exit();
                                                                                                                             
          fprintf(stderr,"Numerical Recipes run-time error...\n");
          fprintf(stderr,"%s\n",error_text);
          fprintf(stderr,"...now exiting to system...\n");
          exit(1);
  }
                                                                                                                             
double *vector(int nl,int nh)
  {
          double *v;
                                                                                                                             
          v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
          if (!v) nrerror("allocation failure in vector()");
          return v-nl;
  }
                                                                                                                             
void free_vector(double *v, int nl, int nh)
  {
          free((char*) (v+nl));
  }
                                        

double **dmatrix(int nrl,int nrh,int ncl,int nch)
{
        int i;
        double **m;
                                                                                                                             
        m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
        if (!m) nrerror("allocation failure 1 in dmatrix()");
        m -= nrl;
                                                                                                                             
        for(i=nrl;i<=nrh;i++) {
                m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
                if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
                m[i] -= ncl;
        }
        return m;
}

void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch)
{
        int i;
                                                                                                                             
        for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
        free((char*) (m+nrl));
}





void ludcmp(double **a,int n,int *indx,double *d)
{ 
        int i,imax,j,k;
        double big,dum,sum,temp;
        double *vv;
        double *vector();
        void nrerror(),free_vector();
 
        vv=vector(1,n);
        *d=1.0;
        for (i=1;i<=n;i++) {
                big=0.0;
                for (j=1;j<=n;j++)
                        if ((temp=fabs(a[i][j])) > big) big=temp;
                if (big == 0.0) nrerror("Singular matrix in routine LUDCMP");
                vv[i]=1.0/big;
        }
        for (j=1;j<=n;j++) {
                for (i=1;i<j;i++) {
                        sum=a[i][j];
                        for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
                        a[i][j]=sum;
                }
                big=0.0;
                for (i=j;i<=n;i++) {
                        sum=a[i][j];
                        for (k=1;k<j;k++)
                                sum -= a[i][k]*a[k][j];
                        a[i][j]=sum;
                        if ( (dum=vv[i]*fabs(sum)) >= big) {
                                big=dum;
                                imax=i;
                        }
                }
                if (j != imax) {
                        for (k=1;k<=n;k++) {
                                dum=a[imax][k];
                                a[imax][k]=a[j][k];
                                a[j][k]=dum;
                        }
                        *d = -(*d);
                        vv[imax]=vv[j];
                }
                indx[j]=imax;
                if (a[j][j] == 0.0) a[j][j]=TINY;
                if (j != n) {
                        dum=1.0/(a[j][j]);
                        for (i=j+1;i<=n;i++) a[i][j] *= dum;
                }
        }
        free_vector(vv,1,n);
}

