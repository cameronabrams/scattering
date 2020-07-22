#include "skdataspace.h"
/*
 *  Structure factor calculation in a periodic system
 *  (c) 2014 cameron f abrams
 */
double max3 ( double a, double b, double c ) {
  if (a>b) {
    if (a>c) return a;
    else return c;
  } else {  // b > a
    if (b>c) return b;
    else return c;
  }
}

SkDataSpace * NewSkDataSpace ( int n, int nframes, int nkx, int nky, int nkz ) {
  SkDataSpace * sds = malloc(sizeof(SkDataSpace));
  int i,j;
  sds->nframes=nframes;
  sds->n=n;
  sds->nk[0]=nkx;
  sds->nk[1]=nky;
  sds->nk[2]=nkz;
  sds->b=(double*)malloc(n*sizeof(double));
  for (i=0;i<n;i++) sds->b[i]=1.0;
  sds->nq=0;
  sds->qmax=0.0;
  sds->dq=0.01; // hardcode
  sds->x=(double**)malloc(3*sizeof(double*));
  for (i=0;i<3;i++) sds->x[i]=(double*)malloc(n*sizeof(double));

  sds->co=(double**)malloc(3*sizeof(double*));
  for (i=0;i<3;i++) sds->co[i]=(double*)malloc(max3(sds->nk[0],sds->nk[1],sds->nk[2])*sizeof(double));
  sds->si=(double**)malloc(3*sizeof(double*));
  for (i=0;i<3;i++) sds->si[i]=(double*)malloc(max3(sds->nk[0],sds->nk[1],sds->nk[2])*sizeof(double));

  sds->C=(double***)malloc(nkx*sizeof(double**));
  sds->S=(double***)malloc(nkx*sizeof(double**));
  for (i=0;i<nkx;i++) {
    sds->C[i]=(double**)malloc(nky*sizeof(double*));
    sds->S[i]=(double**)malloc(nky*sizeof(double*));
    for (j=0;j<nky;j++) {
      sds->C[i][j]=(double*)malloc(nkz*sizeof(double));
      sds->S[i][j]=(double*)malloc(nkz*sizeof(double));
    }
  }

  sds->dump=fopen("dump.sk","w");

  return sds;
}

void FreeSkDataSpace ( SkDataSpace * sds ) {
  int i;
  fclose(sds->dump);
  for (i=0;i<3;i++) {
    free(sds->x[i]);
  }
  free(sds->x);
}

double SkDataSpace_setdq ( SkDataSpace * sds, double dq ) {
  if (sds) {
    return sds->dq=dq;
  } else return 0.0;
}

double * SkDataSpace_getCellSizeAddr ( SkDataSpace * sds ) {
  if (sds) {
    return sds->c;
  } else return NULL;
}

double * SkDataSpace_getScatLengthAddr ( SkDataSpace * sds ) {
  if (sds) {
    return sds->b;
  } else return NULL;
}

double * SkDataSpace_getCoordAddr ( SkDataSpace * sds, int d ) {
  if (sds) {
    return sds->x[d];
  } else return NULL;
}

void SkDataSpace_report ( SkDataSpace * sds ) {
  int i;
  int n=sds->n;
  fprintf(stdout,"SkDataSpace: n %i nframes %i nkx nky nkz %i %i %i\n",sds->n,sds->nframes,sds->nk[0],sds->nk[1],sds->nk[2]);
  fprintf(stdout,"             cell dimensions (angstroms): %.5lf %.5lf %.5lf\n",
	  sds->c[0],sds->c[1],sds->c[2]);
  fprintf(stdout,"             first five atom coordinates:\n");
  for (i=0;i<(n>5?5:n);i++) {
    fprintf(stdout,"              %i % .8lf % .8lf % .8lf\n",i,sds->x[0][i],sds->x[1][i],sds->x[2][i]);
  }
  fprintf(stdout,"----------------------------------------------------------\n");
}

int SkDataSpace_setFrame ( SkDataSpace * sds, int f ) {
  if (sds) {
    return sds->frame=f;
  } else return -1;
}

void SkDataSpace_scaleCoords ( SkDataSpace * sds, int d ) {
  if (sds) {
    int i;
    double f=2*M_PI/sds->c[d];
    for (i=0;i<sds->n;i++) {
      sds->x[d][i]*=f;
    }
  }
}

double qmag ( int j, int k, int l, double * c ) {
  double s=(j/c[0])*(j/c[0])+(k/c[1])*(k/c[0])+(l/c[1])*(l/c[0]);
  return 2*M_PI*sqrt(s);
}

void init_SkBins ( SkDataSpace * sds ) {
  int i;
  sds->qmax=qmag(sds->nk[0],sds->nk[1],sds->nk[2],sds->c);
  sds->nq=(int)(sds->qmax/sds->dq);
  sds->Sk=(SkBin*)malloc(sds->nq*sizeof(SkBin));
  for (i=0;i<sds->nq;i++) {
    sds->Sk[i].i=i;
    sds->Sk[i].s=0.0;
    sds->Sk[i].c=0;
  }
}

void SkDataSpace_updateSk ( SkDataSpace * sds ) {
  if (sds) {
    int i,j,k,l,m;
    double ** c = sds->co;
    double ** s = sds->si;
    double * b = sds->b;
    double ** x = sds->x;
    double *** C = sds->C;
    double *** S = sds->S;
    double q,v,cc,cs,sc,ss,ccss,sccs;
    int bin;
    if (sds->frame==0) {
      init_SkBins(sds);
    }
    
    // blank C and S
    for (j=0;j<sds->nk[0];j++) {
      for (k=0;k<sds->nk[1];k++) {
	for (l=0;l<sds->nk[2];l++) {
	  C[j][k][l]=S[j][k][l]=0.0;
	}
      }
    }

    for (i=0;i<sds->n;i++) {
      // fprintf(stderr,"Update_sk at atom %i\n",i);fflush(stderr);
      for (j=0;j<3;j++) {
	//fprintf(stderr,"   dir %i\n",j);fflush(stderr);
	c[j][0]=1.0;  s[j][0]=0.0;
	c[j][1]=cos(x[j][i]);
	s[j][1]=sin(x[j][i]);
	for (k=2;k<sds->nk[j];k++) {
	  // fprintf(stderr,"      k-n %i out of %i\n",k,sds->nk[j]);fflush(stderr);
	  c[j][k]=c[j][k-1]*c[j][1]-s[j][k-1]*s[j][1];
	  s[j][k]=s[j][k-1]*c[j][1]+c[j][k-1]*s[j][1];
	}
      }
      for (k=0;k<sds->nk[0];k++) {
	for (l=0;l<sds->nk[1];l++) {
	  cc=c[0][k]*c[1][l];
	  cs=c[0][k]*s[1][l];
	  sc=s[0][k]*c[1][l];
	  ss=s[0][k]*s[1][l];
	  ccss=cc-ss;
	  sccs=sc+sc;
	  for (m=0;m<sds->nk[2];m++) {
	    C[k][l][m]+=ccss*c[2][m]-sccs*s[2][m];
	    S[k][l][m]+=sccs*c[2][m]+ccss*s[2][m];
	  }
	}
      }// end triple-loop
    } // end atom-loop
    // update Sk
    for (j=0;j<sds->nk[0];j++) {
      for (k=0;k<sds->nk[1];k++) {
	for (l=0;l<sds->nk[2];l++) {
	  v=b[j]*b[j]*(C[j][k][l]*C[j][k][l]+S[j][k][l]*S[j][k][l])/(sds->n);
	  q=qmag(j,k,l,sds->c);
	  bin=(int)(q/sds->dq);
	  sds->Sk[bin].c++;
	  sds->Sk[bin].s+=v;
	}
      }
    }
  }
}

void SkDataSpace_outputSk ( SkDataSpace * sds ) {
  if (sds) {
    int i;
    for (i=0;i<sds->nq;i++) {
      if (sds->Sk[i].c) fprintf(sds->dump,"%.5lf %.5lf\n",i*sds->dq,sds->Sk[i].s/sds->Sk[i].c);
    }
  }
}
