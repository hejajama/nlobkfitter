#ifndef _NLODIS_SIGMAR_MASSIVE_H
#define _NLODIS_SIGMAR_MASSIVE_H


int integrand_ILdip_massive_LiLogConst(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILdip_massive_Iab(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILdip_massive_Icd(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;

#endif // _NLODIS_SIGMAR_MASSIVE_H
