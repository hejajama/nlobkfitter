#ifndef _DDIS_H
#define _DDIS_H

int integrand_ddis_lo_qqbar_L(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ddis_lo_qqbar_T(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ddis_nlo_qqbarg_T_largeM(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ddis_nlo_qqbarg_T_largeQsq(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ddis_nlo_qqbarg_L(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ddis_nlo_qqbarg_T(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
