#include "opca_binary.h"

#ifdef __cplusplus
extern "C"{
#endif

#include "f2c.h"
#include "blaswrap.h"
#include "clapack.h"

#ifdef __cplusplus
}
#endif

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <cmath>

using namespace::std;

OPCA_BINARY::OPCA_BINARY(int l, double gamma, int Nc, int Nr):
    OPCA(l , gamma, Nc, Nr)
{

}

OPCA_BINARY::~OPCA_BINARY()
{

}


bool OPCA_BINARY::GenerateOPCARealization(double *xi,double *m, double *S){    
	double *a;
    a = new double[Nc_];
	memset(a,0,Nc_*sizeof(double));
    double xS = 0.0;
    double *mu; double *eta;
    mu = new double[Nc_];
	memset(mu,0,Nc_*sizeof(double));
    eta = new double[Nc_];
	memset(eta,0,Nc_*sizeof(double));
    double da_i_dxi_j = 0.0;
	integer M;
	integer N;
	integer K;
	M = Nc_;
	K = l_;
	N = 1;
	doublereal alpha = 1.0;
	doublereal beta = 0.0;
	dgemm_("N","N",&M,&N,&K,&alpha,U_Sig_,&M,xi,&K,&beta,a,&M);
//	OpcaDebug("USig.txt",U_Sig_,Nc_ * l_);	
//a = xm + a; (Nc_,1)
	integer icr = 1;
	beta = 1.0;
	daxpy_(&M, &alpha, xm_, &icr, a, &icr);  
	int index1 = 0;
	int index2 = 0;
    for (int i = 0; i < Nc_ ; i++){
        xS = (a[i] - gamma_/2) / (1 - gamma_);
        if (xS < 1 && xS > 0){
            m[i] = xS;
            mu[i] = 0;
            eta[i] = 0;
        }else{
            if(xS >= 1){
                m[i] = 1.0;
                mu[i] = 0;
                eta[i] = 2 * a[i] + gamma_ - 2;
            }else{
                m[i] = 0.0;
                mu[i] = gamma_ - 2 * a[i];
                eta[i] = 0;
            }
        }
        for(int j=0;j < l_; j++){
			index1 = j * Nc_ + i;
			index2 = j * l_ + j;
            da_i_dxi_j = Sig_[index2] * U_[index1];
			// index of the ith row and jth column element in a column wise matrix data
			// S (Nc_ , l_)
            S[index1] = (2*da_i_dxi_j*m[i]*(1 - m[i])) /
                    (mu[i] + 2*m[i] + eta[i]*m[i] - 2*gamma_*m[i] - mu[i]*m[i]
                     + 2*gamma_*m[i]*m[i] - 2*m[i]*m[i]);
        }
    }
//	OpcaDebug("S1.txt",S, Nc_ * l_);
	delete []a;
	delete []mu;
	delete []eta;
    return true;
}

bool OPCA_BINARY::GenerateOPCARealization(double *xi,double *m){
	double *a;
    a = new double[Nc_];
	memset(a,0,Nc_*sizeof(double));
    double xS = 0.0;
	integer M;
	integer N;
	integer K;
	M = Nc_;
	K = l_;
	N = 1;
	doublereal alpha = 1.0;
	doublereal beta = 0.0;
	dgemm_("N","N",&M,&N,&K,&alpha,U_Sig_,&M,xi,&K,&beta,a,&M);
    	//a = U_Sig_ * xi; (Nc_,l_) * (l_,1);
	integer icr = 1;
	daxpy_(&M, &alpha, xm_, &icr, a, &icr);  
		//a = xm + a; (Nc_,1)
    for (int i = 0; i < Nc_ ; i++){
        xS = (a[i] - gamma_/2) / (1 - gamma_);
        if (xS < 1 && xS > 0){
            m[i] = xS;
        }else{
            if(xS >= 1){
                m[i] = 1.0;
            }else{
                m[i] = 0.0;
            }
        }
    }
    return true;
}
