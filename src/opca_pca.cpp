#include "opca_pca.h"

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

OPCA_PCA::OPCA_PCA(int l, int Nc, int Nr):
    OPCA(l , 0.0, Nc, Nr)
{

}

OPCA_PCA::~OPCA_PCA()
{

}


bool OPCA_PCA::GenerateOPCARealization(double *xi,double *m, double *S){    
	double dm_i_dxi_j = 0.0;
	memset(m,0,Nc_*sizeof(double));

	integer M;
	integer N;
	integer K;
	M = Nc_;
	K = l_;
	N = 1;
	doublereal alpha = 1.0;
	doublereal beta = 0.0;
	// m = U_Sig*xi	
	dgemm_("N","N",&M,&N,&K,&alpha,U_Sig_,&M,xi,&K,&beta,m,&M);
	
//	OpcaDebug("USig.txt",U_Sig_,Nc_ * l_);	
	//m = xm + 1.0*m; (Nc_,1)
	integer icr = 1;
	beta = 1.0;
	daxpy_(&M, &alpha, xm_, &icr, m, &icr);  

	int index1 = 0;
	int index2 = 0;
    for (int i = 0; i < Nc_ ; i++){
        for(int j=0;j < l_; j++){
			index1 = j * Nc_ + i;
			index2 = j * l_ + j;
            dm_i_dxi_j = Sig_[index2] * U_[index1];
			// index of the ith row and jth column element in a column wise matrix data
			// S (Nc_ , l_)
            S[index1] = dm_i_dxi_j;
        }
    }
//	OpcaDebug("S1.txt",S, Nc_ * l_);
    return true;
}

bool OPCA_PCA::GenerateOPCARealization(double *xi,double *m){
	memset(m,0,Nc_*sizeof(double));
	integer M;
	integer N;
	integer K;
	M = Nc_;
	K = l_;
	N = 1;
	doublereal alpha = 1.0;
	doublereal beta = 0.0;
	// m = U_Sig*xi	
	dgemm_("N","N",&M,&N,&K,&alpha,U_Sig_,&M,xi,&K,&beta,m,&M);
	integer icr = 1;
	beta = 1.0;
	daxpy_(&M, &alpha, xm_, &icr, m, &icr);  
    return true;
}
