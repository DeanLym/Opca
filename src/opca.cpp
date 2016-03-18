#include "opca.h"

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
#include <fstream>
#include <algorithm>


using namespace::std;

OPCA::OPCA(int l, double gamma, int Nc, int Nr){
    set_l(l);
    set_gamma(gamma);
	set_Nc(Nc);
	set_Nr(Nr);

	input_U_and_Sig_flag_ = 0;	

	X_ = new double[Nc_ * Nr_];
	xm_ = new double[Nc_];
	U_ = new double[Nc_ * l_];
	Sig_ = new double[l_ * l_];
	U_Sig_ = new double[Nc_ * l_];
	xi_ = new double[l_];

//	m_ = new double[Nc_];
//	S_ = new double[Nc_ * l_];

	memset(X_,0,Nc_ * Nr_ * sizeof(double));

	memset(xm_,0,Nc_ * sizeof(double));

	memset(U_,0,Nc_ * l_ * sizeof(double));
	memset(Sig_,0,l_ * l_ * sizeof(double));
	memset(U_Sig_,0,Nc_ * l_ * sizeof(double));
	memset(xi_,0,l_ * sizeof(double));
//	memset(m_,0,Nc_ * sizeof(double));	
//	memset(S_,0,Nc_ * l_ * sizeof(double));

}

OPCA::~OPCA(){
	delete []X_;
	delete []xm_;
	delete []U_;
	delete []Sig_;
	delete []U_Sig_;
	delete []xi_;
//	delete []m_;
//	delete []S_;
}
bool OPCA::set_Nc(int Nc){
    Nc_ = Nc;
    return true;
}

bool OPCA::set_Nr(int Nr){
    Nr_ = Nr;
    return true;
}

bool OPCA::set_gamma(double gamma){
    gamma_ = gamma;
    return true;
}

bool OPCA::set_l(int l){
    l_ = l;
    return true;
}


bool OPCA::InputDataMatrix(int Nc, int Nr ,double *X){
	integer n = Nc_ * Nr_; 
	integer icr = 1;
	dcopy_(&n, X , &icr, X_ , &icr); 
    	// Copy X to X_;
    return true;
}

bool OPCA::InputDataMatrix(const char* file_name){
    ifstream in;
    in.open(file_name);
	// Data in file should be transpose of X;
	int index=0;
	for (int i = 0; i < Nc_ ; i++){
		for (int j = 0; j < Nr_ ; j++){
			index = j*Nc_ + i;
			in >> X_[index];
		}
	}
	in.close();
//	OpcaDebug("X.dat",X_,Nc_*Nr_);
    return true;
}

bool OPCA::ConstructYAndPerformSVD(){
	// xm = mean(X_);
	double *ones;
	ones = new double[Nr_];
	for (int i=0;i<Nr_;i++){
		ones[i] = 1.0;
	}
	integer M;
	integer N;
	integer K;
	M = Nc_;
	K = Nr_;
	N = 1;
	doublereal alpha = 1.0;
	doublereal beta = 0.0;
	dgemm_("N","N",&M,&N,&K,&alpha,X_,&M,ones,&K,&beta,xm_,&M);
		// xm_ = sum(X_);
	integer icr = 1;
	alpha = 1.0 / (doublereal)(Nr_) - 1.0;
	daxpy_(&M, &alpha, xm_, &icr, xm_, &icr);  
		// xm_ = xm_ / Nr_;
    // X_ = X_ - (xm_,xm_,...,xm_);
	double *XM;
	XM = new double[Nc_ * Nr_];
	memset(XM, 0 , Nc_ * Nr_ * sizeof(double));
	M = Nc_;
	K = 1;
	N = Nr_;
	alpha = 1.0;
	dgemm_("N","T",&M,&N,&K,&alpha,xm_,&M,ones,&N,&beta,XM,&M);	
		// xm -> XM(xm,xm,xm,...,xm);
	delete []ones;  
	alpha = -1.0;
	M = Nc_ * Nr_;
	daxpy_(&M,&alpha, XM, &icr, X_, &icr);
		// X_ = -XM + X_;
	alpha = 1.0 / sqrt((doublereal)Nr_ - 1.0) - 1.0;
	daxpy_(&M, &alpha, X_, &icr, X_, &icr);
		// X_ = X_/sqrt(Nr - 1) -> X_;
    // [U_, SIG, V] = svd(X_);

	M = Nc_;
	N = Nr_;
	integer NCOL;
	NCOL = min(Nc_,Nr_);
	double *U, *Sig,*VT;
	VT = NULL;
	U = new double[Nc_ * NCOL];
	Sig = new double[NCOL];
	integer info = -1;
	doublereal work1;
	integer lwork = -1;
	dgesvd_("S","N",&M,&N,X_,&M,Sig,U,&M,VT,&N,&work1,&lwork,&info);
	if (info != 0)
		cout << "SVD failed" << endl;
	lwork = (integer)(work1);
	doublereal *work;	
	work = new double[lwork];
	dgesvd_("S","N",&M,&N,X_,&M,Sig,U,&M,VT,&N,work,&lwork,&info);
	

	if (info != 0)
		cout << "SVD failed" << endl;
	int index;
	for ( int i = 0; i < l_; i++ ){
		index = i*l_ + i;		
		Sig_[index] = Sig[i];
	}
	integer n = Nc_ * l_; 
	dcopy_(&n, U , &icr, U_ , &icr); 
	
    	// Copy the first l column of U to U_;
	M = Nc_;
	N = l_;
	K = l_;
	alpha = 1.0;
	dgemm_("N","N",&M,&N,&K,&alpha,U_,&M,Sig_,&K,&beta,U_Sig_,&M);	
	delete []U;
	delete []Sig;
	delete []work;
    return true;
}

bool OPCA::InputUSig(double *U_Sig){
    // Copy U_Sig to U_Sig_
	integer n = Nc_ * l_;
	integer icr = 1;
	dcopy_(&n, U_Sig, &icr, U_Sig_ ,&icr);
    return true;
}


bool OPCA::InputUSig(const char* file_name){
    ifstream in;
    in.open(file_name);
	int index=0;
	for (int i = 0; i < Nc_ ; i++){
		for (int j = 0; j < l_ ; j++){
			index = j*Nc_ + i;
			in >> U_Sig_[index];
		}
	}
    in.close();
	return true;
}



bool OPCA::InputXm(double *xm){
    // Copy xm to xm_
	integer n = Nc_;
	integer icr = 1;
	dcopy_(&n, xm, &icr, xm_ ,&icr);
    return true;
}

bool OPCA::InputXm(const char* file_name){
    ifstream in;
    in.open(file_name);
	for (int i = 0; i < Nc_ ; i++){
		in >> xm_[i];
	}
    in.close();
	return true;
}


int OPCA::get_Nc(){
    return Nc_;
}

int OPCA::get_Nr(){
    return Nr_;
}

double OPCA::get_gamma(){
	    return gamma_;
}

bool OPCA::OutputUSig(const char* file_name){
    ofstream out;
    out.open(file_name);
	int index;
	for (int i = 0; i < Nc_ ; i++){
		for (int j = 0; j < l_ ; j++){
			index = j*Nc_ + i;
			out << U_Sig_[index] << "  ";
		}
		out << endl;
	}
    out.close();
	return true;
}


bool OPCA::OutputU(const char* file_name){
	ofstream out;
	out.open(file_name);
	int index;
	for (int i = 0; i < Nc_ ; i++){
		for (int j = 0; j < l_ ; j++){
			index = j*Nc_ + i;
			out << U_[index] << "  ";
		}
		out << endl;
	}
	out.close();
	return true;
}

bool OPCA::OutputSig(const char* file_name){
	ofstream out;
	out.open(file_name);
	int index;	
	for (int i = 0; i < l_ ; i++){
		index = i*l_ + i;
		out << Sig_[index]<< endl;
	}
	out.close();
	return true;
}

bool OPCA::InputU(const char* file_name){
	ifstream in;
	in.open(file_name);
	int index;
	for (int i = 0; i < Nc_ ; i++){
		for (int j = 0; j < l_ ; j++){
			index = j*Nc_ + i;
			in >> U_[index];
		}
	}
	in.close();
	integer M = Nc_, N = l_, K = l_;
	doublereal alpha = 1.0;
	doublereal beta = 0.0;
	if(input_U_and_Sig_flag_==1){
		dgemm_("N","N",&M,&N,&K,&alpha,U_,&M,Sig_,&K,&beta,U_Sig_,&M);	
	}else{
		input_U_and_Sig_flag_ = 1;
	}
	return true;
}

bool OPCA::InputU(double *U){
    // Copy U_Sig to U_Sig_
	integer n = Nc_ * l_;
	integer icr = 1;
	dcopy_(&n, U, &icr, U_ ,&icr);

	integer M = Nc_, N = l_, K = l_;
	doublereal alpha = 1.0;
	doublereal beta = 0.0;
	if(input_U_and_Sig_flag_==1){
		dgemm_("N","N",&M,&N,&K,&alpha,U_,&M,Sig_,&K,&beta,U_Sig_,&M);	
	}else{
		input_U_and_Sig_flag_ = 1;
	}

    return true;
}

bool OPCA::InputSig(double *Sig){
	int index;
	for (int i = 0; i < l_ ; i++){
		index = i*l_ + i;
		Sig_[index] = Sig[i];
	}
	integer M = Nc_, N = l_, K = l_;
	doublereal alpha = 1.0;
	doublereal beta = 0.0;
	if(input_U_and_Sig_flag_==1){
		dgemm_("N","N",&M,&N,&K,&alpha,U_,&M,Sig_,&K,&beta,U_Sig_,&M);	
	}else{
		input_U_and_Sig_flag_ = 1;
	}
	return true;
}

bool OPCA::InputSig(const char* file_name){
	ifstream in;
	in.open(file_name);
	int index;
	for (int i = 0; i < l_ ; i++){
		index = i*l_ + i;
		in >> Sig_[index];
	}
	in.close();
	integer M = Nc_, N = l_, K = l_;
	doublereal alpha = 1.0;
	doublereal beta = 0.0;
	if(input_U_and_Sig_flag_==1){
		dgemm_("N","N",&M,&N,&K,&alpha,U_,&M,Sig_,&K,&beta,U_Sig_,&M);	
	}else{
		input_U_and_Sig_flag_ = 1;
	}
	return true;
}

bool OPCA::OutputXm(const char* file_name){
    ofstream out;
    out.open(file_name);
	for (int i = 0; i < Nc_ ; i++){
		out << xm_[i] << endl;
	}
    out.close();
	return true;
}

void OPCA::OpcaDebug(char* filename,double *data, int size){
	ofstream out;
	out.open(filename);
	for (int i=0;i<size;i++){
		out << data[i] << endl;
	}
	out.close();
}

bool OPCA::GetTrueXi(double* m_true, double *xi_true){
	// m_true = m_true - xm_;
	integer M = Nc_;
	doublereal alpha = -1.0;
	integer icr = 1.0;
	daxpy_(&M, &alpha, xm_ , &icr, m_true , &icr);
	// xi_true = trans(U)*m_true;
	M = l_;
	integer K = Nc_;
	integer N = 1;
	doublereal beta = 0.0;
	alpha = 1.0;
	dgemm_("T","N",&M,&N,&K,&alpha,U_,&K,m_true,&K,&beta,xi_true,&M);	
	// xi_true = xi_true ./ diag(Sig_);
	int index = 0;
	for(int i = 0 ; i < l_ ; i++){
		index = i*l_ + i;
		xi_true[i] = xi_true[i] / Sig_[index];
	}
	return true;
}

