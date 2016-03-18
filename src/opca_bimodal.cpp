#include "opca_bimodal.h"
#include "nlopt.hpp"

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
#include <exception>
#include <stdexcept>

using namespace::std;
using namespace::nlopt;

OPCA_BIMODAL::OPCA_BIMODAL(int l, double gamma, int Nc, int Nr):
    OPCA(l , gamma, Nc, Nr)
{
	mu1_ = 0.0;
	mu2_ = 0.0;
	var1_ = 0.0;
	var2_ = 0.0;
}

OPCA_BIMODAL::~OPCA_BIMODAL()
{

}


bool OPCA_BIMODAL::set_mu(double mu1,double mu2){
	mu1_ = mu1;
	mu2_ = mu2;
	return true;
}

bool OPCA_BIMODAL::set_var(double var1,double var2){
	var1_ = var1;
	var2_ = var2;
	return true;
}

bool OPCA_BIMODAL::set_bounds(double lb, double ub){
	lb_ = lb;
	ub_ = ub;
	return true;
}

double OPCA_BIMODAL::OPCAObjFuncNLOPT(const vector<double> &x,
	vector<double> &grad, void *f_data){
		// f_data[0 - 4] are gamma,mu1,mu2,var1,var2;
		// f_data[5 - (n+4)] are a[0-(n-1)];
	int n = x.size();	
	double *data = (double *) f_data;
	double gamma,mu1,mu2,var1,var2;
	gamma = data[0];
	mu1 = data[1];
	mu2 = data[2];
	var1 = data[3];
	var2 = data[4];
	double pi = 3.141592653589793;
	
	double res = 0.0;
	double ai = 0.0;
	double Ri = 0.0;
	double dRi = 0.0;
	for (int i = 0; i < n; i++){
		ai = data[5+i];
		Ri = 1 - 1*exp(-(x[i]-mu1)*(x[i]-mu1)/(2*var1))/sqrt(2*pi*var1) 
				- 1*exp(-(x[i]-mu2)*(x[i]-mu2)/(2*var2))/sqrt(2*pi*var2);
		dRi = 1*exp(-(x[i]-mu1)*(x[i]-mu1)/(2*var1))/sqrt(2*pi*var1)*(x[i] - mu1)/var1
				+ 1*exp(-(x[i]-mu2)*(x[i]-mu2)/(2*var2))/sqrt(2*pi*var2)*(x[i] - mu2)/var2;
		res += (ai - x[i])*(ai - x[i]) + gamma * Ri;
		grad[i] = 2*x[i] - 2*ai + gamma * dRi;
	}
//	cout << res << endl;
	return res;
}


bool OPCA_BIMODAL::GenerateOPCARealization(double *xi,double *m, double *S){    
	double *a;
    a = new double[Nc_];
	memset(a,0,Nc_*sizeof(double));
//	OpcaDebug("Xi.dat",xi,l_);
	integer M;
	integer N;
	integer K;
	M = Nc_;
	K = l_;
	N = 1;
	doublereal alpha = 1.0;
	doublereal beta = 0.0;
	dgemm_("N","N",&M,&N,&K,&alpha,U_Sig_,&M,xi,&K,&beta,a,&M);
	integer icr = 1;
	beta = 1.0;
	daxpy_(&M, &alpha, xm_, &icr, a, &icr);  
//	OpcaDebug("a.dat",a, Nc_);	
//	OpcaDebug("xm.dat",xm_, Nc_);	
	double *f_data;
	f_data = new double[Nc_ + 5];
	f_data[0] = gamma_;
	f_data[1] = mu1_;
	f_data[2] = mu2_;
	f_data[3] = var1_;
	f_data[4] = var2_;

	dcopy_(&M, a , &icr, f_data + 5 , &icr); 

	// Optimize using nlopt
	int n = Nc_;
	vector<double> x;
	x.assign(a , a+Nc_);
//	cout<< x.size() << endl << x[0] <<" "<< a[0] << endl << x[Nc_-1] << " " << a[Nc_-1]<<endl; 
	opt test_opt(nlopt::LD_LBFGS , n);
	test_opt.set_min_objective( OPCA_BIMODAL::OPCAObjFuncNLOPT , f_data);
	test_opt.set_lower_bounds( lb_ );
	test_opt.set_upper_bounds( ub_ );	
	int max_iter = 400;
	test_opt.set_maxeval(max_iter);
	double opt_f = 0.0;
	BoundInitGuess( x );
	test_opt.optimize( x , opt_f);

	for(int i=0; i < n; i++){
		m[i] = x[i];
	}
	double nu = 0.0;
	double eta = 0.0;
	double da_idxi_j = 0.0;
	double dF_1dx_i = 0.0;
	double *LagrangianGrad = new double[Nc_];
	double temp = 0.0;
	double dmu1 = 0.0, dmu2 = 0.0,exp1 = 0.0, exp2 = 0.0, sqr1 = 0.0, sqr2 = 0.0;
	int index1 = 0, index2 = 0;
	for(int i=0; i < Nc_ ; i++){
		dmu1 = m[i] - mu1_;
		dmu2 = m[i] - mu2_;
		exp1 = exp(-dmu1*dmu1 / (2 * var1_));
		exp2 = exp(-dmu2*dmu2 / (2 * var2_));
		sqr1 = sqrt(2 * OPCA_PI * var1_);
		sqr2 = sqrt(2 * OPCA_PI * var2_);
		temp = -2 * a[i] + 2 * m[i] + gamma_ * (dmu1 * exp1 / (sqr1 * var1_)  + dmu2 * exp2 / (sqr2 * var2_) );
    	if (m[i] > lb_ && m[i] < ub_){
        	nu = 0;
        	eta = 0;
		}
    	else{
			if(m[i] == lb_){
        		eta = 0;
        		nu = temp; // with normalizing
			}
			else{
        	nu = 0;
        	eta = -temp; // with normalizing
    		}
		}
    	dF_1dx_i = 2 + gamma_ * (exp1 / (sqr1 * var1_)  - exp1 * dmu1 * dmu1 / (sqr1  * var1_ * var1_)  + exp2 / (sqr2 * var2_)  - exp2 * dmu2 * dmu2 / (sqr2 * var2_ * var2_) ); // with normalizing
    	for (int j = 0; j < l_ ; j++){
			index1 = j * Nc_ + i;
			index2 = j * l_ + j;
            da_idxi_j = Sig_[index2] * U_[index1];
        	S[index1] =  (2 * da_idxi_j * (m[i] - lb_) * (ub_ - m[i])) / 
            (-lb_ * ub_ * dF_1dx_i - lb_ * eta + ub_ * nu + (ub_ + lb_) * dF_1dx_i * m[i] + eta * m[i] - nu * m[i] - dF_1dx_i * m[i] * m[i]); // sensitivity matrix dx/dxi
    	}
    	LagrangianGrad[i] = temp - nu + eta; // with normalizing
	}

	delete []LagrangianGrad;
	delete []f_data;
	delete []a;
    return true;
}

bool OPCA_BIMODAL::GenerateOPCARealization(double *xi,double *m){
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
//	OpcaDebug("Xi.dat",xi,l_);	
//	OpcaDebug("a.dat",a, Nc_);	
//	OpcaDebug("xm.dat",xm_, Nc_);	
//	OpcaDebug("USig.dat",U_Sig_,Nc_*l_);
	double *f_data;
	f_data = new double[Nc_ + 5];
	f_data[0] = gamma_;
	f_data[1] = mu1_;
	f_data[2] = mu2_;
	f_data[3] = var1_;
	f_data[4] = var2_;

	dcopy_(&M, a , &icr, f_data + 5 , &icr); 

	// Optimize using nlopt
	int n = Nc_;
	vector<double> x;
	x.assign(a , a+Nc_);
//	cout<< x.size() << endl << x[0] <<" "<< a[0] << endl << x[Nc_-1] << " " << a[Nc_-1]<<endl; 
	opt test_opt(nlopt::LD_LBFGS , n);
	test_opt.set_min_objective( OPCA_BIMODAL::OPCAObjFuncNLOPT , f_data);
	test_opt.set_lower_bounds( lb_ );
	test_opt.set_upper_bounds( ub_ );
	int max_iter = 400;
	test_opt.set_maxeval(max_iter);
	double opt_f;
//	cout << "Input For NLopt Completed" << endl;
	BoundInitGuess( x );
	try{
		test_opt.optimize( x , opt_f);
	}
	catch(std::runtime_error& e){
		for(int i=0; i < n; i++){
			m[i] = xm_[i];
		}
		return false;
	}
//	copy(x.begin(), x.end(), m);
	for(int i=0; i < n; i++){
		m[i] = x[i];
	}

	delete []f_data;
	delete []a;
    return true;
}


bool OPCA_BIMODAL::BoundInitGuess(vector<double> &x0){
	for(int i=0; i<Nc_ ; i++){
		if(x0[i] > ub_){
			x0[i] = ub_;
		}else if(x0[i] < lb_){
			x0[i] = lb_;
		}
	}
	return true;


}
