#ifndef OPCA_BIMODAL_H
#define OPCA_BIMODAL_H

#include "opca.h"
#include <vector>


#define OPCA_PI 3.141592653

class OPCA_BIMODAL : public OPCA
{
public:
    OPCA_BIMODAL(int l, double gamma, int Nc, int Nr);
	~OPCA_BIMODAL();
    bool GenerateOPCARealization(double *xi, double *x, double *S);
    bool GenerateOPCARealization(double *xi,double *x);
public:
	bool set_mu(double mu1,double mu2);
	bool set_var(double var1,double var2);
	bool set_bounds(double lb, double ub);
	static double OPCAObjFuncNLOPT(const vector<double> &x,
		vector<double> &grad, void *f_data);
protected:
	bool BoundInitGuess(vector<double> &x0);
private:
	double mu1_, mu2_;
	double var1_, var2_;
	double lb_,ub_;

};


#endif // OPCA_BIMODAL_H

