#include "../include/nlopt.hpp"
#include <iostream>
#include <vector>


using namespace::std;
using namespace::nlopt;

class LYM{
public:
	static double f(const vector<double> &x, vector<double> &grad, void *f_data);
	bool set_c(double c);
private:
	double c_;
};

bool LYM::set_c(double c){
	c_ = c;
	return true;
}

double LYM::f(const vector<double> &x, vector<double> &grad, void *f_data){
	int n = x.size();
	double *data = (double *) f_data;
	double c = data[0];
	double res = 0.0;	
	for (int i = 0; i < n; i++){
		res += (x[i] - c)*(x[i] - c);
		grad[i] = 2*(x[i] - c);
	}
	return res;
}


int main(){
	int n = 5;
	double x0 = 0.0;
	LYM lym;
	lym.set_c(2.0);
	vector<double> x ( n, x0);
	opt test_opt(nlopt::LD_LBFGS , n);
	double *func_data;
	func_data = new double[1];
	func_data[0] = 3.0;
	test_opt.set_min_objective( LYM::f , func_data);
	test_opt.set_maxeval(100);
	double opt_f = 100.0;

	test_opt.optimize( x , opt_f);
	cout<< "Optimization Result:" << endl;
	cout << "x:" << endl;
	cout << x[0] << x[1] << x[2]<< endl;
	cout << "Final ObjF:" << endl;
	cout << opt_f << endl;	

	delete []func_data;

	return 0;
}

