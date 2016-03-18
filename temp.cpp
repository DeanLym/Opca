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
	test_opt.optimize( x , opt_f);

//	copy(x.begin(), x.end(), m);
	for(int i=0; i < n; i++){
		m[i] = x[i];
	}

	delete []f_data;
	delete []a;
    return true;
}

