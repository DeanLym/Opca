#ifndef OPCA_H
#define OPCA_H


using namespace std;

class OPCA{
public:
    OPCA(int l, double gamma, int Nc, int Nr);
    bool InputDataMatrix(int Nc, int Nr ,double * X);
    bool InputDataMatrix(const char* file_name);
    bool InputUSig(double *U_Sig);
    bool InputUSig(const char* file_name);
	bool InputXm(double *xm);
	bool InputXm(const char* file_name);
	bool InputU(const char* file_name);
	bool InputU(double *U);
	bool InputSig(const char* file_name);
	bool InputSig(double *Sig);

    bool ConstructYAndPerformSVD();
	void OpcaDebug(char* filename,double *data, int size);

	bool OutputXm(const char* file_name);
	bool OutputUSig(const char* file_name);	
	bool OutputU(const char* file_name);
	bool OutputSig(const char* file_name);
	bool GetTrueXi(double* m_true, double *xi_true);
    virtual bool GenerateOPCARealization(double *xi,double *m, double *S) = 0;
    // GenerateOPCARealization takes xi as input, output result to m and dm_dxi
    // xi(l,1) is the reduced dimension standard gaussian variabl
    // m(Nc,1) is the resulting opca_generation
    // S(Nc,l) is the resulting sensitivity matrix

    virtual bool GenerateOPCARealization(double *xi,double *m) = 0;
    // Use this form when sensitivity is not required


    ~OPCA();

protected:
    bool set_Nc(int Nc);
    bool set_Nr(int Nr);
    bool set_gamma(double gamma);
    bool set_l(int l);

	int input_U_and_Sig_flag_;	//Used when U and Sig are inputed seperately

    int get_Nc();
    int get_Nr();
    double get_gamma();

    double gamma_;
    int l_;  // l is the reduced dimension
    int Nc_; // Nc is the number of cells
    int Nr_; // Nr is the number of realizations
    double *X_; // X is the uncentered data matrix
    double *xm_; // Xm is the mean X
    double *U_; double *Sig_; //Y_ = U_*Sig_*V_;
    double *U_Sig_;
    double *xi_; // xi is the reduced stanford gaussian variable
//    double *m_; // m is the opca realization
//    double *S_; // sensitivity matrix S = dm/dxi;


};

#endif // OPCA_H
