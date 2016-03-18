#include <iostream>
#include <fstream>

#include "include/opca_bimodal.h"

using namespace std;

#define S_file "S.txt"
#define m_file "m.txt"
#define xm_file "xm.dat"
#define u_sig_file "u_sig.dat"

void LoadMTrue(double *m_, int Nc_);
void SaveMatrix(const char* filename, double *matrix, int row, int col);
void SaveVector(const char* filename, double *v, int n);

int main()
{
	int l = 70;
	double gamma = 16.5;
	int Nc = 3600;
	int Nr = 1000;
    OPCA_BIMODAL opca_test(l,gamma,Nc,Nr); // declare an opca_bimodal object

	string data_file("data_bimodal.dat");  

	opca_test.InputDataMatrix(data_file.c_str());  // Input the sgems data	
	opca_test.ConstructYAndPerformSVD(); // Perform SVD
	opca_test.OutputUSig("USig.dat");
	opca_test.OutputXm("xm.dat");
	opca_test.OutputU("U.dat");
	opca_test.OutputSig("Sig.dat");

/*
	opca_test.InputU("U.dat");
	opca_test.InputSig("Sig.dat");
	opca_test.InputXm(xm_file);
*/


//	opca_test.InputXm(xm_file);
	opca_test.InputUSig(u_sig_file);

	double *xi_true , *m_opca , *m_true_sgems, * S;
	m_true_sgems = new double[Nc];
	xi_true = new double[l];
	m_opca = new double[Nc]; //m_opca is the opca realization corresponding to xi
	S = new double[Nc * l];
	LoadMTrue(m_true_sgems , Nc);  // Input m_true
	opca_test.GetTrueXi(m_true_sgems, xi_true); // Get xi_true corresponding to m_true
	double logkMin = 0, logkMax = 10.0;
	double optMu1 = 2.5, optMu2 = 8.5;
	double optV1 = 2.1, optV2 = 1.6;
	opca_test.set_bounds( logkMin , logkMax );
	opca_test.set_mu( optMu1 , optMu2 );	
	opca_test.set_var( optV1 , optV2 );
	opca_test.GenerateOPCARealization(xi_true , m_opca , S); // Input xi_true, output m_opca and S
//	opca_test.GenerateOPCARealization(xi_true , m_opca); // Input xi_true, output m_opca and S
	SaveMatrix("xi_true.dat",xi_true,l,1);	
	SaveMatrix(S_file,S,Nc,l);
	SaveMatrix(m_file,m_opca,Nc,1);
	SaveVector("S_vector.dat",S,Nc*l);	

	delete []S;
    delete []xi_true;
	delete []m_opca;
	delete []m_true_sgems;

    return 0;
}

void LoadMTrue(double *m_, int Nc_){
	ifstream in;
	in.open("m_true.dat");
	for (int i = 0; i < Nc_ ; i++)
		in >> m_[i];
	in.close();
}

void SaveVector(const char* filename, double *v, int n){
	ofstream out;
	out.open(filename);
	for (int i = 0; i < n; i++){
		out << v[i] << endl;
	}
	out.close();
}

void SaveMatrix(const char* filename, double *matrix, int row, int col){
	int index = 0;
	ofstream out;
	out.open(filename);
	for (int i = 0; i < row; i++){
		for(int j = 0 ; j < col; j++){
			index = j * row + i;
			out << matrix[index]<<" ";
		}
		out << endl;
	}
	out.close();
}

