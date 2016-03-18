#include <iostream>
#include <fstream>

#include "include/opca_binary.h";
using namespace std;

#define S_file "S.txt"
#define m_file "m.txt"

void LoadMTrue(double *m_, int Nc_);
void SaveMatrix(const char* filename, double *matrix, int row, int col);


int main()
{
	int l = 70;
	double gamma = 0.8;
	int Nc = 3600;
	int Nr = 1000;
    OPCA_BINARY opca_test(l,gamma,Nc,Nr); // declare an opca_binary object

	string data_file("data_matrix.dat");  
	opca_test.InputDataMatrix(data_file.c_str());  // Input the sgems data	
	// opca_test.InputDataMatrix(Nc,Nr,X); // Another way to input the sgems data
	opca_test.ConstructYAndPerformSVD(); // Perform SVD
	double *xi_true , *m_opca , *S , *m_true_sgems;
	m_true_sgems = new double[Nc];
	xi_true = new double[l];
	LoadMTrue(m_true_sgems , Nc);  // Input m_true
	opca_test.GetTrueXi(m_true_sgems, xi_true); // Get xi_true corresponding to m_true
	m_opca = new double[Nc]; //m_opca is the opca realization corresponding to xi
	S = new double[Nc * l]; //S is the gridient of dm/dxi
	opca_test.GenerateOPCARealization(xi_true , m_opca , S); // Input xi_true, output m_opca and S
	SaveMatrix(S_file,S,Nc,l);
	SaveMatrix(m_file,m_opca,Nc,1);
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

