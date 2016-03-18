#include <iostream>
#include <fstream>

//#include "opca.h"
#include "opca_binary.h"
using namespace std;

#define S_file "S2.txt"
#define m_file "m2.txt"
#define xi_true_file "xi_true.dat"
#define usig_file "USig.txt"

void LoadVector(const char* filename, double *m_, int Nc_);
void SaveMatrix(const char* filename, double *matrix, int row, int col);


int main()
{
	int l = 70;
	double gamma = 0.8;
	int Nc = 3600;
	int Nr = 1000;
    OPCA_BINARY opca_test(l,gamma,Nc,Nr); // declare an opca_binary object
	opca_test.InputUSig(usig_file); // Input U*Sig from file
	double *xi_true , *m_opca , *S , *m_true_sgems;
	xi_true = new double[l];
	LoadVector(xi_true_file,xi_true , l);
	m_opca = new double[Nc];
	S = new double[Nc * l];
	opca_test.GenerateOPCARealization(xi_true , m_opca , S);
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

