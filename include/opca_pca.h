#ifndef OPCA_PCA_H
#define OPCA_PCA_H

#include "opca.h"

class OPCA_PCA : public OPCA
{
public:
    OPCA_PCA(int l, int Nc, int Nr);
	~OPCA_PCA();
    bool GenerateOPCARealization(double *xi, double *m, double *S);
    bool GenerateOPCARealization(double *xi,double *m);
};



#endif // OPCA_BINARY_H
