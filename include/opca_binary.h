#ifndef OPCA_BINARY_H
#define OPCA_BINARY_H

#include "opca.h"

class OPCA_BINARY : public OPCA
{
public:
    OPCA_BINARY(int l, double gamma, int Nc, int Nr);
	~OPCA_BINARY();
    bool GenerateOPCARealization(double *xi, double *m, double *S);
    bool GenerateOPCARealization(double *xi,double *m);
};



#endif // OPCA_BINARY_H
