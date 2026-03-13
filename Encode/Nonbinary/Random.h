//Writer : Terrence Zhao
//Time   :2009.04.08
//lastupdate : 2009.04.11
//Function:two random number generator....
#ifndef _RANDOM_H
#define _RANDOM_H

#include <complex>

class CLCRandNum  
{
public:
	CLCRandNum();
	virtual ~CLCRandNum();

	void SetSeed(int flag);
	void PrintState(FILE *fp);
	double Uniform();
	void Normal(double *nn, int len_nn);


private:
    long int state;

    static const int A;
    static const long M;
    static const int Q;
    static const int R;
};


/////////////////////////////////////////////////////////////
//The following generator employs the Wichman-Hill algorithm
/////////////////////////////////////////////////////////////

class CWHRandNum  
{
public:
	CWHRandNum();
	virtual ~CWHRandNum();

	void SetSeed(int flag);
	void PrintState(FILE *fp);
	double Uniform();
	void Normal(double *nn, int len_nn);
	void Normal_c(std::complex<double> *nn, int len_nn  , double sigma);

private:
    int X, Y, Z;
};


#endif