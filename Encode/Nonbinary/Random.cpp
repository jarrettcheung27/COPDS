// Writer : Terrence Zhao
// Time   :2009.04.08
// Function:two random number generator...

#include "stdafx.h"
#include <complex>
#include "Random.h"

// Global random generator instance used across translation units
CLCRandNum rndGen0;

CLCRandNum::CLCRandNum()
{
}

CLCRandNum::~CLCRandNum()
{
}

///////////////////////////////////////////////////////////////
const int CLCRandNum::A = 48271;
const long CLCRandNum::M = 2147483647;
const int CLCRandNum::Q = M / A;
const int CLCRandNum::R = M % A;

////////////////////////////////////////////////////////////////
void CLCRandNum::SetSeed(int flag)
{
	if (flag < 0)
		state = 17;
	else if (flag == 0)
	{
		state = 0;
		while (state == 0)
		{
			srand((unsigned)time(NULL));
			state = rand();
		}
	}
	else
	{
		fprintf(stdout, "\nEnter the initial state: ");
		fscanf(stdin, "%ld", &state);
	}

	return;
}

///////////////////////////////////////////////////////////////
void CLCRandNum::PrintState(FILE *fp)
{
	fprintf(fp, "\n***init_state = %ld***\n", state);

	return;
}
///////////////////////////////////////////////////////////////
double CLCRandNum::Uniform()
{
	double u;

	int tmpState = A * (state % Q) - R * (state / Q);
	if (tmpState >= 0)
		state = tmpState;
	else
		state = tmpState + M;

	u = state / (double)M;

	return u;
}

///////////////////////////////////////////////////////////////
void CLCRandNum::Normal(double *nn, int len_nn)
{
	double x1, x2, w;
	int t;

	for (t = 0; 2 * t + 1 < len_nn; t++)
	{
		w = 2.0;
		while (w > 1.0)
		{
			x1 = 2.0 * Uniform() - 1.0;
			x2 = 2.0 * Uniform() - 1.0;

			w = x1 * x1 + x2 * x2;
		}

		w = sqrt(-2.0 * log(w) / w);

		nn[2 * t] = x1 * w;
		nn[2 * t + 1] = x2 * w;
	}

	if (len_nn % 2 == 1)
	{
		w = 2.0;
		while (w > 1.0)
		{
			x1 = 2.0 * Uniform() - 1.0;
			x2 = 2.0 * Uniform() - 1.0;

			w = x1 * x1 + x2 * x2;
		}

		w = sqrt(-2.0 * log(w) / w);

		nn[len_nn - 1] = x1 * w;
	}

	return;
}

/////////////////////////////////////////////////////////////
// The following generator employs the Wichman-Hill algorithm
/////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CWHRandNum::CWHRandNum()
{
}

CWHRandNum::~CWHRandNum()
{
}

//////////////////////////////////////////////////////////////////////
void CWHRandNum::SetSeed(int flag)
{
	if (flag < 0)
	{
		X = 13;
		Y = 37;
		Z = 91;
	}
	else if (flag == 0)
	{
		X = 0;
		Y = 0;
		Z = 0;
		while (X == 0 || Y == 0 || Z == 0)
		{
			srand((unsigned)time(NULL));
			X = rand();
			Y = rand();
			Z = rand();
		}
	}
	else
	{
		fprintf(stdout, "\nEnter the initial state (X Y Z): ");
		fscanf(stdin, "%d %d %d", &X, &Y, &Z);
	}

	return;
}

///////////////////////////////////////////////////////////////
void CWHRandNum::PrintState(FILE *fp)
{
	fprintf(fp, "\n***init_state (X Y Z) = %d %d %d***\n", X, Y, Z);

	return;
}
///////////////////////////////////////////////////////////////
double CWHRandNum::Uniform()
{
	double U;

	X = 171 * X % 30269;
	Y = 172 * Y % 30307;
	Z = 170 * Z % 30323;

	U = X / 30269.0 + Y / 30307.0 + Z / 30323.0;
	U = U - int(U);

	return U;
}

///////////////////////////////////////////////////////////////
void CWHRandNum::Normal(double *nn, int len_nn)
{
	double x1, x2, w;
	int t;

	for (t = 0; 2 * t + 1 < len_nn; t++)
	{
		w = 2.0;
		while (w > 1.0)
		{
			x1 = 2.0 * Uniform() - 1.0;
			x2 = 2.0 * Uniform() - 1.0;

			w = x1 * x1 + x2 * x2;
		}

		w = sqrt(-2.0 * log(w) / w);

		nn[2 * t] = x1 * w;
		nn[2 * t + 1] = x2 * w;
	}

	if (len_nn % 2 == 1)
	{
		w = 2.0;
		while (w > 1.0)
		{
			x1 = 2.0 * Uniform() - 1.0;
			x2 = 2.0 * Uniform() - 1.0;

			w = x1 * x1 + x2 * x2;
		}

		w = sqrt(-2.0 * log(w) / w);

		nn[len_nn - 1] = x1 * w;
	}

	return;
}
/////////////////////////////////////////////////////////////////
void CWHRandNum::Normal_c(complex<double> *nn, int len_nn, double sigma)
{
	double *ss = new double[len_nn];
	this->Normal(ss, len_nn);
	double t = 0.0;
	/*for(int i = 0 ; i < len_nn ; i ++)
	{
		t +=pow(ss[i] , 2.0);
	}
	std::cout<<t/len_nn;*/
	for (int i = 0; i < len_nn; i++)
	{
		nn[i] = complex<double>(nn[i].real(), nn[i].imag() + sigma * ss[i]);
	}
	this->Normal(ss, len_nn);

	/*for(int i = 0 ; i < len_nn ; i ++)
	{
		t +=pow(ss[i] , 2.0);
	}
	std::cout<<t/len_nn;*/

	for (int i = 0; i < len_nn; i++)
	{
		nn[i] = complex<double>(nn[i].real() + sigma * ss[i], nn[i].imag());
	}
	delete[] ss;
}