//Writer : Terrence Zhao
//Time   :2009.04.08
//Function:random interleaver generator..........

#include "stdafx.h"
#include "Interleaver.h"
#include "Random.h"

extern CLCRandNum rndGen0;
extern CWHRandNum rndGen1;

Interleaver::Interleaver()
{}

Interleaver::~Interleaver()
{}

void Interleaver::Init_Interleaver(int *pai , int period)
{
	int t;
	double random_num;
	int temp;
	int position;

	for (t = 0; t < period; t++)
		pai[t] = t;

	for (t = period - 1; t > 0; t--){
		random_num = rndGen0.Uniform();
		position = (int) (random_num * t);
		temp = pai[position];
		pai[position] = pai[t];
		pai[t] = temp;
	}

	return;
}