//Writer : Terrence Zhao
//Time   :2009.04.08
//lastupdate : 2009.04.11
//Function:We use this class to generate a certain long random interleaver..
#ifndef _INTERLEAVER_H
#define _INTERLEAVER_H



class Interleaver
{
public:

	//constructor
	Interleaver();

	//destructor
	virtual ~Interleaver();

	//first random interleaver generator.It generate a random number and change the present number with the number at the generating position.
	/*
	      in : period : the length of the interleaver
		  out : pai   : the interleaver
	*/
	void Init_Interleaver(int *pai , int period);

	void Interleaving(double *inpput, int *pai);
	void DeInterleaving(double *inpput, int *pai);
};

#endif