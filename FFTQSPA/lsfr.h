//----author : terence zhao
//----time   : 2008.11.20
//----lastUpdate : 2008.11.27

#ifndef _LSFR_H
#define _LSFR_H



#include <math.h>
class lsfr
{
public:
	int _conn_polynomial;
	int _degree;
	int _init_state;
	char * _curr_state;
	int count;
	char *_storage;
public:
	/*input :
    //       conn_polynomial : the decimal form of the polynomial from high to low...
	                           example :x^4 + x^2 + x + 1 corresponds  to 10111, and to  23
							   and it is stored converse.
		     init_state      :the init state of the register.
			                    example: state 1000 represents by 1.
			 degree          : the order of the polynomial.

	*/
	lsfr(int conn_polynomial ,int degree ,int init_state = 1);
	lsfr()
	{
		this->_conn_polynomial = this->_degree = this->_init_state = this->count =0;
	}
	~lsfr();
	void setstate(int state);
	int getcurrstate();
	int getcount();
	char step();
	char step(int &state);
	void steps(int steps ,char *outputs);

	static int chartoint(char c);    //make one char to its corresponds int.....
	static int btod(char *bi,int de);//left is the least singnificant bit....
	static void dtob(char *bi,int de , int num);//left is the least singnificant bit....
};


#endif