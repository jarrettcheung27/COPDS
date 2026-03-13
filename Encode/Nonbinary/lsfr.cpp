#include "stdafx.h"
#include "lsfr.h"
int lsfr::chartoint(char c)
{
	char * r = new char[1];
	r[0] = c;
	int re = atoi(r);
	delete []r;
	return re;
}

int lsfr::btod(char *bi,int de)
{
	int re = 0;
	for(int i = 0 ; i < de ; i ++)
	{
		if(bi[i] != '0')
			re += static_cast<int>(pow(static_cast<double>(2),i));
	}
	return re;
}

void lsfr::dtob(char *bi,int de , int num)
{
	int re = de;
	int s = 1;
	for(int i = 0 ; i < num ;i ++)
	{
		bi[i] = de%2 + '0';
		de = de>>1;
	}

	/*while(s <= num)
	{
		bi[num - s] = '0';
		if(re % 2 == 0)
			bi[num - s] = '0';
		else
			bi[num - s] = '1';
		re = re / 2;
		s ++;
	}*/
}

lsfr::lsfr(int conn_polynomial, int degree, int init_state)
{
	this->_conn_polynomial = conn_polynomial;
	this->_degree = degree;
	this->_init_state = init_state;
	this->count = 0;
	this->_storage = new char[_degree+1];
	this->_curr_state = new  char[_degree];
	for(int i = 0 ; i < degree+1 ; i ++)
	{
		if( i < degree)
		{
			this->_curr_state[i] = init_state%2+'0';
		    init_state = init_state>>1;
		}
		this->_storage[i] = conn_polynomial%2+'0';
		conn_polynomial = conn_polynomial>>1;
	}
}

lsfr::~lsfr()
{
	delete [] _storage;
	delete [] _curr_state;
}
void lsfr::setstate(int state)
{
	this->_init_state = state;
	this->count = 0;
}

char lsfr::step()
{
	count++;
	int last = chartoint(this->_curr_state[_degree-1]);
	for(int i = _degree-1 ; i >=1; i--)
	{
		this->_curr_state[i] = (last * chartoint(this->_storage[i])+chartoint(this->_curr_state[i-1])) % 2 + '0';
	}
	this->_curr_state[0] = last*chartoint(this->_storage[0]) + '0';
	return last + '0';
}

char lsfr::step(int &state)
{
	state = btod(this->_curr_state,_degree);
	char re = this->step();
	return re;
}

void lsfr::steps(int steps, char *outputs)
{
	outputs = new char[steps];
	for(int i = 0 ; i < steps ; i ++)
	{
		outputs[i] = lsfr::step();
	}
}

int lsfr::getcurrstate()
{
	return btod(this->_curr_state,this->_degree);
}

int lsfr::getcount()
{
	return this->count;
}