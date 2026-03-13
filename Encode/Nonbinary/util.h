//Writer : Terrence Zhao
//Time   :2009.04.08
//Lastupdate : 2009.04.11
//Function:utility..Basic common part of the whole program.

#ifndef _UTILITY_H
#define _UTILITY_H
/////////
#include "stdafx.h"
//#include "int.cpp"
#define SMALLPROB 1.0e-12





void ProbClip(double *xx, int len_xx);//for stability.

void LogLikeClip(double *xx, int len_xx);//for stability.

void Dec2Bin(int d, int *b, int len_b);//For dec number to binary number array

int Bin2Dec(int *b, int len_b);//from binary number array to dec number.

void Dec2BinLS(int d, int *b, int len_b);

int Bin2DecLS(int *b, int len_b);//from binary number array to dec number.

void BubbleSortTerminate(double *value, int *index, int len, int ter_len);

void FFTTruncLEFT(double *msg, int q_ary, int nt);

void FFTTruncRIGHT(double *msg, int q_ary, int nt);

int btod(char *bi,int de);
//const int operator+(const int& gfa , const int& gfb);
//
//
//const int operator*(const int& gfa , const int& gfb);
//
//
//const int operator-(int& gfa , const int& gfb);

//int powerof(const int& gfa ,int pow);

void Hadamard_Transform(double *data_in , double *data_out , int stage);



int int_pow(int q_ary);
//{
//	int re = 1;
//	for(int i = 0 ; i < q_ary ; i ++)
//	{
//		re = re<<1;
//	}
//	return re;
//}

template<typename data_type>
void fft(data_type *in , data_type *out , int q_ary)
{
	int length = int_pow(q_ary);
	Swap_Via_Bit_Re(in , length , q_ary);
	for(int i = 0 ; i < q_ary ; i ++)
	{
		int L = int_pow(i + 1);
		int r = length / L;
		int L_c = L >>1;
		for(int j = 0 ; j < L_c ; j ++)
		{
			complex<double> womiga = complex<double>(cos(2 * 3.14159 * j / L) , -sin(2 * 3.14159 * j / L));
			for(int k = 0 ; k < r ; k ++)
			{
				complex<double> tao = complex<double>(womiga * in[k * L + j + L_c]);
				in[k * L + j + L_c] = in[k * L + j] - tao;
				in[k * L + j] = in[k * L + j] + tao;
			}
		}
	}
}

template<typename data_type>
void fft_field(data_type *in , data_type *out , int q_ary)
{
	int length = int_pow(q_ary);
	Swap_Via_Bit_Re(in , length , q_ary);
	for(int i = 0 ; i < q_ary ; i ++)
	{
		int L = int_pow(i + 1);
		int r = length / L;
		int L_c = L >>1;
		for(int j = 0 ; j < L_c ; j ++)
		{
			complex<double> womiga = complex<double>(cos(2 * 3.14159 * j / L) , -sin(2 * 3.14159 * j / L));
			for(int k = 0 ; k < r ; k ++)
			{
				complex<double> tao = complex<double>(womiga * in[k * L + j + L_c]);
				in[k * L + j + L_c] = in[k * L + j] - tao;
				in[k * L + j] = in[k * L + j] + tao;
			}
		}
	}
}

template<typename data_type>
void ifft(data_type *in , data_type *out , int q_ary)
{
	int length = int_pow(q_ary);
	Swap_Via_Bit_Re(in , length , q_ary);
	for(int i = 0 ; i < q_ary ; i ++)
	{
		int L = int_pow(i + 1);
		int r = length / L;
		int L_c = L >>1;
		for(int j = 0 ; j < L_c ; j ++)
		{
			complex<double> womiga = complex<double>(cos(2 * 3.14159 * j / L) , sin(2 * 3.14159 * j / L));
			for(int k = 0 ; k < r ; k ++)
			{
				complex<double> tao = complex<double>(womiga * in[k * L + j + L_c]);
				in[k * L + j + L_c] = in[k * L + j] - tao;
				in[k * L + j] = in[k * L + j] + tao;
			}
		}
	}
	for(int i = 0 ; i < length ; i ++)
	{
		in[i] = in[i] / 16.0;
	}
}


int bit_reversing(int k , int q_ary);


template<class data_type>
void Swap_Via_Bit_Re(data_type *k , int length , int q_ary)
{
	data_type temp;
	int curr_location;
	for(int i = 0 ; i < length ; i ++)
	{
		curr_location = bit_reversing(i , q_ary);
		if(curr_location > i)
		{
			temp = k[i];
			k[i] = k[curr_location];
			k[curr_location] = temp;
		}
	}
}



const long long PrimitivePolynomials[17][2] = 
{
	{0, 1},//1
	{1, 3},//1 1//3
	{2, 7},//
	{3, 11},
	{4, 19},//1 1 0 0 1=19
	{5, 37},
	{6, 67},//1 + 2 + 64
	{7, 137},
	{8, 285},//1 0 1 1 1 0 0 0 1=1+4+8+16+256=256+29=285
	{9, 529},
	{10, 1033},//1 + 8 + 1024 = 1033
	{11, 2053},
	{12, 4179},//1 1 0 0 1 0 1 0 0 0 0 0 1 = 1 + 2 + 16 + 64 + 4096 = 4179
	{13, 11011000000001LL},
	{14, 110000100010001LL},
	{15, 110000000000001LL},
	{16, 11010000000010001LL}
};
//{5, 0, 2, 3, 4, 8} meaning: 5 nonzero terms--x^0 + x^2 + x^3 + x^4 + x^8
//////////////////////////////////////////////////////////////////////




////////
//void Count_Mapping(int *or , int *map_result , int up , int max , int q_ary);
#endif