//Writer : Terrence Zhao
//Time   :2009.04.08
//Function:utility...

#include "stdafx.h"
#include "util.h"


double msgTemp[256];
int index[256];

void ProbClip(double *xx, int len_xx)
{
	int i;
	double sum = 0.0;
	for (i = 0; i < len_xx; i++){
		if (xx[i] < SMALLPROB)
			xx[i] = SMALLPROB;
		else if (xx[i] > 1.0 - SMALLPROB)
			xx[i] = 1.0 - SMALLPROB;
		sum += xx[i];
	}
//norm to make sure that the sum of vector [xx] is equal to 1.0.
	for (i = 0; i < len_xx; i++)
	{
		xx[i] /= sum;	
	}
	return;
}


void LogLikeClip(double *xx, int len_xx)
{
	int i = 0; 
	for(i = 0 ; i < len_xx; i ++)
	{
		if(xx[i] > 200)
			xx[i] = 200;
		else if(xx[i] < -200)
			xx[i] = -200;
	}
	return;
}
//From dec to Binary string]
//Here Left(zero place) is the most singnificant
void Dec2Bin(int d, int *b, int len_b)
{
	for (int i = 0; i < len_b; i++)
		b[len_b-i-1] = (d >> i) % 2;
	return;
}

//From dec to Binary string]
//Here right(zero place) is the most singnificant
void Dec2BinLS(int d, int *b, int len_b)
{
	for (int i = 0; i < len_b; i++)
		b[i] = (d >> i) % 2;
	return;
}

//From binary string to its corresponding deci format.
//Here the left(that is the zero place) is the most singnificant
int Bin2Dec(int *b, int len_b)
{
	int d;
	d = b[0];
	for (int i = 1; i < len_b; i++)
		d = (d << 1) + b[i];
	return d;
}

//From binary string to its corresponding deci format.
//Here the right(that is the zero place) is the most singnificant
int Bin2DecLS(int *b, int len_b)
{
	int d;
	d = b[len_b - 1];
	for (int i = 1; i < len_b; i++)
		d = (d << 1) + b[len_b - 1 - i];
	return d;
}


int bit_reversing(int k , int q_ary)
{
	int j = 0;
	int m = k;
	for(int i = 0 ; i < q_ary ; i ++)
	{
		j = 2 * j + m - 2 * (m>>1);
		m = m>>1;
	}
	return j;
}


void Hadamard_Transform(double *data_in , double *data_out , int stage)
{
	int num = 1<<stage;
	int curr_dist = 0;
	int block = 0;
	int block_row_num = 0 ;
	double *curr = new double[num];
	for(int i = 0 ; i < num ; i ++)
	{
		curr[i] = data_in[i];
	}
	for(int i = 0 ; i < stage ; i ++)
	{
		curr_dist = 1<<i;
		block = num>>(i + 1);
		block_row_num = 1<<(i + 1);
		for(int j = 0 ; j < block ; j ++)
		{
			for(int s = 0 ; s < curr_dist ; s ++)
			{
				data_out[j * block_row_num + s] = curr[j * block_row_num + s] + curr[j * block_row_num + curr_dist + s];
				data_out[j * block_row_num + s + curr_dist] = curr[j * block_row_num + s] - curr[j * block_row_num + curr_dist + s];
			}
		}
		for(int j = 0 ; j < num ; j ++)
		{
			curr[j] = data_out[j];
		}
	}
	for(int j = 0 ; j < num ; j ++)
	{
		data_out[j] /= (double)(num>>1);
	}
	delete []curr;
}


int int_pow(int q_ary)
{
	int re = 1;
	for(int i = 0 ; i < q_ary ; i ++)
	{
		re = re<<1;
	}
	return re;
}

//Left  is LSB
//Right is MSB
int btod(char *bi,int de)
{
	int re = 0;
	for(int i = 0 ; i < de ; i ++)
	{
		if(bi[i] != '0')
			re += static_cast<int>(pow(static_cast<double>(2),i));
	}
	return re;
}



void FFTTruncLEFT(double *msg, int q_ary, int nt)
{
	double sum = 0.0;
	for(int i = 0 ; i < q_ary ; i ++)
	{
		msgTemp[i] = msg[i];
	}
	BubbleSortTerminate(msgTemp, index, q_ary, nt);
	for(int i = 0 ; i < q_ary ; i ++)
	{
		if(i>=nt)
			msg[index[i]] = SMALLPROB;
		else
			sum += msg[index[i]];
	}
	sum = 0.0;
	for(int i = 0 ; i < nt ; i ++)
	{
		msg[index[i]] *= 1.0;
		sum += msg[index[i]];
	}
	for(int i = 0 ; i < q_ary ; i ++)
	{
		if(i>=nt)
		{
			msg[index[i]] = (1.0-sum)/(q_ary-nt);
			//msg[index[i]] = SMALLPROB;
		}
		else
		{
			//msg[index[i]] += (1.0-sum)/q_ary;
			//msg[index[i]] /= (sum+nt*SMALLPROB);
		}
	}
	//for(int i = 0 ; i < nt ; i ++)
	//{
	//	msg[index[i]] /= (sum + (q_ary-nt)*SMALLPROB);
	//	//msg[index[i]] /= sum;
	//}
}



void FFTTruncRIGHT(double *msg, int q_ary, int nt)
{
	double sum = 0.0;
	for(int i = 0 ; i < q_ary ; i ++)
	{
		msgTemp[i] = msg[i];
	}
	BubbleSortTerminate(msgTemp, index, q_ary, nt);
	for(int i = 0 ; i < nt ; i ++)
	{
			msg[index[i]] *= 1.0;
			sum += msg[index[i]];
	}
	for(int i = nt ; i < q_ary ; i ++)
	{
		msg[index[i]] = (1.0 - sum)/(q_ary-nt);
	}
	/*for(int i = 0 ; i < q_ary ; i ++)
	{
		if(i>=nt)
			msg[index[i]] = (1.0-sum)/q_ary;
		else
			msg[index[i]] += (1.0-sum)/q_ary;
	}*/
	/*for(int i = 0 ; i < q_ary ; i ++)
	{
		if(i>=nt)
			msg[index[i]] = SMALLPROB;
		else
			msg[index[i]] /= (sum + (q_ary-nt)*SMALLPROB);
	}*/
}


void BubbleSortTerminate(double *value, int *index, int len, int ter_len)
{
	int i, j;
	int cnt;
	int itemp;
	double vtemp;

	for (i = 0; i < len; i++)
		index[i] = i;

	for (i = 0; i < ter_len; i++){
		cnt = 0;
		for (j = len-1; j > i; j--){
			if (value[j] > value[j-1]){
				vtemp = value[j-1];
				value[j-1] = value[j];
				value[j] = vtemp;
				
				itemp = index[j-1];
				index[j-1] = index[j];
				index[j] = itemp;

				cnt++;
			}
		}
		if (cnt == 0)
			break;
	}


	return;
}


