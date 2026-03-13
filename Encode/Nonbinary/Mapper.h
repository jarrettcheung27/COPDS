#pragma once
#include <cstdio>

//created by liangchulong @ 2014-08-04
class CModem
{
public:
	CModem(void);
	~CModem(void);

	int m_num_signal;
	int m_len_signal;
	double **m_signal_set;
	double *m_sym_prob;
	double m_Es;

	void Malloc(const char* file_name);
	void Free();
	void Mapping(int *cc, double *xx, int len_symseq);
	void Demapping(double *yy, double *sym_prob, double sigma2, int len_symseq);
	void Demapping(double *yy, double **sym_prob, double sigma, int len_symseq);

	//void PrintCodeParamter(FILE* fp);//added by liangchulong @ 2013-08-25-19-38
};

