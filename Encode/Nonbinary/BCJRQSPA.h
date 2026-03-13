#if !defined(AFX_BCJRQSPA_H__46CA0220_A3AA_432F_BBA0_5392F4A3E7D2__INCLUDED_)
#define AFX_BCJRQSPA_H__46CA0220_A3AA_432F_BBA0_5392F4A3E7D2__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include "stdafx.h"
#include "util.h"
#include "QaryLDPC.h"
#include "Mapper.h"


class BCJRQSPA
{
public:
	BCJRQSPA();
	BCJRQSPA(string filename , int max_iteration, string mapping_name);
	~BCJRQSPA();

	int get_block_length() const;
	int get_total_length() const;
	void encoder4BiBo(int *uu , int * cc);
	bool is_codeword(int *b_code);
	

private:
	int parity_row;   //row number of the parity check matrix
	int parity_column;//column number of the parity check matrix
	int rank;         //rank of the parity check matrix
	                  //determine the check bit of the encoding process
	int **parity_matrix;
	int **enc_parity_matrix;
	int d;

	int max_iteration;

	int * row_weight_parity;
	int **row_location_parity;
	int * row_weight_enc;
	int **row_location_enc;

public:
	int q_ary;
	void Malloc(string filename, string mapping_name);
	CModem m_modem;

	Tanner_Graph tanner;	
	int m_degree;//GF(2^m)�е�m
	int FFTQSPA4BiBo(double *rr, int *uu);

	double *total_mul;

//for message exchange
	double **channel_for_spa;
};

#endif // _MSC_VER > 1000
