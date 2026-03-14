#ifndef _QARY_GAUSS_E_H
#define _QARY_GAUSS_E_H



#include "stdafx.h"
#include "util.h"

class Qary_Matrix_Mani
{
public:
	int m_num_col;
	int m_num_row;
	int m_rank;
	int m_codechk;

	int **matrix;
	int **enc_matrix;

	int *tempP;//store rearrange information
	Qary_Matrix_Mani();
	Qary_Matrix_Mani(int row , int col , int rank , int **_matrix);
	void SysMatrix();
	void Gene_Enc_Matrix();
	~Qary_Matrix_Mani();
};


#endif