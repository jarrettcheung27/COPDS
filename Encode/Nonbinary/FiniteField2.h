// FiniteField2.h: interface for the CFiniteField2 class.
//
//////////////////////////////////////////////////////////////////////
#include "stdafx.h"

#if !defined(AFX_FINITEFIELD2_H__46CA0220_A3AA_432F_BBA0_5392F4A3E7D2__INCLUDED_)
#define AFX_FINITEFIELD2_H__46CA0220_A3AA_432F_BBA0_5392F4A3E7D2__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Designed by Xiao Ma, Sun Yat-sen University (maxiao@mail.sysu.edu.cn). 
//The program can only be employed for academic research.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class CFiniteField2  
{
public:
	CFiniteField2();
	virtual ~CFiniteField2();
	
	int m_len_element;
	int m_num_element;
	int *m_gen_poly;
	int *Exp;
	int *Log;
	int **Vec;

	int Add(int alpha, int beta);
	int Mult(int alpha, int beta);
	int Div(int alpha, int beta);
	int Pow(int alpha, int exp);
	int VecInd(int alpha, int index);

	void Malloc(int m);
	void Free();

	void VecAdd(int *alpha, int *beta, int *gamma, int len);
	void VecMult(int *alpha, int *beta, int *gamma, int len);
	void ScalarMult(int alpha, int *beta, int *gamam, int len);
	int PolyEvaluate(int x, int *poly, int deg);

	int PolyMult(int *A, int deg_A, int *B, int deg_B, int *C);
	int PolyAdd(int *A, int deg_A, int *B, int deg_B, int *C);

private:
	static const int PrimitivePolynomial[17][6];
};

#endif // !defined(AFX_FINITEFIELD2_H__46CA0220_A3AA_432F_BBA0_5392F4A3E7D2__INCLUDED_)
