#include "stdafx.h"
#include "Qary_Gauss_E.h"


extern CFiniteField2 GF;

Qary_Matrix_Mani::Qary_Matrix_Mani()
{}
Qary_Matrix_Mani::Qary_Matrix_Mani(int row , int col , int rank, int **_matrix)
{
	this->m_num_col = col;
	this->m_num_row = row;
	this->m_rank = rank;
	this->matrix = new int*[this->m_num_row];
	for(int i = 0 ; i < this->m_num_row ; i ++)
	{
		this->matrix[i] = new int[this->m_num_col];
	}

	for(int i = 0 ; i < this->m_num_row ; i ++)
	{
		for(int j = 0 ; j < this->m_num_col ; j ++)
		{
			this->matrix[i][j] = _matrix[i][j];
		}
	}
	this->enc_matrix = new int*[m_num_row];
	for (int i = 0; i < m_num_row; i ++)
	{
		this->enc_matrix[i] = new int[m_num_col];
	}
	tempP = new int[this->m_num_col];
}
//the output matrix is of the form : [I , -P]
//the corresponding generator matrix is of the form [T(P) , I]
//[information bit on the right] when use the elim result for encoding
void Qary_Matrix_Mani::SysMatrix()
{
	int flag;
	int i;
	int ii;
	int jj;
	int j;
	int n;
	int m;
	m_codechk = 0;
	int temp(0);
	int temp_index;
	for(i = 0 ; i < this->m_num_col ; i ++)
	{
		tempP[i] = i;
	}

	int **tempH = new int*[this->m_num_row];
	for(i = 0 ; i < this->m_num_row ; i ++)
	{
		tempH[i] = new int[this->m_num_col];
	}

	for (i = 0; i < m_num_row; i ++)
	{
		for (j = 0; j < m_num_col; j ++)
		{
			this->enc_matrix[i][j] = this->matrix[i][j];
			tempH[i][j] = this->matrix[i][j];
		}
	}

	for (i = 0; i < m_num_row; i ++)
	{
		//从第i行第i个位置开始，找第一个非0元素，从上往下，从左往往右。先找行，再找列。
		flag = 0;
		for (jj = i; jj < m_num_col; jj ++)
		{
			for (ii = i; ii < m_num_row; ii ++)
			{
				if (enc_matrix[ii][jj] != 0)
				{
					flag = 1;
					break;
				}
			}
			if (flag == 1)
			{
				m_codechk++;
				break;
			}
		}

		if (flag == 0)
			break;
		else
		{
//swap i and ii row
//换行
			if (ii != i)
			{
				for (n = 0; n < m_num_col; n++)
				{
					temp = enc_matrix[i][n];
					enc_matrix[i][n] = enc_matrix[ii][n];
					enc_matrix[ii][n] = temp;
				}
			}
//swap i and jj col
//换列
			if (jj != i)
			{
				temp_index = tempP[i];
				tempP[i] = tempP[jj];
				tempP[jj] = temp_index;

				for (m = 0; m < m_num_row; m++)
				{
					temp = enc_matrix[m][i];
					enc_matrix[m][i] = enc_matrix[m][jj];
					enc_matrix[m][jj] = temp;
				}
			}
//elimination
//消元
			for (m = 0; m < m_num_row; m++)
			{
				if (m != i && enc_matrix[m][i] != 0)
				{
					temp = GF.Div(1, enc_matrix[i][i]);
					temp = GF.Mult(temp, enc_matrix[m][i]);
					for (n = 0; n < m_num_col; n++)
					{
						enc_matrix[m][n] = GF.Add(enc_matrix[m][n], GF.Mult(temp, this->enc_matrix[i][n]));
					}
				}
			}
			temp = GF.Div(1, enc_matrix[i][i]);
			for (n = 0; n < m_num_col; n++)
			{
				enc_matrix[i][n] = GF.Mult(temp, this->enc_matrix[i][n]);
			}
		}
	}

	for(i = 0 ; i <  this->m_num_row ; i ++)
	{
		for(j = 0 ; j < this->m_num_col ; j ++)
		{
			matrix[i][j] = tempH[i][tempP[j]];
		}
	}

	string parity_filename = "parity0099.txt";
	std::ofstream oFile(parity_filename.c_str() , ios_base::app); 
	if(!oFile)
	{
		std::cerr<<"Can't open the file ber.txt or the file fer.txt !Please check it !"<<std::endl;
		exit(-1);
	}

	//输出译码矩阵
	oFile<<"row_number"<<std::endl;
	oFile<<this->m_num_row<<std::endl;
	
	oFile<<"col_number"<<std::endl;
	oFile<<this->m_num_col<<std::endl;

	oFile<<"rank"<<std::endl;
	oFile<<this->m_codechk<<std::endl;

	oFile<<"row_number"<<std::endl;
	oFile<<this->m_num_row<<std::endl;

	oFile<<"*****1's_location/per_row****"<<std::endl;
	for(int i = 0 ; i < this->m_num_row ; i ++)
	{
		n = 0;
		//oFile<<i<<"   "<<this->q_ary - 1<<"  ";
		for(int j = 0 ; j < this->m_num_col ; j ++)
		{
			if(matrix[i][j] != 0)
			{
				n ++;
			}
		}
		oFile<<i<<"   "<<n<<"   ";
		for(int j = 0 ; j < m_num_col ; j ++)
		{
			if(matrix[i][j] != 0)
			{
				oFile<<j<<"   "<<matrix[i][j]<<"  ";
			}
		}
		oFile<<std::endl;
	}

	oFile<<"*****1's_location/per_col****"<<std::endl;
	for(int i = 0 ; i < this->m_num_col ; i ++)
	{
		n = 0;
		//oFile<<i<<"   "<<this->q_ary - 1<<"  ";
		for(int j = 0 ; j < this->m_num_row ; j ++)
		{
			if(matrix[j][i] != 0)
			{
				n ++;
			}
		}
		oFile<<i<<"   "<<n<<"   ";
		for(int j = 0 ; j < m_num_row ; j ++)
		{
			if(matrix[j][i] != 0)
			{
				oFile<<j<<"   "<<matrix[j][i]<<"  ";
			}
		}
		oFile<<std::endl;
	}

	oFile<<"*****1's_location/per_row****"<<std::endl;
	for(int i = 0 ; i < this->m_codechk ; i ++)
	{
		n = 0;
		//oFile<<i<<"   "<<this->q_ary - 1<<"  ";
		for(int j = 0 ; j < m_num_col ; j ++)
		{
			if(enc_matrix[i][j] != 0)
			{
				n ++;
			}
		}
		oFile<<i<<"   "<<n<<"   ";
		for(int j = 0 ; j < m_num_col ; j ++)
		{
			if(enc_matrix[i][j] != 0)
			{
				oFile<<j<<"   "<<enc_matrix[i][j]<<"  ";
			}
		}
		oFile<<std::endl;
	}
//输出，检验消元结果
	/*std::cout<<std::endl;
	for(i = 0 ; i < this->m_num_row ; i ++)
	{
		for(j = 0 ; j < this->m_num_col ; j ++)
		{
			std::cout<<this->enc_matrix[i][j].get_vec_num()<<"  ";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;*/
}
void Qary_Matrix_Mani::Gene_Enc_Matrix()
{
	//int temp(0);
	//for(int i = 0 ; i < this->m_codechk ; i ++)//row
	//{
	//	temp.inv(this->enc_matrix[i + this->m_num_col - this->m_codechk][i + this->m_num_col - this->m_codechk]);
	//}
}
Qary_Matrix_Mani::~Qary_Matrix_Mani()
{
	for(int i = 0 ; i < this->m_num_row ; i ++)
	{
		delete []this->matrix[i];
	}
	delete []this->matrix;



	for (int i = 0; i < m_num_row; i ++)
	{
		delete []this->enc_matrix[i];
	}
	delete []this->enc_matrix;


	delete []tempP;
}