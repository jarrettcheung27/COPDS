#include "stdafx.h"
#include "Mapper.h"
#include "util.h"



CModem::CModem(void)
{
}


CModem::~CModem(void)
{
}


//////////////////////////////////////////////////////////////////////
void CModem::Malloc(const char* file_name)
{
	char temp_str[80];
	char mark[80];	
	int i,j;

	FILE *fp;

	//sprintf(file_name, "Modulation_chart.txt");
	if ((fp = fopen(file_name, "r")) == NULL){
		fprintf(stderr, "\nCan't open the %s file!\n", file_name);
		exit(1);
	}
	fscanf(fp, "%s", temp_str);

	sprintf(mark, "***MappingChart***");
	while (strcmp(temp_str, mark) != 0)
		fscanf(fp, "%s", temp_str);

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &m_num_signal);
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &m_len_signal);
	fscanf(fp, "%s", temp_str);

	m_signal_set = new double *[m_num_signal];
	for (i=0; i<m_num_signal; i++)
		m_signal_set[i] = new double[m_len_signal]; 
	
	m_Es = 0.0;
	for (i=0; i<m_num_signal; i++)
	{
		for (j=0; j<m_len_signal; j++){
			fscanf(fp, "%lf", &m_signal_set[i][j]);
			m_Es += m_signal_set[i][j] * m_signal_set[i][j];
		}
	}
	fclose(fp);
	m_Es = m_Es / m_num_signal;

	//for (i=0; i<m_num_signal; i++)
	//	for (j=0; j<m_len_signal; j++)
	//		m_signal_set[i][j] /= sqrt(m_Es);

	//m_Es = 1.0;//归一化
	m_sym_prob = new double [m_num_signal];

	return;
}

//////////////////////////////////////////////////////////////////////
void CModem::Mapping(int *cc, double *xx, int len_symseq)
{
	int i, j;

	for (i=0; i<len_symseq; i++)
	{
		for (j=0; j<m_len_signal; j++)
		{
			xx[i*m_len_signal+j] = m_signal_set[cc[i]][j];
		}
	}

	return;
}

//////////////////////////////////////////////////////////////////////
void CModem::Free()
{
	delete []m_signal_set;
	delete []m_sym_prob;
	return;
}


void CModem::Demapping( double *yy, double *sym_prob, double sigma2, int len_signalseq )
{
	int i,q,j;
	double sum;
	double sqr_sum;

	for (i=0; i<len_signalseq; i++) //计算 prob(x=0 | y); +1
	{
		for (q=0; q<m_num_signal; q++)
		{
			sqr_sum = 0.0;
			for (j=0; j<m_len_signal; j++)
			{
				sqr_sum += (yy[i*m_len_signal+j]-m_signal_set[q][j]) * (yy[i*m_len_signal+j]-m_signal_set[q][j]);
			}
			m_sym_prob[q] = exp(-0.5 * sqr_sum / sigma2);
		}
		ProbClip(m_sym_prob, m_num_signal);
		sum = 0.0;
		for (q=0; q<m_num_signal; q++)
		{
			sum +=  m_sym_prob[q];
		}
		for (q=0; q<m_num_signal; q++)
		{
			sym_prob[i*m_num_signal+q] = m_sym_prob[q]/sum;
		}			
	} //end i

	return;
}




void CModem::Demapping( double *yy, double **sym_prob, double sigma, int len_signalseq)
{
	int i,q,j;
	double sum;
	double sqr_sum;

	for (i=0; i<len_signalseq; i++) //计算 prob(x=0 | y); +1
	{
		for (q=0; q<m_num_signal; q++)
		{
			sqr_sum = 0.0;
			for (j=0; j<m_len_signal; j++)
			{
				sqr_sum += (yy[i*m_len_signal+j]-m_signal_set[q][j]) * (yy[i*m_len_signal+j]-m_signal_set[q][j]);
			}
			m_sym_prob[q] = sqr_sum / (sigma*sigma);
		}
		//find the minimum distance
		sqr_sum = m_sym_prob[0];
		for (q=1; q<m_num_signal; q++)
		{
			if(m_sym_prob[q] < sqr_sum)
			{
				sqr_sum = m_sym_prob[q];
			}
		}
		//norm the probability
		for (q=0; q<m_num_signal; q++)
		{	
			m_sym_prob[q] -= sqr_sum;		
			if(m_sym_prob[q] > 40)
				m_sym_prob[q] = 40;
		}
		//compute the probability and norm the probability
		sqr_sum = 0.0;
		for (q=0; q<m_num_signal; q++)
		{
			m_sym_prob[q] = exp(-0.5 * m_sym_prob[q]);
			sqr_sum += m_sym_prob[q];
		}
		for (q=0; q<m_num_signal; q++)
		{
			sym_prob[i][q] = m_sym_prob[q]/sqr_sum;
		}		
		ProbClip(m_sym_prob, m_num_signal);
		sum = 0.0;
		for (q=0; q<m_num_signal; q++)
		{
			sum +=  m_sym_prob[q];
		}
		for (q=0; q<m_num_signal; q++)
		{
			sym_prob[i][q] = m_sym_prob[q]/sum;
		}			
	} //end i

	return;
}

//void CModem::PrintCodeParamter( FILE* fp )
//{
//	fprintf(fp, "%%#######%s***MappingChart***Parameters#######\n", ToString(CModem));
//	fprintf(fp, "%%%-20s = [%d, %d];\n", "[num_signal, len_signal]", m_num_signal, m_len_signal);
//	fprintf(fp, "%%%-20s = %lf;\n", "m_Es", m_Es);
//	fprintf(fp, "%%signal_set = [...\n");
//	for (int i = 0; i < m_num_signal; i++)
//	{
//		fprintf(fp, "%%");
//		for (int j=0; j<m_len_signal; j++)
//			fprintf(fp, "%.12lf ", m_signal_set[i][j]);
//		fprintf(fp, "\n");
//	}
//	fprintf(fp, "%%];\n");
//}
