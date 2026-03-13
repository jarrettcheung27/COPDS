#include "stdafx.h"
#include "QaryLDPC.h"

extern CFiniteField2 GF;

#define SMALLPROB 1.0e-12

void ldpc_edge::init(int q_ary , int degree)
{
	ldpc_edge::q_ary = q_ary;
	ldpc_edge::m = degree;
	ldpc_edge::curr = new double[q_ary];
	ldpc_edge::ar = new double[q_ary];
	ldpc_edge::curr_array = new double[q_ary];
		//from function
}

ldpc_edge::ldpc_edge()
{
	this->m_c2v = new double[ldpc_edge::q_ary];
	this->m_v2c = new double[ldpc_edge::q_ary];

	this->m_alpha = new double[ldpc_edge::q_ary];
	this->m_beta = new double[ldpc_edge::q_ary];

	this->pai_element = int(1);
}
ldpc_edge::ldpc_edge(int row, int col , const int& pai_e)
{
	this->row_num = row;
	this->col_num = col;

	this->m_c2v = new double[ldpc_edge::q_ary];
	this->m_v2c = new double[ldpc_edge::q_ary];

	this->m_alpha = new double[ldpc_edge::q_ary];
	this->m_beta = new double[ldpc_edge::q_ary];

	this->pai_element = pai_e;
}

void ldpc_edge::set(int row, int col , const int& pai_e)
{
	this->row_num = row;
	this->col_num = col;

	this->pai_element = pai_e;
}

ldpc_edge::~ldpc_edge()
{
	delete []this->m_c2v;
	delete []this->m_v2c;
	delete []m_alpha;
	delete []m_beta;
}

int ldpc_edge::q_ary = 16;
int ldpc_edge::m = 4;
double *ldpc_edge::curr = 0;
double *ldpc_edge::ar = 0;
double *ldpc_edge::curr_array = 0;

Tanner_Graph::Tanner_Graph(std::string filename)
{
	/*this->check_num = check_num;
	this->variable_num = vari_num;
	this->m_row_head = new ldpc_edge[this->check_num];
	this->m_col_head = new ldpc_edge[this->variable_num];*/
	this->Construct_Tanner_Graph(filename);
}
Tanner_Graph::~Tanner_Graph()
{

	ldpc_edge *temp_edge;

	for (int i = 0 ; i < this->check_num ; i ++)
	{
		while (m_row_head[i].right->col_num != -1)
		{
			temp_edge = m_row_head[i].right;
			(m_row_head + i)->right = temp_edge->right;
			delete temp_edge;
		}
	}
	
	delete []this->m_row_head;
	delete []this->m_col_head;
}

void Tanner_Graph::Construct_Tanner_Graph(std::string filename)
{
	ldpc_edge *temp_edge;
	FILE *fp;
	int row_no;
	int col_no;
	int row_deg;
	int field_number;
	int degree = 0;
	int gf = 0;
	if( (fopen_s(&fp,filename.c_str(),"r")) !=0 )
	{   
		fprintf(stderr,"\n cannot open %s",filename.c_str());
		exit(1);
	}
	char temp_string[80];

	fscanf_s(fp,"%s",temp_string,80);
	fscanf_s(fp,"%d",&this->check_num);

	fscanf_s(fp,"%s",temp_string,80);
	fscanf_s(fp,"%d",&this->variable_num);

	fscanf_s(fp,"%s",temp_string,80);
	fscanf_s(fp,"%d",&this->q_ary);

	fscanf_s(fp,"%s",temp_string,80);
	fscanf_s(fp,"%d",&this->q_ary);


	degree = 1;
	while((1<<degree) < this->q_ary)
		degree ++;
	ldpc_edge::init(this->q_ary , degree);

	this->m_row_head = new ldpc_edge[this->check_num];
	this->m_col_head = new ldpc_edge[this->variable_num];

	

//begin init tanner graph
	for(int i = 0 ; i < this->check_num ; i ++)
	{
		this->m_row_head[i].col_num = -1;
		this->m_row_head[i].row_num = i;

		this->m_row_head[i].left = this->m_row_head + i;
		this->m_row_head[i].right = this->m_row_head + i;
		this->m_row_head[i].up = this->m_row_head + i;
		this->m_row_head[i].down = this->m_row_head + i;
	}

	for(int i = 0 ; i < this->variable_num ; i ++)
	{
		this->m_col_head[i].col_num = i;
		this->m_col_head[i].row_num = -1;

		this->m_col_head[i].left = this->m_col_head + i;
		this->m_col_head[i].right = this->m_col_head + i;
		this->m_col_head[i].up = this->m_col_head + i;
		this->m_col_head[i].down = this->m_col_head + i;
	}
//end init of tanner graph
	fscanf_s(fp ,"%s" , temp_string , 80);
//begin construct tanner graph
	for (int i = 0 ; i < this->check_num ; i ++)
	{
		fscanf(fp, "%d %d", &row_no, &row_deg);
		for (int j = 0 ; j < row_deg ; j ++)
		{
			temp_edge = new ldpc_edge();
			fscanf(fp, "%d", &col_no);
			fscanf(fp, "%d", &field_number);

			//temp_edge->col_num = col_no;
			//temp_edge->row_num = row_no;
			temp_edge->set(row_no , col_no , field_number);

			temp_edge->right = m_row_head[i].right;//point the left one....
			m_row_head[i].right = temp_edge;//the ->right of the last elements of a row is always pointing to himself

			temp_edge->left = m_row_head + i;//point to the left most
			(temp_edge->right)->left = temp_edge;//point to the right one



			temp_edge->down = m_col_head[col_no].down;//point to the up one
			m_col_head[col_no].down = temp_edge;//point to the current one

			temp_edge->up = m_col_head + col_no;//point to the upmost
			(temp_edge->down)->up = temp_edge;//point to the lower one
		}
	}
//end construct tanner graph
}
void ldpc_edge::Permutation(int direction)
{
	int gf(0);
	int gf_curr(0);

	if(direction == 0)
		gf = this->pai_element;
	else
		gf = GF.Div(1,this->pai_element);

	for(int i = 0 ; i < this->q_ary ; i ++)
	{
		gf_curr = GF.Mult(gf, i);
		if(direction == 0)
			curr[gf_curr] = this->m_v2c[i];
		else
			curr[gf_curr] = this->m_c2v[i];
	}
	
	for(int i = 0 ; i < this->q_ary ; i ++)
	{
		if(direction == 0)
			this->m_v2c[i] = curr[i];
		else
			this->m_c2v[i] = curr[i];
	}
}
void ldpc_edge::Transform(int direction)
{
	if(direction == 0)
	{
		//for(int i = 0 ; i < this->q_ary ; i ++)
		//{
		//	std::cout<<this->m_v2c[i]<<"  ";
		//}
		//std::cout<<std::endl;
		HadamardTransform(this->m_v2c , ar , this->m);
		for(int i = 0 ; i < this->q_ary ; i ++)
		{
			this->m_v2c[i] = ar[i];
		}
	}
	else
	{
		/*for(int i = 0 ; i < this->q_ary ; i ++)
		{
			std::cout<<this->m_c2v[i]<<"  ";
		}
		std::cout<<std::endl;*/
		HadamardTransform(this->m_c2v , ar , this->m);
		for(int i = 0 ; i < this->q_ary ; i ++)
		{
			this->m_c2v[i] = ar[i];
		}
	}
}

void ldpc_edge::Norm(int direction)
{
	double curr_sum = 0.0;
	if(direction == 0)
	{
		for(int i = 0  ; i < this->q_ary ; i ++)
		{
			curr_sum += this->m_v2c[i];
		}
	}
	else
	{
		for(int i = 0  ; i < this->q_ary ; i ++)
		{
			curr_sum += this->m_c2v[i];
		}
	}
	if(direction == 0)
	{
		for(int i = 0  ; i < this->q_ary ; i ++)
		{
			this->m_v2c[i] /= curr_sum ;
			if(this->m_v2c[i] < SMALLPROB)
				this->m_v2c[i] = SMALLPROB;
			else if(this->m_v2c[i] > 1 - SMALLPROB)
				this->m_v2c[i] = 1 - SMALLPROB;
		}
	}
	else
	{
		for(int i = 0  ; i < this->q_ary ; i ++)
		{
			this->m_c2v[i] /= curr_sum;
			if(this->m_c2v[i] < SMALLPROB)
				this->m_c2v[i] = SMALLPROB;
			else if(this->m_c2v[i] > 1 - SMALLPROB)
				this->m_c2v[i] = 1 - SMALLPROB;
		}
	}
	curr_sum = 0.0;
//To make sure that the sum is euqal to 1.0.
	if(direction == 0)
	{
		for(int i = 0  ; i < this->q_ary ; i ++)
		{
			curr_sum += this->m_v2c[i];
		}
	}
	else
	{
		for(int i = 0  ; i < this->q_ary ; i ++)
		{
			curr_sum += this->m_c2v[i];
		}
	}

	if(direction == 0)
	{
		for(int i = 0  ; i < this->q_ary ; i ++)
		{
			this->m_v2c[i] /= curr_sum ;
		}
	}
	else
	{
		for(int i = 0  ; i < this->q_ary ; i ++)
		{
			this->m_c2v[i] /= curr_sum;
		}
	}

}

void ldpc_edge::HadamardTransform(double *data_in , double *data_out , int stage)
{
	int num = 1<<stage, curr_dist = 0, block = 0, block_row_num = 0;
	for(int i = 0 ; i < num ; i ++)
	{
		curr_array[i] = data_in[i];
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
				data_out[j * block_row_num + s] = curr_array[j * block_row_num + s] + curr_array[j * block_row_num + curr_dist + s];
				data_out[j * block_row_num + s + curr_dist] = curr_array[j * block_row_num + s] - curr_array[j * block_row_num + curr_dist + s];
			}
		}
		for(int j = 0 ; j < num ; j ++)
		{
			curr_array[j] = data_out[j];
		}
	}
	/*for(int j = 0 ; j < num ; j ++)
	{
		  data_out[j] /= (double)(num>>1);
	}*/
}