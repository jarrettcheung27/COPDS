#ifndef _QARYLDPC_H
#define _QARYLDPC_H

#include "stdafx.h"

#include "util.h"


class ldpc_edge
{
public:
	int row_num;//row number
	int col_num;//colomn number

	double *m_alpha;//forward
	double *m_beta;//backward

	double *m_v2c;//metrics from variable node to check node
	double *m_c2v;//metrics from check node to variable node

	int pai_element;

	static int q_ary;
	static int m;

	//from function
	static double *curr;
	static double *ar;

	static double *curr_array;

	void HadamardTransform(double *data_in , double *data_out , int stage);


	ldpc_edge();
	ldpc_edge(int row , int col , const int& pai_e);
	~ldpc_edge();

	void set(int row , int col , const int& pai_e);
	void Permutation(int direction);//0 represent permutation from variable node to check node 
	                                //1 represent permutation from check node to variable node

	void Transform(int direction);//0 represent forward transform
	                              //1 represent inverse transform

	void Norm(int direction);//0 represent normilize the vector v2c
	                         //1 represent normilize the vector c2v

//According to Prof. Ma's Program
	/*
	this form a net of edges
	we can tranverse from left to right 
	                 from right to left
					 from up to bottom
					 from bottom to up
	 it form a circle.
	 the last edges point to the begining of the link list

	 left represent tranverse from left to right
	 right                         right   left
	 up                            up       bottom
	 down                          bottom   up
	*/
	ldpc_edge *up;//point to the lower one
	ldpc_edge *down;//point to the upper one
	ldpc_edge *right;//point to the left one
	ldpc_edge *left;//point to the right one

	static void init(int q_ary , int degree);//must be called first before any use of this class to specify the Order of the Finite Field
};

class Tanner_Graph
{
public:
	typedef ldpc_edge *  p_Edge;
	ldpc_edge *m_row_head;
	ldpc_edge *m_col_head;
	Tanner_Graph(string filename);
	~Tanner_Graph();
	Tanner_Graph(){}
private:
	//typedef ldpc_edge LEdge;
public:
	void Construct_Tanner_Graph(string filename);
public:
	int check_num;
	int variable_num;
	int q_ary;
};


#endif