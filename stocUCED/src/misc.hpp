#ifndef _MISC_H
#define _MISC_H

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include <stdlib.h>
#include <map>
#include <climits>

#define MISCPRECISION 1e-4

using namespace std;

// Commonly used functions

/****************************************************************************
 * numToStr
 * - Converts numbers to string
 * - Source: https://stackoverflow.com/questions/5590381/easiest-way-to-convert-int-to-string-in-c
 *****************************************************************************/
template <typename T>
string numToStr (T Number) {
	ostringstream ss;
	ss << Number;
	return ss.str();
}

istream& safeGetline(istream& is, string& t);

int getDirs (string dir, vector<string> &subdirs);
int getFiles (string dir, vector<string> &files);
int getContents (string dir, vector<string> &contents);

vector<string> splitString(string &line, char delimiter);

bool open_file (ifstream &file, string filename);
void move_cursor (ifstream &file, char character);

void disp_vector (vector<int> &x);
void disp_vector (vector<double> &x);
void disp_vector (vector<bool> &x);

void print_vector (vector <int> &x, ofstream &output);
void print_vector (vector <bool> &x, ofstream &output);
void print_vector (vector <double> &x, ofstream &output);

void disp_matrix (vector< vector<double> > &x);
void disp_matrix (vector< vector<bool> > &x);
void disp_matrix (vector< vector<int> > &x);

void print_matrix (vector < vector<double> > &x, ofstream &output);
void print_matrix (vector < vector<bool> > &x, ofstream &output);

void alloc_mem();

void resize_matrix(vector< vector<int> > &mat, int rows, int cols);
void resize_matrix(vector< vector<bool> > &mat, int rows, int cols);
void resize_matrix(vector< vector<double> > &mat, int rows, int cols);

vector< vector<int> > create_int_matrix(int rows, int cols);

vector< vector<double> > create_dec_matrix(int rows, int cols);

vector< vector< vector<int> > > create_int_matrix(int size1, int size2, int size3);

vector< vector< vector<bool> > > create_boolean_quad (int size1, int size2, int size3);
vector< vector< vector<bool> > > create_boolean_quad (int size1, int size2, int size3, bool default_val);

vector< vector< vector<int> > > create_int_quad(int size1, int size2, int size3);

vector<int> sort_this (vector<int> &x, bool ascend);
vector<int> sort_this (vector<double> &x, bool ascend);

template <class object>
bool isFractional (object x)
{
	return ( abs(x - round(x)) > MISCPRECISION );
}

template <class object>
bool isFractional (object x, double prec)
{
	return ( abs(x - round(x)) > prec );
}

template <class object>
object sum (vector<object> &x)
{
	object res = 0;
	for (unsigned int i=0; i<x.size(); i++)	res += x[i];
	return res;
}

template <class object>
int count (vector<object> &x, object obj)
{
	int count = 0;
	for (unsigned int i=0; i<x.size(); i++)	
		if (x[i] == obj)
			count += x[i];
	return count;
}

template <class object> 
vector<object> sum (vector< vector<object> > &x, bool rows)
{
	vector<object> res ((!rows)*(x.size()) + (rows)*(x[0].size()), 0);
	for (unsigned int i=0; i<x.size(); i++)	
		for (unsigned int j=0; j<x[0].size(); j++)
		{
			if (rows)	res[j]	+= x[i][j];
			else		res[i]  += x[i][j];
		}
	return res;
}

template <class object1, class object2>
double scal_prod (vector<object1> &x, vector<object2> &y)
{
	double res = 0;
	for (unsigned int i=0; i<x.size(); i++)	res += x[i] * y[i];
	return res;
}

template <class object>
object maximum (vector<object> &x)
{
	object temp = -999999;
	for (int i=0; i<x.size(); i++) {
		if (x[i] > temp) temp = x[i];
	}
	return temp;
}

template <class object>
vector<object> Subtract (vector<object> x, vector<object> y)
{
	// x - y
	vector<object> difference (x.size());
	for (int i=0; i<x.size(); i++) 	difference[i] = x[i] - y[i];
	return difference;
}

const std::string currentDateTime();

double get_cpu_time();
double get_wall_time();

#endif
