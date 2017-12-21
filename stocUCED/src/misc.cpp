#include "misc.hpp"

/****************************************************************************
 * safeGetline
 * - Works the same as getline, however, can handle issues where the end of
 * line tokens might be either '\n', '\r', or '\n\r'.
 * - Source: https://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
 *****************************************************************************/
istream& safeGetline(istream& is, string& t)
{
	t.clear();

	// The characters in the stream are read one-by-one using a std::streambuf.
	// That is faster than reading them one-by-one using the std::istream.
	// Code that uses streambuf this way must be guarded by a sentry object.
	// The sentry object performs various tasks,
	// such as thread synchronization and updating the stream state.

	std::istream::sentry se(is, true);
	std::streambuf* sb = is.rdbuf();

	for(;;) {
		int c = sb->sbumpc();
		switch (c) {
		case '\n':
			return is;
		case '\r':
			if(sb->sgetc() == '\n')
				sb->sbumpc();
			return is;
		case std::streambuf::traits_type::eof():
						// Also handle the case when the last line has no line ending
						if(t.empty())
							is.setstate(std::ios::eofbit);
		return is;
		default:
			t += (char)c;
		}
	}
}

/* The subroutines lists all the sub-directories of a directory */
int getDirs (string dir, vector<string> &subdirs) {
	DIR *dp;
	struct dirent *dirp;

	if((dp  = opendir(dir.c_str())) == NULL) {
		cout << "Error opening directory : " << dir << endl;
		return 1;
	}

	while ((dirp = readdir(dp)) != NULL) {
		if ( dirp->d_type == DT_DIR && strcmp(dirp->d_name, ".") && strcmp(dirp->d_name, "..") )
			subdirs.push_back(string(dirp->d_name));
	}

	closedir(dp);
	return 0;
}//END getDirs()

/* The subroutines lists all the files in a directory */
int getFiles (string dir, vector<string> &files) {
	DIR *dp;
	struct dirent *dirp;

	if((dp  = opendir(dir.c_str())) == NULL) {
		cout << "Error opening directory : " << dir << endl;
		return 1;
	}

	while ((dirp = readdir(dp)) != NULL) {
		if ( dirp->d_type != DT_DIR && dirp->d_name[0] != '.')
			files.push_back(string(dirp->d_name));
	}

	closedir(dp);

	return 0;
}//genFiles()

/* The subroutines lists all the contents of a directory */
int getContents (string dir, vector<string> &contents) {
	DIR *dp;
	struct dirent *dirp;

	if((dp  = opendir(dir.c_str())) == NULL) {
		cout << "Error opening directory : " << dir << endl;
		return 1;
	}

	while ((dirp = readdir(dp)) != NULL) {
		if ( strcmp(dirp->d_name, ".") && strcmp(dirp->d_name, "..") )
			contents.push_back(string(dirp->d_name));
	}

	closedir(dp);
	return 0;
}//END getDirs()

/* The subroutine splits the line of type string along the delimiters into a vector of shorter strings */
vector<string> splitString(string &line, char delimiter) {

	stringstream ss(line);
	string item;
	vector<string> tokens;
	while (getline(ss, item, delimiter)) {
		tokens.push_back(item);
	}
	return tokens;
}//END splitString()

bool open_file (ifstream &fptr, string filename) {

	fptr.open( filename.c_str() );
	if (fptr.fail()) {
		cout << "Error opening file: " << filename << endl;
		return false;
	}
	return true;
}

bool open_file (ofstream &fptr, string filename) {
	
	fptr.open( filename.c_str() );
	if (fptr.fail()) {
		cout << "Error opening file: " << filename << endl;
		return false;
	}
	return true;
}

void move_cursor (ifstream &file, char character) {
	char x;
	while (!file.eof())	{
		file.get(x);
		if(x == character)	break;
	}
}

void disp_vector (vector<int> &x)
{
	for (unsigned int i=0; i<x.size(); i++)
		cout << fixed << setprecision(0) << x[i] << " ";
	cout << endl;
}

void disp_vector (vector<double> &x)
{
	if (x.size() == 0)  {
		cout << "Empty\n";
		return;
	}

	for (unsigned int i=0; i<x.size(); i++)
		cout << fixed << setprecision(2) << x[i] << " ";
	cout << endl;
}

void disp_vector (vector<bool> &x)
{
	if (x.size() == 0)  {
		cout << "Empty\n";
		return;
	}

	for (unsigned int i=0; i<x.size(); i++)
		cout << x[i] << " ";
	cout << endl;
}

void disp_matrix (vector< vector<double> > &x)
{
	for (unsigned int i=0; i<x.size(); i++)
	{
		for (unsigned int j=0; j<x[0].size(); j++)
			cout << fixed << setprecision(3) << setw(6) << x[i][j] << " ";
		cout << endl;
	}
	cout << endl;
}

void disp_matrix (vector< vector<bool> > &x)
{
	for (unsigned int i=0; i<x.size(); i++)
	{
		for (unsigned int j=0; j<x[0].size(); j++)
			cout << fixed << setw(1) << x[i][j] << ", ";
		cout << endl;
	}
	cout << endl;
}

void disp_matrix (vector< vector<int> > &x)
{
	for (unsigned int i=0; i<x.size(); i++)
	{
		for (unsigned int j=0; j<x[0].size(); j++)
			cout << fixed << setprecision(0) << setw(3) << x[i][j] << ", ";
		cout << endl;
	}
	cout << endl;
}

void print_matrix (ofstream &output, vector < vector<double> > &x, char delimiter, unsigned short precision)
{
	for (unsigned int i=0; i<x.size(); i++) {
		for (unsigned int j=0; j<x[0].size(); j++) {
			output << setprecision(precision) << fixed << x[i][j] << delimiter;
		}
		output << endl;
	}
}

void print_matrix (ofstream &output, vector < vector<bool> > &x, char delimiter)
{
	for (unsigned int i=0; i<x.size(); i++) {
		for (unsigned int j=0; j<x[0].size(); j++) {
			output << x[i][j] << delimiter;
		}
		output << endl;
	}
}

void print_matrix (ofstream &output, vector < vector<int> > &x, char delimiter)
{
	for (unsigned int i=0; i<x.size(); i++) {
		for (unsigned int j=0; j<x[0].size(); j++) {
			output << x[i][j] << delimiter;
		}
		output << endl;
	}
}

void print_vector (ofstream &output, vector <double> &x, char delimiter, unsigned short precision)
{
	for (unsigned int i=0; i<x.size(); i++)
		output << setprecision(precision) << fixed << x[i] << delimiter;
	output << endl;
}

void print_vector (ofstream &output, vector <bool> &x, char delimiter)
{
	for (unsigned int i=0; i<x.size(); i++)
		output << x[i] << delimiter;
	output << endl;
}

void print_vector (ofstream &output, vector <int> &x, char delimiter)
{
	for (unsigned int i=0; i<x.size(); i++)
		output << x[i] << delimiter;
	output << endl;
}

void resize_matrix(vector< vector<bool> > &mat, int rows, int cols)
{
	mat.resize(rows);
	for (unsigned int j=0; j<mat.size(); j++)
		mat[j].resize(cols);
}

void resize_matrix(vector< vector<int> > &mat, int rows, int cols)
{
	mat.resize(rows);
	for (unsigned int j=0; j<mat.size(); j++)
		mat[j].resize(cols);
}

void resize_matrix(vector< vector<double> > &mat, int rows, int cols)
{
	mat.resize(rows);
	for (unsigned int j=0; j<mat.size(); j++)
		mat[j].resize(cols);
}

vector< vector<int> > create_int_matrix(int rows, int cols)
				{
	vector< vector<int> > new_matrix (rows, vector<int> (cols) );
	return new_matrix;
				}

vector< vector< vector<int> > > create_int_quad(int size1, int size2, int size3)
				{
	vector< vector< vector<int> > > new_matrix (size1, vector< vector<int> > (size2, vector<int> (size3)));
	return new_matrix;
				}

vector< vector<double> > create_dec_matrix(int rows, int cols)
				{
	vector< vector<double> > new_matrix (rows, vector<double> (cols) );
	return new_matrix;
				}

vector< vector< vector<bool> > > create_boolean_quad (int size1, int size2, int size3)
				{
	vector< vector< vector<bool> > > new_matrix (size1, vector< vector<bool> > (size2, vector<bool> (size3)));
	return new_matrix;
				}

vector< vector< vector<bool> > > create_boolean_quad (int size1, int size2, int size3, bool default_val)
				{
	vector< vector< vector<bool> > > new_matrix (size1, vector< vector<bool> > (size2, vector<bool> (size3, default_val)));
	return new_matrix;
				}

struct int_pair {
	int value;
	int index;
};

struct double_pair {
	double value;
	int index;
};

bool compare_int_ascend (int_pair i, int_pair j) { return i.value < j.value; };
bool compare_double_ascend (double_pair i, double_pair j) { return i.value < j.value; };
bool compare_int_descend (int_pair i, int_pair j) { return i.value > j.value; };
bool compare_double_descend (double_pair i, double_pair j) { return i.value > j.value; };

vector<int> sort_this (vector<int> &x, bool ascend)
				{
	//returns the indices of jobs after sorting (indexing starts with 0)
	vector<int_pair> myVec (x.size());
	for (unsigned int i=0; i<x.size(); i++)
	{
		myVec[i].value = x[i];
		myVec[i].index = i;
	}

	if (ascend)
		sort(myVec.begin(), myVec.end(), compare_int_ascend);
	else
		sort(myVec.begin(), myVec.end(), compare_int_descend);

	vector<int> indices (x.size());
	for (unsigned int i=0; i<x.size(); i++)
	{
		x[i] = myVec[i].value;
		indices[i] = myVec[i].index;
	}

	return indices;
				}

vector<int> sort_this (vector<double> &x, bool ascend)
				{
	//returns the indices of jobs after sorting (indexing starts with 0)
	vector<double_pair> myVec (x.size());
	for (unsigned int i=0; i<x.size(); i++)
	{
		myVec[i].value = x[i];
		myVec[i].index = i;
	}

	if (ascend)
		sort(myVec.begin(), myVec.end(), compare_double_ascend);
	else
		sort(myVec.begin(), myVec.end(), compare_double_descend);

	vector<int> indices (x.size());
	for (unsigned int i=0; i<x.size(); i++)
	{
		x[i] = myVec[i].value;
		indices[i] = myVec[i].index;
	}

	return indices;
				}

#include <stdio.h>
#include <string>

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
// Reference: http://stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c
const std::string getCurrentDateTime() {
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);
	// Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
	// for more information about date/time format
	strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);

	return buf;
}


// Get the wall time and cpu time
// Reference: http://stackoverflow.com/questions/17432502/how-can-i-measure-cpu-time-and-wall-clock-time-on-both-linux-windows

#ifdef _WIN32
#include <Windows.h>
double get_wall_time(){
	LARGE_INTEGER time,freq;
	if (!QueryPerformanceFrequency(&freq)){
		//  Handle error
		return 0;
	}
	if (!QueryPerformanceCounter(&time)){
		//  Handle error
		return 0;
	}
	return (double)time.QuadPart / freq.QuadPart;
}
double get_cpu_time(){
	FILETIME a,b,c,d;
	if (GetProcessTimes(GetCurrentProcess(),&a,&b,&c,&d) != 0){
		//  Returns total user time.
		//  Can be tweaked to include kernel times as well.
		return
				(double)(d.dwLowDateTime |
						((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
	}else{
		//  Handle error
		return 0;
	}
}

//  Posix/Linux
#else
#include <time.h>
#include <sys/time.h>
double get_wall_time(){
	struct timeval time;
	if (gettimeofday(&time,NULL)){
		//  Handle error
		return 0;
	}
	return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
	return (double)clock() / CLOCKS_PER_SEC;
}
#endif
