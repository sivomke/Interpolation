#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
//#include <math.h>
#include <vector>
//#define M_PI

using namespace std;

const double pi = 3.14159265358979323846;
const int nodes = 3;
const double left_edge = -1 * pi,
right_edge = pi;


const double a = 2, w = 1, b = 1;

const double max = b + abs(a), min = b - abs(a);

double f(double x)
{
	return a * sin(w * x) + b;
}

void get_nodes_equidist(vector<double>& x_nodes, double lleft, double rright, double num)
{
	x_nodes.resize(0);
	double tmp_node = lleft;
	double step = (rright - lleft) / (num - 1);
	for (int i = 0; i < num; ++i)
	{
		x_nodes.push_back(tmp_node);
		tmp_node += step;
	}
}

void get_nodes_chebyshev(vector<double>& x_nodes, double lleft, double rright, double num)
{
	x_nodes.resize(0);
	for (int i = 0; i < num; ++i)
	{
		x_nodes.push_back((rright + lleft) / 2 + (rright - lleft) / 2 * cos(pi * ((2 * i + 1) / (2 * num))));
	}
}

void get_nodes_f_vals(const vector<double> xnodes, vector<double>& fnodes)
{
	int nums = xnodes.size();
	fnodes.resize(0);
	for (int i = 0; i < nums; ++i)
	{
		fnodes.push_back(f(xnodes[i]));
	}
}

void compute_differences(const vector<double> xnodes, const vector<double> fnodes, vector<double>& diffs)
{
	int nums = fnodes.size();
	if (nums != xnodes.size())
	{
		cout << "\n Arrays have different sizes!\n";
		return;
	}
	diffs.resize(0);
	vector<double> diffs_k_order = fnodes;
	for (int i = 0; i < nums - 1; ++i)
	{
		for (int j = 0; j < nums - 1 - i; ++j)
		{
			diffs_k_order[(nums - 1) - j] = (diffs_k_order[(nums - 1) - j] - diffs_k_order[(nums - 1) - j - 1]) / (xnodes[nums - j - 1] - xnodes[nums - j - 1 - i - 1]);
		}
	}
	diffs = diffs_k_order;
}

double mypol(const vector<double> diffs, const vector<double> xnodes, const double x) //Newton`s pol
{
	double res = 0;
	double tmp = 1;
	int n = xnodes.size();
	for (int i = 0; i < n; ++i)
	{
		res += diffs[i] * tmp;
		tmp *= (x - xnodes[i]);
	}
	return res;
}

void inverse_int() {
	//double left_end = -abs(pi / w), right_end = abs(pi / w);
	double left_end = 0.3, right_end = 0.4;
	double target = min + 2. * (max - min) / 3.;
	vector <double> x_nodes;
	get_nodes_chebyshev(x_nodes, left_end, right_end, nodes);
	vector <double> ch_val;
	get_nodes_f_vals(x_nodes, ch_val);
	vector <double> diffs;
	compute_differences(ch_val, x_nodes, diffs);
	double res = mypol(diffs, ch_val, target);
	cout << "Target: " << target << endl;
	cout << "Approximate solution of f(x)=target: " << res << endl;
	cout << "Approximate target: " << f(res) << endl;
	cout << "Diff between actual and computed value: " << abs(f(res) - target) << endl;
}

void print_to_file(const vector<double> x_nodes, const vector<double> diffs)
{
	int nums = x_nodes.size();
	if (diffs.size() != nums)
	{
		cout << "\n Arrays have different sizes (in file printing method)!\n";
		return;
	}
	ostringstream ss;
	ofstream myfile;
	myfile.open("C:/NumAnalysis/Interpolation/lab.txt");
	ss << nodes << "\n";
	for (int i = 0; i < nums; ++i)
	{
		ss << x_nodes[i] << " ";
	}
	ss << "\n";
	for (int i = 0; i < nums; ++i)
	{
		ss << diffs[i] << " ";
	}
	myfile << ss.str();
	myfile.close();
}

int main()
{
	cout.precision(8);
	/*
	vector<double> xnds, fnds, dffs;
	get_nodes_equidist(xnds, left_edge, right_edge, nodes);
	get_nodes_f_vals(xnds, fnds);
	compute_differences(xnds, fnds, dffs);
	print_to_file(xnds, dffs);
	*/

	vector<double> xnds, fnds, dffs;
	get_nodes_chebyshev(xnds, left_edge, right_edge, nodes);
	get_nodes_f_vals(xnds, fnds);
	compute_differences(xnds, fnds, dffs);
	print_to_file(xnds, dffs);

	inverse_int();
	return 0;
}
