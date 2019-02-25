#pragma once
#include <vector>
#include "graph.h"

using namespace std;

double Jij_sigma_calc(Sample& s, int i, int j, const vector<double>& Ji_vector);

vector<double> Ln_Si_grad_calc(const vector<double>& S_grad, double Si);

void Arg_min_calc(Graph& g, double eps,
                  vector<vector<double>>& Ji_start_vector, double lambda);

pair<double, vector<double>> Si_calc(Graph& g, int i, const vector<double>& Ji);

double Func(double Ln_Si, double lambda, double Ji_norm);

double Ji_norm_calc(const vector<double>& Ji, int i);

vector<vector<double>> Gen_Ji_start_vector();

void Graph_reconst_calc(Graph& g, double eps);

vector<double> Gradient();