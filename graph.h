#pragma once
#include <map>
#include <memory>
#include <set>
#include <vector>
#include <random>
#include <ctime>
#include <iostream>

using namespace std;

class Sample {
public:
    Sample(mt19937& gen, random_device& rd);
    void Print_sample();
    vector<int>& Get_sigma_vector();
private:
    vector<int> sigma_vector;
    const pair<int, int> sigma_value = {-1, 1};
};

class Graph {
public:
    Graph(mt19937& gen, random_device& rd);
    void Print_graph();
    void Print_Ji();
    vector<Sample>& Get_m_samples();
    vector<vector<double>>& Get_Ji_vector();

    static const int M = 5000, N = 25, beta = 1;
    static constexpr double alfa = 0.4;
private:
    vector<vector<double>> Ji_vector;
    vector<Sample> M_samples;
};

ostream& operator << (ostream& os, const vector<double>& v);