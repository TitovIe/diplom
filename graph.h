#pragma once

#include <iomanip>
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
    explicit Sample(mt19937& gen);
    vector<int>& Get_sigma_vector();
    void Print_sample();
private:
    vector<int> sigma_vector;
    const pair<int, int> sigma_value = {-1, 1};
};

class Graph {
public:
    explicit Graph(mt19937& gen);
    vector<Sample>& Get_m_samples();
    vector<vector<double>>& Get_Ji_vector();
    void Gen_Ji_start_vector(mt19937& gen);
    void Print_sigma();
    void Print_Ji();

    static const int M = 5000, N = 25, beta = 1;
    static constexpr double alfa = 0.4;
private:
    vector<vector<double>> Ji_vector;
    vector<Sample> M_samples;
};

ostream& operator << (ostream& os, const vector<double>& v);