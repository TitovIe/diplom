#include "graph.h"

Sample::Sample(mt19937& gen, random_device& rd) {
    uniform_int_distribution<int> urd_false_true(0, 1);
    for (int j = 0; j < Graph::N; j++) {
        gen.seed(rd());
        if (urd_false_true(gen) == 0)
            sigma_vector.push_back(sigma_value.first);
        else
            sigma_vector.push_back(sigma_value.second);
    }
}

vector<int>& Sample::Get_sigma_vector() { return sigma_vector;}

Graph::Graph(mt19937& gen, random_device& rd){
    for (int i = 0; i < M; i++) {
        gen.seed(rd());
        Sample s = Sample(gen, rd);
        M_samples.push_back(s);
    }
}

vector<Sample>& Graph::Get_m_samples() {
    return M_samples;
}

void Graph::Print_graph(){
    for(auto& i : M_samples) {
        i.Print_sample();
        cout << endl;
    }
}

void Graph::Print_Ji() {
    for(const auto& i : Get_Ji_vector()){
        cout << i << endl;
    }
}

vector<vector<double>>& Graph::Get_Ji_vector() {
    return Ji_vector;
}

ostream& operator << (ostream& os, const vector<double>& v) {
    for (const auto& j : v) {
        os << j << " ";
    }
    return os;
}

void Sample::Print_sample() {
    for(int i = 0; i < Graph::N; i++){
        cout << Get_sigma_vector()[i] << " ";
    }
}