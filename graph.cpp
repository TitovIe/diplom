#include "graph.h"

Sample::Sample(mt19937& gen) {
    uniform_int_distribution<int> urd_false_true(0, 1);
    for (int j = 0; j < Graph::N; j++) {
        if (urd_false_true(gen) == 0)
            sigma_vector.push_back(sigma_value.first);
        else
            sigma_vector.push_back(sigma_value.second);
    }
}

vector<int>& Sample::Get_sigma_vector() { return sigma_vector;}

Graph::Graph(mt19937& gen){
    Gen_Ji_start_vector(gen);
    for (int i = 0; i < M; i++) {
        Sample s = Sample(gen);
        M_samples.push_back(s);
    }
}

/*Генерируем стартовый вектор*/
void Graph::Gen_Ji_start_vector(mt19937& gen){
    uniform_real_distribution<double>
            urd_right(Graph::alfa, Graph::beta);

    uniform_real_distribution<double>
            urd_left(-Graph::beta, -Graph::alfa);

    uniform_int_distribution<int> urd_false_true(0, 1);

    for(int i = 0; i < Graph::N; i++) {
        vector<double> Jij_start_vector;
        for (int j = 0; j < Graph::N - 1; j++) {
            if (urd_false_true(gen) == 0) {
                Jij_start_vector.push_back(urd_left(gen));
            } else Jij_start_vector.push_back(urd_right(gen));
        }
        Ji_vector.push_back(Jij_start_vector);
    }
}

vector<Sample>& Graph::Get_m_samples() {
    return M_samples;
}

void Graph::Print_sigma() {
    for(auto& i : M_samples) {
        i.Print_sample();
        cout << endl;
    }
}

void Graph::Print_Ji() {
    for(const auto& i : Get_Ji_vector()){
        cout << i << endl;
    }
    cout << endl;
}

vector<vector<double>>& Graph::Get_Ji_vector() {
    return Ji_vector;
}

ostream& operator << (ostream& os, const vector<double>& v) {
    for (const auto& j : v) {
        os.width(10);
        os << j << " ";
    }
    return os;
}

void Sample::Print_sample() {
    for(int i = 0; i < Graph::N; i++){
        cout << Get_sigma_vector()[i] << " ";
    }
}