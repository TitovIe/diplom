#include "graph.h"
#include <fstream>

Graph::Graph(mt19937& gen){
    Gen_Ji_start_vector(gen);
    Gen_sigma_vector(gen);
}

void Graph::Gen_sigma_vector(mt19937 &gen) {
    //Glauber(sigmai_vector, gen);
    int number_conff, spin, number_conff_all = 0;
    string line;
    stringstream ls;
    ifstream file_spin("/home/titov/CLionProjects/Diplom_ver2/samples.txt");
    if (file_spin.is_open())
    {
        while (getline(file_spin, line))
        {
            ls << line;
            ls >> number_conff;
            number_conff_all += number_conff;
            vector<int> spins;
            while(!ls.eof()){
                ls.ignore(1);
                ls >> spin;
                spins.push_back(spin);
            }
            ls.clear();
            sigmai_vector.push_back({number_conff, spins});
        }
        sigmai_vector.pop_back();
    }
    file_spin.close();
}

/*Генерируем стартовые векторы Jij и pij*/
void Graph::Gen_Ji_start_vector(mt19937& gen){
    uniform_real_distribution<double>
            urd_right(Graph::alfa, Graph::beta);

    uniform_real_distribution<double>
            urd_left(-Graph::beta, -Graph::alfa);

    uniform_int_distribution<int> urd_false_true(0, 1);

    Ji_vector.reserve(Graph::N);
    constraint_vector.reserve(Graph::N);
    
    vector<double> Jij_start_vector(Graph::N - 1);
    vector<double> constraint_start_vector(Graph::N - 1);
    
    for(int k = 0; k < Graph::N; k++){
        Ji_vector.push_back(Jij_start_vector);
        constraint_vector.push_back(constraint_start_vector);
    }

    for(int i = 0; i < Graph::N; i++) {
        for (int j = i; j < Graph::N - 1; j++) {
            if (urd_false_true(gen) == 0) {
                Ji_vector[i][j] = urd_left(gen);
            } else Ji_vector[i][j] = urd_right(gen);
            Ji_vector[j + 1][i] = Ji_vector[i][j];
            constraint_vector[i][j] = 1 / lambda * Ji_vector[i][j];
            constraint_vector[j + 1][i] = constraint_vector[i][j];
        }
    }
}

//Генерация M выборок спинов методом метрополиса
void Glauber(vector<vector<int>>& samples,
                const vector<vector<double>>& Jij_true, mt19937 &gen){

    uniform_int_distribution<int> urd_false_true(0, 1);
    uniform_real_distribution<double> probability(0, 1);
    uniform_int_distribution<int> choise_spin(0, Graph::N - 1);
    
    const int number_step = 100;
    int number_spin;
    double E;
    
    // Начальная конфигурация рандомно заполняется
    vector<int> sigma_vector;
    for (int j = 0; j < Graph::N; j++) {
        if (urd_false_true(gen) == 0)
            sigma_vector.push_back(-1);
        else
            sigma_vector.push_back(1);
    }

    /*Генерируем M выборок. Каждая новая конфигурация - через 100 итераций.
     * Каждая итерацию выбираем какой то узел и с вероятностью 50% делаем следующие операции.
     * Считаем энергию системы с текущим значением спина. Переворачиваем его 
     * и снова считаем. Если deltaЕ < 0, оставляем перевернутым. В противном случае
     * с вероятностью ext(-beta * E_delta) оставляем перевернутым. */
    for(int j = 0; j < Graph::M; j++) {
        for (int i = 0; i < number_step; i++){
            
            number_spin = choise_spin(gen);
            E = Jij_sum_calc(Jij_true, sigma_vector, number_spin);
            
            if(probability(gen) > 1 / (1 + exp(-2 * E)))
                sigma_vector[number_spin] = 1;
            else
                sigma_vector[number_spin] = -1;
        }
        samples.push_back(sigma_vector);
    }
}

//Считаем энергию системы
double Jij_sum_calc(const vector<vector<double>>& Jij_vector,
        const vector<int>& sigma_vector, int i){
    double Jij_sum = 0;
    for(int j = 0; j < Jij_vector.size(); i++){
        if(j == i)
            j++;
        
            Jij_sum += -Jij_vector[i][j]
                        * sigma_vector[j];
    }
    return Jij_sum;
}


vector<pair<int, vector<int>>>& Graph::Get_m_samples() {
    return sigmai_vector;
}

vector<vector<double>>& Graph::Get_Ji_vector() {
    return Ji_vector;
}

void Graph::Print_Ji() {
    for(const auto& i : Get_Ji_vector()){
        cout << i << endl;
    }
    cout << endl;
}

void Graph::Print_sigma_all() {
    for(const auto& sample : sigmai_vector)
        Print_sigma_sample(sample);
}

void Graph::Print_sigma_sample(const pair<int, vector<int>>& sample) {
    cout << sample.first << " ";
    for(const auto& i : sample.second) {
        cout << i << " ";
    }
    cout << endl;
}

vector<vector<double>>& Graph::Get_constraint_vector() {
    return constraint_vector;
}

ostream& operator << (ostream& os, const vector<double>& v) {
    for (const auto &j : v) {
        os.width(10);
        os << j << " ";
    }
    return os;
}

Graph& GetClobalObject()
{
    static mt19937 gen(0);
    static Graph graph(gen);
    return graph;
}

/*Приводим к нулю Jij такие,которые < alfa/2 */
void Jij_to_zero(Graph& g){
    for(int i = 0; i < Graph::N; i++){
        for(int j = i; j < Graph::N - 1; j++){
            g.Get_Ji_vector()[i][j] = (g.Get_Ji_vector()[i][j]
                                       + g.Get_Ji_vector()[j + 1][i]) / 2;

            if(g.Get_Ji_vector()[i][j] > -Graph::alfa / 2
               && g.Get_Ji_vector()[i][j] < Graph::alfa / 2)
                g.Get_Ji_vector()[i][j] = 0;

            g.Get_Ji_vector()[j + 1][i] = g.Get_Ji_vector()[i][j];
        }
    }
}
