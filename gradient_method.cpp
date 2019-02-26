#include "gradient_method.h"

/*Оператор, реализующий сумму элементов вектора1 и вектора2 */
template <typename T>
vector<T> operator+ (const vector<T>& v1,
        const vector<T>& v2){
    vector<T> v_new;
    for(size_t i = 0; i < v1.size(); i++) {
        v_new.push_back(v1[i] + v2[i]);
    }

    return v_new;
}

/*Оператор, реализующий разность элементов вектора1 и вектора2 */
template <typename T>
vector<T> operator- (const vector<T>& v1,
                     const vector<T>& v2){
    vector<T> v_new;
    for(size_t i = 0; i < v1.size(); i++) {
        v_new.push_back(v1[i] - v2[i]);
    }

    return v_new;
}

/*Оператор, реализующий умножение элементов вектора на число */
template <typename T, typename K>
vector<double> operator* (const vector<T>& v, K number){
    vector<double> v_new;
    for(const auto& i : v){
        v_new.push_back(i * number);
    }
    return v_new;
}

/*Оператор, реализующий деление элементов вектора на число */
template <typename T, typename K>
vector<double> operator/ (const vector<T>& v, K number){
    vector<double> v_new;
    for(const auto& i : v){
        v_new.push_back(i / number);
    }
    return v_new;
}

/*Вычсление градиента от натур. логорифма Si для i узла */
vector<double> Ln_Si_grad_calc(const vector<double>& Si_grad,
                              double Si) {
    return Si_grad / Si;
}

/*Вычисление нормы Ji*/
double Ji_norm_calc(const vector<double>& Ji){
    double Ji_norm = 0;
    for(int j = 0; j < Ji.size(); j++){
            Ji_norm += abs(j);
    }
    return Ji_norm;
}

/*Вычисление функции,минимум которой ищем*/
double Func(double Ln_Si, double lambda, double Ji_norm) {
    return Ln_Si + lambda * Ji_norm;
}

/*Вычисление градиента и логарифма Si*/
pair<double, vector<double>>
Si_calc(Graph& g, int i, const vector<double>& Ji_vector) {
    double s_sum = 0;
    double Jij_sigma_sum = 0;
    double Si = 0;
    double Ln_Si = 0;
    vector<double> Si_grad(Graph::N - 1);
    vector<double> Ln_Si_grad(Graph::N - 1);

    for (auto &m : g.Get_m_samples()) {
        for (int j = 0, k = 0; j < Graph::N - 1; j++, k++) {
            if (k == i)
                k++;

            /*Считаем сумму по j: Sum(-Jij * sigmai * sigmaj)
            * для одной выборки */
            Jij_sigma_sum += -Ji_vector[j]
                    * m.Get_sigma_vector()[k]
                    * m.Get_sigma_vector()[i];
        }
        /*Считаем сумму для M выборок, чтобы найти среднее.
         * Также находим градиент от Si,
         * в цикле - сумму градиенов для М выборок */
        s_sum += exp(Jij_sigma_sum);

        for(int j = 0, k = 0; j < Graph::N - 1; j++, k++){
            if(k == i)
                k++;

            Si_grad[j] += m.Get_sigma_vector()[k]
                    * exp(Jij_sigma_sum);
        }

        Jij_sigma_sum = 0;
    }
    /*Находим средний градиент по M выборкам
     * а также среднее Si */
    Si_grad = Si_grad * (-1) / Graph::M;
    Si = s_sum / Graph::M;

    /*Вычисляем градиент от натур. логорифма и
     * находим его значение*/
    Ln_Si = log(Si);
    Ln_Si_grad = Ln_Si_grad_calc(Si_grad, Si);

    return {Ln_Si, Ln_Si_grad};
}

/*Нахождение минимума с помощью градиентного спуска*/
void Arg_min_calc(Graph& g, double eps,
        vector<vector<double>>& Ji_start_vector, double lambda) {
    pair<double, vector<double>> F_prev, F_new;
    vector<double>  Func_grad_prev, Func_grad_new, Ji1, Ji0;

    double h = 1, Func_value_prev = 0, Func_value_new = 0;

    for(int i = 0; i < Graph::N; i++) {
        Ji0 = Ji_start_vector[i];

        do {
            Ji1 = Ji0;
            F_prev = Si_calc(g, i, Ji0);
            Func_value_prev = Func(F_prev.first, lambda,
                                   Ji_norm_calc(Ji0));
            Func_grad_prev = F_prev.second;

            for(int j = 0; j < Graph::N - 1; j++) {
                if(Ji1[j] != 0)
                    Ji1[j] = Ji0[j] - Func_grad_prev[j] * h;
            }

            F_new = Si_calc(g, i, Ji1);
            Func_value_new = Func(F_new.first, lambda,
                                  Ji_norm_calc(Ji1));
            Func_grad_new = F_new.second;

            if (Func_value_new > Func_value_prev)
                h /= 2;
            else {
                h *= 2;
                Ji0 = Ji1;
            }
        } while (abs(Func_value_new - Func_value_prev) > eps);
        Ji_start_vector[i] = Ji0;
        h = 1;
    }

    g.Get_Ji_vector() = Ji_start_vector;
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

/*Строим руконструированный граф*/
void Graph_reconst_calc(Graph& g, double eps){
    const double c_lambda = 0.8;
    double lambda =
            c_lambda * sqrt(log(pow(Graph::N, 2) / eps) / Graph::M);

    Arg_min_calc(g, eps, g.Get_Ji_vector(), lambda);
    Jij_to_zero(g);

    /*Повторяем тоже самое для lambda = 0,
     * чтобы получить искомый граф*/
    lambda = 0;
    Arg_min_calc(g, eps, g.Get_Ji_vector(), lambda);
    Jij_to_zero(g);
}

/*
vector<double> Gradient(){
    vector<double> start = {4,5};
    vector<double> finish = start;
    double Func_value_new, Func_value_prev, eps = 0.00005, h = 0.1;


    do{
        finish = start - start * 2 * h;

        Func_value_prev = pow((start[0] + start[1]), 2);
        Func_value_new = pow((finish[0] + finish[1]), 2);
        if (Func_value_new >= Func_value_prev)
            h /= 2;
        else {
            start = finish;
        }
    } while (abs(Func_value_new - Func_value_prev) > eps);
    return finish;
}
*/











