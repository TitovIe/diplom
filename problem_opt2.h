#include <ifopt/variable_set.h>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>
#include <vector>
#include <graph.h>

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
template <typename T>
vector<double> operator* (const vector<T>& v, int number){
    vector<double> v_new;
    for(const auto& i : v){
        v_new.push_back(i * number);
    }
    return v_new;
}

/*Оператор, реализующий деление элементов вектора на число */
template <typename T>
vector<double> operator/ (const vector<T>& v, int number) {
    vector<double> v_new;
    for(const auto& i : v){
        v_new.push_back(i / number);
    }
    return v_new;
}

namespace ifopt {
    using namespace Eigen;

    //Получаем текущий объект(в данном случае наш граф)
    Graph graph = GetClobalObject();
    
    //Здесь считаем значение функции S для каждого вектора
    pair<double, vector<double>> S;
    pair<double, vector<double>> Si_calc(VectorXd& Ji){
        double s_sum = 0;
        double Jij_sigma_sum = 0;
        double Si = 0;
        double Ln_Si = 0;
        vector<double> Si_grad(Graph::N - 1);
        vector<double> Ln_Si_grad(Graph::N - 1);

        for (const auto& sample : graph.Get_m_samples()) {
            for (int j = 0, k = 0; j < Graph::N - 1; j++, k++) {
                if (k == global_counter)
                    k++;

                /*Считаем сумму по j: Sum(-Jij * sigmai * sigmaj)
                * для одной выборки */
                Jij_sigma_sum += -Ji(j)
                                 * sample.second[k]
                                 * sample.second[global_counter];
            }
            /*Считаем сумму для M выборок, чтобы найти среднее.
             * Также находим градиент от Si,
             * в цикле - сумму градиенов для М выборок */
            s_sum += exp(Jij_sigma_sum) * sample.first;

            for(int j = 0, k = 0; j < Graph::N - 1; j++, k++){
                if(k == global_counter)
                    k++;

                Si_grad[j] += sample.second[k]
                              * exp(Jij_sigma_sum) * sample.first;
            }

            Jij_sigma_sum = 0;
        }
        /*Находим средний градиент по M выборкам
         * а также среднее Si */
        Si_grad = Si_grad * (-1) / Graph::M;
        Si = s_sum / Graph::M;
        Ln_Si = log(Si);
        Ln_Si_grad = Si_grad / Si;

        return {Ln_Si, Ln_Si_grad};
    }
    

    class ExVariables : public VariableSet {
    public:
        ExVariables() : ExVariables("var_set1") {};
        ExVariables(const std::string& name) : VariableSet(2 * (Graph::N - 1), name)
        {
            Ji_vector = graph.Get_Ji_vector()[global_counter];
            pi_vector = graph.Get_constraint_vector()[global_counter];
        }

        void SetVariables(const VectorXd& x) override
        {
            for(size_t i = 0; i < Graph::N - 1; i++){
                Ji_vector[i] = x(i);
            }
            
            for(size_t i = Graph::N - 1, j = 0; i < 2 * (Graph::N - 1); i++, j++){
                pi_vector[j] = x(i);
            }
        };
        
        VectorXd GetValues() const override
        {
            VectorXd v(2 * (Graph::N - 1));
            for(size_t i = 0; i < Graph::N - 1; i++){
                v[i] = Ji_vector[i];
            }

            for(size_t i = 0, j = Graph::N - 1; i < Graph::N - 1; i++, j++){
                v[j] = pi_vector[i];
            }
            
            return v;
        };

        VecBound GetBounds() const override
        {
            VecBound bounds(GetRows());
            for(size_t i = 0; i < 2 * (Graph::N - 1); i++)
                bounds.at(i) = NoBound;
            
            return bounds;
        }

    private:
        vector<double> Ji_vector;
        vector<double> pi_vector;
    };


    class ExConstraint : public ConstraintSet {
    public:
        ExConstraint() : ExConstraint("constraint1") {}
        ExConstraint(const std::string& name) : ConstraintSet(2 * (Graph::N - 1), name) {}
        
        VectorXd GetValues() const override
        {
            VectorXd g(GetRows());
            VectorXd x = GetVariables()->GetComponent("var_set1")->GetValues();
            S = Si_calc(x);
            
            for(size_t i = 0, j = Graph::N - 1; i < Graph::N - 1; j++, i++){
                g(i) = x(i) - x(j);
            }

            for(size_t i = 0, j = Graph::N - 1; i < Graph::N - 1; j++, i++){
                g(j) = x(i) + x(j);
            }
                        
            return g;
        };
        
        VecBound GetBounds() const override
        {
            VecBound b(GetRows());

            for(size_t i = 0; i < Graph::N - 1; i++){
                b.at(i) = Bounds(-inf, 0);
            }

            for(size_t i = Graph::N - 1; i < 2 * (Graph::N - 1); i++){
                b.at(i) = Bounds(0, inf);
            }
            
            return b;
        }
        
        void FillJacobianBlock (std::string var_set, Jacobian& jac_block) const override
        {
            if (var_set == "var_set1") {
                VectorXd x = GetVariables()->GetComponent("var_set1")->GetValues();
                
                for(int j = 0; j < Graph::N - 1; j++) {
                    for (int i = 0; i < 2 * (Graph::N - 1); i++) {
                        int value;
                        
                        if (i == j)
                            value = 1;
                        else if(i == j + Graph::N - 1)
                            value = -1;
                        else value = 0;
                            jac_block.coeffRef(j, i) = value;
                    }
                }
                
                for(int j = Graph::N - 1; j < 2 * (Graph::N - 1); j++) {
                    for (int i = 0; i < 2 * (Graph::N - 1); i++) {
                        int value;

                        if (j == i || j == i + Graph::N - 1)
                            value = 1;
                        else value = 0;
                        jac_block.coeffRef(j, i) = value;
                    }
                }
            }
        }
    };


    class ExCost: public CostTerm {
    public:
        ExCost() : ExCost("cost_term1") {}
        ExCost(const std::string& name) : CostTerm(name) {}

        double GetCost() const override
        {
            VectorXd x = GetVariables()->GetComponent("var_set1")->GetValues();
            double func_value = 0;

            func_value += S.first;
            for(size_t i = Graph::N - 1; i < 2 * (Graph::N - 1); i++)
                func_value += x(i) * graph.lambda;
            
            return func_value;
        };

        void FillJacobianBlock (std::string var_set, Jacobian& jac) const override
        {
            if (var_set == "var_set1") {
                VectorXd x = GetVariables()->GetComponent("var_set1")->GetValues();

                for(size_t i = 0; i < Graph::N - 1; i++)
                    jac.coeffRef(0, i) = S.second[i];            

                for(size_t i = Graph::N - 1; i < 2 * (Graph::N - 1); i++)
                    jac.coeffRef(0, i) = graph.lambda;
            }
        }
    };

} 