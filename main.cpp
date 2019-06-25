#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>
#include "graph.h"
#include "problem_opt2.h"

using namespace ifopt;

int main() {
    Graph graph = GetClobalObject();
    graph.Print_sigma_all();
    graph.Print_Ji();
    
    for(size_t i = 0; i < 2; i++) {
        for (; global_counter < Graph::N; global_counter++) {
       
            Problem nlp;
            nlp.AddVariableSet(make_shared<ExVariables>());
            nlp.AddConstraintSet(make_shared<ExConstraint>());
            nlp.AddCostSet(make_shared<ExCost>());
            nlp.PrintCurrent();

            IpoptSolver ipopt;
            ipopt.SetOption("linear_solver", "mumps");
            ipopt.SetOption("jacobian_approximation", "exact");
            
            ipopt.Solve(nlp);
            Eigen::VectorXd x = nlp.GetOptVariables()->GetValues();
            cout << x.transpose() << endl;
            
        }
        Jij_to_zero(graph);
        graph.lambda = 0;
        global_counter = 0;
    }
    
    return 0;
}

