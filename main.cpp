#include "graph.h"
#include "gradient_method.h"
int main() {
    /*
    ifstream input_bonds(R"(C:\Users\ILYA\CLionProjects\Diplom\alanin.psf)"),
             input_q(R"(C:\Users\ILYA\CLionProjects\Diplom\alanin.psf)"),
             input_xyz(R"(C:\Users\ILYA\CLionProjects\Diplom\alanin.txt)");

    Atom atom;

    atom.Create_bonds(input_bonds);
    //PrintBonds(atom);
    atom.Create_q(input_q);
    //PrintQ(atom);
    atom.Create_xyz(input_xyz);
    PrintAll(atom);
    */
    random_device rd;
    mt19937 gen;
    gen.seed(rd());

    Graph graph(gen, rd);
    Graph_reconst_calc(graph, 0.00005);
    graph.Print_Ji();
    //cout << S_calculate(graph);
    //cout << Gradient();

    return 0;
}

/*
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <tuple>
#include <vector>

//#include <random>
//#include <ctime>
//#include <ip>

using namespace std;

using Type_atom = map<int, tuple<set<int>, double, vector<double>>>;

ostream& operator<<(ostream& os, const set<int>& s) {
    for(const auto& i : s)
        cout << i << " ";

    return os;
}

ostream& operator<<(ostream& os, const vector<double>& v) {
    for(const auto& i : v)
        cout << i << " ";

    return os;
}

ostream& operator<<(ostream& os, const tuple<set<int>, double, vector<double>>& t) {
    cout << get<0>(t) << " ";
    cout << get<1>(t) << " ";
    cout << get<2>(t);

    return os;
}

ostream& operator<<(ostream& os, const Type_atom& t) {
    for(const auto& item : t) {
        cout << item.first << ": " << item.second << endl;
    }
    return os;
}

class Atom {
private:
    Type_atom atoms;
public:
    void Create_bonds(ifstream& input) {
        string line;
        int atom1, atom2;

        while (getline(input, line)) {
            if(line.find("bonds") != string::npos) {
                break;
            }
        }

        getline(input, line);
        while(line.find("angles") == string::npos) {
            stringstream ss(line);

            while (ss) {
                ss >> atom1;

                if(ss.eof())
                    break;

                ss >> atom2;

                get<0>(atoms[atom1]).insert(atom2);
                get<0>(atoms[atom2]).insert(atom1);
            }

            getline(input, line);
        }
    }
    void Create_q(ifstream& input) {
        string line;
        string q;

        while (getline(input, line)) {
            if(line.find("NATOM") != string::npos) {
                break;
            }
        }

        getline(input, line);
        for (auto &atom : atoms) {
            stringstream ss(line);

            for(int i = 0; i < 8; i++) {
                ss >> q;
            }

            get<1>(atom.second) = atof(q.c_str());
            getline(input, line);
        }
    }
    void Create_xyz(ifstream& input) {
        string line, x, y, z;

        getline(input, line);
        while (getline(input, line)) {
            if(line.find("created by user:") != string::npos) {
                break;
            }
        }

        getline(input, line);
        for (auto &atom : atoms) {
            stringstream ss(line);

            for(int i = 0; i < 8; i++) {
                if(i < 6)
                    ss >> x;
                else if(i == 6)
                    ss >> y;
                else if(i == 7 )
                    ss >> z;
            }

            get<2>(atom.second).push_back(atof(x.c_str()));
            get<2>(atom.second).push_back(atof(y.c_str()));
            get<2>(atom.second).push_back(atof(z.c_str()));
            getline(input, line);
        }
    }
    Type_atom GetAtoms() const {
        return atoms;
    }
};

void PrintBonds(const Atom& atom) {
    for(const auto& b : atom.GetAtoms())
        cout << b.first << ": " << get<0>(b.second) << endl;
}

void PrintQ(const Atom& atom) {
    for(const auto& b : atom.GetAtoms())
        cout << b.first << ": " << get<1>(b.second) << endl;
}

void PrintAll(const Atom& atom) {
    cout << atom.GetAtoms();
}
*/
