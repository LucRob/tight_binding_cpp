#include <iostream>
#include "params.hpp"
#include "lattice.hpp"

int main() {
//-------------------------------- Define parameters -----------------------------------
    double r = 10; // QD radius in Angs
    double h = 5; // QD height in Angs
    double x = 0.3; // QD As fraction
    std::vector<double> b = {10, 10, 10}; // Barrier thickness around the QD

    // Construct the full set of material parameters
    Params params(r, h, x, b);
    
//-------------------------------- Print parameters -----------------------------------
    // Print lattice constants
    std::cout << "Lattice constant a (InP):   " << params.inp.a << "\n";
    std::cout << "Lattice constant a (InAs):  " << params.inas.a << "\n";
    std::cout << "Lattice constant a (InAsP): " << params.inasp.a << "\n";

    // Print one rNN vector from each material (atom 0, neighbor 0)
    int atom = 0;
    int neighbor = 0;
    
    std::cout << "\nrNN vector (atom 0, neighbor 0):\n";
    std::cout << "InP   : [" << params.inp.rNN[0][neighbor][atom] << ", "
                           << params.inp.rNN[1][neighbor][atom] << ", "
                           << params.inp.rNN[2][neighbor][atom] << "]\n";

    std::cout << "InAs  : [" << params.inas.rNN[0][neighbor][atom] << ", "
                           << params.inas.rNN[1][neighbor][atom] << ", "
                           << params.inas.rNN[2][neighbor][atom] << "]\n";

    std::cout << "InAsP : [" << params.inasp.rNN[0][neighbor][atom] << ", "
                           << params.inasp.rNN[1][neighbor][atom] << ", "
                           << params.inasp.rNN[2][neighbor][atom] << "]\n";
                           
                           
    // Print system parameters
    std::cout << "QD radius: " << params.sys_params.r << " Angs" << "\n";
    std::cout << "QD height: " << params.sys_params.h << " Angs" << "\n";
    std::cout << "QD composition: " << params.sys_params.x*100 << "% arsenic" << "\n";
    std::cout << "QD barrier: " << "[" << params.sys_params.b[0] << ", " << params.sys_params.b[1] << ", " << params.sys_params.b[2] << "]" << " Angs" << "\n";
    
//--------------------------------- Generate lattice ----------------------------------
    Lattice lattice(params);
    lattice.bravais_lattice_gen();
    std::vector<std::array<double, 3>> Ls = lattice.Ls;
    
    std::cout << "Lattice vectors (Ls):\n";
    for (size_t i = 0; i < Ls.size(); ++i) {
        const auto& vec = Ls[i];
        std::cout << "Ls[" << i << "] = (" 
                  << vec[0] << ", " 
                  << vec[1] << ", " 
                  << vec[2] << ")\n";
    }
    
    return 0;
}

