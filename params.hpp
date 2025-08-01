#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <vector>
#include <cmath>
#include <iostream>

// -------------------------- Utility functions -----------------------
inline double dot_product(const std::vector<double>& A, const std::vector<double>& B) {
    double result = 0.0;
    for (size_t i = 0; i < A.size(); ++i) result += A[i] * B[i];
    return result;
}

inline std::vector<double> cross_product(const std::vector<double>& A, const std::vector<double>& B) {
    return {
        A[1]*B[2] - A[2]*B[1],
        A[2]*B[0] - A[0]*B[2],
        A[0]*B[1] - A[1]*B[0]
    };
}

inline void insert_vector(double rNN[3][4][4], int neighbor, int atom, const std::vector<double>& vec) {
    rNN[0][neighbor][atom] = vec[0];
    rNN[1][neighbor][atom] = vec[1];
    rNN[2][neighbor][atom] = vec[2];
}

void setup_rNN(double rNN[3][4][4], double a, double c) {
    insert_vector(rNN, 0, 0, { a / 2, std::sqrt(3) * a / 6, c / 8 });
    insert_vector(rNN, 1, 0, { 0, -a / std::sqrt(3), c / 8 });
    insert_vector(rNN, 2, 0, { -a / 2, std::sqrt(3) * a / 6, c / 8 });
    insert_vector(rNN, 3, 0, { 0, 0, -3 * c / 8 });

    insert_vector(rNN, 0, 1, { 0, a / std::sqrt(3), -c / 8 });
    insert_vector(rNN, 1, 1, { a / 2, -std::sqrt(3) * a / 6, -c / 8 });
    insert_vector(rNN, 2, 1, { -a / 2, -std::sqrt(3) * a / 6, -c / 8 });
    insert_vector(rNN, 3, 1, { 0, 0, 3 * c / 8 });

    insert_vector(rNN, 0, 2, { 0, a / std::sqrt(3), c / 8 });
    insert_vector(rNN, 1, 2, { a / 2, -std::sqrt(3) * a / 6, c / 8 });
    insert_vector(rNN, 2, 2, { -a / 2, -std::sqrt(3) * a / 6, c / 8 });
    insert_vector(rNN, 3, 2, { 0, 0, -3 * c / 8 });

    insert_vector(rNN, 0, 3, { a / 2, std::sqrt(3) * a / 6, -c / 8 });
    insert_vector(rNN, 1, 3, { 0, -a / std::sqrt(3), -c / 8 });
    insert_vector(rNN, 2, 3, { -a / 2, std::sqrt(3) * a / 6, -c / 8 });
    insert_vector(rNN, 3, 3, { 0, 0, 3 * c / 8 });
}

// -------------------------- InP Material --------------------------
struct InP {
    // For shortness
    double sq3 = std::sqrt(3);

    // Lattice parameters
    double a = 5.856;  // lattice constant (angs)
    double c = sqrt(8.0 / 3.0) * a;  // lattice constant (angs)
    
    std::vector<double> r1 = {0,0,0};  // atom 1 in unit cell (angs)
    std::vector<double> r2 = {a/2, sq3*a/6, c/8};  // atom 2 in unit cell (angs)
    std::vector<double> r3 = {a/2, sq3*a/6, c/2};  // atom 3 in unit cell (angs)
    std::vector<double> r4 = {0, 0, 5*c/8};  // atom 4 in unit cell (angs)
    
    std::vector<double> a1 = {a, 0, 0};  // lattice vector  (angs)
    std::vector<double> a2 = {a/2, a*sq3/2, 0.0};  // lattice vector (angs)
    std::vector<double> a3 = {0, 0, c};  // lattice vector (angs)
    
    double vol = dot_product(cross_product(a1, a2), a3);  // unit cell volume (angs^3)
    
    double rNN[3][4][4]; // Nearest-neighbor relative-positions based on which unit cell atom, indexing is (position coordinates,neighbors,unit cell atom number) (angs)
    
    // strain parameters
    double alpha = 39520;
    double beta  = 6600;
    
    InP() {
        setup_rNN(rNN, a, c);
    }
};

// ---------------------- InAs Material -------------------------
struct InAs {
    // For shortness
    double sq3 = std::sqrt(3);
    
    // Lattice parameters
    double a = 6.044;  // lattice constant (angs)
    double c = sqrt(8.0 / 3.0) * a;  // lattice constant (angs)
    
    std::vector<double> r1 = {0,0,0};  // atom 1 in unit cell (angs)
    std::vector<double> r2 = {a/2, sq3*a/6, c/8};  // atom 2 in unit cell (angs)
    std::vector<double> r3 = {a/2, sq3*a/6, c/2};  // atom 3 in unit cell (angs)
    std::vector<double> r4 = {0, 0, 5*c/8};  // atom 4 in unit cell (angs)
    
    std::vector<double> a1 = {a, 0, 0};  // lattice vector (angs)
    std::vector<double> a2 = {a/2, a*sq3/2, 0.0};  // lattice vector (angs)
    std::vector<double> a3 = {0, 0, c};  // lattice vector (angs)
    
    double vol = dot_product(cross_product(a1, a2), a3);  // unit cell volume  (angs^3)
    
    double rNN[3][4][4]; // Nearest-neighbor relative-positions based on which unit cell atom, indexing is (position coordinates,neighbors,unit cell atom number) (angs)
        
    // Strain parameters
    double alpha = 35180;
    double beta  = 5490;
    
    InAs() {
        setup_rNN(rNN, a, c);
    }
};

// ---------------------------- InAsP Material ----------------------------
struct InAsP {
    double a, c, alpha, beta;
    std::vector<double> r1, r2, r3, r4, a1, a2, a3;
    double vol;
    double rNN[3][4][4];
};

// ---------- Interpolation ----------
inline InAsP interpolate(const InP& inp, const InAs& inas, double x) {
    InAsP inasp;
    
    // For shortness
    double sq3 = std::sqrt(3);
        
    // Lattice constants
    inasp.a = (1 - x) * inp.a + x * inas.a;
    inasp.c = (1 - x) * inp.c + x * inas.c;
    
    // Unit cell atom positions
    inasp.r1 = {0,0,0};  // atom 1 in unit cell (angs)
    inasp.r2 = {inasp.a/2, sq3*inasp.a/6, inasp.c/8};  // atom 2 in unit cell (angs)
    inasp.r3 = {inasp.a/2, sq3*inasp.a/6, inasp.c/2};  // atom 3 in unit cell (angs)
    inasp.r4 = {0, 0, 5*inasp.c/8};
    
    // Lattice vectors
    inasp.a1 = {inasp.a, 0, 0};
    inasp.a2 = {inasp.a/2, inasp.a*sq3/2, 0};
    inasp.a3 = {0, 0, inasp.c};
    
    // Unit cell volume
    inasp.vol = dot_product(cross_product(inasp.a1, inasp.a2), inasp.a3);
    
    // Interpolate rNN
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                inasp.rNN[i][j][k] = (1 - x) * inp.rNN[i][j][k] + x * inas.rNN[i][j][k];
            }
        }
    }
    
    // Strain parameters
    inasp.alpha = (1 - x) * inp.alpha + x * inas.alpha;
    inasp.beta  = (1 - x) * inp.beta  + x * inas.beta;
    
    return inasp;
}


//----------- System parameters ----------
struct Sys_params {
  double r {};
  double h {};
  double x {};
  std::vector<double> b {};
};


// ---------- Params Container ----------
struct Params {
    InP inp;
    InAs inas;
    InAsP inasp;
    Sys_params sys_params;
    
    Params(double r, double h, double x, std::vector<double> b) : inp(), inas(), inasp(interpolate(inp, inas, x)), sys_params{r, h, x, b} {}
};

#endif // PARAMS_HPP

