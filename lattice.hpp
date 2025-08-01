#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include "params.hpp"

struct Lattice {
  Params &params;
  int N_atoms;
  std::vector<std::array<double, 3>> Ls;
  std::vector<std::array<double, 3>> r0s; 
  std::vector<int> atom_types;
  std::vector<int> mat_types;
  std::vector<int> ion_types;
  std::vector<int> dot_domaines;
  std::vector<int> unit_cell_number;
  
  // Constructor
  Lattice(Params &params_in): params(params_in) {}
  
  void bravais_lattice_gen()
  {
      // material parameters
      double a = params.inp.a;
      double c = params.inp.c;
      std::vector<double> a1 = params.inp.a1;
      std::vector<double> a2 = params.inp.a2;
      std::vector<double> a3 = params.inp.a3;
            
      // system parameter
      double r = params.sys_params.r;
      double h = params.sys_params.h;
      std::vector<double> b = params.sys_params.b;
      
      // Find span of basis vectors
      double Lx = 2 * (r + b[0]);
      double Ly = 2 * (r + b[1]);
      double Lz = h + 2*b[2];
      
      int m1 = (int)ceil(Lx/(2*a));
      int m2 = (int)ceil(Ly/(2*a));
      int m3 = (int)ceil(Lz/(2*c));
      int N1 = 2*m1+1;
      int N2 = 2*m2+1;
      int N3 = 2*m3+1;
      
      int cnt = 0;
      for (int i1 = -m1; i1 <= m1; i1++) {
          for (int i2 = -m2; i2 <= m2; i2++) {
              for (int i3 = -m3; i3 <= m3; i3++) {
                  std::array<double, 3> pos = {
                      i1 * a1[0] + i2 * a2[0] + i3 * a3[0],
                      i1 * a1[1] + i2 * a2[1] + i3 * a3[1],
                      i1 * a1[2] + i2 * a2[2] + i3 * a3[2]
                  };
                  Ls.push_back(pos);
                  cnt++;
              }
          }
      }
  }
  
  void atomic_lattice_gen() 
  {
      std::vector<double> r1 = params.inp.r1;
      std::vector<double> r2 = params.inp.r2;
      std::vector<double> r3 = params.inp.r3;
      std::vector<double> r4 = params.inp.r4;
      
      for (size_t i = 0; i < Ls.size(); i++) {
          std::array<double, 3> Li = Ls[i];
          Li
          
          r0s.pushback()
      }
  }
};

#endif // LATTICE_HPP
