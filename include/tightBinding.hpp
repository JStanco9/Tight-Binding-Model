//
//  TightBinding.hpp
//  
//
//  Created by John Stanco on 6/28/17.
//
//

#ifndef TightBinding_hpp
#define TightBinding_hpp

#include "dataInput.hpp"
#include <time.h>

class TightBinding
{
    
private:
  std::string seedname;
  Lattice     lat;
  mat         hr_dat,
              hr_weights,
              wsvec_dat,
              wsvec_weights;
  int         hamSize; //Dimension of Hamiltonian matrix;
    
public:
  TightBinding();
  TightBinding(std::string seed);
  TightBinding(TightBinding &other);
    
  //Methods that access Lattice
  mat latVecs();
  mat kVecs();
  cube kPoints();
  std::vector<std::string> kPt_names_from();
  std::vector<std::string> kPt_names_to();
    
  //Calculations
  cx_mat Ham(vec k);
  cx_cube expandHam_order1(vec k);
  cx_cube expandHam_order2(vec k);
  //cx_cube expandHam(vec k);
  int computeBands(); 
  double bandGap(vec k);
  vec locateWeylNodes(vec k);
  cx_mat fermiVelocity(vec k);
  vec injectionCurrent(double omega, vec A, double T = 0);
  cx_vec shiftCurrent(double omega, vec A, double T = 0);
  int plotGap(int h, int k, int l);
  };

#endif /* TightBinding_hpp */
