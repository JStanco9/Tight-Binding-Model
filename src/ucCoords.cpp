//
//  ucCoords.cpp
//  
//
//  Created by John Stanco on 8/7/17.
//
//

#include <stdio.h>
#include "../include/tightBinding.hpp"

/* Aim is to obtain unit vectors of conventional unit cell,
 * find kVectors, and then convert cartesian coords to lattice basis.
 */

#define BOHR 1.889726125454

int
main()
{
  TightBinding tb("TaAs");
  vec     k(3), G(3), S(3), Z(3);
  mat     B1(3, 3),
          B2(3, 3);
  
  //k << 0.2596415695094647 << -0.2519177560710199 + 1 << 0.2519177560710199;   //conventional uc node W1
  //k << 0.2519177560710199 << 0.2596415695094647 << -0.2596415695094647;       //G-S, G-Z node W1
  //k << 0.4217482327947525 << 0.4416990632397959 << 0.1399152466992413;          //conventional uc node W2
  k << 0.1598660771442794 << 0.4416990632397931 << 0.1399152466992356;        //G-S, G-Z node W2
  
  G << 0 << 0 << 0;
  S << -0.27176 << -0.27176 << 0.27176;
  Z << 0.50000 << 0.50000 << 0.50000;
  
  double  a = 3.437 * BOHR,
          c = 11.656 * BOHR,
          s = norm(tb.kVecs() * (S - G)),
          z = norm(tb.kVecs() * (Z - G));
  
  B1  << a << 0 << 0 << endr
      << 0 << a << 0 << endr
      << 0 << 0 << c << endr;
  
  B2  << s << 0 << 0 << endr
      << 0 << s << 0 << endr
      << 0 << 0 << z << endr;
  
  
  //changeBasis(recipLat(B1), tb.kVecs(), k).print();
  //std::cout << std::endl;
  changeBasis(B2, tb.kVecs(), k).print();
  std::cout << std::endl;
  
  //k << .0185 << .2831 << .6;
  //k << .4827 << .0072 << 1;
  //k << .520  << .037  << .592;
  //changeBasis(tb.kVecs(), recipLat(B1), k).print();
  //std::cout << std::endl;
  //changeBasis(tb.kVecs(), B2, k).print();
   
  std::cout << tb.diff(k) << std::endl;
  
  return 0;
}
