//
//  test_tb.cpp
//  
//
//  Created by John Stanco on 8/7/17.
//
//

#include <stdio.h>
#include "../include/tightBinding.hpp"

#define PI 3.141592653589793

int main()
{
  clock_t t = clock();
  std::string seedname = "TaAs";
  TightBinding tb(seedname);
  try
  {
    
    vec w(3);
    w << 0.0072 << 0.4827 << 1 << endr;
    
    vec a(3);
    a << -0.0995 << 0.5069 << 0.9005;
    
    vec d(3);
    d << .5 << 0 << 0;
    
    vec d_(3);
    d_ << .51 << 0 << 0;
    
    vec e(3);
    e << 0.2596415695094647 << -0.2519177560710199 + 1 << 0.2519177560710199;
    
    vec G(3);
    G << 0 << 0 << 0;
    
    //(tb.fermiVelocity(d) - tb.Ham(d * (1 + .0001))).print();
    //tb.fermiVelocity(G).print();
    //std::cout << "\n\n\n" << std::endl;
    //tb.Ham(G).print();
    //std::cout << tb.kVecs() << std::endl;
    
    //std::cout << tb.fermiVel(e) - tb.buildHam(e);
    tb.fermiVelocity(e);
    std::cout << '\n' << tb.bandGap(e) << std::endl;
    
    /*
     vec v0(3),
     v1(3),
     v2(3);
     v0 << .5 << 0  << 0;
     v1 << 0  << .5 << 0;
     v2 << 0  << 0  << .5;
     
     
     std::cout << tb.bandGap(v0) << std::endl;
     std::cout << tb.bandGap(v1) << std::endl;
     std::cout << tb.bandGap(v2) << std::endl;
     
     vec G(3);
     G  << 0  << 0  << 0;
     
     std::cout << tb.diff(G) << std::endl;;
     */
    //tb.locateWeylNodesß(w);
    std::cout << "Finished" << std::endl;
    t = clock() - t;
    int  time = t / CLOCKS_PER_SEC,
    hrs = time / 3600,
    min = time / 60 - hrs * 60,
    sec = time - (hrs * 3600 + min * 60);
    std::cout << "Time Elapsed:   " <<  hrs << " hours, " << min
    << " minutes, " << sec << " seconds" << std::endl;
  }
  catch(std::string err)
  {
    cerr << err << std::endl;
  }
  catch(...)
  {
    return EXIT_FAILURE;
  }
  
  return 0;
}
