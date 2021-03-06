//
//  findWeyl.cpp
//  
//
//  Created by John Stanco on 6/27/17.
//
//

#include "../include/tightBinding.hpp"

//Object oriented style that creates a w90 object, consisting of a lattice and the data from the input files
//contains method to find Weyl nodes

#define BOHR 1.889726125454

int
main(int argc, char* argv[])
{
  if(argc != 2)
  {
    cerr << "Improper number of command line arguments specified" << std::endl;
  }
  else
  {
    clock_t t = clock();
    std::string seedname = argv[1];
    TightBinding tb(seedname);
    try
    {
      vec     w(3), k(3), G(3), S(3), Z(3);
      mat     B1(3, 3),
              B2(3, 3);
      k << 0.2596415695094647 << -0.2519177560710199 + 1 << 0.2519177560710199;   //conventional uc node W1
      //k << 0.2519177560710199 << 0.2596415695094647 << -0.2596415695094647;       //G-S, G-Z node W1
      //k << 0.4217482327947525 << 0.4416990632397959 << 0.1399152466992413;          //conventional uc node W2
      //k << 0.1598660771442794 << 0.4416990632397931 << 0.1399152466992356;        //G-S, G-Z node W2

      //w << 0.4272 << 0.4473 << 0.1447;
      //w << 0.1677 << 0.4508 << 0.1492;

      w = tb.locateWeylNodes(k);
      
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
      
      std::cout << "\nLattice Coords: ";
      printMat(trans(changeBasis(tb.kVecs(), w)));
      std::cout << "\nCartesian Coords: ";
      printMat(trans(w));
      std::cout << "\nConventional UC units: ";
      printMat(trans(changeBasis(recipLat(B1), w)));
      std::cout << "\nS - G, Z - G units: ";
      printMat(trans(changeBasis(B2, w)));
      std::cout << std::endl;
      //k << .0185 << .2831 << .6;
      //k << .4827 << .0072 << 1;
      //k << .520  << .037  << .592;
      //changeBasis(tb.kVecs(), recipLat(B1), k).print();
      //std::cout << std::endl;
      //changeBasis(tb.kVecs(), B2, k).print();
      
      
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
  }
  return 0;
}
