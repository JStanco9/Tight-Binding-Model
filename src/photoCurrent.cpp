//
//  photoCurrent.cpp
//  
//
//  Created by John Stanco on 8/19/17.
//
//

#include "../include/tightBinding.hpp"

int
main(int argc, char* argv[])
{
  if(argc != 2)
  {
    cerr << "routine tbBands: Improper number of command line arguments specified (1)" << std::endl;
  }
  else
  {
    clock_t t = clock();
    std::string seedname = argv[1];
    TightBinding tb(seedname);
    try
    {
      double  omega = 750 * pow(10, 12);
      vec     A(3);
      A << 3 << 5 << 0;
      //tb.injectionCurrent(omega, A).print();
      tb.shiftCurrent(omega, A).print();
      
      std::cout << "\n\nFinished" << std::endl;
      t = clock() - t;
      int   time = t / CLOCKS_PER_SEC,
            hrs = time / 3600,
            min = time / 60 - hrs * 60,
            sec = time - (hrs * 3600 + min * 60);
      std::cout << "Time Elapsed:   " <<  hrs << " hours, " << min << " minutes, " << sec << " seconds" << std::endl;
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

