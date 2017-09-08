//
//  plotGap.cpp
//  
//
//  Created by John Stanco on 7/12/17.
//
//

#include <stdio.h>
#include "../include/tightBinding.hpp"

int
main(int argc, char* argv[])
{
  if(argc != 5)
  {
    cerr << "routine plotGap: Imporoper number of command line arguements specified (4)" << std::endl;
  }
  else{
    clock_t t = clock();
    std::string seedname = argv[1];
    int         h, k, l;
    sscanf(argv[2], "%d", &h);
    sscanf(argv[3], "%d", &k);
    sscanf(argv[4], "%d", &l);
        
    TightBinding tb(seedname);
    try
    {
      tb.plotGap(h, k, l);
      std::cout << "Finished" << std::endl;
      t = clock() - t;
      int  time = t / CLOCKS_PER_SEC,
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
