//
//  dataInput.hpp
//  TightBindingModel
//
//  Created by John Stanco on 4/18/17.
//  Copyright © 2017 John Stanco. All rights reserved.
//

#ifndef dataInput_hpp
#define dataInput_hpp

#include "lattice.hpp"

using namespace arma;

std::vector<mat> read_hrFile(std::string seedname);
Lattice readWinFile(std::string seedname);
std::vector<mat> read_wsvecFile(std::string seedname);


#endif /* dataInput_hpp */
