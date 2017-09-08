/*
 * lattice.hpp

 *
 *  Created on: Apr 8, 2017
 *      Author: johnstanco
 */

#include <sstream>
#include "tb_help.hpp"

using namespace arma;

//Currently supports any output of wannier90;
class Lattice{

private:
	mat                      latticeVecs,
                           kVecs;
  cube                     kPts;
  std::vector<std::string> kPtsFrom,
                           kPtsTo;
public:
  Lattice();
  Lattice(mat &lat, cube &ks, std::vector<std::string> &kFrom, std::vector<std::string> &kTo);
  Lattice(const Lattice &other);

	mat latticeVectors();
	mat kVectors();
  cube kPoints();
  std::vector<std::string> kPt_names_from();
  std::vector<std::string> kPt_names_to();
  unsigned dim();
};


