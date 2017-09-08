/*
 * lattice.cpp

 *
 *  Created on: Apr 8, 2017
 *      Author: johnstanco
 */

#include "../include/lattice.hpp"

Lattice::Lattice(){}

Lattice::Lattice(mat &lat, cube &ks, std::vector<std::string> &kFrom, std::vector<std::string> &kTo)
{
  latticeVecs = lat;
  kVecs = recipLat(lat);
  kPts = ks;
  kPtsTo = kTo;
  kPtsFrom = kFrom;
}


Lattice::Lattice(const Lattice &other)
{
  latticeVecs = other.latticeVecs;
  kVecs = other.kVecs;
  kPts = other.kPts;
  kPtsTo = other.kPtsTo;
  kPtsFrom = other.kPtsFrom;
}


mat
Lattice::latticeVectors()
{
  return latticeVecs;
}

mat
Lattice::kVectors()
{
  return kVecs;
}

cube
Lattice::kPoints()
{
  return kPts;
}

std::vector<std::string>
Lattice::kPt_names_from()
{
  return kPtsFrom;
}

std::vector<std::string>
Lattice::kPt_names_to()
{
  return kPtsTo;
}

unsigned
Lattice::dim()
{
  return (unsigned)latticeVectors().n_cols;
}


