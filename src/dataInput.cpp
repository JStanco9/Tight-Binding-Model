//
//  dataInput.cpp
//  TightBindingModel
//
//  Created by John Stanco on 4/18/17.
//  Copyright © 2017 John Stanco. All rights reserved.
//

#include "../include/dataInput.hpp"

/***************************** Data Input ******************************/
std::vector<mat>
read_hrFile(std::string seedname)
{
  int           currentRow = 0,
                currentCol,
                numCols = 0,
                numBands = 0,
                numSites = 0;
  double        x;
  std::string 	inputLine,
                path = "../data/" + seedname + "_hr.dat";
  mat           data;
  vec           weights;
  std::ifstream inputStream(path);

  std::cout << "Attempting to read data from \'" << path << "\'\n" << "...\n" << std::endl;
    
  if(inputStream.fail())
  {
    throw "The file you are trying to access cannot be found or opened. (\'\')";
  }
  while(inputStream.good())
  {
    getline(inputStream, inputLine);
    inputStream >> numBands;
    inputStream >> numSites;
    weights.resize(numSites);
    data.resize(numSites * numBands * numBands, 0);
    int index = 0;
    
    while(getline(inputStream, inputLine) && index < numSites - 1)
    {
      while(index < numSites && inputStream >> x)
      {
        weights(index++) = x;
      }
    }
          
    while(getline(inputStream, inputLine))
    {
      std::istringstream ss(inputLine);
      currentCol = 0;
      while(ss >> x)
      {
        if (currentCol >= data.n_cols)
        {
          data.resize(numBands * numBands * numSites, numCols + 1);
        }
        data(currentRow, currentCol++) = x;
        if (currentCol >= numCols)
        {
          numCols = currentCol;
        }
      }
      currentRow++;
    }
  }
    
  //Reading _wsvec file
    
  int     length = numBands * numBands * numSites,
          wslength = 0;
  mat     wsweights(length, 1),
          wsvec(length * 10, 3);
  path = "../data/" + seedname + "_wsvec.dat";
    
  std::cout << "Attempting to read data from \'" << path << "\'\n" << "...\n" << std::endl;
    
  inputStream.close();
  inputStream.open(path);
  getline(inputStream, inputLine);

  for(int i = 0; i < length; i++)
  {
    getline(inputStream, inputLine);
    getline(inputStream, inputLine);
    std::istringstream ss(inputLine);
    ss >> x;
    wslength += x;
    if (wsvec.n_rows < wslength)
    {
      wsvec.resize(wslength, 3);
    }
    wsweights(i, 0) = x;
    for(int j = 0; j < x; j++)
    {
      getline(inputStream, inputLine);
      std::istringstream iss(inputLine);
      int tmp, k = 0;
      while(iss >> tmp)
      {
        wsvec(wslength - x + j, k++) = tmp;
      }
    }
    if(!(i % 50000))
    {
      //std::cout << i << " out of " << length << std::endl;
    }
  }
  wsvec.resize(wslength, 3);
  
  std::vector<mat> rdata(4);
  rdata[0] = data;
  rdata[1] = weights;
  rdata[2] = wsvec;
  rdata[3] = wsweights;
  return rdata;
}

/***************************** Reading .win file *****************************/
Lattice
readWinFile(std::string seedname)
{
  mat                      latticeVectors;
  cube                     kPts;
  std::vector<vec>         kPts_from,
                           kPts_to;
  std::vector<std::string> kPt_names_from,
                           kPt_names_to;
  std::vector<double>      consts;
  std::vector<std::string> currentLine;
  std::string              inputLine,
                           x,
                           path = "../data/" + seedname + ".win";
  std::ifstream            inputStream(path);
  unsigned                 currentRow,
                           currentCol,
                           numRows = 0,
                           kDim = 0;
    
  std::cout << "\nAttempting to read data from \'" << path << "\'\n" << "...\n" << std::endl;
    
  if (inputStream.fail())
  {
    throw "The file you are trying to access cannot be found or opened. (\'" "\')";
  }
  
  while (inputStream.good())
  {
    while (getline(inputStream, inputLine))
    {
      if (!inputLine.empty())
      {
        currentLine = split(inputLine, ' ');
        if(inputLine == "begin unit_cell_cart")
        {
          currentCol = 0;
          getline(inputStream, inputLine);
          currentLine = split(inputLine, ' ');
                    
          while (currentLine[0] != "end")
          {
            std::istringstream ss(inputLine);
            currentRow = 0;
            while(ss >> x)
            {
              if(x == "bohr")
              {
                currentCol--;
              }
              else
              {
                if (latticeVectors.n_cols < currentCol + 1)
                {
                  latticeVectors.resize(numRows, currentCol + 1);
                }
                if (latticeVectors.n_rows < currentRow + 1)
                {
                  latticeVectors.resize(currentRow + 1, currentCol + 1);
                }
                latticeVectors(currentRow, currentCol) = atof(x.c_str());
                currentRow++;
              }
              if (numRows < currentRow)
              {
                numRows = currentRow;
              }
            }
              currentCol++;
              getline(inputStream, inputLine);
              currentLine = split(inputLine, ' ');
          }
        }
      }
    }
          
          
    kDim = (unsigned)latticeVectors.n_rows;
    inputStream.clear();
    inputStream.seekg(0, ios::beg);
              
    while (inputStream.good())
    {
      while (getline(inputStream, inputLine))
      {
        if (inputLine == "begin kpoint_path")
        {
          std::string kPoints;
          while(getline(inputStream, inputLine) && inputLine != "end kpoint_path")
          {
            if(inputLine != "end kpoint_path")
            {
              kPoints += inputLine + " ";
            }
          }
          std::istringstream ss(kPoints);
          int index = 0;
          while(ss >> x)
          {
            if(index % 2 == 0)
            {
              kPt_names_from.push_back(x);
            }
            else
            {
              kPt_names_to.push_back(x);
            }
            vec kPt(kDim);
            for(int i = 0; i < kDim; i++)
            {
              if (ss >> x)
              {
                kPt(i) = atof(x.c_str());
              }
            }
            if(index % 2 == 0)
            {
              kPts_from.push_back(kPt);
            }
            else
            {
              kPts_to.push_back(kPt);
            }
            index++;
          }
        }
      }
    }
  }
  if (kPts_from.size() != kPt_names_from.size()
      || kPts_from.size() != kPts_to.size()
      || kPts_to.size() != kPt_names_to.size())
  {
    throw "Error in input file '.win': kPoints not properly specified.";
  }
    
  unsigned int numPts = (unsigned)kPts_from.size() + 1;
  kPts.resize(kDim, numPts, 2);
  kPts.zeros();
    
  for(int i = 0; i < kDim; i++)
  {
    for(int j = 0; j < numPts - 1; j++)
    {
      kPts(i, j, 0) = kPts_from[j](i);
      }
  }
    
  for(int i = 0; i < kDim; i++)
  {
    for(int j = 0; j < numPts - 1; j++)
    {
      kPts(i, j + 1, 1) = kPts_to[j](i);
    }
  }
  Lattice rLattice(latticeVectors, kPts, kPt_names_from, kPt_names_to);
  return rLattice;
}

