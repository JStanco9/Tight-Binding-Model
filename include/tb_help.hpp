//
//  round.hpp
//  
//
//  Created by John Stanco on 6/27/17.
//
//

#ifndef round_hpp
#define round_hpp

//#include <cmath>
#include <stdio.h>
#include "./armadillo-7.960.1/include/armadillo"
#include <tgmath.h>
#include <iomanip>

using namespace arma;

template<typename Out>
void split(const std::string &s, char delim, Out result);
std::vector<std::string> split(const std::string &s, char delim);
bool isDouble(const std::string& s);
std::complex<double> myRound(std::complex<double> x, double n = .0005);
cx_double cexp(double x);
mat mergeSort(mat m);
int factorial(int n);
mat buildSimplex(vec v, double r);
double simplexSize(mat m);
vec f_min(double (*f)(vec), mat x, int n);
field<vec> grid(mat vertices, int ld);
mat recipLat(mat lat);
vec changeBasis(mat A, mat y);
vec changeBasis(mat A, mat B, mat y);
cx_mat downFold(cx_mat H, int val, int cond);
void printMat(mat x);
void print_cx_mat(cx_mat x);
vec gradient(double (*f)(vec), vec k, double h = .00001);
double fermiDirac(double E, double Ef, double T);
double delta(double a, double alpha = 25);

#endif /* round_hpp */
