//
//
//
//  Created by John Stanco on 6/27/17.
//
//

#include "../include/tb_help.hpp"

#define kB 8.61733035e-5  //Boltzmann's constant in units of eV / K

//'dataInput.cpp' help functions

template<typename Out>
void
split(const std::string &s, char delim, Out result)
{
  std::stringstream ss;
  ss.str(s);
  std::string item;
  while (std::getline(ss, item, delim))
  {
    *(result++) = item;
  }
}

std::vector<std::string>
split(const std::string &s, char delim)
{
  std::vector<std::string> words;
  split(s, delim, std::back_inserter(words));
  return words;
}

bool
isDouble(const std::string& s)
{
  std::istringstream iss(s);
  double d;
  char c;
  return iss >> d && !(iss >> c);
}

std::complex<double>
myRound(std::complex<double> x, double n)
{
  double  rtmp = x.real() / n,
          itmp = x.imag() / n;
  return  std::complex<double>(round(rtmp), round(itmp)) * n;
}

//'TightBinding' help functions

std::complex<double>
cexp(double x)
{
  cx_double rVal(cos(x), sin(x));
  return rVal;
}

mat
merge(mat left, mat right)
{
  if(left.n_cols != right.n_cols)
  {
    throw "cannot merge matrices of different column dimensions";
  }
  int     size = left.n_rows + right.n_rows,
          Ind = 0,
          lInd = 0,
          rInd = 0;
  mat     rMat(size, left.n_cols);
    
  for(int i = 0; i < size; i++)
  {
    if(lInd < left.n_rows && rInd < right.n_rows)
    {
      if(left(lInd, 0) < right(rInd, 0))
      {
        rMat.row(i) = left.row(lInd++);
      }
      else
      {
        rMat.row(i) = right.row(rInd++);
      }
    }
    else if(rInd == right.n_rows && lInd < left.n_rows)
    {
      rMat.row(i) = left.row(lInd++);
    }
    else if(lInd == left.n_rows && rInd < right.n_rows)
    {
      rMat.row(i) = right.row(rInd++);
    }
  }
  return rMat;
}

mat
mergeSort(mat m)
{
  int     size = m.n_rows,
          mid = size / 2;
  mat     rMat;
  if(size == 1)
  {
    rMat = m;
  }
  else
  {
    mat     left(mid, m.n_cols),
            right(size - mid, m.n_cols);
    for(int i = 0; i < mid; i++)
    {
      left.row(i) = m.row(i);
    }
    for(int i = mid; i < size; i++)
    {
      right.row(i - mid) = m.row(i);
    }
    rMat = merge(mergeSort(left), mergeSort(right));
  }
  return rMat;
}

/*
vec
merge(vec left, vec right)
{
  if(left.size() != right.size())
  {
    throw "cannot merge matrices of different column dimensions";
  }
    
  int     size = left.size() + right.size(),
          Ind = 0,
          lInd = 0,
          rInd = 0;
  vec     rVec(size);
    
  for(int i = 0; i < size; i++)
  {
    if(lInd < left.size() && rInd < right.size())
    {
      if(left(lInd, 0) < right(rInd, 0))
      {
        rVec(i) = left(lInd++);
      }
      else
      {
        rVec(i) = right(rInd++);
      }
    }
    else if(rInd == right.n_rows && lInd < left.n_rows)
    {
      rVec(i) = left(lInd++);
    }
    else if(lInd == left.n_rows && rInd < right.n_rows)
    {
      rVec(i) = right(rInd++);
    }
  }
  return rVec;
}

vec
mergeSort(vec v)
{
  int     size = v.size(),
          mid = size / 2;
  vec     rVec;
  if(size == 1)
  {
    rVec = v;
  }
  else
  {
    vec      left(mid),
             right(size - mid);
    for(int i = 0; i < mid; i++)
    {
      left(i) = v(i);
    }
    for(int i = mid; i < size; i++)
    {
      right(i - mid) = v(i);
    }
    rVec = merge(mergeSort(left), mergeSort(right));
  }
  return rVec;
}*/


int factorial(int n)
{
  if (n < 0)
  {
    throw "Cannot compute factorial of negative integer";
  }
  else if (n == 0)
  {
    return 1;
  }
  else{
    return n * factorial(n - 1);
  }
}

//constructs simplex with centroid v and distance from centroid to vertex r;
mat buildSimplex(vec v, double r)  //v in any basis
{
  int     n = v.size();
  double  tmp = r;
  mat     x(n, n + 1);
    
  x.zeros();
  for(int i = n; i > 0; i--)
  {
    x(i - 1, i) = tmp;
    for(int j = i; j > 0; j--)
    {
      x(i - 1, j - 1) += -tmp / i;
    }
    tmp *= sqrt(i * i - 1) / i;
  }
  for(int i = 0; i < n; i++)
  {
    x.col(i) += v;
  }
  return x;
}

double
simplexSize(mat x)
{
  int n = x.n_rows;
  if(x.n_cols != (n + 1))
  {
    throw "Simplex must have dimensions n x (n + 1)";
  }
  mat D(n, n);
  for(int i = 0; i < n; i++)
  {
    D.col(i) = x.col(i + 1) - x.col(0);
  }
  return det(D) / (double)factorial(n);
}

field<vec>
grid(mat vertices, int ld)
{
  if(vertices.n_cols != 3)
  {
    throw "function \"grid\": must specify 3 points in order to uniquely specify plane in R^3";
  }
    
  mat kVecs(3, 3);
  kVecs << 0.9674 << 0 << -0.9674 << endr
        << -0.9674 << 0.9674 << 0 << endr
        << 0 << 0.2853 << 0.2853 << endr;
    
  vec     a = vertices.col(1), //defining starting corner
          ab = vertices.col(0) - a,
          ad = vertices.col(2) - a;
  double  l1 = norm(kVecs * ab),
          l2 = norm(kVecs * ad),
          ratio = l1 / l2;
  if(ratio > 1)
  {
    ad *= ratio;
  }
  else if(ratio < 1)
  {
    ab /= ratio;
  }
  
  double  xmult = 1, ymult = 1;
    
  field<vec> pts(ld, ld, 1);
  for(int i = 0; i < ld; i++)
  {
    ymult = (double)i / (ld - 1);
    for(int j = 0; j < ld; j++)
    {
      xmult = (double)j / (ld - 1);
      pts(i, j, 0) = a + xmult * ab + ymult * ad;
    }
  }
  return pts;
}

/*
vec
f_min(double (*f)(vec), mat x, int n)
{
  if((x.n_cols != n + 1) || (x.n_rows != n))
  {
    throw "wrong input guess dimensions";
  }
  mat         fs(n + 1, 2),
              P(n + 1, n + 1);
  vec         x_new(n),
              x_mean(n),
              x_r(n),
              x_e(n), x_c(n);
  double      f_new, f_r, f_e, f_c;

  //evaluate function at initial vertices
  for(int i = 0; i < n + 1; i++)
  {
    fs(i, 0) = f(x.col(i));
    fs(i, 1) = i; //To keep track of initial indices
  }
    
  fs = mergeSort(fs); //Sort by function value
  P.zeros();
    
  for(int i = 0; i < n + 1; i++)
  {
    P(fs(i, 1), i) = 1;
  }
  x *= P;
    
  while(simplexSize(x) > .0001 || (fs(n, 0) - fs(0, 0) > .0001))
  {
    x.print();
    std::cout << std::endl;
    fs.print();
    std::cout << std::endl;
    // 2) Reflection Step
    x_mean.zeros();
    for(int i = 0; i < n; i++)
    {
      x_mean += x.col(i);
    }
        
    x_mean /= n;
    x_r = x_mean + (x_mean - x.col(n));
    f_r = f(x_r);
        
    if(fs(0, 0) <= f_r < fs(n - 1, 0))
    {
      x.col(n) = x_r;
      fs(n, 0) = f_r;
    }
    else if(f_r < fs(0, 0)) // 3) Expansion step
    {
      vec x_e = x_mean + 2 * (x_r - x_mean);
      f_e = f(x_e);
      if(f_e < f_r)
      {
        x.col(n) = x_e;
        fs(n, 0) = f_e;
      }
      else
      {
        x.col(n) = x_r;
        fs(n, 0) = f_r;
      }
    }
    else if(f_r >= fs(n - 1, 0)) // 4) Contraction Step
    {
      if(f_r < fs(n, 0))
      {
        x_c = x_mean + 1 / 2 * (x_r - x_mean);
        f_c = f(x_c);
        if(f_c <= f_r)
        {
          x.col(n) = x_c;
          fs(n, 0) = f_c;
        }
      }
      else if(f_r >= fs(n, 0))
      {
        x_c = x_mean + 1 / 2 * (x.col(n) - x_mean);
        f_c = f(x_c);
        if(f_c < fs(n, 0))
        {
          x.col(n) = x_c;
          fs(n, 0) = f_c;
        }
      }
    }
    else  // 5) Shrink step
    {
      for(int i = 1; i < n + 1; i++)
      {
        x.col(i) = x.col(0) + 1 / 2 * (x.col(i) - x.col(0));
        fs(i, 0) = f(x.col(i));
      }
    }
    for(int i = 0; i < n + 1; i++)
    {
      fs(i, 1) = i;
    }
    fs = mergeSort(fs);
    P.zeros();
    for(int i = 0; i < n + 1; i++)
    {
      P(fs(i, 1), i) = 1;
    }
    x *= P;
  }
  x.print();
  std::cout << std::endl;
  fs.print();
  std::cout << std::endl;
  return x.col(0);
}



double
testFunc(vec x)
{
  return dot(x, x) + 10;
}

int
main()  //Testing purposes
{
  //test simplexSize
 
  //test f_min
  mat simplex(2, 3);
  simplex << .11 << -.11 << 0 << endr
          << .11 << .11 << .22 << endr;
 
  simplex.row(0).print();
  std::cout << std::endl;
  vec x = trans(simplex.row(0));
  x.print();
  std::cout << std::endl;
  trans(x);
  x.print();
   
  vec v = f_min(*testFunc, simplex, 2);
  v.print();
}*/

vec changeBasis(mat A, mat B, mat y)
{
  return inv(A) * B * y;
}

vec changeBasis(mat A, mat y)
{
  return inv(A) * y;
}

mat recipLat(mat lat)
{
  unsigned  latticeDim  = (unsigned)lat.n_cols;
  mat				recip(latticeDim, latticeDim);
  double    pi = std::acos(-1);
  
  //3D case, constructing b1, b2, b3, from a1, a2, a3.
  if (latticeDim == 3)
  {
    double volume = dot(lat.col(0), cross(lat.col(1), lat.col(2)));
    recip.col(0) = 2 * pi * cross(lat.col(1), lat.col(2)) / volume;
    recip.col(1) = 2 * pi * cross(lat.col(2), lat.col(0)) / volume;
    recip.col(2) = 2 * pi * cross(lat.col(0), lat.col(1)) / volume;
  }
  //2D case, constructing b1, b2 from a1, a2.
  else if (latticeDim == 2)
  {
    recip.col(0) = vec(2);
    recip.col(1) = vec(2);
    recip(0, 0) = -lat(1, 1);
    recip(0, 1) = lat(1, 0);
    recip(1, 0) = lat(0, 1);
    recip(1, 1) = -lat(0, 0);
    recip.col(0) *= 2 * pi / dot(lat.col(0), recip.col(0));
    recip.col(1) *= 2 * pi / dot(lat.col(1), recip.col(1));
  }
  else
  {
    throw "function recipLat : Improper Lattice Dimension. Lattice must be of dimension 2 or 3.";
  }
  return recip;
}


cx_mat
downFold(cx_mat H, int val, int cond) //
{
  if(H.n_rows != H.n_cols)
  {
    throw "function downFold : Hamiltonian matrix must be square";
  }
 
  int     dif = cond - val + 1,
          hamSize = H.n_cols;
  cx_mat  H00(dif, dif),
          H11(hamSize - dif, hamSize - dif),
          T01(dif, hamSize - dif),
          T10(hamSize - dif, dif);
  
  //Filling H00
  for(int i = 0; i < dif; i++)
  {
    for(int j = 0; j < dif; j++)
    {
      H00(i, j) = H(i + val, j + val);
    }
  }
  
  //Filling H11
  for(int i = 0; i < val; i++)
  {
    for(int j = 0; j < val; j++)
    {
      H11(i, j) = H(i, j);
    }
  }
  
  for(int i = cond + 1; i < hamSize; i++)
  {
    for(int j = 0; j < val; j++)
    {
      H11(i - dif, j) = H(i, j);
    }
  }
  
  for(int i = 0; i < val; i++)
  {
    for(int j = cond + 1; j < hamSize; j++)
    {
      H11(i, j - dif) = H(i, j);
    }
  }
  
  for(int i = cond + 1; i < hamSize; i++)
  {
    for(int j = cond + 1; j < hamSize; j++)
    {
      H11(i - dif, j - dif) = H(i, j);
    }
  }
  
  
  //Filling T01
  for(int i = 0; i < dif; i++)
  {
    for(int j = 0; j < val; j++)
    {
      T01(i, j) = H(i + val, j);
    }
    for(int j = cond + 1; j < hamSize; j++)
    {
      T01(i, j - dif) = H(i + val, j);
    }
  }
  
  //Filling T10
  for(int i = 0; i < val; i++)
  {
    for(int j = 0; j < dif; j++)
    {
      T10(i, j) = H(i, j + val);
    }
  }
  for(int i = cond + 1; i < hamSize; i++)
  {
    for(int j = 0; j < dif; j++)
    {
      T10(i - dif, j) = H(i, j + val);
    }
  }

  return H00 + T01 * inv(H11) * T10;
}

void
printMat(mat x)
{
  for(uword i = 0; i < x.n_rows; i++)
  {
    for(uword j = 0; j < x.n_cols; j++)
    {
      std::cout << std::setprecision(16) << x(i, j) << "  ";
    }
    std::cout << std::endl;
  }
}

void print_cx_mat(cx_mat x)
{
  for(uword i = 0; i < x.n_rows; i++)
  {
    for(uword j = 0; j < x.n_cols; j++)
    {
      std::cout << std::setprecision(16) << x(i, j) << "  ";
    }
    std::cout << std::endl;
  }
}

vec gradient(double (*f)(vec), vec k, double h)
{
  double  f_l, f_s;
  vec     tmp(k),
          grad(k.size());
  //Second order error
  for(int i = 0; i < k.size(); i++)
  {
    tmp(i) += h;
    f_l = f(tmp);
    tmp(i) -= 2 * h;
    f_s = f(tmp);
    grad(i) = (f_l - f_s);
  }
  
  return grad / h;
}

double fermiDirac(double E, double Ef, double T) //Energies in eV
{
  double ex;
  if(T == 0)
  {
    if(E > Ef)
    {
      return 0;
    }
    else if(E < Ef)
    {
      ex = 0;
    }
    else
    {
      ex = 1;
    }
  }
  else
  {
    ex = exp((E - Ef) / (kB * T));
  }
  return pow(1 + ex , -1);
}

double delta(double a, double alpha)
{
  return exp(-pow(alpha * a, 2));
}
