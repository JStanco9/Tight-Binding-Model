//
//  BasisChange.cpp
//  
//
//  Created by John Stanco on 5/5/17.
//
//

#include <stdio.h>
#include <sstream>
#include <iomanip>
#include <armadillo>

/* This program is designed to look at lattice information in the Springer Materials website format
 * and print the input required for a PWscf calculation.  Note, the iBrav lattice index can be obtained
 * by cross-referencing the space group with pw.x input description site:
 *   http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PW.html#idm140629872010784
 * This allows one to specify only the space-group and primitive atom coordinates in the PWscf input
 * file (as opposed to every atom in the unit cell)
 */

using namespace arma;

#define PI 3.1415926
#define BOHR 1.88973

double
min(double a, double b)
{
    if(a < b)
    {
        return a;
    }
    return b;
}

void
printMat(mat M)    //For printing with 16-digit precision
{
    for(int i = 0; i < M.n_rows; i++)
    {
        for(int j = 0; j < M.n_cols; j++)
        {
            std::cout << std::setprecision(16) << M(i, j) << "  ";
        }
        std::cout << std::endl;
    }
}

/*
void
swap(mat A, int i1, int i2)
{
    assert(i1 <= A.n_rows & i2 <= A.n_rows & i1 >= 0 & i2 >= 0);
    vec temp = A.row(i1);
    A.row(i1) = A.row(i2);
    A.row(i2) = temp;
    //std::cout << "swap operator" << std::endl;
}


mat
solve(mat A, vec b)
{
    assert((b.size() == A.n_cols) & (A.n_rows == A.n_cols));
    
    vec         pivots(A.n_rows);
    int         maxRow;
    double      currentPivot;
    
    //Perform LU factorization
    for(int index = 0; index < min(A.n_rows, A.n_cols); index++)
    {
        //Pivoting
        maxRow = index;
        for(int currentRow = index; currentRow < A.n_rows; currentRow++)
        {
            if (std::abs(A(currentRow, index)) > std::abs(A(maxRow,index)))
            {
                maxRow = currentRow;
            }
        }
        
        currentPivot = A(maxRow, index);
        swap(A, maxRow, index);
        
        //Save each pivot, allowing reversal of permutation after
        pivots(index) = maxRow;
        
        
        //Multiply first column by a11
        for(int currentRow = index + 1; currentRow < A.n_rows; currentRow++)
        {
            A(currentRow, index) /= currentPivot;
        }
        
        //Act on each interior matrix entry:
        for(int i = index + 1; i < A.n_rows; i++)
        {
            for(int j = index + 1; j < A.n_cols; j++)
            {
                A(i, j) -= A(i, index) * A(index, j);
            }
        }
    }
    
    double temp;
    //Must permute vector b according to how we permuted matrix
    for(int index = 0; index < pivots.size(); index++)
    {
        if(pivots(index) != index)
        {
            temp = b(index);
            b(index) = b(pivots(index));
            b(pivots(index)) = temp;
        }
    }
    
    //Now have LU factorization, must solve Ly = b via forward substitution
    vec y(b.size());
    vec x(y.size());
    double sum;
    
    for(int i = 0; i < y.size(); i++)
    {
        sum = 0;
        for(int j = 0; j < i; j++){
            if(i == j) sum += y(j);
            else sum += A(i, j) * y(j);
        }
        y(i) = (b(i) - sum);  //L[i][i] = 1, so no need to divide by it
        //std::cout << y[i] << std::endl;
    }
    
    //Now have y, must solve Ux = y via backwards substitution
    
    for(int i = x.size() - 1; i >= 0; i--)
    {
        sum = 0;
        for(int j = i + 1; j < x.size(); j++)
        {
            sum += A(i, j) * x(j);
        }
        x(i) = ( y(i) - sum ) / A(i, i);
    }
    
    for(int i = 0; i < x.size(); i++)
    {
        //std::cout << x[i] << std::endl;
    }
    return x;
}*/

void
pwInput(std::string seedname)
{
    std::string     inputLine,
                    path = "../data/" + seedname + "_coords.dat";
    double          a,
                    b,
                    c,
                    alpha,
                    beta,
                    gamma,
                    x;
    int             iBrav = 0,
                    dataFlag = 0;
    mat             vecs(3, 3),
                    atoms_cart,
                    atoms_frac;
    std::ifstream   stream;
    
    alpha *= PI / 180;
    beta  *= PI / 180;
    gamma *= PI / 180;
    
    std::cout << "Enter the lattice dimension \'a\'. (nm)" << std::endl;
    std::cin >> a;
    std::cout << "Enter the lattice dimension \'b\'. (nm)" << std::endl;
    std::cin >> b;
    std::cout << "Enter the lattice dimension \'c\'. (nm)" << std::endl;
    std::cin >> c;
    std::cout << "Enter the lattice dimension \'alpha\'. (degrees)" << std::endl;
    std::cin >> alpha;
    std::cout << "Enter the lattice dimension \'beta\'. (degrees)" << std::endl;
    std::cin >> beta;
    std::cout << "Enter the lattice dimension \'gamma\'. (degrees)" << std::endl;
    std::cin >> gamma;
    
    while(!dataFlag)
    {
        std::cout << "Enter the Bravais Lattice Index of the Lattice. (0-14 & -5, -9, -12)" << std::endl;
        std::cin >> iBrav;
        dataFlag = 1;
        double      c4,
                    tx,
                    ty,
                    tz,
                    u,
                    v, v2y, v2z;
        switch(iBrav)
        {
            case 1 :                            //cubic P (sc)
                vecs << 1 << 0 << 0 << endr
                << 0 << 1 << 0 << endr
                << 0 << 0 << 1 << endr;
                vecs *= a;
                break;
            case 2 :                            //cubic F (fcc)
                vecs << -1 << 0 << 1 << endr
                << 0 << 1 << 1 << endr
                << -1 << 1 << 0 << endr;
                vecs *= (a / 2);
                break;
            case 3 :                            //cubic I (bcc)
                vecs << 1 << 1 << 1 << endr
                << -1 << 1 << 1 << endr
                << -1 << -1 << 1 << endr;
                vecs *= (a / 2);
                break;
            case 4 :                            //Hex & Trigonal P
                vecs << 1 << 0 << 0 << endr
                << -.5 << sqrt(3)/2 << 0 << endr
                << 0 << 0 << c / a << endr;
                vecs *= a;
                break;
            case 5 :                            //Trigonal R (3 - fold axis around <111>)
                c4 = cos(alpha);
                tx = sqrt((1 - c4) / 2);
                ty = sqrt((1 - c4) / 6);
                tz = sqrt((1 + 2 * c4) / 3);
                vecs << tx << -ty << tz << endr
                << 0 << 2 * ty << tz << endr
                << -tx << -ty << tz << endr;
                vecs *= a;
                break;
            case -5 :                           //Trigonal R (3 - fold axis around z)
                c4 = cos(alpha);
                tx = sqrt((1 - c4) / 2);
                ty = sqrt((1 - c4) / 6);
                tz = sqrt((1 + 2 * c4) / 3);
                u = tz - 2 * sqrt(2) * ty;
                v = tz + sqrt(2) * ty;
                vecs << u << v << v << endr
                << v << u << v << endr
                << v << v << u << endr;
                vecs *= (a / sqrt(3));
                break;
            case 6 :                            //Tetragonal P (st)
                vecs << 1 << 0 << 0 << endr
                << 0 << 1 << 0 << endr
                << 0 << 0 << c / a << endr;
                vecs *= a;
                break;
            case 7 :                            //Tetragonal I (bct)
                vecs << 1 << -1 << c / a << endr
                << 1 << 1 << c / a << endr
                << -1 << -1 << c / a << endr;
                vecs *= (a / 2);
                break;
            case 8 :                            //Orthorhombic P
                vecs << a << 0 << 0 << endr
                << 0 << b << 0 << endr
                << 0 << 0 << c << endr;
                break;
            case 9 :                            //Orthorhombic base-centered (bco)
                vecs << a / 2 << b / 2 << 0 << endr
                << -a / 2 << b / 2 << 0 << endr
                << 0 << 0 << c << endr;
                break;
            case -9 :                           //Orthorhombic base-centered (alternate descr.)
                vecs << a / 2 << -b / 2 << 0 << endr
                << a / 2 << b / 2 << 0 << endr
                << 0 << 0 << c << endr;
                break;
            case 10 :                           //Orthorhombic face-centered
                vecs << a / 2 << 0 << c / 2 << endr
                << a / 2 << b / 2 << 0 << endr
                << 0 << b / 2 << c / 2 << endr;
                break;
            case 11 :                           //Orthorhombic body-centered
                vecs << a / 2 << b / 2 << c / 2 << endr
                << -a / 2 << b / 2 << c / 2 << endr
                << -a / 2 << -b / 2 << c / 2 << endr;
                break;
            case 12 :                           //Monoclinic P unique axis c
                vecs << a << 0 << 0 << endr
                << b * sin(gamma) << b * cos(gamma) << 0 << endr
                << 0 << 0 << c << endr;
                break;
            case -12 :                          //Monoclinic P unique axis b
                vecs << a << 0 << 0 << endr
                << 0 << b << 0 << endr
                << c * cos(beta) << 0 << c * sin(beta) << endr;
                break;
            case 13 :                           //Monoclinic base-centered
                vecs << a / 2 << 0 << -c / 2 << endr
                << b * cos(gamma) << b * sin(gamma)<< 0 << endr
                << a / 2 << 0 << c / 2 << endr;
                break;
            case 14 :                           //Triclinic
                v2y = c * (cos(alpha) - cos(beta) * cos(gamma)) / sin(alpha);
                v2z = c * sqrt(1 + 2 * cos(alpha) * cos(beta) * cos(gamma)
                    - pow(cos(alpha), 2) - pow(cos(beta), 2) - pow(cos(gamma), 2))
                    / sin(gamma);
                vecs << a << 0 << 0 << endr
                << b * cos(gamma) << b * sin(gamma)<< 0 << endr
                << c * cos(beta) << v2y << v2z << endr;
                
                break;
            default : std::cout << iBrav << " is not a valid Bravais lattice index" << std::endl;
                dataFlag = 0;
                break;
        }
    }
    
    std::cout << "\nAttempting to open file \'" << path << "\'.\n...\n" << std::endl;
    stream.open(path);
    if(stream.fail())
    {
        throw "The file you are trying to access cannot be found or opened.";
    }

    int     rowIndex = 0,
            numCols = 0;
    while(getline(stream, inputLine))
    {
        std::istringstream ss(inputLine);
        int colIndex = 0;
        if(atoms_cart.n_rows <= rowIndex)
        {
            atoms_cart.resize(rowIndex + 1, numCols);
        }

        while(ss >> x)
        {
            if(atoms_cart.n_cols <= colIndex)
            {
                atoms_cart.resize(rowIndex + 1, colIndex + 1);
            }
            atoms_cart(rowIndex, colIndex++) = x;
            if(colIndex > numCols)
            {
                numCols = colIndex;
            }
            //std::cout << x;
        }
        rowIndex++;
    }
    
    atoms_cart *= 10 * BOHR * a;
    vecs *= 10 * BOHR;
    
    atoms_frac = trans(inv(trans(vecs)) * trans(atoms_cart));
    //atoms_frac = atoms_cart * inv(vecs);

    std::cout << "Atomic Positions in Cartesian Coordinates (bohr): " << std::endl;
    printMat(atoms_cart);
    std::cout << "\nLattice Vectors in Cartesian Coordinates (bohr): " << std::endl;
    printMat(vecs);
    std::cout << "\nAtomic Positions in Lattice Coordinates (crystal_sg): " << std::endl;
    printMat(atoms_frac);
}

/***************************** Main Method **********************************/

int
main(int argc, char* argv[])
{
    if(argc != 2)
    {
        cerr << "Improper number of command line arguments specified" << std::endl;
    }
    else{
        std::string seedname = argv[1];
        try
        {
            pwInput(seedname);
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
