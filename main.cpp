/*! \file main.cpp
 *  \brief This is a demonstration on how to use the LinearOpticalTransform object.
 *
 *     The LinearOpticalTransform object constructs a quantum evolution operator from the unitary matrix
 *     describing a linear optical quantum circuit.
 *
 *
 *     For more information, refer to:
 *      https://arxiv.org/abs/1711.01319
 *
 *
 *     Here, we construct the evolution operator A(U) for a simple Mach-Zehnder interferometer; we apply a phase shift
 *     of pi/3 to the first optical mode, and then apply a 50-50 beam splitter between modes 1 and 2.
 *
 *
 * \author Jake Smith <jsmith74@tulane.edu>
 */

#include "LinearOpticalTransform.h"
#define PI 3.141592653589793
#include "omp.h"
#include <unsupported/Eigen/MatrixFunctions>


void printState(Eigen::VectorXcd& vec,Eigen::MatrixXi& basis);
void setToFullHilbertSpace(const int& subPhotons, const int& subModes,Eigen::MatrixXi& nv);
void setToRandomBasisStates(Eigen::MatrixXi& basis,int photons,int modes,int basisDim);


int main(){

    /** Establish number of photons and modes and input and output Fock basis  */

        int photons = 2;
        int modes = 4;

        Eigen::MatrixXi inBasis(4,modes);
        Eigen::MatrixXi outBasis(10,modes);

        inBasis << 2,0,0,0,
                   0,2,0,0,
                   0,0,2,0,
                   0,0,0,2;

        outBasis << 2,0,0,0,
                    1,1,0,0,
                    1,0,1,0,
                    1,0,0,1,
                    0,2,0,0,
                    0,1,1,0,
                    0,1,0,1,
                    0,0,2,0,
                    0,0,1,1,
                    0,0,0,2;

    /** Initialize the LinearOpticalTransform object */

        LinearOpticalTransform LOCircuit;

        LOCircuit.initializeCircuit(inBasis,outBasis);

    /** Set the Mach-Zehnder interferometer */

        Eigen::MatrixXcd U1(modes,modes);
        Eigen::MatrixXcd U2(modes,modes);

        std::complex<double> I(0.0,1.0);

        U1 << exp( I * PI /3.0),          0,      0,      0,
                            0,            1,      0,      0,
                            0,            0,      1,      0,
                            0,            0,      0,      1;

        U2 << cos( PI/4.0 ),  sin( PI/4.0 ),      0,      0,
             -sin( PI/4.0 ),  cos( PI/4.0 ),      0,      0,
                        0,            0,          1,      0,
                        0,            0,          0,      1;


        Eigen::MatrixXcd U;

        U = U1 * U2;

    /** Construct A(U) */

        LOCircuit.setA(U);

    /** Simulate the the evolution of the quantum state 0.5 * |2,0,0,0> + 0.5 * |0,2,0,0> + 0.5 * |0,0,2,0> + 0.5 * |0,0,0,2>
        through the Mach Zehnder interferometer */

        Eigen::VectorXcd psi(4);

        psi << 0.5,
               0.5,
               0.5,
               0.5;

        Eigen::VectorXcd psiPrime(10);

        psiPrime = LOCircuit.A * psi;

    /** Print Results */

        inBasis.resize(4,modes);
        outBasis.resize(10,modes);

        inBasis << 2,0,0,0,
                   0,2,0,0,
                   0,0,2,0,
                   0,0,0,2;

        outBasis << 2,0,0,0,
                    1,1,0,0,
                    1,0,1,0,
                    1,0,0,1,
                    0,2,0,0,
                    0,1,1,0,
                    0,1,0,1,
                    0,0,2,0,
                    0,0,1,1,
                    0,0,0,2;

        std::cout << "Input State:\n\n";

        printState(psi,inBasis);

        std::cout << "Output State:\n\n";

        printState(psiPrime,outBasis);

        std::cout << "A(U):\n\n";

        std::cout << LOCircuit.A << std::endl << std::endl;


    return 0;

}

void setToRandomBasisStates(Eigen::MatrixXi& basis,int photons,int modes,int basisDim){

    basis = Eigen::MatrixXi::Zero(basisDim,modes);

    for(int i=0;i<basisDim;i++) for(int j=0;j<modes;j++){

        basis(i,j) = rand() % ( 1 + photons - basis.row(i).sum() );

	if(j==modes-1) basis(i,j) += photons - basis.row(i).sum();

    }

    return;

}

inline double doublefactorial(int x){

    assert(x < 171);

    double total=1.0;
    if (x>=0){
        for(int i=x;i>0;i--){
            total=i*total;
        }
    }
    else{
        std::cout << "invalid factorial" << std::endl;
        total=-1;
    }
    return total;
}

inline int g(const int& n,const int& m){
    if(n==0 && m==0){
        return 0;
    }
    else if(n==0 && m>0){
        return 1;
    }

    else{
        return (int)(doublefactorial(n+m-1)/(doublefactorial(n)*doublefactorial(m-1))+0.5);
    }
}

void setToFullHilbertSpace(const int& subPhotons, const int& subModes,Eigen::MatrixXi& nv){

    if(subPhotons==0 && subModes == 0){

        nv.resize(0,0);

        return;

    }

    int markers = subPhotons + subModes - 1;
    int myints[markers];
    int i = 0;
    while(i<subPhotons){
        myints[i]=1;
        i++;
    }
    while(i<markers){
        myints[i]=0;
        i++;
    }
    nv = Eigen::MatrixXi::Zero(g(subPhotons,subModes),subModes);
    i = 0;
    int j,k = 0;
    do {
        j = 0;
        k = 0;
        while(k<markers){
        if(myints[k]==1){
            nv(i,j)=nv(i,j)+1;
        }
        else if(myints[k]==0){
            j++;
        }

        k++;
        }
        i++;
    } while ( std::prev_permutation(myints,myints+markers) );
    return;;
}


void printState(Eigen::VectorXcd& vec,Eigen::MatrixXi& basis){

    for(int i=0;i<vec.size();i++){

        std::cout << vec(i) << " * |" << basis.row(i) << ">\n";

    }

    std::cout << std::endl;

    return;

}
