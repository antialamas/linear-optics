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
#include "omp.h"
#include <unsupported/Eigen/MatrixFunctions>

void printState(Eigen::VectorXcd& vec,Eigen::MatrixXi& basis);
void setToFullHilbertSpace(const int& subPhotons, const int& subModes,Eigen::MatrixXi& nv);
void setToRandomBasisStates(Eigen::MatrixXi& basis,int photons,int modes,int basisDim);


int main(){

    /** Establish number of photons and modes and input and output Fock basis  */

        int photons = 7;
        int modes = 13;

        Eigen::MatrixXi inBasis,outBasis;

        setToRandomBasisStates(inBasis,photons,modes,4);
        setToRandomBasisStates(outBasis,photons,modes,5);

        //for(int i=0;i<photons;i++) inBasis(0,i) = 1;
        //for(int i=0;i<photons;i++) outBasis(0,i) = 1;

        std::cout << inBasis << std::endl << std::endl;
        std::cout << outBasis << std::endl << std::endl;

    /** Initialize the LinearOpticalTransform object */

        LinearOpticalTransform LOCircuit;

        LOCircuit.initializeCircuit(inBasis,outBasis);

    /** Set the Mach-Zehnder interferometer */

        Eigen::MatrixXcd H = Eigen::MatrixXcd::Random(modes,modes);

        H += H.conjugate().transpose().eval();

        H *= std::complex<double>(0.0,1.0);

        Eigen::MatrixXcd U = H.exp();

    /** Construct A(U) */

        double startTime = omp_get_wtime();

        LOCircuit.setA(U);

        double endTime = omp_get_wtime();

        std::cout << LOCircuit.A << std::endl << std::endl;

        std::cout << "Running time: " << endTime - startTime << std::endl << std::endl;

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
