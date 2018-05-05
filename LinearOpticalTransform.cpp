/*! \file LinearOpticalTransform.cpp
 *  \brief This the general algorithm for simulating a linear optical quantum circuit.
 *
 *     For more information, refer to:
 *      https://arxiv.org/abs/1711.01319
 *
 *
 *
 * \author Jake Smith <jsmith74@tulane.edu>
 */


#include "LinearOpticalTransform.h"

LinearOpticalTransform::LinearOpticalTransform(){


}

void LinearOpticalTransform::setA(Eigen::MatrixXcd& U){

    for(int i=0;i<A.rows();i++){

        if( useRysers[i] == true ) rysersAlgorithm(U,i);

        else permutationAlgorithm(U,i);

    }

    return;

}

void LinearOpticalTransform::rysersAlgorithm(Eigen::MatrixXcd& U,int& i){

    double bosonOutput = 1.0;

    for(int p=0;p<nPrime[i].size();p++) bosonOutput *= factorial[ nPrime[i][p] ];

    for(int j=0;j<A.cols();j++){

        bool even = true;

        A(i,j) = 0;

        Eigen::ArrayXcd weights = Eigen::ArrayXcd::Zero( photons );

        while( graycode.iterate() ){

            for(int l=0;l<photons;l++){

                weights(l) +=
                -std::pow( -1,graycode.sign ) * U( m[j][l],mPrime[i][graycode.j] );

            }

            A(i,j) += std::pow( -1,even ) * weights.prod();

            even = !even;

        }

        double bosonInput = 1.0;

        for(int p=0;p<n[j].size();p++) bosonInput *= factorial[ n[j][p] ];

        A(i,j) /= sqrt( bosonInput * bosonOutput );

    }

    return;

}

void LinearOpticalTransform::permutationAlgorithm(Eigen::MatrixXcd& U,int& i){

    A.row(i) = Eigen::VectorXcd::Zero(A.cols());

    do{

        for(int j=0;j<A.cols();j++){

            std::complex<double> Uprod(1.0,0.0);

            for(int k=0;k<m[j].size();k++){

                Uprod *= U( m[j][k],mPrime[i][k] );

            }

            A(i,j) += Uprod;

        }

    } while( std::next_permutation( mPrime[i].begin(), mPrime[i].end() ) );

    double bosonNum = 1.0;

    for(int p=0;p<U.rows();p++) bosonNum *= factorial[ nPrime[i][p] ];

    for(int j=0;j<A.cols();j++){

        double bosonDen = 1.0;

        for(int p=0;p<U.rows();p++) bosonDen *= factorial[ n[j][p] ];

        A(i,j) *= sqrt( bosonNum/bosonDen );

    }

    return;

}

void LinearOpticalTransform::initializeCircuit(Eigen::MatrixXi& inBasis, Eigen::MatrixXi& outBasis){

    A.resize( outBasis.rows(), inBasis.rows() );

    n.resize( inBasis.rows() );
    m.resize( inBasis.rows() );

    nPrime.resize( outBasis.rows() );
    mPrime.resize( outBasis.rows() );

    photons = inBasis.row(0).sum();

    for(int i=0;i<inBasis.rows();i++){

        assert( inBasis.row(i).sum() == photons && "Error: Photon number must be preserved you have included some input basis states that do not have the correct number of photons." );

        n.at(i).resize( inBasis.cols() );

        for(int j=0;j<inBasis.cols();j++) n.at(i).at(j) = inBasis(i,j);

        m.at(i).resize( photons );

        setmVec( m.at(i), n.at(i) );

    }

    for(int i=0;i<outBasis.rows();i++){

        assert( outBasis.row(i).sum() == photons && "Error: Photon number must be preserved you have included some output basis states that do not have the correct number of photons." );

        nPrime.at(i).resize( outBasis.cols() );

        for(int j=0;j<outBasis.cols();j++) nPrime.at(i).at(j) = outBasis(i,j);

        mPrime.at(i).resize( photons );

        setmVec( mPrime.at(i), nPrime.at(i) );

    }

    factorial.resize( photons + 1 );

    for(int i=0;i<factorial.size();i++) factorial[i] = doublefactorial(i);

    useRysers.resize( outBasis.rows() );

    for(int i=0;i<outBasis.rows();i++){

        if( std::pow(2,photons) < numbPermutations(i) ) useRysers[i] = true;
        else useRysers[i] = false;

    }

    inBasis.resize(0,0);
    outBasis.resize(0,0);

    graycode.initialize( photons );

    return;

}

int LinearOpticalTransform::numbPermutations(int& i){

    int output = factorial[photons] + 0.5;

    int intFactorial[factorial.size()];

    for(int j=0;j<factorial.size();j++) intFactorial[j] = factorial[j] + 0.5;

    for(int j=0;j<nPrime[i].size();j++) output /= intFactorial[ nPrime[i][j] ];

    int counter = 0;

    do{

        counter++;

    } while( std::next_permutation( mPrime[i].begin(),mPrime[i].end() ) );

    assert( counter == output );

    return output;

}

void LinearOpticalTransform::setmVec(std::vector<int>& m, std::vector<int>& n){

    int k=0;

    for(int i=0;i<n.size();i++){

        for(int j=0;j<n.at(i);j++){

            m.at(k) = i;

            k++;

        }

    }

    return;
}

template <typename T>
void LinearOpticalTransform::printVec(std::vector<T>& a){

    for(int i=0;i<a.size();i++) std::cout << a[i] << " ";

    std::cout << std::endl;

    return;

}

double LinearOpticalTransform::doublefactorial(int x){

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
