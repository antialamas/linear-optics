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

        if( std::pow(2.0,photons) < numbPermutations(i) ) useRysers[i] = true;
        else useRysers[i] = false;

    }

    inBasis.resize(0,0);
    outBasis.resize(0,0);

    graycode.initialize( photons );

    return;

}

void LinearOpticalTransform::initializeCircuit(int& ancillaP,int& ancillaM){

    ancillaPhotons = ancillaP;
    ancillaModes = ancillaM;
    int HSDimension = g( ancillaPhotons + 2,ancillaModes + 4 );

    photons = ancillaPhotons + 2;

    assert( ancillaModes >= ancillaPhotons );

    factorial.resize( ancillaPhotons + 2 + 1 );

    for(int i=0;i<factorial.size();i++) factorial[i] = doublefactorial(i);

    useRysers.resize( HSDimension );

    nPrime.resize( HSDimension );
    mPrime.resize( HSDimension );

    for(int i=0;i<HSDimension;i++){ nPrime.at(i).resize(ancillaModes + 4); mPrime.at(i).resize(ancillaPhotons + 2); }

    setNPrimeAndMPrime(nPrime,mPrime);

    for(int i=0;i<HSDimension;i++){

        if( std::pow(2.0,ancillaPhotons + 2) < numbPermutations(i) ) useRysers[i] = true;
        else useRysers[i] = false;

    }

    graycode.initialize( ancillaPhotons + 2 );

    return;

}

void LinearOpticalTransform::setMutualInformation(Eigen::MatrixXcd& U){

    mutualInformation = 0;

    for(int y=0;y<useRysers.size();y++){

        if( useRysers[y] == true ) rysersHXY(U,y);

        else permutationHXY(U,y);

    }

    return;

}

void LinearOpticalTransform::rysersHXY(Eigen::MatrixXcd& U,int& y){

    assert( false && "UP TO HERE - WRITE CONTRIBUTION TO MUTUAL INFORMATION WITH RYSERS");

    return;

}

void LinearOpticalTransform::permutationHXY(Eigen::MatrixXcd& U,int& y){

    std::complex<double> stateAmplitude[4];

    stateAmplitude[0] = 0.0;
    stateAmplitude[1] = 0.0;
    stateAmplitude[2] = 0.0;
    stateAmplitude[3] = 0.0;

    double pyx[4];

    do{

        setStateAmplitude(stateAmplitude,U,y);

    } while( std::next_permutation( mPrime[y].begin(), mPrime[y].end() ) );

    normalizeStateAmplitude(stateAmplitude,y);

    pyx[0] = std::norm( stateAmplitude[0] );
    pyx[1] = std::norm( stateAmplitude[1] );
    pyx[2] = std::norm( stateAmplitude[2] );
    pyx[3] = std::norm( stateAmplitude[3] );

    if(pyx[0] != 0.0) mutualInformation += pyx[0] * log2( ( pyx[0] + pyx[1] + pyx[2] + pyx[3] ) / pyx[0] );
    if(pyx[1] != 0.0) mutualInformation += pyx[1] * log2( ( pyx[0] + pyx[1] + pyx[2] + pyx[3] ) / pyx[1] );
    if(pyx[2] != 0.0) mutualInformation += pyx[2] * log2( ( pyx[0] + pyx[1] + pyx[2] + pyx[3] ) / pyx[2] );
    if(pyx[3] != 0.0) mutualInformation += pyx[3] * log2( ( pyx[0] + pyx[1] + pyx[2] + pyx[3] ) / pyx[3] );

    return;

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

        bool even = ( photons + 1 ) % 2; // true;

        A(i,j) = 0;

        Eigen::ArrayXcd weights = Eigen::ArrayXcd::Zero( photons );

        while( graycode.iterate() ){

            /** ========= UNROLL THIS LOOP MANUALLY FOR SPEED ===================== */

            for(int l=0;l<photons;l++){

                weights(l) += boolPow( graycode.sign ) * U( m[j][l],mPrime[i][graycode.j] );

            }

            A(i,j) -= boolPow( even ) * weights.prod();

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

double LinearOpticalTransform::numbPermutations(int& i){

    double output = factorial[photons];

    for(int j=0;j<nPrime[i].size();j++) output /= factorial[ nPrime[i][j] ];

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

void LinearOpticalTransform::setNPrimeAndMPrime(std::vector< std::vector<int> >& nPrime,std::vector< std::vector<int> >& mPrime){

    Eigen::MatrixXi outBasis;

    setToFullHilbertSpace( ancillaPhotons + 2,ancillaModes + 4,outBasis );

    for(int i=0;i<outBasis.rows();i++){

        for(int j=0;j<nPrime.at(i).size();j++) nPrime.at(i).at(j) = outBasis(i,j);

        setmVec(mPrime.at(i),nPrime.at(i));

    }

    return;

}


void LinearOpticalTransform::setToFullHilbertSpace(int subPhotons, int subModes,Eigen::MatrixXi& nv){

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

template <typename T>
void LinearOpticalTransform::printVec(std::vector<T>& a){

    for(int i=0;i<a.size();i++) std::cout << a[i] << " ";

    std::cout << std::endl;

    return;

}

int LinearOpticalTransform::g(const int& n,const int& m){
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
