#ifndef LINEAROPTICALTRANSFORM_H_INCLUDED
#define LINEAROPTICALTRANSFORM_H_INCLUDED

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "GrayCode.h"

class LinearOpticalTransform{

    public:

        LinearOpticalTransform();
        void initializeCircuit(Eigen::MatrixXi& inBasis, Eigen::MatrixXi& outBasis);
        void initializeCircuit(int& ancillaP,int& ancillaM);
        void setA(Eigen::MatrixXcd& U);
        void setMutualInformation(Eigen::MatrixXcd& U);

        double mutualInformation;

        Eigen::MatrixXcd A;

    private:

        std::vector< std::vector<int> > n,m,nPrime,mPrime;
        std::vector<double> factorial;
        std::vector<bool> useRysers;
        int photons, ancillaPhotons, ancillaModes;
        GrayCode graycode;

        template <typename T>
        void printVec(std::vector<T>& a);

        double numbPermutations(int& i);
        double doublefactorial(int x);
        void permutationAlgorithm(Eigen::MatrixXcd& U,int& i);
        void setmVec(std::vector<int>& m, std::vector<int>& n);
        void rysersAlgorithm(Eigen::MatrixXcd& U,int& i);

        void rysersHXY(Eigen::MatrixXcd& U,int& y);
        void permutationHXY(Eigen::MatrixXcd& U,int& y);

        inline void setStateAmplitude(std::complex<double> stateAmplitude[],Eigen::MatrixXcd& U,int y);
        inline void normalizeStateAmplitude(std::complex<double> stateAmplitude[],int y);

        inline double boolPow(bool& x);
        int g(const int& n,const int& m);
        void setNPrimeAndMPrime(std::vector< std::vector<int> >& nPrime,std::vector< std::vector<int> >& mPrime);
        void setToFullHilbertSpace(int subPhotons, int subModes,Eigen::MatrixXi& nv);

};

inline double LinearOpticalTransform::boolPow(bool& x){

    return -1 + 2 * x;

}

inline void LinearOpticalTransform::normalizeStateAmplitude(std::complex<double> stateAmplitude[],int y){

    stateAmplitude[0] *= 0.7071067811865475;
    stateAmplitude[1] *= 0.7071067811865475;
    stateAmplitude[2] *= 0.7071067811865475;
    stateAmplitude[3] *= 0.7071067811865475;

    for(int p=0;p<ancillaModes+4;p++){

        stateAmplitude[0] *= sqrt( factorial[ nPrime[y][p] ] );
        stateAmplitude[1] *= sqrt( factorial[ nPrime[y][p] ] );
        stateAmplitude[2] *= sqrt( factorial[ nPrime[y][p] ] );
        stateAmplitude[3] *= sqrt( factorial[ nPrime[y][p] ] );

    }

    return;

}

inline void LinearOpticalTransform::setStateAmplitude(std::complex<double> stateAmplitude[],Eigen::MatrixXcd& U,int y){

    std::complex<double> UProdTemp(1.0,0.0);

    for(int i=0;i<ancillaPhotons;i++) UProdTemp *= U( i,mPrime[y][i] );

    stateAmplitude[0] += UProdTemp * ( U(ancillaModes,mPrime[y][ancillaPhotons]) * U(ancillaModes+2,mPrime[y][ancillaPhotons+1])
                                    + U(ancillaModes + 1,mPrime[y][ancillaPhotons]) * U(ancillaModes + 3,mPrime[y][ancillaPhotons+1]) );

    stateAmplitude[1] += UProdTemp * ( U(ancillaModes,mPrime[y][ancillaPhotons]) * U(ancillaModes+3,mPrime[y][ancillaPhotons+1])
                                    + U(ancillaModes + 1,mPrime[y][ancillaPhotons]) * U(ancillaModes + 2,mPrime[y][ancillaPhotons+1]) );

    stateAmplitude[2] += UProdTemp * ( U(ancillaModes,mPrime[y][ancillaPhotons]) * U(ancillaModes+2,mPrime[y][ancillaPhotons+1])
                                    - U(ancillaModes + 1,mPrime[y][ancillaPhotons]) * U(ancillaModes + 3,mPrime[y][ancillaPhotons+1]) );

    stateAmplitude[3] += UProdTemp * ( U(ancillaModes,mPrime[y][ancillaPhotons]) * U(ancillaModes+3,mPrime[y][ancillaPhotons+1])
                                    - U(ancillaModes + 1,mPrime[y][ancillaPhotons]) * U(ancillaModes + 2,mPrime[y][ancillaPhotons+1]) );

    return;

}

#endif // LINEAROPTICALTRANSFORM_H_INCLUDED

