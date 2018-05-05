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
        void setA(Eigen::MatrixXcd& U);

        Eigen::MatrixXcd A;

    private:

        std::vector< std::vector<int> > n,m,nPrime,mPrime;
        std::vector<double> factorial;
        std::vector<bool> useRysers;
        int photons;
        GrayCode graycode;

        template <typename T>
        void printVec(std::vector<T>& a);

        int numbPermutations(int& i);
        double doublefactorial(int x);
        void permutationAlgorithm(Eigen::MatrixXcd& U,int& i);
        void setmVec(std::vector<int>& m, std::vector<int>& n);
        void rysersAlgorithm(Eigen::MatrixXcd& U,int& i);

};

#endif // LINEAROPTICALTRANSFORM_H_INCLUDED

