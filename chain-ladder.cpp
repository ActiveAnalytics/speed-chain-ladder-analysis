/**
 * @title Speed Chain Ladder Analysis with Rcpp
 * @author Chibisi Chima-Okereke
 * @license GPL (>= 2)
 * @tags modeling armadillo
 * @summary Speed-up of Chain Ladder analysis by calling C++ routines from R
 */

/**
 * The Chain Ladder method is an actuarial technique used for
 * projecting incurred insurance claims to their ultimate loss
 * values. The data exists as claims triangles where the claims for
 * each accounting year increments down the rows and the claims for
 * each development period increments along the columns.  This claims
 * triangle can be represented in a triangular upper matrix (along the
 * anti-diagonal) and the Chain Ladder technique works by filling in
 * the lower part of the matrix using ratios of claims in previous
 * accounting years and development periods.
 * 
 * In this example, we show how an implementation in R is sped up by
 * calling an equivaluent implementation in C++ from R using the Rcpp
 * interface.
 * 
 * We start with the C++ code for carrying out the Chain Ladder calculation.
 *
 */

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>  // most of the algorithms 
#include <numeric> // some numeric algorithm

using namespace std;
using namespace arma;
using namespace Rcpp;

// Code for the age-to-age factor when the column index is given 
double GetFactor(int index, mat mTri) {
    int nRow = mTri.n_rows;
    mat subMat = mTri.submat(0, index, nRow - (index + 2), index + 1);
    rowvec colSums = arma::sum(subMat, 0);
    double inFact = colSums(1)/colSums(0);
    return inFact;
}

// Code for getting all the factors from the triangle 
vec GetFactors(mat mTri) {
    int nCol = mTri.n_cols;
    vec dFactors(nCol - 1);
    for (int i=0; i < nCol - 1; ++i) {
        dFactors(i) = GetFactor(i, mTri);
    }
    return dFactors;
}

// This is code for the cumulative product of a vector 
vec cumprod(vec mvec) {
    int nElem = mvec.n_elem;
    double cprod = mvec(0);
    for (int i = 1; i < nElem; ++i) {
        cprod *= mvec(i);
        mvec(i) = cprod;
    }
    return mvec;
}

/**
 * The following function returns the  fully projected triangle 
 */

// [[Rcpp::export]]
SEXP GetChainSquareCpp(SEXP mClaimTri) {
    NumericMatrix nMat(mClaimTri);
    int nRow = nMat.nrow(), nCol = nMat.ncol();
    mat armMat(nMat.begin(), nRow, nCol, FALSE);

    vec dFactors = GetFactors(armMat);
    mat revMat = fliplr(armMat);
    vec dAntiDiag = diagvec(revMat);
    dAntiDiag = dAntiDiag.subvec(1, nCol - 1);
    double dMult;
    vec prodVec;
    for (unsigned int index = 0; index < dAntiDiag.n_elem; ++index) {
        dMult = dAntiDiag(index);
        prodVec = dFactors.subvec(nCol - index - 2, nCol - 2);
        prodVec = cumprod(prodVec);
        armMat(index + 1, span(nCol - index - 1, nCol - 1)) = dMult*prodVec.st();
    }
    return wrap(armMat);
}
