/* % QuantiSNP v1.0
% --------------
% An Objective Bayesian Hidden Markov Model approach for detecting copy
% number variation from SNP genotyping data.
%
% Copyright (C) 2007  The University of Oxford.
%
% For all enquiries please consult the contacts page on the program website: 
%   http://www.well.ox.ac.uk/quantisnp 
% */

#include "mex.h"
#include "matrix.h"
    
/* If you are using a compiler that equates NaN to zero, you must
* compile this example using the flag -DNAN_EQUALS_ZERO. For 
* example:
*
*     mex -DNAN_EQUALS_ZERO findnz.c  
*
* This will correctly define the IsNonZero macro for your
compiler. */

#if NAN_EQUALS_ZERO
    #define IsNonZero(d) ((d) != 0.0 || mxIsNaN(d))
#else
    #define IsNonZero(d) ((d) != 0.0)
#endif

double normalise(double *array, int elements);
double* normaliseVec(double *array, int elements);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    /* Declare variables. */ 
    int i, j, k;
    int count, nStates, nSnps, ndims;
    int dims[3]; 
	int* dimsPt;
    double *gammaVecPt, *beta, *betaOld, *betaMatrix, *betaScale, *alpha, *alphaOld, *alphaMatrix, *alphaScale, *gammaMatrix;
    double *etaMat2DPt, *etaMatrix;
    double *priorMatrix, *transitionMatrix, *obsMatrix, *tMatrix;  
    mxArray *tMat, *alphaVec, *alphaVecOld, *betaVecOld, *betaVec, *gammaVec, *etaMat2D;
    double Z;      
    
    /* Get the data. */
    priorMatrix         = (double *)mxGetPr(prhs[0]);
    transitionMatrix    = (double *)mxGetPr(prhs[1]);
    obsMatrix           = (double *)mxGetPr(prhs[2]);
    nStates             = (int )mxGetScalar(prhs[3]);
    nSnps               = (int )mxGetScalar(prhs[4]);
    
    plhs[0] = mxCreateDoubleMatrix(nStates, nSnps, mxREAL);
    alphaMatrix = mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(1, nSnps, mxREAL);
    alphaScale = mxGetPr(plhs[1]);
    
    tMat = mxCreateDoubleMatrix(nStates, nStates, mxREAL);
    tMatrix = mxGetPr(tMat);
    
    alphaVec = mxCreateDoubleMatrix(nStates, 1, mxREAL);
    alpha = mxGetPr(alphaVec);

    alphaVecOld = mxCreateDoubleMatrix(nStates, 1, mxREAL);
    alphaOld = mxGetPr(alphaVecOld);
    
    plhs[2] = mxCreateDoubleMatrix(nStates, nSnps, mxREAL);
    betaMatrix = mxGetPr(plhs[2]);
    
    plhs[3] = mxCreateDoubleMatrix(1, nSnps, mxREAL);
    betaScale = mxGetPr(plhs[3]);
    
    plhs[4] = mxCreateDoubleMatrix(nStates, nSnps, mxREAL);
    gammaMatrix = mxGetPr(plhs[4]);
    
    tMat = mxCreateDoubleMatrix(nStates, nStates, mxREAL);
    tMatrix = mxGetPr(tMat);
    
    alphaVec = mxCreateDoubleMatrix(nStates, 1, mxREAL);
    alpha = mxGetPr(alphaVec);

    alphaVecOld = mxCreateDoubleMatrix(nStates, 1, mxREAL);
    alphaOld = mxGetPr(alphaVecOld);
    
    betaVec = mxCreateDoubleMatrix(nStates, 1, mxREAL);
    beta = mxGetPr(betaVec);

    betaVecOld = mxCreateDoubleMatrix(nStates, 1, mxREAL);
    betaOld = mxGetPr(betaVecOld);  
    
    gammaVec = mxCreateDoubleMatrix(nStates, 1, mxREAL);
    gammaVecPt = mxGetPr(gammaVec);     
    
    /* forward algorithm  */
    for ( i = 0; i < nStates; i++ ) {
        alpha[i] = priorMatrix[i]*obsMatrix[i];
    }
    alphaScale[0] = normalise(alpha, nStates);
    for ( i = 0; i < nStates; i++ ) {
        alphaMatrix[i] = alpha[i]/alphaScale[0];
    }
    
    for ( k = 1; k < nSnps; k++ ) {
        
        for ( i = 0; i < nStates; i++ ) {
            for ( j = 0; j < nStates; j++ ) {
                tMatrix[nStates*i + j] = transitionMatrix[nStates*nStates*(k-1) + nStates*i + j];
            }

            alphaOld[i] = alphaMatrix[nStates*(k-1) + i];
        }
        
        Z = 0.0;
        for ( i = 0; i < nStates; i++ ) {
            alpha[i] = 0.0;   
        }
        for ( j = 0; j < nStates; j++ ) {
            for ( i = 0; i < nStates; i++ ) {
                alpha[j] = alpha[j] + tMatrix[nStates*j + i]*alphaOld[i];
            }    
        }
        for ( j = 0; j < nStates; j++ ){
            alpha[j] = alpha[j]*obsMatrix[nStates*k + j]; 
        }
        alphaScale[k] = normalise(alpha, nStates);
        for ( j = 0; j < nStates; j++ ) {
            alphaMatrix[nStates*k + j] = alpha[j]/alphaScale[k];
        }
      
    }
 
 
    /* backward algorithm */
    for ( i = 0; i < nStates; i++ ) {
        beta[i] = 1;
    }
    betaScale[nSnps-1] = 1;
    for ( i = 0; i < nStates; i++ ) {
        betaMatrix[nStates*(nSnps-1) + i] = beta[i]/betaScale[nSnps-1];
    }
    
    for ( k = nSnps-2; k > -1; k-- ) {
        
        for ( i = 0; i < nStates; i++ ) {
            for ( j = 0; j < nStates; j++ ) {
                tMatrix[nStates*j + i] = transitionMatrix[nStates*nStates*k + nStates*j + i];
            }
            betaOld[i] = betaMatrix[nStates*(k+1) + i];
        }
        
        for ( i = 0; i < nStates; i++ ) {
            beta[i] = 0.0;   
        }
        for ( j = 0; j < nStates; j++ ){
            betaOld[j] = betaOld[j]*obsMatrix[nStates*(k+1) + j]; 
        }
        for ( j = 0; j < nStates; j++ ) {
            for ( i = 0; i < nStates; i++ ) {
                beta[i] = beta[i] + tMatrix[nStates*j + i]*betaOld[j];
            }    
        }
        
        betaScale[k] = normalise(beta, nStates);
        for ( j = 0; j < nStates; j++ ) {
            betaMatrix[nStates*k + j] = beta[j]/betaScale[k];
        }
      
    }
    
    /* compute marginal posterior probabilities */
    for ( k = 0; k < nSnps; k++ ) {
	    
	    for ( i = 0; i < nStates; i++ ) {
	       gammaVecPt[i] = alphaMatrix[nStates*k + i]*betaMatrix[nStates*k + i];
	 	}
	 	
	 	gammaVecPt = normaliseVec(gammaVecPt, nStates);
	 	
	 	for ( i = 0; i < nStates; i++ ) {
	       gammaMatrix[nStates*k + i] = gammaVecPt[i];
	 	}
	 	
	}    
	
	/* compute marginal switching probability */
	ndims = 3;
    dims[0] = nStates;
    dims[1] = nStates;
    dims[2] = nSnps-1;
	dimsPt = dims;
    
    plhs[5] = mxCreateNumericArray(ndims, dimsPt, mxDOUBLE_CLASS, mxREAL);
    etaMatrix = mxGetPr(plhs[5]);
    etaMat2D = mxCreateDoubleMatrix(nStates, nStates, mxREAL);
    etaMat2DPt = mxGetPr(etaMat2D); 
 
    for ( k = 0; k < nSnps-1; k++ ) {
        
        Z = 0;
        for ( i = 0; i < nStates; i++ ) {
            for ( j = 0; j < nStates; j++ ) {
                
                tMatrix[nStates*j + i] = transitionMatrix[nStates*nStates*k + nStates*j + i];
                
                etaMat2DPt[nStates*j + i] = 
                    alphaMatrix[nStates*k + i]*tMatrix[nStates*j + i]*obsMatrix[nStates*(k+1) + j]*betaMatrix[nStates*(k+1) + j];
                
                Z = Z + etaMat2DPt[nStates*j + i];
                
            }
  
      }
        
        for ( i = 0; i < nStates; i++ ) {
            for ( j = 0; j < nStates; j++ ) {
                etaMatrix[nStates*nStates*k + nStates*i + j] = etaMat2DPt[nStates*i + j]/Z;
            }
        } 
       
        
    }  
 
}



double normalise( double *array, int elements ) {
    
    int i;
    double sum = 0.0;
    
    for ( i = 0; i < elements; i++ ) {
        sum = sum + array[i];   
    }
    
    return(sum);
    
}

double* normaliseVec( double *array, int elements ) {
    
    int i;
    double sum = 0.0;
    
    for ( i = 0; i < elements; i++ ) {
        sum = sum + array[i];   
    }
    
    for ( i = 0; i < elements; i++ ) {
        array[i] = array[i]/sum;   
    }
    
    return(array);
    
}
