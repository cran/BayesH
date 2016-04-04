///////////////////////////////////////////////////////////////////////// 
//
// bayesH.c 
//
// copyright (c) 2016-2016, Renato Rodrigues Silva
// last modified Mar, 2016
// first written Jan, 2016
//
//     This program is free software; you can redistribute it and/or
//     modify it under the terms of the GNU General Public License,
//     version 3, as published by the Free Software Foundation.
//
//     This program is distributed in the hope that it will be useful,
//     but without any warranty; without even the implied warranty of
//     merchantability or fitness for a particular purpose.  See the GNU
//     General Public License, version 3, for more details.
//
//     A copy of the GNU General Public License, version 3, is available
//     at http://www.r-project.org/Licenses/GPL-3
//
// Part of the BayesH package
// Contains: bayesH.c
// 
/////////////////////////////////////////////////////////////////////////// 


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include "util.h"
#include <unistd.h>
#include <errno.h>

//Gibbs sampling for Bayesian Hierarchical Regression Model with mixture of two t scaled
//distribution as prior of coefficients

SEXP bayesH(SEXP X, SEXP y, SEXP mu, SEXP beta, SEXP sigma, SEXP group, 
            SEXP nu0, SEXP s0, SEXP niter, SEXP nburnin){


   //////////////////////////////////////////////////////////////////////////////
   //  Declare the variables 
   //////////////////////////////////////////////////////////////////////////////

   //Iterators    
    int nprot = 0, i, j, iter;
        
   //Genotypic and Phenotypic data
    int *dimxptr, *dimyptr;
    double *xptr, *yptr; 
   //Number of iterations
    int *niterptr, *nburninptr;     
   //Others Parameters 
    int *dimbetaptr;
    double  *muptr, *betaptr, *sigmaptr, *nu0ptr, *s0ptr;
   
    
    double *resptr, *predptr;
    double *betaGibbsptr, *sigmaGibbsptr, *muGibbsptr;     
    double *predGibbsptr; 
    double *sxsqptr; 

    int *groupptr; 

    double sigmu=0.0;
    double sxytilde=0.0, sse=0;
    double betatilde=0.0, sigmatilde=0.0;
    double betabar=0.0, sigmabar=0.0;
    double mudot=0.0, sigmadot=0.0;   
    double df=0.0, scale=0.0;
    double *sigbeta1ptr;

    SEXP betaGibbs, sigmaGibbs, muGibbs;
    SEXP predGibbs; 

     
    double mutmp = 0.0;    
    double dfp1 = 0.0;
    double s0p1 = 0.0;  
    double dfp2 = 0.0;
    double s0p2 = 0.0; 
    double ssbstand = 0.0;


    double sigbetatmp=0.0, sigmatmp=0.0;
    //double epsilon = 0.0000001;
    
    //pi 
    double pi=0.5; 
    double alphaprior= 1.0;
    double gammaprior= 1.0;

    //z
    double r = 1.0000001;
    int *z; 
    double *probz;   
    double nz = 0.0;
    double a = 0.0;
    double b = 0.0;
    double ab = 0.0;  
  
    //double tmp = 0.0; 

    
    //Output
    SEXP output;

    /////////////////////////////////////////////////
    // Defining the output
    /////////////////////////////////////////////////
    PROTECT(output = allocVector(VECSXP, 4)); nprot ++;    

    ////////////////////////////////////////////////////////////////////////////////
    //  Read variables
    ////////////////////////////////////////////////////////////////////////////////

    //Read the genotypic data  
    PROTECT(X = AS_NUMERIC(X)); nprot ++;
    xptr = REAL(X);
    dimxptr = INTEGER(GET_DIM(X)); 
    //Read the phenotypic data  
    PROTECT(y = AS_NUMERIC(y)); nprot ++;
    yptr = REAL(y);
    dimyptr = INTEGER(GET_DIM(y));     
    //Read the mu parameter
    PROTECT(mu = AS_NUMERIC(mu)); nprot ++;
    muptr = REAL(mu);       
    //Read the beta parameter
    PROTECT(beta = AS_NUMERIC(beta)); nprot ++;
    betaptr = REAL(beta);
    dimbetaptr = INTEGER(GET_DIM(beta)); 
    //Read the sigma parameter 
    PROTECT(sigma = AS_NUMERIC(sigma)); nprot ++;
    sigmaptr = REAL(sigma);       
    //Read the nu0 parameter
    PROTECT(nu0 = AS_NUMERIC(nu0)); nprot ++;
    nu0ptr = REAL(nu0);    
    //Read the s0 parameter
    PROTECT(s0 = AS_NUMERIC(s0)); nprot ++;
    s0ptr = REAL(s0);    
    //Read the niter parameter
    PROTECT(niter = AS_INTEGER(niter)); nprot ++;
    niterptr = INTEGER(niter);
    //Read the nburnin parameter
    PROTECT(nburnin = AS_INTEGER(nburnin)); nprot ++;
    nburninptr = INTEGER(nburnin);

    //Read the group
    PROTECT(group = AS_INTEGER(group)); nprot ++;
    groupptr =INTEGER(group);
   

    ///////////////////////////////////////////////////////////////////////////////////
    //  Load and Initialize Variables
    ///////////////////////////////////////////////////////////////////////////////////

    double G = (double) niterptr[0] + 0.0;     
    double N = (double) dimyptr[0] + 0.0; 
    double p = (double) dimbetaptr[0] + 0.0;
 
    double tmpRate = 0.0;
    double tmpShape = 0.0; 
    double scaleUsed = s0ptr[1];

  
   //Load residuals  
    resptr =(double*) malloc( dimyptr[0]*sizeof(double));    
    for(i = 0; i < dimyptr[0]; i++){       
       resptr[i] = yptr[i] - muptr[0];
    }  
     
    //Load pred
    predptr =(double*) malloc( dimyptr[0]*sizeof(double));    
    for(i = 0; i < dimyptr[0]; i++){
       predptr[i] = muptr[0];
    } 

   
    //Load the sumxsq
    sxsqptr = (double*)malloc(dimbetaptr[0]*sizeof(double));   
    for(j=0; j < dimbetaptr[0]; j++){
        sxsqptr[j] = 0.0;
    }

    //Define sumxsq
    for(j=0; j < dimbetaptr[0]; j++){
       for(i=0; i < dimyptr[0]; i++){ 
          sxsqptr[j] = sxsqptr[j] + xptr[(dimyptr[0]*j+i)]*xptr[(dimyptr[0]*j+i)];                           
       } 
    } 

    //Load betaGibbs
    PROTECT(betaGibbs = allocMatrix(REALSXP, dimbetaptr[0],1)); nprot ++; 
    betaGibbsptr = REAL(betaGibbs);  
    for(j = 0; j < dimbetaptr[0]; j++){
       betaGibbsptr[j] = 0.0;
    } 

    //Load sigmaGibbs
    PROTECT(sigmaGibbs = allocMatrix(REALSXP, niterptr[0],1)); nprot ++; 
    sigmaGibbsptr = REAL(sigmaGibbs);
    for(iter = 0; iter < niterptr[0]; iter++){
       sigmaGibbsptr[iter] = 0.0;
    } 

    //Load muGibbs
    PROTECT(muGibbs = allocMatrix(REALSXP, niterptr[0],1)); nprot ++; 
    muGibbsptr = REAL(muGibbs);
    for(iter = 0; iter < niterptr[0]; iter++){
       muGibbsptr[iter] = 0.0;
    }    
   
    
   //Load predGibbs
    PROTECT(predGibbs = allocMatrix(REALSXP, dimyptr[0],1)); nprot ++; 
    predGibbsptr = REAL(predGibbs);
    for(i = 0; i < dimyptr[0]; i++){
       predGibbsptr[i] = 0.0;
    }   

   //Load z
    z = (int*) malloc(dimbetaptr[0]*sizeof(int));
    for(j = 0; j < dimbetaptr[0]; j++){
       GetRNGstate(); 
       z[j] = rbinom(1, pi);
       //Rprintf("%d\n", z[j]); 
       PutRNGstate();

    }  

    //Load probz
    probz = (double*) malloc(dimbetaptr[0]*sizeof(double));
    for(j = 0; j < dimbetaptr[0]; j++){
       probz[j] = 0.0;
    }  
  
   //Load the sigmabeta1
    sigbeta1ptr = (double*)malloc(dimbetaptr[0]*sizeof(double));   
    for(j=0; j < dimbetaptr[0]; j++){
        sigbeta1ptr[j] = 1; 
    }

   

    ///////////////////////////////////////////////////////////////////////////////
    //Gibbs Sampling
    ///////////////////////////////////////////////////////////////////////////////
   
     Rprintf("Gibbs Sampling running, please wait!\n"); 
      
     for(iter = 0; iter < (niterptr[0]+nburninptr[0]); iter++){
         /*if(iter >= nburninptr[0]){              
             printf("%d\n", iter-nburninptr[0]);            
         }*/  
         
         mutmp = 0.0;         
         sse = 0.0; 
         ssbstand = 0.0; 
         sigbetatmp = sigmatmp = 0.0;

         //z
         nz = 0.0;
         a = 0.0;
         b = 0.0;
         ab = 0.0;
          
       ///////////////////////////////////////////////////////////////////////
       //Step 1- Generate intercept from Gaussian distribution 
       //////////////////////////////////////////////////////////////////////
       //Residuals
         for(i=0; i < dimyptr[0]; i++){  
             resptr[i] = resptr[i] + muptr[0];
         }       
         
         for(i=0; i < dimyptr[0]; i++){ 
            mutmp = mutmp + resptr[i]; 
         }                                 
         
         mudot =  (mutmp / N  ); 
         //Rprintf("%lf\n",  mudot);  
         sigmadot = sigmaptr[0] / N  ;               
         GetRNGstate();  
            muptr[0] =   rnorm(  mudot, sqrt( sigmadot));  
            //Rprintf("%lf\n",  muptr[0]);                         
            if(iter >= nburninptr[0]){              
               muGibbsptr[(iter-nburninptr[0])] = muptr[0];             
            }       
         PutRNGstate();
         for(i=0; i < dimyptr[0]; i++){
            resptr[i] = resptr[i]  - muptr[0]; 
         }
                 
                 
       ///////////////////////////////////////////////////////////////////////
       //Step 2 - For each marker, generate beja_j from Gaussian distribution
       //mean= (1/varE x_j' e)/ (1/varE x_j' x_j + 1/(varE tau_j^2))
       //variance= 1/ (1/varE x_j' x_j + 1/(varE tau_j^2)) 
       //////////////////////////////////////////////////////////////////////      

           for(j=0; j < dimbetaptr[0]; j++){  
            sxytilde = 0.0;
            for(i=0; i < dimyptr[0]; i++){ 
               resptr[i] =  resptr[i] + betaptr[j]*xptr[(dimyptr[0]*j+i)];
            }
                                
          //Define sumxytilde
            for(i=0; i < dimyptr[0]; i++){
               sxytilde = sxytilde+ xptr[(dimyptr[0]*j+i)]*resptr[i]; 
            }                           
               
           //Define betatilde               
            betatilde = (sxytilde/sigmaptr[0])/((sxsqptr[j]/sigmaptr[0])+(1/(sigbeta1ptr[j]*sigmaptr[0]*r)));                             
           //Define sigmatilde                 
            sigmatilde = 1/((sxsqptr[j]/sigmaptr[0])+(1/(sigbeta1ptr[j]*sigmaptr[0]*r))); 
            GetRNGstate(); 
            //Generate beta_j                    
            betaptr[j] = rnorm(betatilde,sqrt(sigmatilde)); 
            PutRNGstate();       
            //ssb  = ssb + (betaptr[j])*(betaptr[j]);  
            ssbstand  = ssbstand +  (betaptr[j]/sqrt(sigbeta1ptr[j]))*(betaptr[j]/sqrt(sigbeta1ptr[j]));  
            
           if(iter > nburninptr[0]){                         
              betaGibbsptr[j] = betaGibbsptr[j] + betaptr[j];                                     
           }               
           for(i=0; i < dimyptr[0]; i++){   
               resptr[i] =  resptr[i] - betaptr[j]*xptr[(dimyptr[0]*j+i)];
           }             
                                    
        }       
                 
        
        
       ///////////////////////////////////////////////////////////////////////
       //Step 3 - For each marker, generate the variance of marker
       //////////////////////////////////////////////////////////////////////         
         
         for(j=0; j < dimbetaptr[0]; j++){                               
               
            if(z[j] == 1){        
               //Generate marker variance
               dfp1 = nu0ptr[0] + 1;
               s0p1 = (pow(betaptr[j]/sqrt(sigmaptr[0]), 2.0)  + s0ptr[0]*nu0ptr[0]); 
               GetRNGstate(); 
               sigbeta1ptr[j] = s0p1 / rchisq(dfp1); 
               PutRNGstate();                         
           } else {             
               //Generate marker variance
               dfp2 = nu0ptr[1] + 1;
               s0p2 = (pow(betaptr[j]/sqrt(sigmaptr[0]), 2.0)  + s0ptr[1]*nu0ptr[1]);                
               GetRNGstate(); 
               sigbeta1ptr[j] = s0p2 / rchisq(dfp2); 
               PutRNGstate();                   
                                          
          }                         
             // Rprintf("%lf\n",  sigbeta1ptr[j]);                        
        }  
         
        //////////////////////////////////////////////////////////////////////////
        //Step 4 - For each marker, simulate Z        
        /////////////////////////////////////////////////////////////////////////// 

        //Rprintf("%lf\n", pi);
        GetRNGstate(); 
        for(j = 0; j < dimbetaptr[0]; j++){           
           a = dbinom(1, 1.0, pi, 1);
           //a = a + dnorm( betaptr[j], 0,  sqrt(sigbeta1ptr[j]), 1);
           a = a + dgamma( (1 / sigbeta1ptr[j]), (nu0ptr[0] / 2.0), (2.0/(nu0ptr[0] * s0ptr[0])) ,1) + log(1 / pow(sigbeta1ptr[j], 2.0) );
           a = exp(a);
           b = dbinom(0, 1.0, pi, 1);
           //b = b + dnorm( betaptr[j], 0,  sqrt(sigbeta1ptr[j]), 1);
           b = b + dgamma( (1 / sigbeta1ptr[j]), (nu0ptr[1] / 2.0), (2.0/(nu0ptr[1] * s0ptr[1])) ,1) + log(1 / pow(sigbeta1ptr[j], 2.0) );  
           b = exp(b);
           //probz[j] = a/ (a + b);
           //Rprintf("%lf\n", probz[j]);  
           z[j] = rbinom(1, probz[j]);
           if(z[j] == 1){
               nz = nz +  1.0;                  
           }                
        }
       PutRNGstate(); 

        ////////////////////////////////////////////////////////////////////////
        //Step 5 - Estimate pi 
        // pi = mixing probability of variables included in the model
        ///////////////////////////////////////////////////////////////////////

          GetRNGstate();  
          pi = rbeta((alphaprior + nz), (gammaprior + (p - nz)) ); 
          //Rprintf("%lf\n", pi); 
          PutRNGstate();  
   
         //Predicted values
         for(i=0; i < dimyptr[0]; i++){
            predptr[i] = yptr[i] - resptr[i]; 
         }    
    
            
         /////////////////////////////////////////////////////////////////////////
         //Step 6 - Generate sigma from scaled inverse chi-squared distribution
         //////////////////////////////////////////////////////////////////////// 
        
         //Predicted values
         //Residuals
         for(i=0; i < dimyptr[0]; i++){
            sse = sse + resptr[i]*resptr[i];              
         }
         
         //Inverse Scaled Inverse           
         df = (N + nu0ptr[2] +p );    
         scale = ( sse  + ssbstand + s0ptr[0]*nu0ptr[0]);  
         //scale = ( sse  + s0ptr[2]*nu0ptr[2]);                           
         //Generate values from scaled inverse chi-squared distribution
         GetRNGstate();       
         sigmaptr[0]  =  (scale) /  rchisq(df);                               
         if(iter >= nburninptr[0]){             
            sigmaGibbsptr[(iter-nburninptr[0])] = sigmaptr[0];               
         }      
         PutRNGstate();
    
        
         for(i=0; i < dimyptr[0]; i++){              
            if(groupptr[i] == 1){                                   
                yptr[i] =  predptr[i] + rnorm(0, sqrt(sigmaptr[0])); 
                resptr[i] = yptr[i] - predptr[i];               
            }              
         } 

       
         for(i=0; i < dimyptr[0]; i++){
            if(iter >= nburninptr[0]){             
               predGibbsptr[i] = predGibbsptr[i] + predptr[i];                
            }
         } 
               
     } //end of Gibbs Sampling  
     for(j=0; j < dimbetaptr[0]; j++){
       betaGibbsptr[j] = betaGibbsptr[j]*(1 / G + 0.0);
     }
     for(i=0; i < dimyptr[0]; i++){
        predGibbsptr[i] = predGibbsptr[i]*(1 / G + 0.0);
     }
     //Output
     SET_VECTOR_ELT(output, 0, sigmaGibbs);
     SET_VECTOR_ELT(output, 1, muGibbs);
     SET_VECTOR_ELT(output, 2, betaGibbs); 
     SET_VECTOR_ELT(output, 3, predGibbs);     
     free(resptr);
     free(predptr); 
     free(z);
     free(probz);
     free(sxsqptr);
     free(sigbeta1ptr);          
     UNPROTECT(nprot);         
     return output;      
  
}

