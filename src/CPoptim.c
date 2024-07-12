// File CPoptim.R
// Part of the CPoptim R package
// Copyright 2013 Erick G.G. de Paz, Humberto Vaquera Huerta, Francisco J Albores Velasco, John R. Bauer Mengelberg, Juan Manuel Romero Padilla
// Distributed under GPL 2

/* Reviewed 09/09/2023 Notes for RStudio Users
Compile this code with the following command:
R CMD SHLIB CPoptim.c

Run this code in RStudio
 tools::package_native_routine_registration_skeleton(".")
*/

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rmath.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP CPoptim(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"CPoptim", (DL_FUNC) &CPoptim, 15},
  {NULL, NULL, 0}
};

void R_init_CPoptim(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}


/* Auxiliary routine for evaluating an objective function coded in R */
int increaseSampleXY(int maxFE, int dimension, double X[maxFE][dimension],
                     double Y[maxFE], int ID[maxFE],
                    double minX[][dimension], double maxX[][dimension],
                    int curLenY, int curPart, int sampleSz,
                    SEXP rFunction, SEXP env);

const int minimalSampleInAGroup=10;

/**
* Implementation of Algorithm 1
*
* @param argInt    An R vector that contains:
*                    1) number of real parameters of the obj. function
*                    2) the sample size by partition (default 1000)
*                    3) number of function evaluations to compute
*                    4) number of trees to compute (default 1)
* @param minimalX  minimal value for each real parameter
* @param maximalX  maximal value for each real parameter
* @param rFunction an objective function coded in R
* @param env       R environment to run the objective function
* @return          An R vector that contains all fun. evaluations
*/
SEXP CPoptim(SEXP argInt, SEXP minimalX,
             SEXP maximalX, SEXP rFunction, SEXP env,
SEXP    aPriori_arg,
SEXP    minCondP_arg,
SEXP    hatMu_arg,
SEXP    hatSD_arg,
SEXP    minX_arg,
SEXP    maxX_arg,
SEXP    X_arg,
SEXP    Y_arg,
SEXP    ID_arg,
SEXP    N_arg)
{
    //////////////////////////////////////////////////////////////
    //       Decomposing R-readable vectors into C variables    //

	int *auxInt=INTEGER(argInt);
    double *auxDouble;

	int dimension=auxInt[0];
	int sampleSize=auxInt[1];
	int maxFE=auxInt[2];
	int forestSz=auxInt[3];

	double minXorg[dimension];
	double maxXorg[dimension];

	int curPart,i,j;

    GetRNGstate();
	auxDouble=REAL(minimalX);
	for(i=0;i<dimension;i++)
		minXorg[i]=auxDouble[i];

	auxDouble=REAL(maximalX);
	for(i=0;i<dimension;i++)
		maxXorg[i]=auxDouble[i];

    int partitions=maxFE/minimalSampleInAGroup+1;

    //           END OF Decomposing R-readable vectors          //
    //////////////////////////////////////////////////////////////


	////////////////////////////////////////////////////////////
    ///                    Memory Allocation                ////
    double *aPriori=(double*)REAL(aPriori_arg);
    double *minCondP=(double*)REAL(minCondP_arg);
    double *hatMu=(double*)REAL(hatMu_arg);
    double *hatSD=(double*)REAL(hatSD_arg);
    double (*minX)[dimension]=(double (*)[dimension])REAL(minX_arg);
    double (*maxX)[dimension]=(double (*)[dimension])REAL(maxX_arg);
    double (*X)[dimension]=(double (*)[dimension])REAL(X_arg);
    double *Y=(double*)REAL(Y_arg);
    int *ID=(int*)INTEGER(ID_arg);
    int *N=(int*)INTEGER(N_arg);

    for(curPart=0;curPart<partitions;curPart++)
    {
        aPriori[curPart]=0;
        hatMu[curPart]=0;
        hatSD[curPart]=1;
        minCondP[curPart]=0;
    }

    int curLenY=0;
    for(i=0;i<maxFE;i++)
        ID[i]=-1;
    ///                  END of Memory Allocation           ////
    ////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////
    ///               Algorithm 1 (lines 1-2)               ////
    /*
        By default forestSz=1, so the following loop iterates
        only one time. For more details, see Section 4 Future
        Works (Point 1).
    */
    for(curPart=0;curPart<forestSz;curPart++)
    {
        // Algorithm 1 (line 1):
        // The first partition is the original domain
        for(i=0;i<dimension;i++)
        {
            minX[curPart][i]=minXorg[i];
            maxX[curPart][i]=maxXorg[i];
        }
        aPriori[curPart]=1.0/forestSz;

        // Algorithm 1 (line 2):
        // An uniform sample is computed over the original domain
        curLenY=curLenY+increaseSampleXY(maxFE,dimension,X,Y,ID,
                                         minX,maxX,curLenY,curPart,
                                         sampleSize/2,rFunction,env);
    }
    ///               END OF Algorithm 1 (lines 1-2)         ////
    /////////////////////////////////////////////////////////////



    ///////////////////////////////////////////////////////////////
    ////// Algorithm 1 (line 3): Computing initial statistics /////
    double mean=0;
    for(i=0;i<curLenY;i++)
        mean=mean+Y[i];
    mean=mean/curLenY;

    double var=0;
    for(i=0;i<curLenY;i++)
        var=var+R_pow(Y[i]-mean,2.0);
    var=var/curLenY;

    double minMean=mean;
    for(i=0;i<forestSz;i++)
    {
        hatMu[i]=mean;
        hatSD[i]=sqrt(var/sampleSize);
        minCondP[i]=.5;
    }
    ////////////// END of Computing initial statistics ///////////
    //////////////////////////////////////////////////////////////

    int selPart,lenSubID;
    double bestCriterion;

    int subID[sampleSize];
    double subX[sampleSize];

    double score, bestScore, tot, area;
    double sum[2];
    double P[2];
    double bestMu[2];
    double bestP[2];
    double bestX;
    int bestDim;
    int memSize,autre;

    ///////////////////////////////////////////////////////////////
    //// This For-loop corresponds to Algorithm 1 (lines 3-16) ////
    for(curPart=forestSz;curPart<partitions;curPart++)
    {

        // Algorithm 1 (line 5):
        // The most promessing subset to split is selected by
        // finding the subset with the greatest criterion.
        //
        // See Conjecture 3 for details about the selection
        // criterion
        bestCriterion=-1;
        for(i=0;i<curPart;i++)
            if(minCondP[i]*aPriori[i]>bestCriterion)
            {
                bestCriterion=minCondP[i]*aPriori[i];
                selPart=i;
            }


        // Algorithm 1 (line 6):
        // Select previous observations belonging to the selected subset.
        // The vector subID contains references to these observations
        lenSubID=0;
        for(i=0;i<curLenY;i++)
            if(ID[i]==selPart)
                subID[lenSubID++]=i;


        // Algorithm 1 (line 7-10):
        // Compute new observations and merge them with the previous
        // observations to complete a standard sample size 'sampleSize'.
        curLenY=curLenY+increaseSampleXY(maxFE,dimension,X,Y,ID,
                                         minX,maxX,curLenY,selPart,
                                         sampleSize-lenSubID,rFunction,env);


        // Algorithm 1 (line 11-13):
        // If the limit of function evaluations maxFE is reach,
        // break the loop
        if(curLenY>=maxFE)
            break;


        // References to observations in the selected subset are updated
        // to include the observations just created by increaseSampleXY.
        for(   ;i<curLenY;i++)
            if(ID[i]==selPart)
                subID[lenSubID++]=i;


        ////////////////////////////////////////////////////////////////
        ////// Algorithm 1 (line 14): Spliting the selected subset /////
        tot=0;
        for(j=0; j<sampleSize; j++)
            tot=tot+Y[subID[j]];

        // Find the best axis and interval to create a new partition
        bestScore=-1;
        for(i=0; i<dimension; i++)
        {
            // Sort subID according to the current dim
            for(j=0; j<sampleSize; j++)
                subX[j]=X[subID[j]][i];
            rsort_with_index(subX,subID,sampleSize);

            // Compute statistics for the minimal partition
            sum[0]=0;
            for(j=0; j<minimalSampleInAGroup; j++)
                sum[0]=sum[0]+Y[subID[j]];

            // Update statistics progresively for different partitions
            area=maxX[selPart][i]-minX[selPart][i];
            for(;j<sampleSize-minimalSampleInAGroup;j++)
            {
                sum[1]=tot-sum[0];
                P[0]=(X[subID[j]][i]-minX[selPart][i])/area;
                P[1]=(maxX[selPart][i]-X[subID[j]][i])/area;

                // See Definition 4 (Comparative score)
                score=P[0]*R_pow_di(sum[0]/j,2)+
                      P[1]*R_pow_di(sum[1]/(sampleSize-j),2);

                if(score>bestScore)
                {
                    bestScore=score;

                    bestMu[0]=sum[0]/j;
                    bestMu[1]=sum[1]/(sampleSize-j);

                    bestP[0]=P[0];  bestP[1]=P[1];

                    bestX=X[subID[j]][i];
                    bestDim=i;
                    memSize=j;
                }

                sum[0]=sum[0]+Y[subID[j]];
            }
        }

        // The best partition is a Convex partition (See Definition 5)
        // with proportion bestX/area in the index(dimension) bestDim.

        // This partition is performed (See Formulae 2-3)
        for(i=0;i<dimension;i++)
        {
            maxX[curPart][i]=maxX[selPart][i];
            minX[curPart][i]=minX[selPart][i];
        }
        maxX[selPart][bestDim]=bestX;
        minX[curPart][bestDim]=bestX;

        //////////////////// END of Spliting method //////////////////
        //////////////////////////////////////////////////////////////


        //////////////////////////////////////////////////////////////////
        ////// Algorithm 1 (line 15): Replace the selected subset by /////
        ////// the resultant subsets after splitting.               //////

        hatMu[selPart]=bestMu[0];
        hatMu[curPart]=bestMu[1];

        lenSubID=0;
        autre=0;
        hatSD[curPart]=0;
        hatSD[selPart]=0;

        // Sort subID according to the BEST dim
        for(j=0; j<sampleSize; j++)
            subX[j]=X[subID[j]][bestDim];
        rsort_with_index(subX,subID,sampleSize);

        // Distribute the inner sample to the resultant subsets
        // after splitting
        for(i=0; i<sampleSize; i++)
            if(i<memSize)
            {
               hatSD[selPart]=hatSD[selPart]+R_pow(Y[subID[i]]-hatMu[selPart],2.0);
               lenSubID=lenSubID+1;
            }
            else
            {
              ID[subID[i]]=curPart;
              hatSD[curPart]=hatSD[curPart]+R_pow(Y[subID[i]]-hatMu[curPart],2.0);
            }

        N[curPart]=sampleSize-lenSubID;
        N[selPart]=lenSubID;
        ////////////////////  END of Replacing method  //////////////////
        /////////////////////////////////////////////////////////////////


        //////////////////////////////////////////////////////////////
        ////// Algorithm 1 (line 3): Updating initial statistics /////

        // See Formula 8
        hatSD[curPart]=R_pow(hatSD[curPart]/((N[curPart]-1)*(sampleSize)),.5);
        hatSD[selPart]=R_pow(hatSD[selPart]/((N[selPart]-1)*(sampleSize)),.5);

        aPriori[curPart]=aPriori[selPart]*bestP[1];
        aPriori[selPart]=aPriori[selPart]*bestP[0];

        //// Finding the minimal mu
        minMean=hatMu[0];
        for(i=1;i<=curPart;i++)
            if(minMean>hatMu[i])
                minMean=hatMu[i];

        // See Formula 9
        for(i=0;i<=curPart;i++)
            minCondP[i]=pnorm(minMean,hatMu[i],hatSD[i],1,0);

        //////////////////// END of Updating method ///////////////////
        //////////////////////////////////////////////////////////////


    }
    ////            END OF LOOP Algorithm 1 (lines 3-16)             ////
    /////////////////////////////////////////////////////////////////////



    /////////////////////////////////////////////////////////////
    ////    Creating a R-readable matrix as final result     ////
	SEXP ans=allocVector(REALSXP,1);
    auxDouble=REAL(ans);
	auxDouble[0]=curPart;
    ////           END of Creating a R-readable matrix       ////
    /////////////////////////////////////////////////////////////

    //free(aPriori);
    //free(minCondP);
    //free(hatMu);
    //free(hatSD);
    //free(minX);
    //free(maxX);
    //free(X);
    //free(Y);
    //free(ID);
    //free(N);

	return(ans);
}


/* Auxiliary routine for evaluating an objective function coded in R */
int increaseSampleXY(int maxFE, int dimension, double X[maxFE][dimension],
                 double Y[maxFE], int ID[maxFE],
                 double minX[][dimension], double maxX[][dimension],
                 int curLenY, int curPart, int sampleSz,
                 SEXP rFunction, SEXP env)
{
    int i,j,auxError;
    double x;
    double *aux;

    if((curLenY+sampleSz)>maxFE)
        sampleSz=maxFE-curLenY;

    SEXP arg,call,returned;

    for(i=0;i<sampleSz;i++)
    {
        arg=allocVector(REALSXP,dimension);
        PROTECT(arg);  aux=REAL(arg);
        for(j=0;j<dimension;j++)
        {
            x=runif(minX[curPart][j],maxX[curPart][j]);
            aux[j]=x; X[curLenY+i][j]=x;
        }
        ID[curLenY+i]=curPart;

        call=lang2(rFunction,arg);
        PROTECT(call);
        returned=R_tryEval(call,env,&auxError);
        Y[curLenY+i]=REAL(returned)[0];
        UNPROTECT(1);
        UNPROTECT(1);
    }
    return(sampleSz);
}
