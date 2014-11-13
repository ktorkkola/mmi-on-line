#include "mex.h"
#include "math.h"
#ifdef SUNOPT
	#include "sunmath.h"
#endif



/*=========================================================
  pfMex
  
  Computes the information potentials and information forces
  without storing the mutual distance matrix.

  The MATLAB call is 
    I                      = pfMex(A,sigma,J); or
    [I,dC]                 = pfMex(A,sigma,J); or
    [I,dC,maxdist,avedist] = pfMex(A,sigma,J);
 
  ver 1.00  2-10-00  
  ver 1.50  2-10-00   returns also the max and ave distances between two output datapoints
  ver 1.60  4-10-01   possibility to only compute I

=========================================================*/



int pf(            // Matrix sizes below are expressed as (rows x columns) 
  double *I,       // total information potential (mutual information) 
  double *dC,      // total information forces (od x Ntot)  
  int od,          // dimension of the output space 
  int Ntot,        // number of items of training data 
  int Nc,          // number of classes 
  int *J,          // indices to the 1st and the last item of each class (Nc x 2) 
  double *A,       // training data in the output space (od x Ntot), class after another 
  double sigma,    // kernel width 
  double *maxdist, // returns max distance between two points 
  double *avedist  // returns ave distance between two points 
)
{
  double *pdC;  // pointer to forces for current item 
  double iI;    // potential for current item 
  double *idC;  // intermediate forces for current sample 
  
  double largestdistance = 0.0; // distance between two farthest points in the output space 
  double sumdistance = 0.0;     // sum of distances between two points in the output space 
  int distcounter = 0;
  
  int cc;       // column class index 
  int ci;       // column index within class
  int rc;       // row class index 
  int ri;       // row index within class

  // constants
  double sigmasqr, zeropotential, pi;
  double *N;
  double const_Vcy2,  const_Vc2y2,  const_Vcy,  const_V;
  double const_dVcy2, const_dVc2y2, const_dVcy, const_dV;

  // temporaries
  int i;
  double *diff, diffsq, v;
  
  // allocate temporary memory & check
  N    = malloc(Nc*sizeof(double));
  idC  = malloc(od*sizeof(double));
  diff = malloc(od*sizeof(double));
  if (!(N && idC && diff)) return(0);
    
  // Pre-compute several constants
  const_Vc2y2 = 0.0;
  for (i=0; i<Nc; i++) {
        N[i] = J[i+Nc] - J[i] +1;
        const_Vc2y2 += N[i]*N[i];
  }
  const_Vcy2   = 1.0 / ((double)Ntot*(double)Ntot);
  const_Vc2y2 *= (const_Vcy2 * const_Vcy2);
  const_Vcy    = -2.0 * const_Vcy2 / Ntot;

  const_dVcy2  = 1.0 / ((sigma*Ntot)*(sigma*Ntot));
  const_dVc2y2 =  const_Vc2y2 / (sigma*sigma);
  const_dVcy   = -const_dVcy2 / Ntot;

  sigmasqr = 1.0 / (-2.0 * sigma*sigma);
  pi = acos( -1.0 );
  zeropotential = 1.0 / ( sigma * pow( 2.0*pi, 0.5*od) );
     
  // Initialize MI and information forces
  // Matlab is supposed to initialize its allocated variables to zeroes, 
  //  but this is here mainly for use outside Matlab environment.
  *I = 0.0;
  if (dC != NULL) {
     for (i=0; i<od*Ntot; i++) dC[i]=0.0;
  }
  
  // for each item in the training set compute the mutual information and its derivative
  for (cc=0; cc<Nc; cc++) {                     // each class
     for (ci=J[cc]; ci<=J[Nc+cc]; ci++) {       // each item in the class
        
        pdC = dC + od*ci;                       // pointer to the resulting force for item ci
                                                // compute anyway even if forces not requested
        // by computing the distance to all other items   
        for (rc=0; rc<Nc; rc++) {                       // each class
           iI = 0.0;
           for (i=0;i<od;i++) idC[i]=0.0;
           
           for (ri=J[rc]; ri<=J[Nc+rc]; ri++) {         // each item in the class
              diffsq = 0.0;
              for (i=0;i<od;i++) {
                     v = A[od*ri+i] - A[od*ci+i];       // distance between items ri and ci
                     diff[i] = v;
                     diffsq += v*v;
              }
#ifdef SUNOPT
              v = expf( sigmasqr*diffsq ) * zeropotential;// potential between items ri and ci
#else
              v = exp( sigmasqr*diffsq ) * zeropotential;// potential between items ri and ci
#endif
              iI += v;

              if (dC!=NULL)   // don't compute if forces not requested 
                 for (i=0;i<od;i++) idC[i] += v*diff[i];   // force between items ri and ci
              
              // keeping books about distances 
              if (diffsq>largestdistance) largestdistance=diffsq;
              sumdistance += diffsq;
              distcounter++;
           }


           // now we have potential and sum of forces between classes cc and rc
           // sum potential to total potential 
           const_V = const_Vc2y2 + N[cc]*const_Vcy;
           if (cc==rc) const_V += const_Vcy2;
           *I += const_V * iI;


           if (dC!=NULL) { // don't compute if forces not requested 
              // sum force to total force 
              const_dV = const_dVc2y2 + (N[rc]+N[cc])*const_dVcy;
              if (cc==rc) const_dV += const_dVcy2;
              for (i=0;i<od;i++) pdC[i] += const_dV * idC[i];
           } // end for dC!=NULL

           
        } // end for rc
        
     } // end for ci 
     
  } // end for cc
  
  *avedist = sqrt(sumdistance/distcounter);
  *maxdist = sqrt(largestdistance);
    
  free(N); free(idC); free(diff);
  return(1);
}





/*=========================================================
        the gateway function 
=========================================================*/
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *A, *I, *dC;
  double  sigma;
  int     od, Ntot, Nc;
  double  *Jreal; int *Jint;
  int i, ok;
  
  double maxdist, avedist;
  double *pmaxdist, *pavedist;
  
  //  check for proper number of arguments 
  if(nrhs!=3) mexErrMsgTxt("Three inputs required.");
  if((nlhs!=1) & (nlhs!=2) & (nlhs!=4)) mexErrMsgTxt("One, two or four outputs required.");
  
  //--------------------------- input args ---------------------------
  
  // get 1st input argument - A 
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetN(prhs[0])*mxGetM(prhs[0])<=1 ) {
    mexErrMsgTxt("Input A must be a double array.");
  }
  A = mxGetPr(prhs[0]);   // create a pointer to the input data matrix A 
  od   = mxGetM(prhs[0]); // get the dimensions of A 
  Ntot = mxGetN(prhs[0]);  
  
  
  // get 2nd input argument - sigma 
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetN(prhs[1])*mxGetM(prhs[1])!=1 ) {
    mexErrMsgTxt("Input sigma must be a scalar.");
  }
  sigma = mxGetScalar(prhs[1]);  // get the scalar input sigma 
    
    
  // get 3rd input argument - J 
  if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetN(prhs[2])*mxGetM(prhs[2])<=1 ) {
    mexErrMsgTxt("Input J must be a double array.");
  }
  Jreal = mxGetPr(prhs[2]);   // create a pointer to the input class index matrix J 
  Nc = mxGetM(prhs[2]);       // get the number of classes 
  
  // create an integer copy of the (double) index array 
  // change to zero-based indexing, too 
  Jint = mxCalloc(2*Nc, sizeof(int));
  for (i=0; i<2*Nc; i++) Jint[i] = Jreal[i]-1;
    
  //--------------------------- output args ---------------------------
    
  plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL); // create output data I as 1x1 matrix 
  I = mxGetPr(plhs[0]);         // create a C pointer to I 

  if (nlhs>1) {
    plhs[1] = mxCreateDoubleMatrix(od,Ntot, mxREAL); // create output data dC as od x Ntot matrix 
    dC = mxGetPr(plhs[1]);                // create a C pointer to dC 
  } else {
    dC = NULL; // dC (=forces, derivatives) is not requested 
  }
  
  // optional output arguments (these are calculated anyway, but not returned if not requested) 
  if (nlhs==4) {
  	plhs[2] = mxCreateDoubleMatrix(1,1, mxREAL); // create output data maxdist as 1x1 matrix 
  	pmaxdist = mxGetPr(plhs[2]);         // create a C pointer to maxdist 
  	plhs[3] = mxCreateDoubleMatrix(1,1, mxREAL); // create output data avedist as 1x1 matrix 
  	pavedist = mxGetPr(plhs[3]);         // create a C pointer to avedist 
  } else {
  	pmaxdist = &maxdist;
  	pavedist = &avedist;
  }
    
  //  call the C subroutine 
  ok = pf(I,dC, od,Ntot,Nc,Jint,A,sigma,pmaxdist,pavedist);

  // free whatever needs to be freed 
  mxFree(Jint);
  
  if (!ok) mexErrMsgTxt("Potentials and forces failed!");
    
}
