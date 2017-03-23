#include "mex.h"
#include <math.h>
#include <algorithm>
//#include <sstream>
//#include <string>
/* MI
 *Jon Chappelow
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *mi, *nmi, *ecc, eoh;
    double *I1, *I2, *N;
    double **H, *H1, *H2, S1, S2, S12, I1min, I1max, I2min, I2max, pfraction;
    //mwSize hndims=2; mwSize hdims[2];
    size_t i, j, npixI1, npixI2;
    //std::ostringstream sstream;
    //std::string *errstring;
    //const char *errstring;
    
    if (nlhs>7)
        mexErrMsgTxt("Wrong number of output parameters, usage:\n\t[mi,nmi,ecc,s12,s1,s2,eoh] = mimex(Image1, Image2, N)");
    if (nrhs!=3)
        mexErrMsgTxt("Wrong number of input parameters, usage:\n\t[mi,nmi,ecc,s12,s1,s2,eoh] = mimex(Image1, Image2, N)");
    
    if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]))
        mexErrMsgTxt("Inputs not doubles.");
    
    npixI1=mxGetNumberOfElements(prhs[0]);
    npixI2=mxGetNumberOfElements(prhs[1]);
    if (npixI1 != npixI2) mexErrMsgTxt("Data sets must have the same number of elements.");
    
    I1=mxGetPr(prhs[0]);
    I2=mxGetPr(prhs[1]);
    N=mxGetPr(prhs[2]);
    
    I1min=*std::min_element(I1,I1+npixI1);
    I1max=*std::max_element(I1,I1+npixI1);
    I2min=*std::min_element(I2,I2+npixI2);
    I2max=*std::max_element(I2,I2+npixI2);
    
    /*
    mexPrintf("Highest value in I1: %g\n",I1max);
    mexPrintf("Highest value in I2: %g\n",I2max);
    mexPrintf("Lowest value in I1: %g\n",I1min);
    mexPrintf("Lowest value in I2: %g\n",I2min);
    */
    
    if (I1min<0 || I2min<0 || I1max>*N-1 || I2max>*N-1) {
        //sstream << "Data out of range for N of " << *N << std::endl;
        //const char *errstring = sstream.str().c_str();
        mexPrintf("ERROR: Data out of range for N of %g.\n",*N);
        mexErrMsgTxt("Rescale your data.\n");
    }
    
    //hdims[0]=(mwSize) *N; // dereference to get double and cast to mwSize
    //hdims[1]=(mwSize) *N;
    //plhs[0]=mxCreateNumericArray(hndims, hdims, mxDOUBLE_CLASS, mxREAL);
    //H=mxGetPr(plhs[0]);
    
    plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    mi=mxGetPr(plhs[0]);
    if (nlhs>1) {
        plhs[1]=mxCreateDoubleMatrix(1,1,mxREAL);
        nmi=mxGetPr(plhs[1]);
        if (nlhs>2) {
            plhs[2]=mxCreateDoubleMatrix(1,1,mxREAL);
            ecc=mxGetPr(plhs[2]);
        }
    }
    
    H = (double **) mxMalloc((size_t) *N*sizeof(double *));
    for (i=0; i<*N; i++) H[i]=(double *) mxCalloc((size_t) *N, sizeof(double));
    //H = (double **) malloc((size_t) *N*sizeof(double *));
    //for (i=0; i<*N; i++) H[i]=(double *) calloc((size_t) *N, sizeof(double));
    
    pfraction=1./npixI1;
    
    for (i=0; i<npixI1; i++) {
        H[(size_t) I1[i]][(size_t) I2[i]]+=pfraction;
        //*(*(H+I1[i])+I2[i])
        //H[round(I1[i])][round(I2[i])]+=pfraction;
        //H[I1[i]][I2[i]]++;
        //H[*(I1+i)][*(I2+i)]++;
        
        //H1[I1[i]]+=pfraction;
        //H2[I2[i]]+=pfraction;
    }
    
    H1 = (double *) mxCalloc((size_t) *N, sizeof(double));
    H2 = (double *) mxCalloc((size_t) *N, sizeof(double));
    S1=S2=S12=eoh=0;
    for (i=0; i<*N; i++) {
        for (j=0; j<*N; j++) {
            H1[i]+=H[i][j];
            //H2[j]+=H[i][j];
            H2[i]+=H[j][i];
            S12-=(H[i][j]>0?H[i][j]*log(H[i][j]):0);
            eoh+=H[i][j]*H[i][j];
        }
        S1-=(H1[i]>0?H1[i]*log(H1[i]):0);
        S2-=(H2[i]>0?H2[i]*log(H2[i]):0);
    }

    /*
    S1=S2=S12=0;
    for (i=0; i<*N; i++) {
        S1-=(H1[i]>0?H1[i]*log(H1[i]):0);
        S2-=(H2[i]>0?H2[i]*log(H2[i]):0);
        for (j=0; j<*N; j++)
            S12-=(H[i][j]>0?H[i][j]*log(H[i][j]):0);
    }
     */
    //for (i=0; i<*N; i++)
    //    for (j=0; j<*N; j++) S12-=(H[i][j]>0?H[i][j]*log(H[i][j]):0);
    
    S1/=log(2.);
    S2/=log(2.);
    S12/=log(2.);
    
    *mi=S1+S2-S12;
    if (nlhs>1) {
        *nmi=(S1+S2)/S12;
        if (nlhs>2) {
            *ecc=2-2*S12/(S1+S2);
            if (nlhs>3) {
                plhs[3]=mxCreateDoubleMatrix(1,1,mxREAL);
                *mxGetPr(plhs[3])=S12;
                if (nlhs>4) {
                    plhs[4]=mxCreateDoubleMatrix(1,1,mxREAL);
                    *mxGetPr(plhs[4])=S1;
                    if (nlhs>5) {
                        plhs[5]=mxCreateDoubleMatrix(1,1,mxREAL);
                        *mxGetPr(plhs[5])=S2;
                        if (nlhs>6) {
                            plhs[6]=mxCreateDoubleMatrix(1,1,mxREAL);
                            *mxGetPr(plhs[6])=eoh;
                        }
                    }
                }
            }
        }
    }
    
    //free(H);
    
}
