/* 2D linear interpolation
 * based on code by Andriy Myronenko, myron@csee.ogi.edu
 *
 * DI=linterp2mex(D,XI,YI,[BG]) interpolates the 2D data at XI,YI.
 * D is defined on a regular grid 1:size(D,1), 1:size(D,2).
 * Points outside the boundary are set to BG (optional), or NaN if not specified.
 * Using MATLAB, this is: DI = interp2(D,XI,YI,'linear',BG)
 *
 * D can be a stack of many 3D volumes (4D).
 */

#include <math.h>
#include "mex.h"

// Check windows arch
#if _WIN32 || _WIN64
#if _WIN64
#define ENVIRONMENT64
#pragma message("64-bit Windows")
#else
#define ENVIRONMENT32
#pragma message("32-bit Windows")
#endif
#endif

// Check GCC arch
#if __GNUC__
#if __x86_64__ || __ppc64__
#define ENVIRONMENT64
#pragma message("64-bit GCC")
#else
#define ENVIRONMENT32
#pragma message("32-bit GCC")
#endif
#endif

#ifdef ENVIRONMENT64
typedef signed long long int signedsize_t;
#else
typedef signed long int signedsize_t;
#endif

void linterp(double* D, double* X, double* Y, double* I,
        size_t ninterps, size_t nrows, size_t ncols, size_t narrays, double bgval) {
    
    // unsigned long long int (__int64) for the array indexes:
    size_t n, i, nrc, zskip;
    // signed long long ints for the relative indexes:
    signedsize_t ind1, ind2, ind3, ind4;
    signedsize_t ndx, fy, fx;
    // doubles for the interpolated data, fractional coefficients, and potentially oob locations
    double y, x, dx, dy;
    double w1, w2, w3, w4;
    
    nrc=nrows*ncols;
    
    for (n=0; n<ninterps; n++) {
        
        y=*(Y+n);
        x=*(X+n);
        fx=(size_t) floor(x);
        fy=(size_t) floor(y);
        
        if (fx<1 || x>ncols || fy<1 || y>nrows){
            // Set bgval when outside volume
            for (i=0; i<narrays; i++) I[n+i*ninterps]=bgval;
        } else {
            // linear index into I (lower left corner of interpolation neighborhood)
            ndx = fy + (fx-1)*nrows;
           
            if (x==ncols) {
                //x=x+1;
                ndx=ndx-nrows; // move left
            }
            if (y==nrows){
                //y=y+1;
                ndx=ndx-1; // move up
            }
            // decimal value component of coords
            dx=x-fx;
            dy=y-fy;
            
            ind1=ndx-1;     // upper left
            ind2=ndx;       // lower left
            ind4=ndx+nrows; // lower right
            ind3=ind4-1;    // upper right
            
            w4=dy*dx;       // lower right
            w1=1+w4-dy-dx;  // upper left
            w2=dy-w4;       // lower left
            w3=dx-w4;       // upper right
            
            for (i=0; i<narrays; i++){
                zskip=i*nrc;
                I[n+i*ninterps]=D[ind1+zskip]*w1 + D[ind2+zskip]*w2 + 
                        D[ind3+zskip]*w3 + D[ind4+zskip]*w4;
            }
        }
    }
}

// ------- Main function definitions ----------
/* Input arguments */
#define IN_D		prhs[0]
#define IN_X		prhs[1]
#define IN_Y		prhs[2]

/* Output arguments */
#define OUT_I		plhs[0]

/* Gateway routine */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    double *D, *X, *Y, *I, bgval;
    size_t ninterps, nrows, ncols, npages;
    mwSize ndim, *newdims, Xndim;
    const mwSize *dims, *Xdims;
    
    /* Check for input errors */
    if (nlhs>1)
        mexErrMsgTxt("Wrong number of output parameters, usage:  DI = linterp2mex(D, X, Y, [bgval])");
    if (nrhs>4 || nrhs<3)
        mexErrMsgTxt("Wrong number of input parameters, usage:  DI = linterp2mex(D, X, Y, [bgval])");
    
    if (!mxIsDouble(IN_D) || !mxIsDouble(IN_X) || !mxIsDouble(IN_Y))
        mexErrMsgTxt("linterp2mex: Input arguments must be double.");
    
    if ((mxGetNumberOfDimensions(IN_X) != mxGetNumberOfDimensions(IN_Y)) ||
            (mxGetNumberOfElements(IN_X) != mxGetNumberOfElements(IN_Y))) mexErrMsgTxt("Input parameters X, Y must have the same size");
    
    if (nrhs==4) {
        if (mxGetNumberOfElements(prhs[3])!=1) mexErrMsgTxt("linterp2mex: bgval must be a scalar.");
        bgval=mxGetScalar(prhs[3]);
    } else {
        bgval=mxGetNaN();
    }
    
    /* Get the sizes of each input argument */
    ndim = mxGetNumberOfDimensions(IN_D);
    dims = mxGetDimensions(IN_D);
    nrows = dims[0];
    ncols = dims[1];
    if ((nrows==1) || (ncols==1)) mexErrMsgTxt("linterp2mex: Data must not be 1D.");
    
    Xndim = mxGetNumberOfDimensions(IN_X);
    Xdims = mxGetDimensions(IN_X);
    ninterps=mxGetNumberOfElements(IN_X);
    
    newdims=(mwSize*) mxMalloc(ndim*sizeof(mwSize));
    newdims[0]=Xdims[0];
    newdims[1]=Xdims[1];
    if (ndim>2) newdims[2] = npages=dims[2];
    else npages=1;
    
    // Input arguments pointers
    D = mxGetPr(IN_D);
    X = mxGetPr(IN_X);
    Y = mxGetPr(IN_Y);
    
    // Output arguments
    OUT_I = mxCreateNumericArray(ndim, newdims, mxDOUBLE_CLASS, mxREAL);
    I = mxGetPr(OUT_I);
    
    /* Do the actual computations in a subroutine */
    linterp(D, X, Y, I, ninterps, nrows, ncols, npages, bgval);
}
