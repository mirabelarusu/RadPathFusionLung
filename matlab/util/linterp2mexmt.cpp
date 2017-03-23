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
#include "../general/maxNumCompThreads.h"

/*   undef needed for LCC compiler  */
#undef EXTERN_C
#ifdef _WIN32
	#include <windows.h>
	#include <process.h>
#else
	#include <pthread.h>
#endif

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

struct linterp_threadargs{
    double* D;
    double* X;
    double* Y;
    double* I;
    size_t ninterps;
    size_t nrows;
    size_t ncols;
    size_t narrays;
    double bgval;
    size_t ThreadID;
    size_t Nthreads;
};

/* GLOBALS */
size_t THREADLIMIT=2;

#ifdef _WIN32
unsigned __stdcall linterp(void *ThreadArgsV) {
#else
void *linterp(void *ThreadArgsV) {
#endif
    // unsigned long long int (__int64) for the array indexes:
    size_t n, i, nrc, zskip;
    // signed long long ints for the relative indexes:
    signedsize_t ind1, ind2, ind3, ind4;
    signedsize_t ndx, fy, fx;
    // doubles for the interpolated data, fractional coefficients, and potentially oob locations
    double y, x, dx, dy;
    double w1, w2, w3, w4;
    
    linterp_threadargs *ThreadArgs=(linterp_threadargs *) ThreadArgsV;
    
    double *D, *X, *Y, *I, bgval;
    size_t ninterps, nrows, ncols, narrays, ThreadOffset, Nthreads;
    
    D=ThreadArgs->D;
    X=ThreadArgs->X;
    Y=ThreadArgs->Y;
    I=ThreadArgs->I;
    ninterps=ThreadArgs->ninterps;
    nrows=ThreadArgs->nrows;
    ncols=ThreadArgs->ncols;
    narrays=ThreadArgs->narrays;
    bgval=ThreadArgs->bgval;
    ThreadOffset=ThreadArgs->ThreadID;
    Nthreads=ThreadArgs->Nthreads;
    
    nrc=nrows*ncols;
    
    for (n=ThreadOffset; n<ninterps; n+=Nthreads) {
        
        y=*(Y+n); // advance pointer and dereference instead of array indexing
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
            // fractional component of coords
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
    
    /* end thread */
    #ifdef _WIN32
	_endthreadex( 0 );
    return 0;
	#else
	pthread_exit(NULL);
	#endif
	
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
    size_t ninterps, nrows, ncols, npages, i;
    mwSize ndim, *newdims, Xndim;
    const mwSize *dims, *Xdims;
    
    /* Check for input errors */
    if (nlhs>1)
        mexErrMsgTxt("Wrong number of output parameters, usage:  DI = linterp2mexmt(D, X, Y, [bgval])");
    if (nrhs>4 || nrhs<3)
        mexErrMsgTxt("Wrong number of input parameters, usage:  DI = linterp2mexmt(D, X, Y, [bgval])");
    
    if (!mxIsDouble(IN_D) || !mxIsDouble(IN_X) || !mxIsDouble(IN_Y))
        mexErrMsgTxt("linterp2mexmt: Input arguments must be double.");
    
    if ((mxGetNumberOfDimensions(IN_X) != mxGetNumberOfDimensions(IN_Y)) ||
            (mxGetNumberOfElements(IN_X) != mxGetNumberOfElements(IN_Y))) mexErrMsgTxt("Input parameters X, Y must have the same size");
    
    if (nrhs>=4  && !mxIsEmpty(prhs[3])) {
        if (mxGetNumberOfElements(prhs[3])!=1) mexErrMsgTxt("linterp2mexmt: bgval must be a scalar.");
        bgval=mxGetScalar(prhs[3]);
    } else {
        bgval=mxGetNaN();
    }
    
    /* Get the sizes of each input argument */
    ndim = mxGetNumberOfDimensions(IN_D);
    dims = mxGetDimensions(IN_D);
    nrows = dims[0];
    ncols = dims[1];
    if ((nrows==1) || (ncols==1)) mexErrMsgTxt("linterp2mexmt: Data must not be 1D.");
    
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
    
    // Threading
    size_t *ThreadID;
    linterp_threadargs **ThreadArgs;
    linterp_threadargs *ThreadArgsi;
    
    /* Handles to the worker threads */
	#ifdef _WIN32
		HANDLE *ThreadList; 
    #else
		pthread_t *ThreadList;
	#endif
    
    size_t Nthreads;
    size_t maxthreads=(size_t) getNumCores();
    Nthreads=(maxthreads>THREADLIMIT ? THREADLIMIT : maxthreads); // hard thread limit
//     if (Nthreads>threadlim) Nthreads=threadlim; // input thread limit
    
    /* Reserve room for handles of threads in ThreadList */
    #ifdef _WIN32
		ThreadList = (HANDLE *) malloc(Nthreads*sizeof(HANDLE));
    #else
		ThreadList = (pthread_t *) malloc(Nthreads*sizeof(pthread_t));
	#endif
	
	ThreadID = (size_t *) malloc(Nthreads*sizeof(size_t));
	ThreadArgs = (linterp_threadargs **) malloc(Nthreads*sizeof(linterp_threadargs *));
    
    // Find kNNs for each point by traversing kd tree
    for (i=0; i<Nthreads; i++) {
        
        ThreadID[i]=i;
        
        ThreadArgsi = (linterp_threadargs *) malloc(sizeof(linterp_threadargs));
        ThreadArgsi->D=D;
        ThreadArgsi->X=X;
        ThreadArgsi->Y=Y;
        ThreadArgsi->I=I;
        ThreadArgsi->ninterps=ninterps;
        ThreadArgsi->nrows=nrows;
        ThreadArgsi->ncols=ncols;
        ThreadArgsi->narrays=npages;
        ThreadArgsi->bgval=bgval;
        ThreadArgsi->ThreadID=ThreadID[i];
        ThreadArgsi->Nthreads=Nthreads;
        
        /* Start thread  */
        ThreadArgs[i]=ThreadArgsi; // now we can overwrite ThreadArgsi for the next thread
        
        #ifdef _WIN32
            ThreadList[i] = (HANDLE)_beginthreadex( NULL, 0, &linterp, ThreadArgs[i] , 0, NULL );
        #else
            pthread_create((pthread_t*)&ThreadList[i], NULL, &linterp, ThreadArgs[i]);
        #endif
    }
    
    #ifdef _WIN32
            for (i=0; i<Nthreads; i++) WaitForSingleObject(ThreadList[i], INFINITE);
            for (i=0; i<Nthreads; i++) CloseHandle( ThreadList[i] );
    #else
            for (i=0; i<Nthreads; i++) pthread_join(ThreadList[i], NULL);
    #endif
    
	for (i=0; i<Nthreads; i++) free(ThreadArgs[i]);
    free(ThreadArgs);
    free(ThreadID);
    free(ThreadList);
}
