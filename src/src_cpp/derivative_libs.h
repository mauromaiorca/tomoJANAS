/*
 * File: derivative_libs.h
 * (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
 */

#ifndef ___DERIVATIVE_LIBS___H__
#define ___DERIVATIVE_LIBS___H__

#include "recline.h"
#include "recline.c"


// **********************************
//  largest eigen value of symmetric matrix in the form
//    | a  b |
//    | b  c |
//
double _largestEigenvalue(double a, double b, double c){
    double delta = pow(a-c,2.0)+4.0*b*b;
    if (delta<0.0) return 0;
    delta=pow(delta,0.5);
//    return ((a+c+delta)/2.0);
    return (a+c+delta);
}

//reverse array
template<typename T>
void flipImage( T* I, 
             T* outI, 
            unsigned long int nx, unsigned long int ny, unsigned long int nz, 
            int direction){
  unsigned long int nxy=nx*ny;
  int bufferSize=nx;
  if (direction==1) bufferSize=ny;
  if (direction==2) bufferSize=nz;
  T * bufferLine = new T [ bufferSize ];
  
	if (direction==0){			 //along X
        for (unsigned long int kk=0; kk<nz; kk++){
            unsigned long int kk0=nxy*kk;
            for (unsigned long int jj=0; jj<ny; jj++){
                unsigned long int jjnxK0 = jj*nx+kk0;
                for (unsigned long int ii=0; ii<nx; ii++){
                    bufferLine[nx-ii-1]=I[ii+jjnxK0];
                }
                for (unsigned long int ii=0; ii<nx; ii++){
                    outI[ii+jjnxK0]=bufferLine[ii];
                }
            }
        }
    }
    else if (direction==1) {			 //along Y
        for (unsigned long int kk=0; kk<nz; kk++){
            unsigned long int kknxy=kk*nxy;
            for (unsigned long int ii=0; ii<nx; ii++){
                //unsigned long int iikk = ii+kknxy;
                for (unsigned long int jj=0; jj<ny; jj++){
                    bufferLine[ny-jj-1]=I[ii+kknxy+jj*nx];
                }
                for (unsigned long int jj=0; jj<ny; jj++){
                    outI[ii+kknxy+jj*nx]=bufferLine[jj];
                }
            }
        }
    }
    else if (direction==2) {			 //along Z
        for (unsigned long int jj=0; jj<ny; jj++){                
            unsigned long int jjnx = jj*nx;
            for (unsigned long int ii=0; ii<nx; ii++){
                for (unsigned long int kk=0; kk<nz; kk++){
                    bufferLine[nz-kk-1]=I[ii+jjnx+kk*nxy];
                }
                for (unsigned long int kk=0; kk<nz; kk++){
                    outI[ii+jjnx+kk*nxy]=bufferLine[kk];
                }
            }
        }
    }

  delete [] bufferLine;

}


/**
 * gaussRecursiveDerivatives1D is an accessory function for computing 1D zero order and first order gaussian derivatives of a 2D image.
 *
 */
template<typename T>
void gaussRecursiveDerivatives1D (
double    sigma,				 //!< standard deviation of Gaussian
unsigned int     nx,			 //!< image dimension in x direction
unsigned int     ny,			 //!< image dimension in y direction
unsigned int     nz,			 //!< image dimension in z direction
double    angpix,					 //!< pixel size
unsigned long int padding,					 //!< padding
unsigned int direction,			 //!< 0 for X, 1 for Y
unsigned int DerivativeOrder,	 //!< 0 smoothing, 1 first derivative, 2 second derivative
T    *I,					 //!< input: original image ;
T    *I_o = NULL			 //!<output: smoothed (if null, input overwritten)
)								 // input: original image ;  output: smoothed
{
    unsigned long int nxy = nx * ny;
    //unsigned long int nxyz = nx * ny * nz;

	//recursiveFilterType recFilter = ALPHA_DERICHE;
	recursiveFilterType recFilter = GAUSSIAN_DERICHE;
	derivativeOrder derivOrder = NODERIVATIVE;
	if (DerivativeOrder == 0)
	{
		derivOrder = SMOOTHING;
	}
	else if (DerivativeOrder == 1)
	{
		derivOrder = DERIVATIVE_1;
	}
	else if (DerivativeOrder == 2)
	{
		derivOrder = DERIVATIVE_2;
	}
	else if (DerivativeOrder == 3)
	{
		derivOrder = DERIVATIVE_3;
	}
    
	if ( I_o == NULL ){
		I_o = I;
	}

	//put everything in a line
	double * bufferIn = NULL;
	double * bufferOut = NULL;
	double * bufferTmp0 = NULL;
	double * bufferTmp1 = NULL;
	//float filterCoefs;

	RFcoefficientType *rfc = InitRecursiveCoefficients( sigma*angpix, recFilter, derivOrder );
  if (!rfc)
    exit(1);

	//fill the buffer
	//fill the buffer
	if (direction==0)			 //along X
	{

        bufferIn = new double [nx+2*padding];
        bufferOut = new double [nx+2*padding];
        bufferTmp0 = new double [nx+2*padding];
        bufferTmp1 = new double [nx+2*padding];

		if (bufferIn == NULL || bufferOut == NULL || bufferTmp0 == NULL || bufferTmp1 == NULL)
			exit(1);

        for (unsigned long int kk=0; kk<nz; kk++){
        unsigned long int kk0=nxy*kk;
        for (unsigned long int jj=0; jj<ny; jj++){
            unsigned long int jjnx = jj*nx+kk0;
            unsigned long int jjnxLast = nx-1+jj*nx+kk0;
            for (unsigned long int ii=0; ii<padding; ii++){                
				bufferIn[ii]=(double) I[jjnx];
				bufferOut[ii]=(double) I[jjnx];
				bufferTmp0[ii]=(double) I[jjnx];
				bufferTmp1[ii]=(double) I[jjnx];
            }
            for (unsigned long int ii=nx+padding; ii<nx+2*padding; ii++){
				bufferIn[ii]=(double) I[jjnxLast];
				bufferOut[ii]=(double) I[jjnxLast];
				bufferTmp0[ii]=(double) I[jjnxLast];
				bufferTmp1[ii]=(double) I[jjnxLast];
            }
            for (unsigned long int ii=0; ii<nx; ii++){
				bufferIn[ii+padding]=(double) I[ii+jjnx];
				bufferOut[ii+padding]=(double) I[ii+jjnx];
				bufferTmp0[ii+padding]=(double) I[ii+jjnx];
				bufferTmp1[ii+padding]=(double) I[ii+jjnx];
            }
            RecursiveFilter1D( rfc, bufferIn, bufferOut, bufferTmp0, bufferTmp1, nx + 2*padding);
			for (unsigned long int ii=0; ii<nx; ii++) {
				I_o[ii+jjnx]=(float)bufferOut[ii+padding];
			}            
        }
        }
        delete [] bufferIn;
        delete [] bufferOut;
        delete [] bufferTmp0;
        delete [] bufferTmp1;


	}	
    else if (direction==1) {

        bufferIn = new double [ny+2*padding];
        bufferOut = new double [ny+2*padding];
        bufferTmp0 = new double [ny+2*padding];
        bufferTmp1 = new double [ny+2*padding];

		if (bufferIn == NULL || bufferOut == NULL || bufferTmp0 == NULL || bufferTmp1 == NULL)
			exit(1);

        for (unsigned long int kk=0; kk<nz; kk++){
            unsigned long int kknxy=kk*nxy;
            //unsigned long int lastJJ=(ny-1)*ny+kknxy;
            for (unsigned long int ii=0; ii<nx; ii++){
                unsigned long int iikk = ii+kknxy;
                for (unsigned long int jj=0; jj<padding; jj++){                
                    bufferIn[jj]=(double) I[iikk];
                    bufferOut[jj]=(double) I[iikk];
                    bufferTmp0[jj]=(double) I[iikk];
                    bufferTmp1[jj]=(double) I[iikk];
                }
                for (unsigned long int jj=ny+padding; jj<ny+2*padding; jj++){
                    unsigned long int iiLastJJ = ii+(ny-1)*nx+kknxy;
                    bufferIn[jj]=(double) I[iiLastJJ];
                    bufferOut[jj]=(double) I[iiLastJJ];
                    bufferTmp0[jj]=(double) I[iiLastJJ];
                    bufferTmp1[jj]=(double) I[iiLastJJ];
                }
                for (unsigned long int jj=0; jj<ny; jj++){
                    unsigned long int iijjnx = ii+jj*nx+kknxy;
                    bufferIn[jj+padding]=(double) I[iijjnx];
                    bufferOut[jj+padding]=(double) I[iijjnx];
                    bufferTmp0[jj+padding]=(double) I[iijjnx];
                    bufferTmp1[jj+padding]=(double) I[iijjnx];
                }
                RecursiveFilter1D( rfc, bufferIn, bufferOut, bufferTmp0, bufferTmp1, ny + 2*padding);
                for (unsigned long int jj=0; jj<ny; jj++) {
                    I_o[ii+jj*nx+kknxy]=(float)bufferOut[jj+padding];
                }            
            }
        }
        delete [] bufferIn;
        delete [] bufferOut;
        delete [] bufferTmp0;
        delete [] bufferTmp1;


	}	
    else if (direction==2) {

        bufferIn = new double [nz+2*padding];
        bufferOut = new double [nz+2*padding];
        bufferTmp0 = new double [nz+2*padding];
        bufferTmp1 = new double [nz+2*padding];

		if (bufferIn == NULL || bufferOut == NULL || bufferTmp0 == NULL || bufferTmp1 == NULL)
			exit(1);


        for (unsigned long int jj=0; jj<ny; jj++){                
            unsigned long int jjnx = jj*nx;
            for (unsigned long int ii=0; ii<nx; ii++){
                unsigned long int iijjnx = ii+jjnx;
                unsigned long int iijjnxLast = iijjnx+(nz-1)*nxy;
                for (unsigned long int kk=0; kk<nz; kk++){
                    bufferIn[kk]=(double) I[iijjnx];
                    bufferOut[kk]=(double) I[iijjnx];
                    bufferTmp0[kk]=(double) I[iijjnx];
                    bufferTmp1[kk]=(double) I[iijjnx];
                }
                for (unsigned long int kk=nz+padding; kk<nz+2*padding; kk++){
                    bufferIn[kk]=(double) I[iijjnxLast];
                    bufferOut[kk]=(double) I[iijjnxLast];
                    bufferTmp0[kk]=(double) I[iijjnxLast];
                    bufferTmp1[kk]=(double) I[iijjnxLast];
                }
                for (unsigned long int kk=0; kk<nz; kk++){
                    unsigned long int iijjnxkk = iijjnx+kk*nxy;
                    bufferIn[kk+padding]=(double) I[iijjnxkk];
                    bufferOut[kk+padding]=(double) I[iijjnxkk];
                    bufferTmp0[kk+padding]=(double) I[iijjnxkk];
                    bufferTmp1[kk+padding]=(double) I[iijjnxkk];
                }
                RecursiveFilter1D( rfc, bufferIn, bufferOut, bufferTmp0, bufferTmp1, nz + 2*padding);
                for (unsigned long int kk=0; kk<nz; kk++) {
                    I_o[iijjnx+kk*nxy]=(float)bufferOut[kk+padding];
                }            
            }
        }
        delete [] bufferIn;
        delete [] bufferOut;
        delete [] bufferTmp0;
        delete [] bufferTmp1;


	}	

  free(rfc);

}


template<typename T>
void gaussRecursiveDerivatives1D_reversed (
double    sigma,				 //!< standard deviation of Gaussian
unsigned int     nx,			 //!< image dimension in x direction
unsigned int     ny,			 //!< image dimension in y direction
unsigned int     nz,			 //!< image dimension in z direction
double    angpix,					 //!< pixel size
int padding,					 //!< padding
unsigned int direction,			 //!< 0 for X, 1 for Y
unsigned int DerivativeOrder,	 //!< 0 smoothing, 1 first derivative, 2 second derivative
T    *I,					 //!< input: original image ;
T    *I_o = NULL			 //!<output: smoothed (if null, input overwritten)
)	
{

    if ( I_o == NULL ){
		I_o = I;
	}

    flipImage( I_o, I_o,  nx,  ny,  nz, direction);
    gaussRecursiveDerivatives1D (sigma, nx, ny, nz, angpix, padding, direction, DerivativeOrder, I_o);	
    flipImage( I_o, I_o,  nx,  ny,  nz, direction);
}


template<typename T>
void gaussRecursiveDerivatives2D (
double    sigma,                 //!< standard deviation of Gaussian
unsigned int     nx,             //!< image dimension in x direction
unsigned int     ny,             //!< image dimension in y direction
unsigned int     nz,             //!< image dimension in z direction
double    angpix,                     //!< pixel size
int padding,                     //!< padding
unsigned int DerivativeOrder,     //!< 0 smoothing, 1 first derivative, 2 second derivative
T    *I,                     //!< input: original image ;
T    *I_o = NULL             //!<output: smoothed (if null, input overwritten)
)                                 // input: original image ;  output: smoothed
{
	if ( I_o == NULL ){
		I_o = I;
	}
    gaussRecursiveDerivatives1D (sigma,nx,ny,nz,angpix,padding,0,DerivativeOrder,I, I_o);
    gaussRecursiveDerivatives1D (sigma,nx,ny,nz,angpix,padding,1,DerivativeOrder,I_o);
    
}




template<typename T>
void computeEigenImage (T* I1, unsigned long int nx, unsigned long int ny, int depth, double sigma){
    //const double _MIN_VAL=-1.0;//0.0f;//0.000001;
    const unsigned long int nxy = nx * ny;
    //double CC1=0;//crossCorrelationDistance(I1, IReproj, nxy, IMaskReproj);

    double angpix = 1;
    int padding = 3;

    float * Ixx = new float[nxy];
    float * Iyy = new float[nxy];
    float * Ixy = new float[nxy];

    if (depth>1){
      gaussRecursiveDerivatives1D (sigma,  nx,  ny,  1, angpix, padding, 0, depth, I1, Ixx);
      gaussRecursiveDerivatives1D (sigma,  nx,  ny,  1, angpix, padding, 1, depth, I1, Iyy);
      gaussRecursiveDerivatives1D (sigma,  nx,  ny,  1, angpix, padding, 0, depth-1, I1, Ixy);
      gaussRecursiveDerivatives1D (sigma,  nx,  ny,  1, angpix, padding, 2, depth-1, Ixy);

      for (unsigned long int ii=0; ii<nxy; ii++){
          double tmpEigen_I1=_largestEigenvalue(Ixx[ii],Ixy[ii],Iyy[ii]);
          I1[ii]=tmpEigen_I1;
      }
    }
    delete [] Ixx;
    delete [] Iyy;
    delete [] Ixy;
}


#endif

