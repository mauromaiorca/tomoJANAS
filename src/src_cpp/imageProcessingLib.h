/*
 * File: imageProcessingLib.h
 * (C) 2022 Mauro Maiorca 
 */

#ifndef __IMAGEPROCESSING_LIB__
#define __IMAGEPROCESSING_LIB__

#include <cmath>
#include <iostream>
//#include <iomanip>      // std::setprecision
//#include <complex>
#include <cstdio>
#include <cstdlib>
#include "scores.h"


#include "mrcIO.h"

// ***********************************
//  
// 
template <typename T>
void blurMaskEdgeOnImage(T* I, T* maskI, T*imageOut, double sigma, int kernelLength , const unsigned long int nx, const unsigned long int ny){

//std::cerr<<"sigma="<< sigma << "   kernelLengthint="<< kernelLength<<"\n";

const int nxy = nx * ny;
const int length=kernelLength;
const int nxK =  (2*length+1);
const int nyK = (2*length+1);
const int kernelSize = nxK*nyK;
float * kernel = new float [ kernelSize ];
float * tmpI = new float [ nxy ];
//double angpix = 1;


double sum = 0;
double counter=0;

for (int jjT=-length; jjT<length+1; jjT++){
  for(int iiT=-length; iiT<length+1; iiT++){
          counter++;
        double distanceSquared =  pow(jjT,2.0)+pow(iiT,2.0);
        int X=iiT+length;
        int Y=jjT+length;
        kernel [  X + (Y * nxK) ] = 1 / (sigma * sigma * 2*3.1415927) 
              * exp (- (distanceSquared) / (2.0 * sigma * sigma));
        sum+=kernel [  X + (Y * nxK) ];
  }
}

//correct compensate for discretization (sum needs to be 1):
double normSum = (double((1.0-sum)))/double(kernelSize);
for (int ii=0; ii<kernelSize; ii++){
        kernel [ii] += normSum;
}


sum = 0;
for (int jjT=-length; jjT<length+1; jjT++){
  for(int iiT=-length; iiT<length+1; iiT++){
        //double distanceSquared =  pow(jjT,2.0)+pow(iiT,2.0);
        int X=iiT+length;
        int Y=jjT+length;
        sum+=kernel [ X +Y * nxK ];
 //       std::cerr << "kernel ["<< X +Y * nxK << "]=";
 //       std::cerr << kernel [ X +Y * nxK ] << " ";
  }
 // std::cerr<< " \n";
}
//std::cerr << " \n";
//std::cerr << "sum=" << sum <<" \n";



for (unsigned long int YY=length; YY<ny-length; YY++){
 for(unsigned long int XX=length; XX<nx-length; XX++){
   int numBlack=0;
   int numWhite=0;

   for (int jjT=-length; jjT<length; jjT++){
     for(int iiT=-length; iiT<length; iiT++){
        unsigned long int idx = (XX+iiT)+(YY+jjT)*nx;
        if ( maskI[idx] > 0.5f){
                numWhite++;  
        }else{
                numBlack++; 
        }
     }
   }


   //
   if (numBlack>0 && numWhite>0){
        double sum = 0;
        for (int jjT=-length; jjT<length; jjT++){
          for(int iiT=-length; iiT<length; iiT++){
                unsigned long int idx = (XX+iiT)+(YY+jjT)*nx;
                unsigned long int kernelIdx=iiT+length+(jjT+length)*nxK;
                sum+= kernel[kernelIdx]*I[idx];
          }
        }
        tmpI [ XX + YY * nx ] = sum;
   }else{
     tmpI [ XX + YY * nx ] = I [ XX + YY * nx ];
   }
 }
 

}

 for (unsigned long int ii=0;ii<(unsigned long int)nxy;ii++){
   imageOut[ii]= tmpI[ii];
 }

delete [] kernel;
delete [] tmpI;
}


// ************************
//
//      blurMaskedImageEdge2D
//
template <typename T, typename U>
void blurMaskedImageEdge2D(T* imageIn, U* IMask, T* imageOut, unsigned long int nx, unsigned long int ny,  std::vector<std::vector<double> > listImageEdgesSigmaWidthList, double maskThreshold = 0.9){
  unsigned long int nxy = nx * ny; 
  if (listImageEdgesSigmaWidthList.size()==0){
    //do nothing
    for (unsigned long int ij=0; ij<nxy; ij++){
      imageOut[ij]=imageIn[ij];
    }
  }else{
    for (unsigned int hh=0; hh<listImageEdgesSigmaWidthList.size(); hh++){
      if (hh==0 ){
        if ( listImageEdgesSigmaWidthList[0][0]==0 ){
            //do nothing
            for (unsigned long int ij=0; ij<nxy; ij++){
              imageOut[ij]=imageIn[ij];
            }
          }else{
            blurMaskEdgeOnImage(imageIn, IMask, imageOut, listImageEdgesSigmaWidthList[0][0], listImageEdgesSigmaWidthList[0][1], nx, ny, maskThreshold);
          }
        }else{
            blurMaskEdgeOnImage(imageOut, IMask, imageOut, listImageEdgesSigmaWidthList[hh][0], listImageEdgesSigmaWidthList[hh][1], nx, ny, maskThreshold);
        }
    }
  }
}




template<typename T>
std::vector<double> rotAverage ( T * image, //!< input image
    const unsigned long int nx, //!< size of the image
    const unsigned long int ny, //!< size of the image
    const unsigned long int nz //!< size of the image
){
  std::vector<double> returnVector;
  unsigned long int nxyz=nx*ny*nz;
  //T * idxI=new T [nxyz];

  //  get the size of the buffer
  unsigned long int sizeBuf=0;
  double nx0=nx/2.0;
  double ny0=ny/2.0;
  double nz0=ny/2.0;
  for (unsigned long int kk=0; kk<nz; kk++){
    double distanceZ=pow((double)kk-nz0,2.0);
    for (unsigned long int jj=0; jj<ny; jj++){
      double distanceY=pow((double)jj-ny0,2.0);
      for (unsigned long int ii=0; ii<nx; ii++){
        double distanceX=pow((double)ii-nz0,2.0);
        unsigned long int distance= floor(pow( distanceX+distanceY+distanceZ,0.5 ));
        if ( distance > sizeBuf ) sizeBuf = distance;

      }
    }
  }
  sizeBuf++;
  double * bufferI= new double [sizeBuf];
  double * bufferCounterI= new double [sizeBuf];
  for (unsigned long i=0; i<sizeBuf; i++){
    bufferI[i]=bufferCounterI[i]=0;
  }


  for (unsigned long int kk=0; kk<nz; kk++){
    double distanceZ=pow((double)kk-nz0,2.0);
    for (unsigned long int jj=0; jj<ny; jj++){
      double distanceY=pow((double)jj-ny0,2.0);
      for (unsigned long int ii=0; ii<nx; ii++){
        double distanceX=pow((double)ii-nz0,2.0);
        unsigned long int distance= floor(pow( distanceX+distanceY+distanceZ,0.5 ));
        unsigned long int idx =ii+jj*nx+kk*(nx*ny);
        //idxI[idx] = distance;
        if (distance<sizeBuf ){
          bufferI[distance]+=distance;
          bufferCounterI[distance]++;
        }
      }
    }
  }
  for (long int i=0; i<sizeBuf; i++){
    if ( bufferCounterI[i] > 0 ){
      returnVector.push_back(bufferI[i]/bufferCounterI[i]);
    }
  }

  //writeMrcImage("idx.mrc", idxI,  nx,  ny,  nz,  1);
  //delete [] idxI;
  delete [] bufferI;
  delete [] bufferCounterI;
  return returnVector;
}








template<typename T>
void amplitudeNormalization ( T * AmplitudesImage, //!< input image
    T * amplitudesToReplace, //!< input image
    const unsigned long int nx, //!< size of the image
    const unsigned long int ny, //!< size of the image
    const unsigned long int nz //!< size of the image
){

  std::vector<double> meanAmplitudesIn=rotAverage(AmplitudesImage, nx, ny, nz);
  std::vector<double> meanAmplitudesOut=rotAverage(amplitudesToReplace, nx, ny, nz);
  std::vector<double> normalizingVector;
  const long double thresholdLowAmplitude = 0.0000000001;
  for (unsigned long int ii=0; ii<meanAmplitudesIn.size(); ii++){
    double normValue=meanAmplitudesIn[ii];
    if ( meanAmplitudesIn[ii] > thresholdLowAmplitude ){
      normValue = meanAmplitudesOut[ii]/meanAmplitudesIn[ii];
    }
    normalizingVector.push_back(normValue);
  }

  unsigned long int nxyz=nx*ny*nz;
  unsigned long int sizeBuf=meanAmplitudesIn.size();
  double nx0=nx/2.0;
  double ny0=ny/2.0;
  double nz0=ny/2.0;


  for (unsigned long int kk=0; kk<nz; kk++){
    double distanceZ=pow((double)kk-nz0,2.0);
    for (unsigned long int jj=0; jj<ny; jj++){
      double distanceY=pow((double)jj-ny0,2.0);
      for (unsigned long int ii=0; ii<nx; ii++){
        double distanceX=pow((double)ii-nz0,2.0);
        unsigned long int distance= floor(pow( distanceX+distanceY+distanceZ,0.5 ));
        unsigned long int idx =ii+jj*nx+kk*(nx*ny);
        //idxI[idx] = distance;
        
        if ( distance < sizeBuf ){
          AmplitudesImage[idx] = AmplitudesImage[idx] * normalizingVector[distance];
        }
      }
    }
  }
  //writeMrcImage("idx.mrc", idxI,  nx,  ny,  nz,  1);
  //delete [] idxI;
}



#endif
