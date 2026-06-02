/*
 * File: scores.h
 * (C) 2022 Mauro Maiorca
 */

#ifndef __SCORING_FUNCTIONS_LIB___
#define __SCORING_FUNCTIONS_LIB___

#include <cmath>
#include <iostream>
#include <vector>
#include "derivative_libs.h"

//Normalize
std::vector<double> normalize(float * image1, //!< input image
        float * image2, //!< input image
        float * mask,
        const unsigned long int nxyz,  //!< size (in pixel) of the image
        double maskThreshold = 0.1
){

  std::vector<double> result;

  float max1=0;
  float min1=0;
  float max2=0;
  float min2=0;
  double mean1=0;
  double mean2=0;
  long int counter = 0;
  for (unsigned long int ijk = 0; ijk< nxyz; ijk++){
  if (mask[ijk]>maskThreshold){
    if (counter==0){
      max1=image1[ijk];
      max2=image2[ijk];
      min1=image1[ijk];
      min2=image2[ijk];
    }else{
      if (max1<image1[ijk]) max1=image1[ijk];
      if (max2<image2[ijk]) max2=image2[ijk];
      if (min1>image1[ijk]) min1=image1[ijk];
      if (min2>image2[ijk]) min2=image2[ijk];
    }
    mean1+=image1[ijk];
    mean2+=image2[ijk];
    counter++;
   }
  }
  if (counter>0){
    mean1/=counter;
    mean2/=counter;
  }
  double diffI1 = (double)(max1-min1);
  double diffI2 = (double)(max2-min2);
  
  result.push_back(mean1);
  result.push_back(diffI1);
  result.push_back(mean2);
  result.push_back(diffI2);
  return result;

}


// *************************************
/** build the grayscale histogram from an image */
// buildHistogram
// **************************************
template<typename T>
void buildHistogram(const T* image, //!< input image
        double * histogram, //!< output histogram
        const unsigned long int bins, //!< number of bins for the histogram
        const unsigned long int nxyz,  //!< size (in pixel) of the image
        const T * mask = NULL,
        double minMaskValue = 0.0
){
  if (mask==NULL){ //NO MASK
	double min=image[0], max=image[0];
	for (unsigned long int ii = 1; ii< nxyz; ii++){
		if (min>image[ii]) min=image[ii];
		if (max<image[ii]) max=image[ii];
	}
    for (unsigned long int ii = 0; ii< bins; ii++){
		histogram[ii]=0.0f;
	}
	for (unsigned long int ii = 0; ii< nxyz; ii++){
		double val = 0;
		if (max-min>0){
		 val = (image[ii] - min)*(bins-1)/(max-min);
		}
		int index = floor(val);
		histogram[index]++;
	}
   }else{ //MASK
	double min=image[0], max=image[0];
	for (unsigned long int ii = 1; ii< nxyz; ii++){
	   if(mask[ii]>minMaskValue){
		if (min>image[ii]) min=image[ii];
		if (max<image[ii]) max=image[ii];
           }
	}
        for (unsigned long int ii = 0; ii< bins; ii++){
		histogram[ii]=0.0f;
	}
	for (unsigned long int ii = 0; ii< nxyz; ii++){
	   if(mask[ii]>minMaskValue){
		double val = 0;
		if (max-min>0){
		 val = (image[ii] - min)*(bins-1)/(max-min);
		}
		int index = floor(val);
		histogram[index]++;
           }
	}
   }
 }



// *************************************
//  SSIM
// **************************************
template<typename T>
double SSIM(const T* GT, const T* I2, unsigned long int nxyz, const T* mask=NULL, double minMaskValue = 0.1){

  const T* I1 = GT;
  double meanI1=0;
  double meanI2=0;
  double sdI1=0;
  double sdI2=0;
  double sigmaI12=0;
  unsigned long int counter  = 0;
  double min = 0;
  double max = 1;
  double dynamicRange=1;

  if (mask){
    for (unsigned long int ijk=0; ijk<nxyz;ijk++){
        if ( mask[ijk] > minMaskValue ){
            if (counter==0){
                meanI1=I1[ijk];
                meanI2=I2[ijk];
            }else{
                meanI1+=I1[ijk];
                meanI2+=I2[ijk];
                if (max<GT[ijk]) max=GT[ijk];
                if (min>GT[ijk]) min=GT[ijk];

            }
            counter++;
        }
    }
  }else{
    counter = nxyz;
    for (unsigned long int ijk=0; ijk<nxyz;ijk++){
        meanI1+=I1[ijk];
        meanI2+=I2[ijk];
        if (max<GT[ijk]) max=GT[ijk];
        if (min>GT[ijk]) min=GT[ijk];
    }
  }
  if (max-min>1)
    dynamicRange=max-min;


  if (counter > 0){
      meanI1/=counter;
      meanI2/=counter;
  }
  if (mask){
    for (unsigned long int ijk=0; ijk<nxyz;ijk++){
      if ( mask[ijk] > minMaskValue ){
        sdI1+=pow(I1[ijk]-meanI1,2.0);
        sdI2+=pow(I2[ijk]-meanI2,2.0);
        sigmaI12+=(I1[ijk]-meanI1)*(I2[ijk]-meanI2);
      }
    }
  }else{
    for (unsigned long int ijk=0; ijk<nxyz;ijk++){
      sdI1+=pow(I1[ijk]-meanI1,2.0);
      sdI2+=pow(I2[ijk]-meanI2,2.0);
      sigmaI12+=(I1[ijk]-meanI1)*(I2[ijk]-meanI2);
    } 
  }

  if (counter > 1.0){
      sdI1=pow((1.0/(counter-1.0))*counter,0.5);
      sdI2=pow((1.0/(counter-1.0))*counter,0.5);
      sigmaI12=pow((1.0/(counter-1.0))*sigmaI12,0.5);
  }

  //C1=(K1*L)^2, where L is the dynamic range of the pixel values, and K1 << 1
  dynamicRange=1.0;
  double C1 = pow(0.00000001*dynamicRange,2);
  double C2 = pow(0.00000001*dynamicRange,2);
  double C3 = C2/2.0;
  double luminance=(2*meanI1*meanI2+C1)/(meanI1*meanI1+meanI2*meanI2+C1);
  double contrast=(2*sdI1*sdI2+C2)/(sdI1*sdI1+sdI2*sdI2+C2);
  double structure=(sigmaI12+C3)/(sdI1*sdI2+C3);
  double ssim=luminance*luminance + contrast*contrast + structure*structure;
  return ssim;
}



// *************************************
//  PSNR
// **************************************
template<typename T>
double PSNR(const T* GT, const T* I2, unsigned long int nxyz, const T* mask=NULL, double minMaskValue = 0.1){

  double maxGT=0;
  double MSE=0;
  unsigned long int counter  = 0;

  if (mask){
    for (unsigned long int ijk=0; ijk<nxyz;ijk++){
        if ( mask[ijk] > minMaskValue ){
            if (counter==0){
                maxGT=GT[ijk];
            }else{
                if ( maxGT<GT[ijk] ){
                  maxGT=GT[ijk];
                }
                MSE+=pow(GT[ijk]-I2[ijk], 2.0);
            }
            counter++;
        }
    }
  }else{
    counter = nxyz;
    for (unsigned long int ijk=0; ijk<nxyz;ijk++){
      if ( maxGT<GT[ijk] ){
        maxGT=GT[ijk];
      }
      counter++;
      MSE+=pow(GT[ijk]-I2[ijk], 2.0);
    }
  }
//  std::cerr<<"counter="<<counter<<"\n";
//  std::cerr<<"maxGT="<<maxGT<<"\n";
//  std::cerr<<"MSE_unnormalized="<<MSE<<"\n";
  MSE=(1.0/counter)*MSE;
//  std::cerr<<"MSE="<<MSE<<"\n";
  double PSNR_val=20.0*log10(maxGT/pow(MSE,0.5));
//  std::cerr<<"PSNR_val="<<PSNR_val<<"\n";
  if (PSNR_val != PSNR_val){
    PSNR_val = 0;
  }

  return PSNR_val;
}



// *************************************
//  crossCorrelationDistance
// **************************************
template<typename T>
double crossCorrelationDistance(const T* I1, const T* I2, unsigned long int nxyz, const T* mask=NULL, double minMaskValue = 0){
 double mean1 = 0;
 double mean2 = 0;
 double totalValues = 0;
 if (mask==NULL){
  for (unsigned long int ii = 0; ii<nxyz; ii++){
   mean1+=I1[ii];
   mean2+=I2[ii];
  }
  totalValues=nxyz;
 }else{
  for (unsigned long int ii = 0; ii<nxyz; ii++){
   if(mask[ii]>minMaskValue){
     mean1+=I1[ii];
     mean2+=I2[ii];
     totalValues++;
   }
  }
 }
 if ( totalValues > 0 ){
    mean1/=totalValues;
    mean2/=totalValues;
 }
 double numer=0;
 double denom1=0;
 double denom2=0; 
 if (mask==NULL){
         for (unsigned long int ii = 0; ii<nxyz; ii++){
          numer+=(I1[ii]-mean1)*(I2[ii]-mean2);
          denom1+=pow(I1[ii]-mean1,2);
          denom2+=pow(I2[ii]-mean2,2);
         }
 }else{
         for (unsigned long int ii = 0; ii<nxyz; ii++){
          if(mask[ii]>minMaskValue){
                  numer+=(I1[ii]-mean1)*(I2[ii]-mean2);
                  denom1+=pow(I1[ii]-mean1,2);
                  denom2+=pow(I2[ii]-mean2,2);
          }
         }
 }
 
 double CCC=0;
 double denomin=pow(denom1,0.5)*pow(denom2,0.5);
 
 if (numer==0.f){
   CCC=0;
 }else if (pow(denomin,2)<0.0000000001f){
  CCC=1;
 }else{
   CCC=numer/denomin;
 }
 return CCC;
 
}




// *************************************
//  SCI: structural Cross Correlation Index
// **************************************
template<typename T>
double SCI( T* I1, 
             T* I2, 
            unsigned long int nx, unsigned long int ny, unsigned long int nz, 
            double sigma, 
            const T* mask=NULL, 
            double minMaskValue = 0.1){
 
  unsigned long int nxyz = nx * ny * nz;
  double angpix = 1;
  int padding = 2;
  //compute first and second derivatives
  T * I1_proc= new T [ nxyz ];
  T * I2_proc= new T [ nxyz ];

double CC = crossCorrelationDistance(I1, I2,  nxyz, mask, minMaskValue);
if (CC != CC){
  CC=0;
} else if ( CC < 0 ){
  CC=0;
}
//std::cerr<<"CC="<<CC<<"\n";

  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	1,	 (T*) I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	1,	 (T*) I2,	I2_proc);
  double CC_tmp = crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  if (CC_tmp != CC_tmp){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp;
  }
  //  std::cerr<<"first derivative X="<<CC_tmp<<"\n";


  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	1,	(T*)  I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	1,	(T*)  I2,	I2_proc);
  CC_tmp = crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  if (CC_tmp != CC_tmp){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp;
  }
  //std::cerr<<"first derivative Y="<<CC_tmp<<"\n";

  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	2,	(T*)  I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	2,	(T*)  I2,	I2_proc);
  CC_tmp = crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  if (CC_tmp != CC_tmp){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp;
  }
  //std::cerr<<"second derivative X="<<CC_tmp<<"\n";


  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	2,	(T*)  I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	2,	(T*)  I2,	I2_proc);
  CC_tmp = crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  if (CC_tmp != CC_tmp){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp;
  }
  //std::cerr<<"score full="<<CC<<"\n";
  //std::cerr<<"second derivative Y="<<CC_tmp<<"\n";

  delete [] I1_proc;
  delete [] I2_proc;

  return CC;
}




// *************************************
//  SCI: Structural Difference Image 
// **************************************

template<typename T>
double SD_norm ( T* I1, 
              T* I2, 
              T* OutI, unsigned long int nxyz){

  long double mean1=0;
  long double mean2=0;
  for (unsigned long int ijk=0; ijk<nxyz; ijk++){
    mean1+=I1[ijk];
    mean2+=I2[ijk];
  }
  mean1/=nxyz;
  mean2/=nxyz;

  double mean = 0;
  for (unsigned long int ijk=0; ijk<nxyz; ijk++){
    OutI[ijk]=exp(-0.1*pow( I1[ijk]-mean1 -I2[ijk]+mean2, 2.0));
  }
//  mean /= (double)nxyz;
  for (unsigned long int ijk=0; ijk<nxyz; ijk++){
    OutI[ijk]/=5.0;
  }
  return 0;
}


// *************************************
//  SCI: Structural Difference Image 
// **************************************
template<typename T>
double SDIM ( T* I1, 
              T* I2, 
              T* OutI,
              unsigned long int nx, unsigned long int ny, unsigned long int nz, 
              double sigma, 
              const T* mask=NULL, 
              double minMaskValue = 0.1 ) {
 
  unsigned long int nxyz = nx * ny * nz;
  double angpix = 1;
  int padding = 2;
  //compute first and second derivatives
  T * I1_proc= new T [ nxyz ];
  T * I2_proc= new T [ nxyz ];
  T * tmpI= new T [ nxyz ];


  double CC = crossCorrelationDistance(I1, I2,  nxyz, mask, minMaskValue);
  if (CC != CC){
    CC=0;
  } else if ( CC < 0 ){
    CC=0;
  }
  SD_norm(I1,I2,tmpI,nxyz);  
  for (unsigned long int ijk=0; ijk<nxyz; ijk++){
    OutI[ijk]=tmpI[ijk];
  }
  writeMrcImage("I0.mrc", tmpI, nx,ny,nz);
  //std::cerr<<"CC="<<CC<<"\n";

  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	1,	 (T*) I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	1,	 (T*) I2,	I2_proc);
  double CC_tmp = crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  SD_norm(I1_proc,I2_proc,tmpI,nxyz);
  for (unsigned long int ijk=0; ijk<nxyz; ijk++){
    //OutI[ijk]+=tmpI[ijk];
    OutI[ijk]+= tmpI[ijk]; //OK
  }
  writeMrcImage("I1x.mrc", tmpI, nx,ny,nz);
  if (CC_tmp != CC_tmp){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp;
  }
  //  std::cerr<<"first derivative X="<<CC_tmp<<"\n";


  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	1,	(T*)  I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	1,	(T*)  I2,	I2_proc);
  CC_tmp = crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  SD_norm(I1_proc,I2_proc,tmpI,nxyz);
  for (unsigned long int ijk=0; ijk<nxyz; ijk++){
    OutI[ijk]+=tmpI[ijk];
  }
  writeMrcImage("I1y.mrc", tmpI, nx,ny,nz);
  if (CC_tmp != CC_tmp){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp;
  }
  //std::cerr<<"first derivative Y="<<CC_tmp<<"\n";

  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	2,	(T*)  I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	2,	(T*)  I2,	I2_proc);
  CC_tmp = crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  SD_norm(I1_proc,I2_proc,tmpI,nxyz);
  for (unsigned long int ijk=0; ijk<nxyz; ijk++){
    OutI[ijk]+=tmpI[ijk];
  }
  writeMrcImage("I1xx.mrc", tmpI, nx,ny,nz);
  if (CC_tmp != CC_tmp){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp;
  }
  //std::cerr<<"second derivative X="<<CC_tmp<<"\n";


  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	2,	(T*)  I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	2,	(T*)  I2,	I2_proc);
  CC_tmp = crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  SD_norm(I1_proc,I2_proc,tmpI,nxyz);
  for (unsigned long int ijk=0; ijk<nxyz; ijk++){
    OutI[ijk]+=tmpI[ijk];
  }
  writeMrcImage("I1yy.mrc", tmpI, nx,ny,nz);
  if (CC_tmp != CC_tmp){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp;
  }


  delete [] I1_proc;
  delete [] I2_proc;
  delete [] tmpI;
  return CC_tmp;
}





// *************************************
//  SCI: structural Cross Correlation Index
// **************************************
template<typename T>
double SCI_sqr( T* I1, 
             T* I2, 
            unsigned long int nx, unsigned long int ny, unsigned long int nz, 
            double sigma, 
            const T* mask=NULL, 
            double minMaskValue = 0.1){
 
  double minTolerance = 0.0;
  unsigned long int nxyz = nx * ny * nz;
  double angpix = 1;
  int padding = 2;
  //compute first and second derivatives
  T * I1_proc= new T [ nxyz ];
  T * I2_proc= new T [ nxyz ];

  T * I1_procReverse= new T [ nxyz ];
  T * I2_procReverse= new T [ nxyz ];


double CC = crossCorrelationDistance(I1, I2,  nxyz, mask, minMaskValue);
if (CC != CC){
  CC=0;
} else if ( CC < 0 ){
  CC=0;
}
//std::cerr<<"CC="<<CC<<"\n";

  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	1,	 (T*) I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	1,	 (T*) I2,	I2_proc);
  for (unsigned long int ijk=0;ijk<nxyz;ijk++){
    I1_proc[ijk]=pow(I1_proc[ijk],2.0);
    I2_proc[ijk]=pow(I2_proc[ijk],2.0);
  }
  //writeMrcImage("I2_proc.mrc", I2_proc, nx,ny,nz);
  //writeMrcImage("I2_procReverse.mrc", I2_procReverse, nx,ny,nz);
  double CC_tmp1 = minTolerance+crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  //std::cerr<<"dx tmp1="<<CC_tmp1<<"\n";
  if (CC_tmp1 != CC_tmp1){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp1;
  }


  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	1,	(T*)  I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	1,	(T*)  I2,	I2_proc);
  for (unsigned long int ijk=0;ijk<nxyz;ijk++){
    I1_proc[ijk]=pow(I1_proc[ijk],2.0);
    I2_proc[ijk]=pow(I2_proc[ijk],2.0);
  }
  CC_tmp1 = minTolerance+crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  //std::cerr<<"dy tmp1="<<CC_tmp1<<"\n";
  if (CC_tmp1 != CC_tmp1){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp1;
  }


  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	2,	(T*)  I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	2,	(T*)  I2,	I2_proc);
  for (unsigned long int ijk=0;ijk<nxyz;ijk++){
    I1_proc[ijk]=pow(I1_proc[ijk],2.0);
    I2_proc[ijk]=pow(I2_proc[ijk],2.0);
  }
  CC_tmp1 = minTolerance+crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  //std::cerr<<"dxx tmp1="<<CC_tmp1<<"\n";
  if (CC_tmp1 != CC_tmp1){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp1;
  }


  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	2,	(T*)  I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	2,	(T*)  I2,	I2_proc);
  for (unsigned long int ijk=0;ijk<nxyz;ijk++){
    I1_proc[ijk]=pow(I1_proc[ijk],2.0);
    I2_proc[ijk]=pow(I2_proc[ijk],2.0);
  }
  CC_tmp1 = minTolerance+crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  //std::cerr<<"dyy tmp1="<<CC_tmp1<<"\n";
  if (CC_tmp1 != CC_tmp1){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp1;
  }

  delete [] I1_proc;
  delete [] I2_proc;
  delete [] I1_procReverse;
  delete [] I2_procReverse;

  return CC;
}




template<typename T, typename U>
void squareDistanceImage(const T* I1, const T* I2, U* I12, unsigned long int nxyz, const T* mask=NULL){
 //double SD=0;
 if (mask==NULL){
   for (unsigned long int ii = 0; ii<nxyz; ii++){
     I12[ii]=pow(I1[ii]-I2[ii],2.0);
   }
 }else{
   for (unsigned long int ii = 0; ii<nxyz; ii++){
       if (mask[ii]>0){  
         I12[ii]=pow(I1[ii]-I2[ii],2.0);
       }else{
         I12[ii]=0;
       }
   }
 }
}



// *************************************
//  buildJointHistogram
// **************************************
template<typename T>
void buildJointHistogram(const T* image1, //!< input image1
        const T* image2, //!< input image2
        double * jointHistogram, //!< output histogram
        const unsigned long int bins, //!< number of bins for the histogram
        const unsigned long int nxyz,  //!< size (in pixel) of the image
        const T* mask=NULL,
        double minMaskValue=0
){
  //std::cerr<<"starting HIST\n";
  if (mask==NULL){
	double min1=image1[0], max1=image1[0];
	for (unsigned long int ii = 1; ii< nxyz; ii++){
		if (min1>image1[ii]) min1=image1[ii];
		if (max1<image1[ii]) max1=image1[ii];
	}
	double min2=image2[0], max2=image2[0];
	for (unsigned long int ii = 1; ii< nxyz; ii++){
		if (min2>image2[ii]) min2=image2[ii];
		if (max2<image2[ii]) max2=image2[ii];
	}
        for (unsigned long int ii = 0; ii< bins*bins; ii++){
		jointHistogram[ii]=0.0f;
	}
	for (unsigned long int ii = 0; ii< nxyz; ii++){
               	double val1 = 0;
               	double val2 = 0;
	        if (max1>min1)
        		val1 = (image1[ii] - min1)*(bins-1)/(max1-min1);
	        if (max2>min2)
        		val2 = (image2[ii] - min2)*(bins-1)/(max2-min2);
	        int index1 = floor(val1);
	        int index2 = floor(val2);
		jointHistogram[index1+bins*index2]++;
	}

  } else { // MASK
	double min1=image1[0], max1=image1[0];
	for (unsigned long int ii = 1; ii< nxyz; ii++){
	   if (  mask[ii] > minMaskValue ){
		if (min1>image1[ii] ) min1=image1[ii];
		if (max1<image1[ii] ) max1=image1[ii];
	   }
	}
	double min2=image2[0], max2=image2[0];
	for (unsigned long int ii = 1; ii< nxyz; ii++){
	   if ( ((long int) mask[ii]) > (long int)minMaskValue ){
		if ( min2>image2[ii] ) min2=image2[ii];
		if ( max2<image2[ii] ) max2=image2[ii];
           }
	}
  for (unsigned long int ii = 0; ii< bins*bins; ii++){
		jointHistogram[ii]=0.0f;
	}
	for (unsigned long int ii = 0; ii< nxyz; ii++){
	   if ( ((long int) mask[ii]) > (long int)minMaskValue){
               	double val1 = 0;
               	double val2 = 0;
	        if (max1>min1)
        		val1 = (image1[ii] - min1)*(bins-1)/(max1-min1);
	        if (max2>min2)
        		val2 = (image2[ii] - min2)*(bins-1)/(max2-min2);
	        int index1 = floor(val1);
	        int index2 = floor(val2);
		jointHistogram[index1+bins*index2]++;
           }
	}
  }
  //std::cerr<<"ending HIST\n";
}





// *************************************
//  computeEntropy
// **************************************
template<typename T>
double computeEntropy(const T* data1, unsigned long int nxyz, const T* mask=NULL){
 double entropy = 0;
 double sum=0;
 if (mask==NULL){
         for (unsigned long int ii=0; ii<nxyz; ii++){
          sum+=data1[ii];
         }
         if (sum == 0)
          return 0.0;
         for (unsigned long int ii=0; ii<nxyz; ii++){
          const double probability = data1[ii] / sum;
          if (probability > 0)
           entropy += - probability * log( probability ) / log( 2.0 );
         }
         return entropy;
 }else{ //MASK
         for (unsigned long int ii=0; ii<nxyz; ii++){
          if (mask[ii]>0){
            sum+=data1[ii];
          }
         }
         if (sum == 0)
          return 0.0;
         for (unsigned long int ii=0; ii<nxyz; ii++){
            if (mask[ii]>0){
                  const double probability = data1[ii] / sum;
                  if (probability > 0){
                   entropy += - probability * log( probability ) / log( 2.0 );
                  }              
            }
         }
         return entropy;
 }
}


// *************************************
//  NormalisedMutualInformation
// **************************************
template<typename T>
double NormalisedMutualInformation(const T* I1Target, 
                            const T* I2,
                            unsigned long int nxyz, 
                            unsigned long int numBins,
                            const T* mask=NULL,
                            double minMaskValue=0){

 double * histogram1 = new double [numBins];
 double * histogram2 = new double [numBins];
 double * histogram12 = new double [numBins*numBins];
 buildHistogram(I1Target, histogram1, numBins, nxyz,mask,minMaskValue);
 buildHistogram(I2, histogram2, numBins, nxyz,mask,minMaskValue);
 buildJointHistogram(I1Target, I2, histogram12, numBins, nxyz,mask,minMaskValue);
  
 double H1 = computeEntropy(histogram1, numBins);
 double H2 = computeEntropy(histogram2, numBins);
 double H12 = computeEntropy(histogram12, numBins*numBins);
 //std::cerr<<"H1="<<H1<<"     H2="<<H2<<"    H12="<<H12<<"\n";
 //double MutualInformation = H1+H2-H12;
 double NormalizedMutualInformation2 = 0;
 if (H12!=0) NormalizedMutualInformation2 = ( H1 + H2 ) / (2*H12);
 
 delete [] histogram1;
 delete [] histogram2;
 delete [] histogram12;
 return NormalizedMutualInformation2;

}





// **********************************
//  largest eigen value of symmetric matrix in the form
//    | a  b |
//    | b  c |
//
double largestEigenvalue(double a, double b, double c){
    double delta = pow(a-c,2.0)+4.0*b*b;
    if (delta<0) return 0;
    delta=pow(delta,0.5);
//    return ((a+c+delta)/2.0);
    return (a+c+delta);
}




#endif
