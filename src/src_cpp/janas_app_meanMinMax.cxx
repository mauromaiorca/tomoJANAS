/*
 * File: janas_app_meanMinMax.h
 * (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
 */

#include <iostream>
#include <complex>
#include <fstream>
#include <vector>
#include <iterator>     // std::back_inserter
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <sstream>
#include <string>
#include <cstdio>
#include <cmath>
#include <iomanip>      // std::setprecision
#include <arpa/inet.h> //used for htonl (big/little endian check)




#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef TWO_PI
#define TWO_PI 6.28318530717958647
#endif

#include "mrcIO.h"
typedef float WorkingPixelType;
typedef float OutputPixelType;
const unsigned int Dimension = 3;



/* ******************************************
 *  USAGE
 ***************************************** */
void usage(  char ** argv ){
    std::cerr<<"\n";
    std::cerr<<"Usage: " << argv[0] << " locresMapFileName.mrc maskFileName.mrc \n";
    std::cerr<<"\n";
    exit(1);
}


// *************************************
//
// retrieveInputParameters
// *************************************
typedef struct inputParametersType {
    char * locresMapFileName;
    char * maskFileName;
}inputParametersType;


void retrieveInputParameters(inputParametersType * parameters, int argc, char** argv){
    parameters->locresMapFileName= NULL;
    parameters->maskFileName= NULL;

    if (argc >= 2){
     parameters->locresMapFileName=argv[1];
    }

    if (argc >= 3){
     parameters->maskFileName=argv[2];
    }
    
    if (!parameters->locresMapFileName && !parameters->maskFileName){
      usage(argv);
    }

    
  }


// **********************************
//
//    INT MAIN
//
// **********************************
int main( int argc, char **argv ){
	inputParametersType parameters;
	retrieveInputParameters(&parameters, argc, argv);	
        //std::cout<<"DEBUG: start the function\n";
	
	
	
           MRCHeader imageHeader;
           readHeaderMrc(parameters.locresMapFileName, imageHeader);
           unsigned long int nx=imageHeader.nx, ny=imageHeader.ny, nz=imageHeader.nz;
           unsigned long int nxyz=nx*ny*nz;
           WorkingPixelType * I = new WorkingPixelType [nxyz];
           WorkingPixelType * maskI = new WorkingPixelType [nxyz];
           readMrcImage(parameters.locresMapFileName, I, imageHeader);
           readMrcImage(parameters.maskFileName, maskI, imageHeader);
           long int firstIdx = 0;
           double mean=0;
           double min=0;
           double max=0;
           double counter=1;
           std::vector<double> values;
           for (unsigned long int ii=0, found = 0; ii<nxyz && found==0;ii++){
             if (maskI[ii]>0.2){
               found = 1;
               firstIdx=ii;
               mean=I[ii];
               min=I[ii];
               max=I[ii];
               values.push_back(I[ii]);
             }
           }
           for (unsigned long int ii=firstIdx+1; ii<nxyz; ii++){
             if (maskI[ii]>0.2){
               firstIdx=ii;
               mean+=I[ii];
               if (min>I[ii]) min=I[ii];
               if (max<I[ii]) max=I[ii];
               counter++;
               values.push_back(I[ii]);
             }
           }
           mean/=counter;
           
        double sd = 0;
        for (unsigned long int ii=0; ii<nxyz;ii++){
                    if (maskI[ii]>0.2){
                        sd+=pow(I[ii]-mean,2.0);
                    }
        }
        if ( counter - 1.0 > 0 ){
            sd = sd / (counter - 1.0 );
        }
        std::sort (values.begin(), values.end());


        //std::cout<<"DEBUG:"<<min<<","<< mean <<","<<max<<"\n";
        int firstQuartileIdx = values.size()/4.0;
        int lastQuartileIdx = values.size()*3.0/4.0;
        double minFinal = 0;
        double maxFinal = 0;
        if (values.size()>0){
            if (values.size()>4){
                minFinal=values[2];
                maxFinal=values[values.size()-3];
            }else{
                minFinal=values[0];
                maxFinal=values[values.size()-1];
            }

        }
        maxFinal = values[(unsigned long int)values.size()*(3.70/4.0)];
        std::cout<<minFinal<<","<<values[firstQuartileIdx]<<","<<mean<<","<<values[lastQuartileIdx]<<","<<maxFinal<<"\n";

        delete [] I;
        delete [] maskI;

}






