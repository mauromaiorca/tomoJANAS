/*
 * File: randomLibs.h
 * (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
 */

#ifndef ___RANDOM__LIBS___
#define ___RANDOM__LIBS___



#include <cmath>
#include <iostream>
#include <iomanip>      // std::setprecision
#include <complex>
#include <cstdio>
#include <cstdlib>
//#include <time.h>
#include <sys/time.h>
#include <cstdlib>



/*
//TO BE USED THIS WAY:
                FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(F2D)
                {
                    DIRECT_A3D_ELEM(F2D, k, i, j).real += rnd_gaus(0., stddev_white_noise);
                    DIRECT_A3D_ELEM(F2D, k, i, j).imag += rnd_gaus(0., stddev_white_noise);
                }
*/
// FROM RELION
// Gaussian distribution ...................................................
float rnd_gaus(float mu, float sigma, bool initX=false)
{
  float U1, U2, W, mult;
  static float X1, X2;
  static int call = 0;
  static int FirstCall = 0;

  if (initX){
    X1=1;
    X2=1;
    call = 0;
    FirstCall = 0;    
  }
  //std::cerr<<"      ->"<<FirstCall<<"\n";
  
  if (sigma == 0)
	  return mu;
	  
  if (FirstCall==0){
     //srand(time(NULL)+rand()%100);
     struct timeval time; 
     gettimeofday(&time,NULL);
     // microsecond has 1 000 000
     // Assuming you did not need quite that accuracy
     // Also do not assume the system clock has that accuracy.
     srand((time.tv_sec * 1000) + (time.tv_usec / 1000));

   FirstCall=1;
  }


  if (call == 1)
  {
      call = !call;
      return (mu + sigma * (float) X2);
  }

  do
  {
      U1 = -1 + ((float) rand () / RAND_MAX) * 2;
      U2 = -1 + ((float) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
  }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (float) X1);

}


int rnd_01()
{
  static int FirstCall0 = 0;
      
  if (FirstCall0==0){
     //srand(time(NULL)+rand()%100);
     struct timeval time;
     gettimeofday(&time,NULL);
     // microsecond has 1 000 000
     // Assuming you did not need quite that accuracy
     // Also do not assume the system clock has that accuracy.
     srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
     FirstCall0=1;
  }

  return rand() % 2;

}

int rnd_Int()
{
  static int FirstCall0 = 0;
      
  if (FirstCall0==0){
     //srand(time(NULL)+rand()%100);
     struct timeval time;
     gettimeofday(&time,NULL);
     // microsecond has 1 000 000
     // Assuming you did not need quite that accuracy
     // Also do not assume the system clock has that accuracy.
     srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
     FirstCall0=1;
  }

  return rand() % 10000;

}

double rnd_double()
{
  static int FirstCall0 = 0;
      
  if (FirstCall0==0){
     //srand(time(NULL)+rand()%100);
     struct timeval time;
     gettimeofday(&time,NULL);
     // microsecond has 1 000 000
     // Assuming you did not need quite that accuracy
     // Also do not assume the system clock has that accuracy.
     srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
     FirstCall0=1;
  }

  return rand(); //10501 is a prime number

}


#endif
