/*
 * File: ctfLibs.h
 * (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
 */
#ifndef __CTF_LIBS__H___
#define __CTF_LIBS__H___

/**
* The CTF (c) and the envelope function (e) are both spatially invariant
* The noise (n) is independent and additive.
* \f$ i(x, y) = c(x, y) \ast e(x, y) \ast f(x, y) + n(x, y) \f$. Where i is the image, c the the ctf image, e the envelope function, 
* and f the projection of the particle being imaged.
* \cite Mallick
*/
#define PI 3.14159265359

typedef struct CTFParameters {
        double local_Cs;// = cs_vector[ii] * 1e7;
        double local_kV; // = ctfVoltage_vector[ii] * 1e3;
        double rad_azimuth; //= azimuthalAngle_vector[ii] * PI/180.0;
        double defocus_average;
        double defocus_deviation;
        double lambda;
        double K1;
        double K2;
        double K3;
        double K4;
        double K5;
        double scale;
        CTFParameters(double SphericalAberration, double voltage, double DefocusAngle, double DefocusU, double DefocusV, double AmplitudeContrast, double Bfac, double phase_shift, double _scale=1){
            local_Cs = SphericalAberration * 1e7;
            local_kV = voltage * 1e3;
            rad_azimuth = DefocusAngle * PI/180.0;
            // Average focus and deviation
            defocus_average   = -(DefocusU + DefocusV) * 0.5;
            defocus_deviation = -(DefocusU - DefocusV) * 0.5;        
            lambda=12.2643247 / sqrt(local_kV * (1. + local_kV * 0.978466e-6));
            K1 = PI / 2.0 * 2.0 * lambda;
            K2 = PI / 2.0 * local_Cs * pow( lambda , 3.0 );
            K3 = atan(AmplitudeContrast/sqrt(1-AmplitudeContrast*AmplitudeContrast));
            K4 = -Bfac / 4.0;
            K5 = phase_shift* PI/180.0;     // Phase shift in radians
            scale=_scale;
        };
} CTFParameters;


    /// Compute Deltaf at a given coordinate
    double getDeltaF(double X, double Y, 
                    const double & rad_azimuth, 
                    const double & defocus_average, 
                    const double & defocus_deviation)  
    {
        //if (ABS(X) < XMIPP_EQUAL_ACCURACY &&
        //    ABS(Y) < XMIPP_EQUAL_ACCURACY)
        //    return 0;

        double ellipsoid_ang = atan2(Y, X) - rad_azimuth;
        
        // * For a derivation of this formulae confer
        // * Principles of Electron Optics page 1380
        // * in particular term defocus and twofold axial astigmatism
        // * take into account that a1 and a2 are the coefficient
        // * of the zernike polynomials difference of defocus at 0
        // * and at 45 degrees. In this case a2=0
        double cos_ellipsoid_ang_2 = cos(2*ellipsoid_ang);
        return (defocus_average + defocus_deviation*cos_ellipsoid_ang_2);

    }



inline double getCTF(double X, double Y,
                    const double & rad_azimuth, 
                    const double & defocus_average, 
                    const double & defocus_deviation,
                    const double & K1,
                    const double & K2,
                    const double & K3,
                    const double & K4,
                    const double & K5,
                    const double scale=1.0,
                    bool do_abs = false, bool do_only_flip_phases = false,
                    bool do_intact_until_first_peak = false, bool do_damping = true) {
        double u2 = X * X + Y * Y;
        double u4 = u2 * u2;

        // if (u2>=ua2) return 0;
        double deltaf = getDeltaF(X, Y, rad_azimuth, defocus_average, defocus_deviation);
        
        double argument = K1 * deltaf * u2 + K2 * u4 - K5 - K3;
        double returnValue;
        if (do_intact_until_first_peak && abs(argument) < PI/2.0) {
        	returnValue = 1.0;
        } else {
            returnValue = -sin(argument);
        }

        if (do_damping)
        {
        	double E = exp(K4 * u2); // B-factor decay (K4 = -Bfac/4);
        	returnValue *= E;
        }
        if (do_abs) {
        	returnValue = abs(returnValue);
        } else if (do_only_flip_phases)
        {
        	returnValue = (returnValue < 0.0) ? -1.0 : 1.0;
        }
        return scale * returnValue;
    }

inline double getCTF(double X, double Y,
                    const CTFParameters & ctf_parameters,
                    bool do_abs = false, bool do_only_flip_phases = false,
                    bool do_intact_until_first_peak = false, bool do_damping = true){
        double scale = 1;
        return getCTF( X, Y, ctf_parameters.rad_azimuth, 
                    ctf_parameters.defocus_average, 
                    ctf_parameters.defocus_deviation,
                    ctf_parameters.K1,
                    ctf_parameters.K2,
                    ctf_parameters.K3,
                    ctf_parameters.K4,
                    ctf_parameters.K5,
                    scale, do_abs, do_only_flip_phases, do_intact_until_first_peak, do_damping);
}


/*
template<typename T>
void multiplyImageCftCentered( T * mapIn,  T * ctfCenteredImage,  T * mapOut,const unsigned long int nx,  const unsigned long int ny){
    const double small_epsilon = 0.0001;
    const unsigned long int nxy = nx * ny;

    //compute mean
    double meanIn = 0;
    for (unsigned long int jj = 0, ij=0; jj < ny; jj++){
        for (unsigned long int ii = 0; ii < nx; ii++, ij++){
            mapIn[ij]*=pow(-1,ii+jj);
            meanIn+=mapIn[ij];
        }
    }
    meanIn/=nxy;

	fftw_complex *I, *F;
	I = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxy);
    F = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxy);
    


    //MAP_IN
    //for (unsigned long int ij = 0; ij < nxy; ij++){
	//				I[ij][0] = (mapIn[nxy-ij-1]-meanIn);
                    //I[ij][0] *= (-1)^(i + j);
	//				I[ij][1] = 0.0; //no imaginary part
    //}

   for (unsigned long int jj = 0, ij=0; jj < ny; jj++){
       for (unsigned long int ii = 0; ii < nx; ii++, ij++){
					I[ij][0] = (mapIn[nxy-ij-1]-meanIn);
					I[ij][1] = 0.0; //no imaginary part
        }
   }



    fftw_plan planForwardMap = fftw_plan_dft_2d(ny, nx, I, F, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(planForwardMap);
    fftw_destroy_plan(planForwardMap);
    for (unsigned long int yxIndex=0; yxIndex < nxy; yxIndex++){
        F[ yxIndex ][0] = ctfCenteredImage[ yxIndex ] * F[ yxIndex ][0]/nxy;
        F[ yxIndex ][1] = ctfCenteredImage[ yxIndex ] * F[ yxIndex ][1]/nxy;
    }

    
    //fftw_plan planForwardOriginalImage;
    fftw_plan planForwardOriginalImage = fftw_plan_dft_2d(ny, nx, F, I, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(planForwardOriginalImage);
    fftw_destroy_plan(planForwardOriginalImage);

   for (unsigned long int jj = 0, ij=0; jj < ny; jj++){
       for (unsigned long int ii = 0; ii < nx; ii++, ij++){
            mapOut[ij] = (meanIn+I[nxy-ij-1][0])*pow(-1,ii+jj);
        }
   }

	fftw_free(I);
	fftw_free(F);

}
*/

template<typename T>
int computeCtfCenteredImage(T * ctfCenteredImage2D, CTFParameters & ctf_parameters, unsigned long int nx, unsigned long int ny, double angpix){
    
    double xs = nx * angpix;
    double ys = ny * angpix;
    //CTFParameters(double SphericalAberration, double voltage, double DefocusAngle, double DefocusU, double DefocusV, double AmplitudeContrast, double Bfac, double phase_shift, double _scale=1){
    unsigned long int X0 = floor(nx/2.0);
    unsigned long int Y0 = floor(ny/2.0);
    for (unsigned long int yy=0, zyxIndex=0; yy<ny; yy++){
        signed long int yyy = yy - Y0;
        for (unsigned long int xx=0; xx<nx; xx++, zyxIndex++){
            signed long int xxx = xx - X0;
            double x = xxx / xs;
            double y = yyy / ys;
            double value=getCTF(x, y, 
                ctf_parameters.rad_azimuth, ctf_parameters.defocus_average, ctf_parameters.defocus_deviation, 
                ctf_parameters.K1, ctf_parameters.K2, ctf_parameters.K3, ctf_parameters.K4, ctf_parameters.K5, 
                ctf_parameters.scale);
            ctfCenteredImage2D[zyxIndex]=value;
        }
    }
    return 0;
}



// *****************************
// *****************************
//
//  computeBlurredProfileImage
//
// *****************************
template<typename T>
int computeBlurredProfileImage(T * ctfCenteredImage2D, unsigned long int nx, unsigned long int ny, double sigma, double angpix){
    
    double xs = nx * angpix;
    double ys = ny * angpix;
    //CTFParameters(double SphericalAberration, double voltage, double DefocusAngle, double DefocusU, double DefocusV, double AmplitudeContrast, double Bfac, double phase_shift, double _scale=1){
    unsigned long int X0 = floor(nx/2.0);
    unsigned long int Y0 = floor(ny/2.0);
    for (unsigned long int yy=0, zyxIndex=0; yy<ny; yy++){
        signed long int yyy = yy - Y0;
        for (unsigned long int xx=0; xx<nx; xx++, zyxIndex++){
            signed long int xxx = xx - X0;
            double x = xxx / xs;
            double y = yyy / ys;
            //double value=1.0;
        double u2 = x * x + y * y;
        double u4 = u2 * u2;

            //value = exp(-0.5*pow(2*3.14*u4*angpix,2)*pow(sigma,2)); //correct 
            double value = exp(-0.5*u4*sigma*sigma);  

            /*double value=getCTF(x, y, 
                ctf_parameters.rad_azimuth, ctf_parameters.defocus_average, ctf_parameters.defocus_deviation, 
                ctf_parameters.K1, ctf_parameters.K2, ctf_parameters.K3, ctf_parameters.K4, ctf_parameters.K5, 
                ctf_parameters.scale);*/
            ctfCenteredImage2D[zyxIndex]=value;
        }
    }
    return 0;
}



#endif
