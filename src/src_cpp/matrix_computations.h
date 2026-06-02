/*
 * File: matrix_computations.h
 * (C) 2022 Mauro Maiorca 
 */

#ifndef ___MATRIX_COMPUTE___H___
#define ___MATRIX_COMPUTE___H___
#include <vector>
#include "mrcIO.h"
// ****************************
// matrixMultiplication (useful for computing the transform matrix)
//   space of matrixes already allocated
//   size of outM=x2*y1
//
void matrixMultiplication(double M1[], double M2[], double outM[],  unsigned int x1,  unsigned int y1, unsigned int x2,  unsigned int y2, bool verbose=false)
{
    unsigned int xOut = x2;
    unsigned int yOut = y1;
    unsigned int c, d, i;
    double sum =0;

   if(verbose){
        printf("M1 is:\n");
        for ( c = 0 ; c < y1 ; c++ ){
            for ( d = 0 ; d < x1 ; d++ )
                printf("%f\t", M1[d+x1*c]);
            printf("\n");
        }
        printf("M2 is:\n");
        for ( c = 0 ; c < y2 ; c++ ){
            for ( d = 0 ; d < x2 ; d++ )
                printf("%f\t", M2[d+x2*c]);
            printf("\n");
        }
    }


   if (x1!=y2){
     if(verbose)
      std::cerr<<"not possible to multiply\n";
     return;
   }
   for (i=0;i<y1*x2;i++){
     outM[i]=0;
   }    

    // matrix multiplication operation
    for ( unsigned int jj = 0 ; jj < y1 ; jj++ ){
           unsigned int jjx1=jj*x1;
           unsigned int jjx2=jj*x2;
           for ( unsigned int ii = 0 ; ii< x2 ; ii++ ){
                sum = 0;
                for ( unsigned int kk = 0 ; kk < x1 ; kk++ ){
                    sum = sum + M1[kk+jjx1]*M2[ii+kk*x2];
                    if(verbose) std::cerr<<M1[kk+jjx1]<<"*"<<M2[ii+kk*x2]<<", ";
                }
                if(verbose) std::cerr<<"="<<sum<<"\n";
                outM[ii+jjx2] = sum;
                
            }
    }


        //Printing the final product matrix
    if(verbose){
        printf("\nThe product of entered matrices is:\n");
        for (unsigned int jj = 0 ; jj < yOut ; jj++ ){
            for ( unsigned int ii = 0 ; ii < xOut ; ii++ ) {
                printf("%f\t", outM[ii+jj*xOut]);
            }
            printf("\n");
        }
    }
 
}


// ****************************
// transformMatrix
// returns a 3x3 transform matrix, rigid 
void transformMatrixRelion(double M[], double Phi, double Theta, double Psi, double tx=0, double ty=0, double tz=0){
 
    double ca, sa, cb, sb, cg, sg;
    double cc, cs, sc, ss;
    //Phi=fmod(Phi,360.0);
    //Theta=fmod(Theta,180.0);
    //Psi=fmod(Psi,360.0);
    double pi_180=M_PI/double(180.0);
    double alpha = Phi*pi_180;
    double beta  = Theta*pi_180;
    double gamma = Psi*pi_180;

    ca = cos(alpha);
    cb = cos(beta);
    cg = cos(gamma);
    sa = sin(alpha);
    sb = sin(beta);
    sg = sin(gamma);
    cc = cb * ca;
    cs = cb * sa;
    sc = sb * ca;
    ss = sb * sa;

  M[0]=cg * cc - sg * sa;    M[4]=cg * cs + sg * ca;    M[8]=-cg * sb;    M[3]=tx;
  M[1]=-sg * cc - cg * sa;   M[5]=-sg * cs + cg * ca;   M[9]=sg * sb;     M[7]=ty;
  M[2]=sc;                   M[6]=ss;                   M[10]=cb;         M[11]=tz;
  M[12]=0;                   M[13]=0;                   M[14]=0;          M[15]=1;

}





//    [a11 a12 a13 a14]
//    [a21 a22 a23 a24]
// M= [a31 a32 a33 a34]
//    [a41 a42 a43 a44]
void inverseMatrix(double M[]){
/*
    double A11=M[5]*M[10]*M[15] + M[9]*M[14]*M[7] + M[13]*M[6]*M[11]
                -M[1+4*3]*M[2+4*2]*M[3+4*1] - M[1+4*2]*M[2+4*1]*M[3+4*3] - M[1+4*1]*M[2+4*3]*M[3+4*2];
    double A12= -M[0+4*1]*M[2+4*2]*M[3+4*3] - M[0+4*2]*M[2+4*3]*M[3+4*1] - M[0+4*3]*M[2+4*1]*M[3+4*2]
                +M[0+4*3]*M[2+4*2]*M[3+4*1] + M[0+4*2]*M[2+4*1]*M[3+4*3] + M[0+4*1]*M[2+4*3]*M[3+4*2];
    double A13=M[0+4*1]*M[1+4*2]*M[3+4*3] + M[0+4*2]*M[1+4*3]*M[3+4*1] + M[0+4*3]*M[1+4*1]*M[3+4*2]
                -M[0+4*3]*M[1+4*2]*M[3+4*1] - M[0+4*2]*M[1+4*1]*M[3+4*3] - M[0+4*1]*M[1+4*3]*M[3+4*2];
    double A14= -M[0+4*1]*M[1+4*2]*M[2+4*3] - M[0+4*2]*M[1+4*3]*M[2+4*1] - M[0+4*3]*M[1+4*1]*M[2+4*2]
                +M[0+4*3]*M[1+4*2]*M[2+4*1] + M[0+4*2]*M[1+4*1]*M[2+4*3] + M[0+4*1]*M[1+4*3]*M[2+4*2];

    double A21= -M[1+4*0]*M[2+4*2]*M[3+4*3] - M[1+4*2]*M[2+4*3]*M[3+4*0] - M[1+4*3]*M[2+4*0]*M[3+4*2]
                +M[1+4*3]*M[2+4*2]*M[3+4*0] + M[1+4*2]*M[2+4*0]*M[3+4*3] + M[1+4*0]*M[2+4*3]*M[3+4*2];
    double A22=M[0+4*0]*M[2+4*2]*M[3+4*3] + M[0+4*2]*M[2+4*3]*M[3+4*0] + M[0+4*3]*M[2+4*0]*M[3+4*2]
                -M[0+4*3]*M[2+4*2]*M[3+4*0] - M[0+4*2]*M[2+4*0]*M[3+4*3] - M[0+4*0]*M[2+4*3]*M[3+4*2];
    double A23= -M[0+4*0]*M[1+4*2]*M[3+4*3] - M[0+4*2]*M[1+4*3]*M[3+4*0] - M[0+4*3]*M[1+4*0]*M[3+4*2]
                +M[0+4*3]*M[1+4*2]*M[3+4*0] + M[0+4*2]*M[1+4*0]*M[3+4*3] + M[0+4*0]*M[1+4*3]*M[3+4*2];
    double A24=M[0+4*0]*M[1+4*2]*M[2+4*3] + M[0+4*2]*M[1+4*3]*M[2+4*0] + M[0+4*3]*M[1+4*0]*M[2+4*2]
                -M[0+4*3]*M[1+4*2]*M[2+4*0] - M[0+4*2]*M[1+4*0]*M[2+4*3] - M[0+4*0]*M[1+4*3]*M[2+4*2];


    double A31=M[1+4*0]*M[2+4*1]*M[3+4*3] + M[1+4*1]*M[2+4*3]*M[3+4*0] + M[1+4*3]*M[2+4*0]*M[3+4*1]
                -M[1+4*3]*M[2+4*1]*M[3+4*0] - M[1+4*1]*M[2+4*0]*M[3+4*3] - M[1+4*0]*M[2+4*3]*M[3+4*1];
    double A32= -M[0+4*0]*M[2+4*1]*M[3+4*3] - M[0+4*1]*M[2+4*3]*M[3+4*0] - M[0+4*3]*M[2+4*0]*M[3+4*1]
                +M[0+4*3]*M[2+4*1]*M[3+4*0] + M[0+4*1]*M[2+4*0]*M[3+4*3] + M[0+4*0]*M[2+4*3]*M[3+4*1];
    double A33=M[0+4*0]*M[1+4*1]*M[3+4*3] + M[0+4*1]*M[1+4*3]*M[3+4*0] + M[0+4*3]*M[1+4*0]*M[3+4*1]
                -M[0+4*3]*M[1+4*1]*M[3+4*0] - M[0+4*1]*M[1+4*0]*M[3+4*3] - M[0+4*0]*M[1+4*3]*M[3+4*1];
    double A34= -M[0+4*0]*M[1+4*1]*M[2+4*3] - M[0+4*1]*M[1+4*3]*M[2+4*0] - M[0+4*3]*M[1+4*0]*M[2+4*1]
                +M[0+4*3]*M[1+4*1]*M[2+4*0] + M[0+4*1]*M[1+4*0]*M[2+4*3] + M[0+4*0]*M[1+4*3]*M[2+4*1];


    double A41= -M[1+4*0]*M[2+4*1]*M[3+4*2] - M[1+4*1]*M[2+4*2]*M[3+4*0] - M[1+4*2]*M[2+4*0]*M[3+4*1]
                +M[1+4*2]*M[2+4*1]*M[3+4*0] + M[1+4*1]*M[2+4*0]*M[3+4*2] + M[1+4*0]*M[2+4*2]*M[3+4*1];
    double A42=M[0+4*0]*M[2+4*1]*M[3+4*2] + M[0+4*1]*M[2+4*2]*M[3+4*0] + M[0+4*2]*M[2+4*0]*M[3+4*1]
                -M[0+4*2]*M[2+4*1]*M[3+4*0] - M[0+4*1]*M[2+4*0]*M[3+4*2] - M[0+4*0]*M[2+4*2]*M[3+4*1];
    double A43= -M[0+4*0]*M[1+4*1]*M[3+4*2] - M[0+4*1]*M[1+4*2]*M[3+4*0] - M[0+4*2]*M[1+4*0]*M[3+4*1]
                +M[0+4*2]*M[1+4*1]*M[3+4*0] + M[0+4*1]*M[1+4*0]*M[3+4*2] + M[0+4*0]*M[1+4*2]*M[3+4*1];
    double A44=M[0+4*0]*M[1+4*1]*M[2+4*2] + M[0+4*1]*M[1+4*2]*M[2+4*0] + M[0+4*2]*M[1+4*0]*M[2+4*1]
                -M[0+4*2]*M[1+4*1]*M[2+4*0] - M[0+4*1]*M[1+4*0]*M[2+4*2] - M[0+4*0]*M[1+4*2]*M[2+4*1];

*/
    double A11=M[5]*M[10]*M[15] + M[9]*M[14]*M[7] + M[13]*M[6]*M[11]
                -M[13]*M[10]*M[7] - M[9]*M[6]*M[15] - M[5]*M[14]*M[11];
    double A12= -M[4]*M[10]*M[15] - M[8]*M[14]*M[7] - M[12]*M[6]*M[11]
                +M[12]*M[10]*M[7] + M[8]*M[6]*M[15] + M[4]*M[14]*M[11];
    double A13=M[4]*M[9]*M[15] + M[8]*M[13]*M[7] + M[12]*M[5]*M[11]
                -M[12]*M[9]*M[7] - M[8]*M[5]*M[15] - M[4]*M[13]*M[11];
    double A14= -M[4]*M[9]*M[14] - M[8]*M[13]*M[6] - M[12]*M[5]*M[10]
                +M[12]*M[9]*M[6] + M[8]*M[5]*M[14] + M[4]*M[13]*M[10];

    double A21= -M[1]*M[10]*M[15] - M[9]*M[14]*M[3] - M[13]*M[2]*M[11]
                +M[13]*M[10]*M[3] + M[9]*M[2]*M[15] + M[1]*M[14]*M[11];
    double A22=M[0]*M[10]*M[15] + M[8]*M[14]*M[3] + M[12]*M[2]*M[11]
                -M[12]*M[10]*M[3] - M[8]*M[2]*M[15] - M[0]*M[14]*M[11];
    double A23= -M[0]*M[9]*M[15] - M[8]*M[13]*M[3] - M[12]*M[1]*M[11]
                +M[12]*M[9]*M[3] + M[8]*M[1]*M[15] + M[0]*M[13]*M[11];
    double A24=M[0]*M[9]*M[14] + M[8]*M[13]*M[2] + M[12]*M[1]*M[10]
                -M[12]*M[9]*M[2] - M[8]*M[1]*M[14] - M[0]*M[13]*M[10];


    double A31=M[1]*M[6]*M[15] + M[5]*M[14]*M[3] + M[13]*M[2]*M[7]
                -M[13]*M[6]*M[3] - M[5]*M[2]*M[15] - M[1]*M[14]*M[7];
    double A32= -M[0]*M[6]*M[15] - M[4]*M[14]*M[3] - M[12]*M[2]*M[7]
                +M[12]*M[6]*M[3] + M[4]*M[2]*M[15] + M[0]*M[14]*M[7];
    double A33=M[0]*M[5]*M[15] + M[4]*M[13]*M[3] + M[12]*M[1]*M[7]
                -M[12]*M[5]*M[3] - M[4]*M[1]*M[15] - M[0]*M[13]*M[7];
    double A34= -M[0]*M[5]*M[14] - M[4]*M[13]*M[2] - M[12]*M[1]*M[6]
                +M[12]*M[5]*M[2] + M[4]*M[1]*M[14] + M[0]*M[13]*M[6];


    double A41= -M[1]*M[6]*M[11] - M[5]*M[10]*M[3] - M[9]*M[2]*M[7]
                +M[9]*M[6]*M[3] + M[5]*M[2]*M[11] + M[1]*M[10]*M[7];
    double A42=M[0]*M[6]*M[11] + M[4]*M[10]*M[3] + M[8]*M[2]*M[7]
                -M[8]*M[6]*M[3] - M[4]*M[2]*M[11] - M[0]*M[10]*M[7];
    double A43= -M[0]*M[5]*M[11] - M[4]*M[9]*M[3] - M[8]*M[1]*M[7]
                +M[8]*M[5]*M[3] + M[4]*M[1]*M[11] + M[0]*M[9]*M[7];
    double A44=M[0]*M[5]*M[10] + M[4]*M[9]*M[2] + M[8]*M[1]*M[6]
                -M[8]*M[5]*M[2] - M[4]*M[1]*M[10] - M[0]*M[9]*M[6];

    M[0]=A11;M[1]=A12;M[2]=A13;M[3]=A14;
    M[4]=A21;M[5]=A22;M[6]=A23;M[7]=A24;
    M[8]=A31;M[9]=A32;M[10]=A33;M[11]=A34;
    M[12]=A41;M[13]=A42;M[14]=A43;M[15]=A44;
}






//SORT
template<typename T>
std::vector<T> bubbleSortAscendingValues(const std::vector<T> valuesIn){
  std::vector<T> valuesOut = valuesIn;
  long int size = valuesOut.size();
  for (long int i = (size - 1); i > 0; i--)
    {
      for (long int j = 1; j <= i; j++)
	{
	  if (valuesOut[j - 1] > valuesOut[j])
	    {
	      T temp = valuesOut[j - 1];
	      valuesOut[j - 1] = valuesOut[j];
	      valuesOut[j] = temp;
	    }
	}
    }
  return valuesOut;
}


void invertRigidTransformInPlace(double M[16]) {
    // Temporary buffers for rotation and translation
    double R[3][3] = {
        {M[0], M[1], M[2]},
        {M[4], M[5], M[6]},
        {M[8], M[9], M[10]}
    };
    double T[3] = {M[12], M[13], M[14]};
    double tempT[3];

    // Compute the transposed rotation matrix (in-place)
    for (int i = 0; i < 3; i++) {
        for (int j = i+1; j < 3; j++) {
            double temp = R[i][j];
            R[i][j] = R[j][i];
            R[j][i] = temp;
        }
    }

    // Compute the translated vector
    for (int i = 0; i < 3; i++) {
        tempT[i] = -(R[0][i] * T[0] + R[1][i] * T[1] + R[2][i] * T[2]);
    }

    // Update the input matrix M
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            M[i*4 + j] = R[i][j];
        }
    }

    M[12] = tempT[0];
    M[13] = tempT[1];
    M[14] = tempT[2];

    // The bottom row remains unchanged ([0 0 0 1])
}


//Matrix this shape:
//[a11 a12 a13 a14]
//[a21 a22 a23 a24]
//[a31 a32 a33 a34]
//[0   0   0   1  ]
void matrixMultiplyRigid(double M1[], double M2[], double outM[], unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2, bool verbose=false) {
    // Validate matrix dimensions
    if (x1 != y2) {
        if (verbose) {
            std::cerr << "Not possible to multiply\n";
        }
        return;
    }

    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < x2; j++) {
            outM[j + x2 * i] = M1[4 * i] * M2[j] + M1[4 * i + 1] * M2[j + x2] + M1[4 * i + 2] * M2[j + 2 * x2] + M1[4 * i + 3] * M2[j + 3 * x2];
        }
    }

    for (unsigned int j = 0; j < x2; j++) {
        outM[j + 3 * x2] = M2[j + 3 * x2];
    }

    // Printing for verbose mode
    if (verbose) {
        printf("\nThe product of entered matrices is:\n");
        for (unsigned int i = 0; i < y1; i++) {
            for (unsigned int j = 0; j < x2; j++) {
                printf("%f\t", outM[j + x2 * i]);
            }
            printf("\n");
        }
    }
}






template<typename T,typename U>
void backprojectToVolumeRealSpace (T* inProjection, U* outI, unsigned long int nx, unsigned long int ny, unsigned long int nz, double Phi, double Theta, double Psi, double tx, double ty, T* mask3D=NULL, T* externalReferenceMap=NULL, T* locresMap=NULL, double maskThreshold = 0.9) {

    //Phi=0; Theta=0; Psi=0;

    double M[16];
    double C[4];
    double CT[4];
    long int nxy = (long int)nx * ny;
    long int nxyz = nxy*nz;
    transformMatrixRelion(M, Psi, Theta, Phi, 0, 0, 0);
    inverseMatrix(M);
    float * tmpI =new float[nxyz];
    float * newX =new float[nxy];
    float * newY =new float[nxy];

    for (long int ij=0; ij<nxy; ij++) {
        //outI[ij] = 0;
        tmpI[ij] = 0;
        newX[ij]=0;
        newY[ij]=0;
    }
    double sum=0;
    const double avgPixNumber=(nx+ny+nz)/3;
    double nz2 = ((double)nz) / 2.0;
    double ny2 = ((double)ny) / 2.0;
    double nx2 = ((double)nx) / 2.0;
    C[3] = 1;
    C[2]=nz/2.0;
    for ( long int jj=0;jj<(long int)ny;jj++){
              C[1]=jj-(ny)/2.0;
              long int nyjj = ny*jj;
              for ( long int ii=0;ii<(long int)nx;ii++){
                        long int ij=ii+nyjj;
                        C[0]=ii-(nx)/2.0;
                        matrixMultiplication(M,C,CT,4,4,1,4);
                        double X = (nx2+CT[0])+tx;
                        double Y = (ny2+CT[1])+ty;


                        /*
                        //bilinear interpolation
                        double X = (nx2+CT[0])+tx;
                        double Y = (ny2+CT[1])+ty;
                        double a = X-floor(X);
                        double b = Y-floor(Y);
                        double c = 1.0-a;
                        double d = 1.0-b;                
                        long int newX0 = (long int)floor(X);
                        long int newX1 = (long int)ceil(X);
                        long int newY0 = (long int)floor(Y);
                        long int newY1 = (long int)ceil(Y);
                        
                        if(newX0>=0 && newX0<(long int)nx && newY0>=0 && newY0<(long int)ny){
                          outI[newX0+newY0*nx]+=c*d*inI[ijk]; //
                          tmpI[newX0+newY0*nx]+=c*d; //

                        }
                        if(newX1>=0 && newX1<(long int)nx && newY0>=0 && newY0<(long int)ny){
                          outI[newX1+newY0*nx]+=a*d*inI[ijk];//
                          tmpI[newX1+newY0*nx]+=a*d;
                        }
                        if(newX0>=0 && newX0<(long int)nx && newY1>=0 && newY1<(long int)ny){
                          outI[newX0+newY1*nx]+=b*c*inI[ijk];//
                          tmpI[newX0+newY1*nx]+=b*c;
                        }
                        if(newX1>=0 && newX1<(long int)nx && newY1>=0 && newY1<(long int)ny){
                          outI[newX1+newY1*nx]+=a*b*inI[ijk];//
                          tmpI[newX1+newY1*nx]+=a*b;
                        }*/
              }
      }

/*
      for (int ij=0; ij<nxy;ij++){
        if (tmpI[ij]>0.00001)
          outI[ij]=outI[ij]*avgPixNumber/tmpI[ij];
      }
*/

/*
    for ( long int kk=0;kk<(long int)nz;kk++){
      long int kknxy = nxy*kk;
      for ( long int jj=0;jj<(long int)ny;jj++){
        long int jjnx = jj*nx;
        long int jjnx_kknxy = jjnx+kknxy;
        for ( long int ii=0;ii<(long int)nx;ii++){
          long int ii_jjnx = ii+jjnx;
          long int ii_jjnx_kknxy = ii+jjnx_kknxy;
          long int tmpX=floor(newX[ii_jjnx]);
          long int tmpY=floor(newY[ii_jjnx]);
          if (tmpX>=0&&tmpX<nx && tmpY>=0&&tmpY<ny){
            outI[ii_jjnx_kknxy]= inProjection[tmpX+tmpY*nx];
          }
        }
      }
    }

*/

    //writeMrcImage("counter.mrc", tmpI, nx,ny,nz);
    delete []tmpI; 
    delete []newX; 
    delete []newY; 


}

// 2D soft-edge of a mask in-place, using a separable triangular kernel
template<typename T>
void softEdge2D_slow(T* __restrict mask2D,
                unsigned long int nx,
                unsigned long int ny,
                int edgePixels)
{
    if (!mask2D || edgePixels <= 0) return;

    const long nxy = (long)nx * (long)ny;
    const int  r   = edgePixels;
    const int  klen = 2 * r + 1;

    // Build 1D triangular kernel of length 2*r+1, normalised
    std::vector<float> k(klen);
    float sum = 0.0f;
    for (int i = 0; i < klen; ++i) {
        const int x = i - r;
        float val = float(r + 1 - std::abs(x)) / float(r + 1);  // peak at centre
        k[i] = val;
        sum += val;
    }
    if (sum > 0.0f) {
        const float inv = 1.0f / sum;
        for (int i = 0; i < klen; ++i)
            k[i] *= inv;
    }

    std::vector<float> tmp(nxy, 0.0f);

    // First pass: convolve rows (x direction) into tmp
    for (unsigned long y = 0; y < ny; ++y) {
        const long row_off = (long)y * (long)nx;
        for (unsigned long x = 0; x < nx; ++x) {
            double acc = 0.0;
            for (int t = -r; t <= r; ++t) {
                long xx = (long)x + (long)t;
                if (xx < 0) xx = 0;
                if (xx >= (long)nx) xx = (long)nx - 1;
                acc += (double)mask2D[row_off + xx] * (double)k[t + r];
            }
            tmp[row_off + (long)x] = (float)acc;
        }
    }

    // Second pass: convolve columns (y direction), writing back into mask2D
    for (unsigned long x = 0; x < nx; ++x) {
        for (unsigned long y = 0; y < ny; ++y) {
            double acc = 0.0;
            for (int t = -r; t <= r; ++t) {
                long yy = (long)y + (long)t;
                if (yy < 0) yy = 0;
                if (yy >= (long)ny) yy = (long)ny - 1;
                acc += (double)tmp[(long)yy * (long)nx + (long)x] * (double)k[t + r];
            }
            mask2D[(long)y * (long)nx + (long)x] = (T)acc;
        }
    }

    // Normalise mask to [0,1] by dividing by its maximum
    float maxval = 0.0f;
    for (long i = 0; i < nxy; ++i)
        if (mask2D[i] > maxval) maxval = mask2D[i];

    if (maxval > 0.0f) {
        const float inv = 1.0f / maxval;
        for (long i = 0; i < nxy; ++i) {
            float v = mask2D[i] * inv;
            if (v < 0.0f) v = 0.0f;
            if (v > 1.0f) v = 1.0f;
            mask2D[i] = (T)v;
        }
    }
}
// 2D soft-edge of a mask in-place, using two fast 1D box blurs (separable)
// Approximates a triangular kernel; much faster than explicit O(r) convolution.
template<typename T>
void softEdge2D(T* __restrict mask2D,
                unsigned long int nx,
                unsigned long int ny,
                int edgePixels)
{
    if (!mask2D || edgePixels <= 0) return;

    const long nxy = (long)nx * (long)ny;
    const int  r   = edgePixels;

    // Work in float internally
    std::vector<float> buf(nxy);
    std::vector<float> tmp(nxy);
    for (long i = 0; i < nxy; ++i)
        buf[i] = (float)mask2D[i];

    // Precompute clamped window bounds for each x and y to avoid branches in inner loops
    std::vector<int> leftX(nx), rightX(nx);
    for (unsigned long x = 0; x < nx; ++x) {
        int xi = (int)x;
        int lx = xi - r;
        int rx = xi + r;
        if (lx < 0)       lx = 0;
        if (rx >= (int)nx) rx = (int)nx - 1;
        leftX[x]  = lx;
        rightX[x] = rx;
    }

    std::vector<int> topY(ny), bottomY(ny);
    for (unsigned long y = 0; y < ny; ++y) {
        int yi = (int)y;
        int ty = yi - r;
        int by = yi + r;
        if (ty < 0)       ty = 0;
        if (by >= (int)ny) by = (int)ny - 1;
        topY[y]    = ty;
        bottomY[y] = by;
    }

    // Prefix buffer large enough for max(nx, ny) + 1
    const unsigned long maxLen = (nx > ny ? nx : ny) + 1;
    std::vector<double> prefix(maxLen);

    // First pass: horizontal box blur (rows) -> tmp
    for (unsigned long y = 0; y < ny; ++y) {
        const long row_off = (long)y * (long)nx;

        // Build prefix sums for this row: prefix[i] = sum_{0..i-1} buf[row_off + j]
        prefix[0] = 0.0;
        for (unsigned long x = 0; x < nx; ++x)
            prefix[x + 1] = prefix[x] + (double)buf[row_off + (long)x];

        // Sliding window using prefix sums
        for (unsigned long x = 0; x < nx; ++x) {
            const int lx   = leftX[x];
            const int rx   = rightX[x];
            const int len  = rx - lx + 1;
            const double s = prefix[rx + 1] - prefix[lx];
            tmp[row_off + (long)x] = (float)(s / (double)len);
        }
    }

    // Second pass: vertical box blur (columns) -> buf
    for (unsigned long x = 0; x < nx; ++x) {
        // Build prefix sums down the column
        prefix[0] = 0.0;
        for (unsigned long y = 0; y < ny; ++y) {
            const long idx = (long)y * (long)nx + (long)x;
            prefix[y + 1] = prefix[y] + (double)tmp[idx];
        }

        for (unsigned long y = 0; y < ny; ++y) {
            const int ty   = topY[y];
            const int by   = bottomY[y];
            const int len  = by - ty + 1;
            const double s = prefix[by + 1] - prefix[ty];

            const long idx = (long)y * (long)nx + (long)x;
            buf[idx] = (float)(s / (double)len);
        }
    }

    // Normalise buf to [0,1] and write back into mask2D
    float maxval = 0.0f;
    for (long i = 0; i < nxy; ++i)
        if (buf[i] > maxval) maxval = buf[i];

    if (maxval > 0.0f) {
        const float inv = 1.0f / maxval;
        for (long i = 0; i < nxy; ++i) {
            float v = buf[i] * inv;
            if (v < 0.0f) v = 0.0f;
            if (v > 1.0f) v = 1.0f;
            mask2D[i] = (T)v;
        }
    } else {
        // All-zero mask: just copy back zeros
        for (long i = 0; i < nxy; ++i)
            mask2D[i] = (T)0;
    }
}


template<typename T, typename U>
void ProjectVolumeRealSpace_with2DMask(T* __restrict inI, U* __restrict outI,
                                       unsigned long int nx, unsigned long int ny, unsigned long int nz,
                                       double Phi, double Theta, double Psi, double tx, double ty,
                                       T* __restrict mask2D = NULL,
                                       int softEdgePixels = 0,
                                       double maskEps = 1e-2)
{
    const long nxy = (long)nx * (long)ny;
    std::vector<float> tmpI(nxy, 0.0f);
    for (long i = 0; i < nxy; ++i)
        outI[i] = 0;

    // Prepare soft mask once
    if (mask2D && softEdgePixels > 0)
        softEdge2D(mask2D, nx, ny, softEdgePixels);

    // Rotation matrix as in the original projector
    double M[16];
    transformMatrixRelion(M, /*Psi*/Psi, /*Theta*/Theta, /*Phi*/Phi, 0, 0, 0);
    inverseMatrix(M);
    const double a00 = M[0], a01 = M[1], a02 = M[2];
    const double a10 = M[4], a11 = M[5], a12 = M[6];

    const double nx2 = 0.5 * (double)nx;
    const double ny2 = 0.5 * (double)ny;
    const double nz2 = 0.5 * (double)nz;
    const double avgPixNumber = ((double)nx + (double)ny + (double)nz) / 3.0;

    const bool useMask = (mask2D != NULL);

    // Parts of X,Y that do not depend on ii
    const double baseXShift = nx2 + tx;
    const double baseYShift = ny2 + ty;

    for (long kk = 0; kk < (long)nz; ++kk) {
        const double z = (double)kk - nz2;
        const long sliceOffset = kk * nxy;

        for (long jj = 0; jj < (long)ny; ++jj) {
            const double y = (double)jj - ny2;
            const long rowOffset = jj * (long)nx;

            // Xc, Yc for x = 0 (ii = 0) and their contribution to X,Y
            const double baseXc = a01 * y + a02 * z;
            const double baseYc = a11 * y + a12 * z;

            // X, Y at ii = 0
            const double baseX = baseXShift + baseXc - a00 * nx2;
            const double baseY = baseYShift + baseYc - a10 * nx2;

            for (long ii = 0; ii < (long)nx; ++ii) {
                const long ij  = ii + rowOffset;
                const long ijk = ij + sliceOffset;

                // CT = inv(M) * [x,y,z,1], take XY
                // but x now folded into ii via a00/a10
                const double ii_d = (double)ii;
                const double X = baseX + a00 * ii_d;
                const double Y = baseY + a10 * ii_d;

                const double fx = std::floor(X);
                const double fy = std::floor(Y);
                const long   x0 = (long)fx, y0 = (long)fy;
                const long   x1 = x0 + 1,   y1 = y0 + 1;
                const double a  = X - fx,   b = Y - fy;
                const double c  = 1.0 - a,  d = 1.0 - b;

                if (useMask) {
                    double m00 = 0.0, m10 = 0.0, m01 = 0.0, m11 = 0.0;

                    if (x0>=0 && x0<(long)nx && y0>=0 && y0<(long)ny) {
                        long idx = x0 + y0*(long)nx;
                        m00 = (double)mask2D[idx];
                        if (m00 < maskEps) continue;
                    }
                    if (x1>=0 && x1<(long)nx && y0>=0 && y0<(long)ny) {
                        long idx = x1 + y0*(long)nx;
                        m10 = (double)mask2D[idx];
                        if (m10 < maskEps) continue;
                    }
                    if (x0>=0 && x0<(long)nx && y1>=0 && y1<(long)ny) {
                        long idx = x0 + y1*(long)nx;
                        m01 = (double)mask2D[idx];
                        if (m01 < maskEps) continue;
                    }
                    if (x1>=0 && x1<(long)nx && y1>=0 && y1<(long)ny) {
                        long idx = x1 + y1*(long)nx;
                        m11 = (double)mask2D[idx];
                        if (m11 < maskEps) continue;
                    }

                    // Now we really need the map value
                    const double v = (double)inI[ijk];

                    // Accumulate, weighting by both bilinear weight and mask value
                    if (x0>=0 && x0<(long)nx && y0>=0 && y0<(long)ny && m00 > maskEps) {
                        const long idx = x0 + y0*(long)nx;
                        const float w  = (float)(c*d);
                        const float wm = (float)(w * m00);
                        outI[idx] += (U)(wm * v);
                        tmpI[idx] += wm;
                    }
                    if (x1>=0 && x1<(long)nx && y0>=0 && y0<(long)ny && m10 > maskEps) {
                        const long idx = x1 + y0*(long)nx;
                        const float w  = (float)(a*d);
                        const float wm = (float)(w * m10);
                        outI[idx] += (U)(wm * v);
                        tmpI[idx] += wm;
                    }
                    if (x0>=0 && x0<(long)nx && y1>=0 && y1<(long)ny && m01 > maskEps) {
                        const long idx = x0 + y1*(long)nx;
                        const float w  = (float)(b*c);
                        const float wm = (float)(w * m01);
                        outI[idx] += (U)(wm * v);
                        tmpI[idx] += wm;
                    }
                    if (x1>=0 && x1<(long)nx && y1>=0 && y1<(long)ny && m11 > maskEps) {
                        const long idx = x1 + y1*(long)nx;
                        const float w  = (float)(a*b);
                        const float wm = (float)(w * m11);
                        outI[idx] += (U)(wm * v);
                        tmpI[idx] += wm;
                    }
                } else {
                    // Original unmasked path, same maths but cheaper X,Y
                    const double v = (double)inI[ijk];

                    if (x0>=0 && x0<(long)nx && y0>=0 && y0<(long)ny) {
                        const long idx = x0 + y0*(long)nx;
                        const float w = (float)(c*d);
                        outI[idx] += (U)(w * v);
                        tmpI[idx] += w;
                    }
                    if (x1>=0 && x1<(long)nx && y0>=0 && y0<(long)ny) {
                        const long idx = x1 + y0*(long)nx;
                        const float w = (float)(a*d);
                        outI[idx] += (U)(w * v);
                        tmpI[idx] += w;
                    }
                    if (x0>=0 && x0<(long)nx && y1>=0 && y1<(long)ny) {
                        const long idx = x0 + y1*(long)nx;
                        const float w = (float)(b*c);
                        outI[idx] += (U)(w * v);
                        tmpI[idx] += w;
                    }
                    if (x1>=0 && x1<(long)nx && y1>=0 && y1<(long)ny) {
                        const long idx = x1 + y1*(long)nx;
                        const float w = (float)(a*b);
                        outI[idx] += (U)(w * v);
                        tmpI[idx] += w;
                    }
                }
            }
        }
    }

    // Final normalisation (same idea as ProjectVolumeRealSpace) + apply 2D mask
    for (long i = 0; i < nxy; ++i) {
        const float w = tmpI[i];

        if (w > 1e-5f) {
            double val = (double)outI[i] * (avgPixNumber / (double)w);

            if (useMask)
                val *= (double)mask2D[i];

            outI[i] = (U)val;
        } else {
            outI[i] = (U)0;
        }
    }
}


template<typename T,typename U>
void ProjectVolumeRealSpace(T* __restrict inI, U* __restrict outI,
                            unsigned long int nx, unsigned long int ny, unsigned long int nz,
                            double Phi, double Theta, double Psi, double tx, double ty,
                            T* __restrict mask3D = NULL, double maskThreshold = 0.9)
{
    const long nxy = (long)nx * (long)ny;
    std::vector<float> tmpI(nxy, 0.0f);
    for (long i=0;i<nxy;++i) outI[i] = 0;

    // === Build exactly the same matrix as legacy path (angles swizzled) ===
    double M[16];
    transformMatrixRelion(M, /*Psi*/Psi, /*Theta*/Theta, /*Phi*/Phi, 0,0,0);
    inverseMatrix(M);
    const double a00 = M[0], a01 = M[1], a02 = M[2];  // row 0
    const double a10 = M[4], a11 = M[5], a12 = M[6];  // row 1

    const double nx2 = 0.5 * (double)nx;
    const double ny2 = 0.5 * (double)ny;
    const double nz2 = 0.5 * (double)nz;
    const double avgPixNumber = ((double)nx + (double)ny + (double)nz) / 3.0;

    for (long kk = 0; kk < (long)nz; ++kk) {
        for (long jj = 0; jj < (long)ny; ++jj) {
            for (long ii = 0; ii < (long)nx; ++ii) {
                const long ijk = ii + jj*(long)nx + kk*nxy;
                if (mask3D && !((double)mask3D[ijk] > maskThreshold)) continue;

                const double x = (double)ii - nx2;
                const double y = (double)jj - ny2;
                const double z = (double)kk - nz2;

                // EXACT legacy projection: CT = inv(M) * [x,y,z,1], take XY
                const double Xc = a00*x + a01*y + a02*z;
                const double Yc = a10*x + a11*y + a12*z;
                // ---------------------------------------------------------------------

                const double X  = nx2 + Xc + tx;
                const double Y  = ny2 + Yc + ty;

                const double fx = std::floor(X);
                const double fy = std::floor(Y);
                const long   x0 = (long)fx, y0 = (long)fy;
                const long   x1 = x0 + 1,    y1 = y0 + 1;
                const double a  = X - fx,    b = Y - fy;
                const double c  = 1.0 - a,   d = 1.0 - b;

                const double v = (double)inI[ijk];

                if (x0>=0 && x0<(long)nx && y0>=0 && y0<(long)ny) {
                    const long idx = x0 + y0*(long)nx; const float w = (float)(c*d);
                    outI[idx] += (U)(w * v);
                    tmpI[idx] += w;
                }
                if (x1>=0 && x1<(long)nx && y0>=0 && y0<(long)ny) {
                    const long idx = x1 + y0*(long)nx; const float w = (float)(a*d);
                    outI[idx] += (U)(w * v);
                    tmpI[idx] += w;
                }
                if (x0>=0 && x0<(long)nx && y1>=0 && y1<(long)ny) {
                    const long idx = x0 + y1*(long)nx; const float w = (float)(b*c);
                    outI[idx] += (U)(w * v);
                    tmpI[idx] += w;
                }
                if (x1>=0 && x1<(long)nx && y1>=0 && y1<(long)ny) {
                    const long idx = x1 + y1*(long)nx; const float w = (float)(a*b);
                    outI[idx] += (U)(w * v);
                    tmpI[idx] += w;
                }
            }
        }
    }

    for (long i=0;i<nxy;++i) {
        const float w = tmpI[i];
        if (w > 1e-5f) outI[i] = (U)((double)outI[i] * (avgPixNumber / (double)w));
        else           outI[i] = (U)0;
    }
}




template<typename T,typename U>
void ProjectVolumeRealSpace_ORIGINALSLOW(T* inI, U* outI, unsigned long int nx, unsigned long int ny, unsigned long int nz, 
                            double Phi, double Theta, double Psi, double tx, double ty, 
                            T* mask3D=NULL, double maskThreshold = 0.9) {


//Phi=0; Theta=0; Psi=0;

    double M[16];
    double C[4];
    double CT[4];
    long int nxy = (long int)nx * ny;
    long int nxyz = nxy*nz;
    transformMatrixRelion(M, Psi, Theta, Phi, 0, 0, 0);
    inverseMatrix(M);
    float * tmpI =new float[nxy];

    for (long int ij=0; ij<nxy; ij++) {
        outI[ij] = 0;
        tmpI[ij] = 0;
    }
    double sum=0;
    const double avgPixNumber=(nx+ny+nz)/3;


    double nz2 = ((double)nz) / 2.0;
    double ny2 = ((double)ny) / 2.0;
    double nx2 = ((double)nx) / 2.0;
    C[3] = 1;


          for ( long int kk=0, ijk=0;kk<(long int)nz;kk++){
            C[2]=kk-(nz)/2.0;
            for ( long int jj=0;jj<(long int)ny;jj++){
              C[1]=jj-(ny)/2.0;
              for ( long int ii=0;ii<(long int)nx;ii++,ijk++){
                if (!mask3D || (mask3D && mask3D[ijk] > maskThreshold)) {
                        C[0]=ii-(nx)/2.0;
                        matrixMultiplication(M,C,CT,4,4,1,4);
                        
                        //bilinear interpolation
                        double X = (nx2+CT[0])+tx;
                        double Y = (ny2+CT[1])+ty;
                        double a = X-floor(X);
                        double b = Y-floor(Y);
                        double c = 1.0-a;
                        double d = 1.0-b;                
                        long int newX0 = (long int)floor(X);
                        long int newX1 = (long int)ceil(X);
                        long int newY0 = (long int)floor(Y);
                        long int newY1 = (long int)ceil(Y);
                        
                        if(newX0>=0 && newX0<(long int)nx && newY0>=0 && newY0<(long int)ny){
                          outI[newX0+newY0*nx]+=c*d*inI[ijk]; //
                          tmpI[newX0+newY0*nx]+=c*d; //

                        }
                        if(newX1>=0 && newX1<(long int)nx && newY0>=0 && newY0<(long int)ny){
                          outI[newX1+newY0*nx]+=a*d*inI[ijk];//
                          tmpI[newX1+newY0*nx]+=a*d;
                        }
                        if(newX0>=0 && newX0<(long int)nx && newY1>=0 && newY1<(long int)ny){
                          outI[newX0+newY1*nx]+=b*c*inI[ijk];//
                          tmpI[newX0+newY1*nx]+=b*c;
                        }
                        if(newX1>=0 && newX1<(long int)nx && newY1>=0 && newY1<(long int)ny){
                          outI[newX1+newY1*nx]+=a*b*inI[ijk];//
                          tmpI[newX1+newY1*nx]+=a*b;
                        }
               }
              }
            }
          }
          for (int ij=0; ij<nxy;ij++){
            if (tmpI[ij]>0.00001)
              outI[ij]=outI[ij]*avgPixNumber/tmpI[ij];
          }

    //writeMrcImage("counter.mrc", tmpI, nx,ny,nz);
    delete []tmpI; 

}








template<typename T,typename U>
void ProjectVolumeRealSpace_original (T* inI, U* outI, unsigned long int nx, unsigned long int ny, unsigned long int nz, double Phi, double Theta, double Psi, double tx, double ty, T* mask3D=NULL, double maskThreshold = 0.9){

  double M[16];
  double C[4];
  double CT[4];
  //long int nxyz=(long int)nx*ny*nz;
  long int nxy=(long int)nx*ny;
  transformMatrixRelion(M, Psi, Theta, Phi, tx, ty, 0);
  inverseMatrix(M);


  //initialize
  for ( long int ij=0; ij<nxy; ij++){
   outI[ij]=0;
  }

  // ********
  //  3D MASK
  if (mask3D){
          double nz2=  (double)nz/2.0;
          double ny2=  (double)ny/2.0;
          double nx2=  (double)nx/2.0;
          C[3]=1;
          for ( long int kk=0, ijk=0;kk<(long int)nz;kk++){
            C[2]=kk-(nz)/2.0;
            for ( long int jj=0;jj<(long int)ny;jj++){
              C[1]=jj-(ny)/2.0;
              for ( long int ii=0;ii<(long int)nx;ii++,ijk++){
                if ( mask3D[ijk] > maskThreshold ){
                        C[0]=ii-(nx)/2.0;
                        matrixMultiplication(M,C,CT,4,4,1,4);
                        
                        //bilinear interpolation
                        double X = (nx2+CT[0])+tx;
                        double Y = (ny2+CT[1])+ty;
                        double a = X-floor(X);
                        double b = Y-floor(Y);
                        double c = 1.0-a;
                        double d = 1.0-b;                
                        long int newX0 = (long int)floor(X);
                        long int newX1 = (long int)ceil(X);
                        long int newY0 = (long int)floor(Y);
                        long int newY1 = (long int)ceil(Y);
                        
                        if(newX0>=0 && newX0<(long int)nx && newY0>=0 && newY0<(long int)ny){
                          outI[newX0+newY0*nx]+=c*d*inI[ijk]; //
                        }
                        if(newX1>=0 && newX1<(long int)nx && newY0>=0 && newY0<(long int)ny){
                          outI[newX1+newY0*nx]+=a*d*inI[ijk];//
                        }
                        if(newX0>=0 && newX0<(long int)nx && newY1>=0 && newY1<(long int)ny){
                          outI[newX0+newY1*nx]+=b*c*inI[ijk];//
                        }
                        if(newX1>=0 && newX1<(long int)nx && newY1>=0 && newY1<(long int)ny){
                          outI[newX1+newY1*nx]+=a*b*inI[ijk];//
                        }
               }
              }
            }
          }
  }
  // *******************
  // NO MASK
  else if(!mask3D){
          double nz2=  (double)nz/2.0;
          double ny2=  (double)ny/2.0;
          double nx2=  (double)nx/2.0;
          C[3]=1;
          for ( long int kk=0, ijk=0;kk<(long int)nz;kk++){
            C[2]=kk-(nz)/2.0;
            for ( long int jj=0;jj<(long int)ny;jj++){
              C[1]=jj-(ny)/2.0;
              for ( long int ii=0;ii<(long int)nx;ii++,ijk++){
                        C[0]=ii-(nx)/2.0;
                        matrixMultiplication(M,C,CT,4,4,1,4);
                        
                        //bilinear interpolation
                        double X = (nx2+CT[0])+tx;
                        double Y = (ny2+CT[1])+ty;
                        double a = X-floor(X);
                        double b = Y-floor(Y);
                        double c = 1.0-a;
                        double d = 1.0-b;                
                        long int newX0 = (long int)floor(X);
                        long int newX1 = (long int)ceil(X);
                        long int newY0 = (long int)floor(Y);
                        long int newY1 = (long int)ceil(Y);
                        
                        if(newX0>=0 && newX0<(long int)nx && newY0>=0 && newY0<(long int)ny){
                          outI[newX0+newY0*nx]+=c*d*inI[ijk]; //
                        }
                        if(newX1>=0 && newX1<(long int)nx && newY0>=0 && newY0<(long int)ny){
                          outI[newX1+newY0*nx]+=a*d*inI[ijk];//
                        }
                        if(newX0>=0 && newX0<(long int)nx && newY1>=0 && newY1<(long int)ny){
                          outI[newX0+newY1*nx]+=b*c*inI[ijk];//
                        }
                        if(newX1>=0 && newX1<(long int)nx && newY1>=0 && newY1<(long int)ny){
                          outI[newX1+newY1*nx]+=a*b*inI[ijk];//
                        }
              }
            }
          }
  }
}





template<typename T, typename U>
void ProjectMaskRealSpace(
    T* __restrict inI, U* __restrict outI,
    unsigned long nx, unsigned long ny, unsigned long nz,
    double Phi, double Theta, double Psi, double tx, double ty,
    double maskThreshold = 1.0)
{
    const long NX = (long)nx;
    const long NY = (long)ny;
    const long NZ = (long)nz;
    const long NXY = NX * NY;

    // Zero output
    std::fill(outI, outI + NXY, (U)0);

    // Build transform matrix once
    double M[16];
    transformMatrixRelion(M, Psi, Theta, Phi, 0, 0, 0);
    inverseMatrix(M);

    // First two rows of inverse
    const double a00 = M[0], a01 = M[1], a02 = M[2];
    const double a10 = M[4], a11 = M[5], a12 = M[6];

    const double nx2 = 0.5 * (double)NX;
    const double ny2 = 0.5 * (double)NY;
    const double nz2 = 0.5 * (double)NZ;

    // Precompute offsets for incremental update
    const double Xc_i = a00;              // ∂Xc/∂i
    const double Yc_i = a10;
    const double Xc_j = a01;              // ∂Xc/∂j
    const double Yc_j = a11;
    const double Xc_k = a02;              // ∂Xc/∂k
    const double Yc_k = a12;

    std::vector<float> cov(NXY, 0.0f);
    T* __restrict mask3D = inI;

    for (long kk = 0; kk < NZ; ++kk) {
        const double z = (double)kk - nz2;
        const double Xc_z = Xc_k * z;
        const double Yc_z = Yc_k * z;

        for (long jj = 0; jj < NY; ++jj) {
            const double y = (double)jj - ny2;

            // Xc_base/Yc_base are the coordinates at i=0 for this (j,k)
            double Xc_base = a00*(-nx2) + a01*y + Xc_z;
            double Yc_base = a10*(-nx2) + a11*y + Yc_z;

            // Shifted to pixel coordinates and translation applied
            double X = nx2 + Xc_base + tx;
            double Y = ny2 + Yc_base + ty;

            const long base_jk = jj * NX + kk * NX * NY;
            const T* __restrict slice = mask3D + base_jk;

            for (long ii = 0; ii < NX; ++ii) {
                if (slice[ii] > (T)maskThreshold) {
                    const long x0 = (long)std::floor(X);
                    const long y0 = (long)std::floor(Y);

                    if ((unsigned)x0 < (unsigned)NX && (unsigned)y0 < (unsigned)NY) {
                        const long idx = x0 + y0 * NX;
                        cov[idx] += 1.0f;
                    }

                    const long x1 = x0 + 1;
                    const long y1 = y0 + 1;

                    if ((unsigned)x1 < (unsigned)NX && (unsigned)y0 < (unsigned)NY) {
                        const long idx = x1 + y0 * NX;
                        cov[idx] += 1.0f;
                    }
                    if ((unsigned)x0 < (unsigned)NX && (unsigned)y1 < (unsigned)NY) {
                        const long idx = x0 + y1 * NX;
                        cov[idx] += 1.0f;
                    }
                    if ((unsigned)x1 < (unsigned)NX && (unsigned)y1 < (unsigned)NY) {
                        const long idx = x1 + y1 * NX;
                        cov[idx] += 1.0f;
                    }
                }

                // Incremental updates instead of recomputing from scratch
                X += Xc_i;
                Y += Yc_i;
            }
        }
    }

    // Threshold coverage → binary mask
    for (long i = 0; i < NXY; ++i)
        outI[i] = (U)((cov[i] > 0.0f) ? 1 : 0);
}




template<typename T,typename U>
void ProjectMaskRealSpace_SLOWER(T* inI, U* outI, unsigned long int nx, unsigned long int ny, unsigned long int nz, double Phi, double Theta, double Psi, double tx, double ty, double maskThreshold = 1){

  double M[16];
  double C[4];
  double CT[4];
  //long int nxyz=(long int)nx*ny*nz;
  long int nxy=(long int)nx*ny;
  transformMatrixRelion(M, Psi, Theta, Phi, tx, ty, 0);
  inverseMatrix(M);

  //initialize
  for (long int ij=0; ij<nxy; ij++){
   outI[ij]=0;
  }
  T* mask3D=inI;

  // ********
  //  3D MASK
          double nz2=  (double)nz/2.0;
          double ny2=  (double)ny/2.0;
          double nx2=  (double)nx/2.0;
          C[3]=1;
          for (unsigned long int kk=0, ijk=0;kk<nz;kk++){
            C[2]=kk-nz2;
            for ( unsigned long int jj=0;jj<ny;jj++){
              C[1]=jj-ny2;
              for (unsigned long int ii=0;ii<nx;ii++,ijk++){
                if (mask3D[ijk] > maskThreshold ){
                        C[0]=ii-nx2;
                        matrixMultiplication(M,C,CT,4,4,1,4);
                        
                        //bilinear interpolation
                        double X = (nx2+CT[0])+tx;
                        double Y = (ny2+CT[1])+ty;
                        long int newX0 = (long int)floor(X);
                        long int newX1 = (long int)ceil(X);
                        long int newY0 = (long int)floor(Y);
                        long int newY1 = (long int)ceil(Y);
                        
                        if(newX0>=0 && newX0<(long int)nx && newY0>=0 && newY0<(long int)ny){
                          outI[newX0+newY0*nx]=1; //
                        }
                        if(newX1>=0 && newX1<(long int)nx && newY0>=0 && newY0<(long int)ny){
                          outI[newX1+newY0*nx]=1;//
                        }
                        if(newX0>=0 && newX0<(long int)nx && newY1>=0 && newY1<(long int)ny){
                          outI[newX0+newY1*nx]=1;//
                        }
                        if(newX1>=0 && newX1<(long int)nx && newY1>=0 && newY1<(long int)ny){
                          outI[newX1+newY1*nx]=1;//
                        }
               }
              }
            }
          }
}





#endif
