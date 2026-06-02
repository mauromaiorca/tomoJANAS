/*
 * File: mrcIO.h
 * (C) 2022 Mauro Maiorca 
 */

#ifndef __MRC_IO__H___
#define __MRC_IO__H___



#include <cmath>
#include <iostream>
#include <iomanip>      // std::setprecision
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstdlib>
#include <fstream>

#include <arpa/inet.h> //used for htonl (big/little endian check)


enum  DataType {UChar, SChar, Short, Int, UShort, UInt, Float,  RFloat, Double, Bool};


//https://www.sciencedirect.com/science/article/pii/S104784771500074X
struct MRCHeader
{             // file header for MRC data
    int nx;              //  0   0       image size
    int ny;              //  1   4
    int nz;              //  2   8
    int mode;            //  3           0=char,1=short,2=float
    int nxstart;         //  4           unit cell offset
    int nystart;         //  5
    int nzstart;         //  6
    int mx;              //  7           Number of intervals along X of the “unit cell” in voxels
    int my;              //  8		 Number of intervals along Y of the “unit cell” in voxels
    int mz;              //  9           Number of intervals along Z of the “unit cell” in voxels
    float a;             // 10   40      cell dimensions in A
    float b;             // 11
    float c;             // 12
    float alpha;         // 13           cell angles in degrees
    float beta;          // 14
    float gamma;         // 15
    int mapc;            // 16           column axis
    int mapr;            // 17           row axis
    int maps;            // 18           section axis
    float dmin;          // 19           minimum density value
    float dmax;          // 20   80      maximum density value
    float dmean;         // 21           average density value
    int ispg;            // 22           space group number
    int nsymbt;          // 23           bytes used for sym. ops. table
    float extra[25];     // 24           user-defined info
    float xOrigin;       // 49           phase origin in pixels
    float yOrigin;       // 50
    float zOrigin;       // 51
    char map[4];         // 52       identifier for map file ("MAP ")
    char machst[4];      // 53           machine stamp
    float drms;          // 54       RMS deviation
    int nlabl;           // 55           number of labels used
    char labels[800];    // 56-255       10 80-character labels
    MRCHeader(){
        nx=0;
        ny=0;
        nz=0;
        mode=2;
        nxstart=0;
        nystart=0;
        nzstart=0;
        mx=0;
        my=0;
        mz=0;
        a=0;
        b=0;
        c=0;
        alpha=90;
        beta=90;
        gamma=90;
        mapc=1;
        mapr=2;
        maps=3;
        dmin=0;
        dmax=0;
        dmean=0;
        ispg=0;
        nsymbt=0;
        xOrigin=0;
        yOrigin=0;
        zOrigin=0;
        drms=0;
        nlabl=1;
        }
    MRCHeader(int __nx, int __ny, int __nz):nx(__nx),ny(__ny),nz(__nz){
        mode=2;
        nxstart=0;
        nystart=0;
        nzstart=0;
        mx=__nx;
        my=__ny;
        mz=__nz;
        a=__nx;
        b=__ny;
        c=__nz;
        alpha=90;
        beta=90;
        gamma=90;
        mapc=1;
        mapr=2;
        maps=3;
        dmin=0;
        dmax=0;
        dmean=0;
        ispg=0;
        nsymbt=0;
        xOrigin=0;
        yOrigin=0;
        zOrigin=0;
        drms=0;
        nlabl=1;
        }

} ;


void printHeaderInfo(MRCHeader & header){
    std::cerr<<"nx="<<header.nx<<"\n";              //  0   0       image size
    std::cerr<<"ny="<< header.ny<<"\n";              //  1   4
    std::cerr<<"nz="<< header.nz<<"\n";              //  2   8
    std::cerr<<"mode="<< header.mode<<"\n";            //  3           0=char,1=short,2=float
    std::cerr<<"nxstart="<< header.nxstart<<"\n";         //  4           unit cell offset
    std::cerr<<"nystart="<< header.nystart<<"\n";         //  5
    std::cerr<<"nzstart="<< header.nzstart<<"\n";         //  6
    std::cerr<<"mx="<< header.mx<<"\n";              //  7           Number of intervals along X of the “unit cell” in voxels
    std::cerr<<"my="<< header.my<<"\n";              //  8		 Number of intervals along Y of the “unit cell” in voxels
    std::cerr<<"mz="<< header.mz<<"\n";              //  9           Number of intervals along Z of the “unit cell” in voxels
    std::cerr<<"a="<< header.a<<"\n";             // 10   40      cell dimensions in A
    std::cerr<<"b="<< header.b<<"\n";             // 11
    std::cerr<<"c="<< header.c<<"\n";             // 12
    std::cerr<<"alpha="<< header.alpha<<"\n";         // 13           cell angles in degrees
    std::cerr<<"beta="<< header.beta<<"\n";          // 14
    std::cerr<<"gamma="<< header.gamma<<"\n";         // 15
    std::cerr<<"mapc="<< header.mapc<<"\n";            // 16           column axis
    std::cerr<<"mapr="<< header.mapr<<"\n";            // 17           row axis
    std::cerr<<"maps="<< header.maps<<"\n";            // 18           section axis
    std::cerr<<"dmin="<< header.dmin<<"\n";          // 19           minimum density value
    std::cerr<<"dmax="<< header.dmax<<"\n";          // 20   80      maximum density value
    std::cerr<<"dmean="<< header.dmean<<"\n";         // 21           average density value
    std::cerr<<"ispg="<< header.ispg<<"\n";            // 22           space group number
    std::cerr<<"nsymbt="<< header.nsymbt<<"\n";          // 23           bytes used for sym. ops. table
    std::cerr<<"xOrigin="<< header.xOrigin<<"\n";       // 49           phase origin in pixels
    std::cerr<<"yOrigin="<< header.yOrigin<<"\n";       // 50
    std::cerr<<"zOrigin="<< header.zOrigin<<"\n";       // 51
    std::cerr<<"drms="<< header.drms<<"\n";          // 54       RMS deviation
    std::cerr<<"nlabl="<< header.nlabl<<"\n";           // 55           number of labels used
}


bool check_file_exist(const char *fileName){
    std::ifstream infile(fileName);
    return infile.good();
}

// **********************
// **********************
//  I/O
//
// copy headear
//
// **********************
// **********************
void copyHeaderMrcImage(MRCHeader & header, const char * filenameMRCoutput){

  //write output header
  FILE *output_fp = NULL;
	output_fp=fopen(filenameMRCoutput,"rb+");
	if(output_fp==NULL) return;
	fseek(output_fp, 0, SEEK_SET);
	fwrite(&header,1,1024,output_fp);
  fclose(output_fp);

}


void headerToMrcImage(const char * filenameMRCinput, const char * filenameMRCoutput){

  MRCHeader header;

  //read input header
  FILE *input_fp = NULL;
	input_fp=fopen(filenameMRCinput,"r");
	if(input_fp==NULL) return;
	rewind(input_fp);
	if(fread(&header,1,1024,input_fp)<1024) return;
  fclose(input_fp);

  //write output header
  FILE *output_fp = NULL;
	output_fp=fopen(filenameMRCoutput,"rb+");
	if(output_fp==NULL) return;
	fseek(output_fp, 0, SEEK_SET);
	fwrite(&header,1,1024,output_fp);
  fclose(output_fp);

}


// **********************
// **********************
//  I/O
//
// read header
//
// **********************
// **********************
void readHeaderMrc(const char * filenameMRCinput, MRCHeader & header){

  //read input header
  FILE *input_fp = NULL;
	input_fp=fopen(filenameMRCinput,"r");
	if(input_fp==NULL) return;
	rewind(input_fp);
	if(fread(&header,1,1024,input_fp)<1024) return;
  fclose(input_fp);

}



// **********************
// **********************
//  I/O
//
// READ A VECTOR on a image
//
// **********************
// **********************
template<typename T>
int readMrcImage(const char * filenameMRC, T * I, MRCHeader & header){
  //std::cerr<<"MRC read image\n";

  //constant
  const int headerSizeMRC = 1024;
  const int _MAX_UNSIGNED_SHORT_ = 65535;

  //header
  FILE *m_fp = NULL;
  m_fp=fopen(filenameMRC,"r+");
  if(m_fp==NULL) return 1;

  //read file header
  rewind(m_fp);
  if(fread(&header,1,headerSizeMRC,m_fp)<headerSizeMRC)
	return 1;


  //READ THE FILE
  int         i;
  if ( ( abs( header.mode ) > _MAX_UNSIGNED_SHORT_ ) || ( abs(header.nx) > _MAX_UNSIGNED_SHORT_ ) ){

        std::cerr<<"Warning: Swapping header byte order for 4-byte types\n";
        int     extent = headerSizeMRC - 800; // exclude labels from swapping
        char _tmp;
	for ( i=0; i<extent; i+=4 ){
	    //swaps bytes
	    char* v = (char *) (&(header)+i);
	    for ( int iTmp=0; iTmp<2; iTmp++ )
	    {
		_tmp = v[iTmp];
		v[iTmp] = v[3-iTmp];
		v[3-iTmp] = _tmp;
	    }
	}
    }



    const unsigned long int nx =(int) header.nx;
    const unsigned long int ny =(int) header.ny;
    const unsigned long int nz =(int) header.nz;
    const unsigned long int nxyz = nx * ny * nz;

    DataType datatype;
    size_t datatypesize;
    unsigned long int blockMemory=268435456; //256 Mb
    void * buf;


    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
	//std::cerr<<"File mode 0: UChar\n";
	datatype=UChar;
	datatypesize = 1;
        buf = new unsigned char [blockMemory/datatypesize];
    }else if (header.mode==1){ //image : 16-bit halfwords
	//std::cerr<<"File mode 1: Short\n";
	datatype=Short;
	datatypesize = 2;
	buf = new short [blockMemory/datatypesize];
    }else if (header.mode==2){ //image : 32-bit reals
	//std::cerr<<"File mode 2: Float\n";
	datatype=Float;
        datatypesize = 4;
	buf = new float [blockMemory/datatypesize];
    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
	//std::cerr<<"File mode 5: Char\n";
	datatype=SChar;
	datatypesize = 1;
	buf = new char [blockMemory/datatypesize];
    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
	//std::cerr<<"File mode 6: UShort\n";
	datatype=UShort;
	datatypesize = 2;
	buf = new unsigned short [blockMemory/datatypesize];
    }else{
	//std::cerr<<"Unsupported type\n";
	datatype=UChar;
	datatypesize = 1;
	buf = new unsigned char [blockMemory/datatypesize];
    }

    unsigned long int blockSize=(unsigned long int)blockMemory/(double)datatypesize; 

    //std::cerr<<"read header\n";
    int error_fseek = fseek( m_fp, headerSizeMRC, SEEK_SET );
    if (error_fseek != 0)
           return -1;
    
    //std::cerr<<"read data\n";
    for ( unsigned long int tt=0; tt<nxyz; tt+=blockSize ){
	unsigned long int endRead=tt+blockSize;
	if ( endRead >= nxyz ){
		endRead=nxyz;
	}
	unsigned long int toRead=endRead-tt;

	// IT MAY BE SLIGHTLY MORE EFFICIENT TO PUT IT IN DIRECTLY
	// BUT NEED TO MANAGE POLYMORPHISM in EXECUTION
	/*if (sizeof(T)==datatypesize ){
		std::cerr<<"same type t="<<tt<<"->"<<endRead<<"   (max="<<nxyz<<")\n";
		size_t result = fread( &I[tt], datatypesize, toRead, m_fp );
		if (result < 0 ){
			std::cerr<<"Error: wrong read\n";
        	        return -1;
		}
	}*/	
		size_t result = fread( buf, datatypesize, toRead, m_fp );
		if (result < 0 ){
			std::cerr<<"Error: wrong read\n";
        	        return -1;
		}
	
		//std::cerr<<" block ="<<tt<<"->"<<endRead<<"   (max="<<nxyz<<")\n";

		    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
			for (unsigned long int ii=0;ii<toRead;ii++)
				I[tt+ii]=(T)(((unsigned char*) buf)[ii]);
		    }else if (header.mode==1){ //image : 16-bit halfwords
			for (unsigned long int ii=0;ii<toRead;ii++)
				I[tt+ii]=(T)(((short*) buf)[ii]);
		    }else if (header.mode==2){ //image : 32-bit reals
			for (unsigned long int ii=0;ii<toRead;ii++)
				I[tt+ii]=(T)(((float*) buf)[ii]);
		    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
			for (unsigned long int ii=0;ii<toRead;ii++)
				I[tt+ii]=(T)(((char*) buf)[ii]);
		    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
			for (unsigned long int ii=0;ii<toRead;ii++)			
				I[tt+ii]=(T)(((unsigned short*) buf)[ii]);
		    }else{
			//std::cerr<<"Unsupported type\n";
			for (unsigned long int ii=0;ii<toRead;ii++)
				I[tt+ii]=(T)(((unsigned char*) buf)[ii]);
		    }
	}

    
    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
        delete [] (unsigned char *) buf;
    }else if (header.mode==1){ //image : 16-bit halfwords
	delete [] (short *) buf;
    }else if (header.mode==2){ //image : 32-bit reals
	delete [](float *) buf;
    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
	delete [](char *) buf;
    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
	delete [](unsigned short *) buf;
    }else{
	delete [] (unsigned char *) buf;
    }

  //close the file
  fclose(m_fp);
  return 0;
}


template<typename T>
void readMrcImage(const char * filenameMRC, T * I, const unsigned long int nx, const unsigned long int ny, const unsigned long int nz, double APix=1.0, DataType datatype= Float){

   MRCHeader header;
    if (datatype==UChar){
	header.mode=0;
    }else if (datatype==Short){
	header.mode=1;
    }else if (datatype==Float){
	header.mode=2;
    }else if (datatype==SChar){
	header.mode=5;
    }else if (datatype==UShort){
	header.mode=6;
    }else{
	//std::cerr<<"Unsupported type\n";
	datatype=UChar;
    }

   header.nx=nx;
   header.ny=ny;
   header.nz=nz;
   header.mx=nx;
   header.my=ny;
   header.mz=nz;
   header.a=nx*APix;
   header.b=ny*APix;
   header.c=nz*APix;


    // fix some portions of the header
    header.mapc = 1;
    header.mapr = 2;
    header.maps = 3;
    header.alpha = (float)90.;
    header.beta = (float)90.;
    header.gamma = (float)90.;
    header.xOrigin = (float)0.;
    header.yOrigin = (float)0.;
    header.zOrigin = (float)0.;
    header.nxstart = (int)0;
    header.nystart = (int)0;
    header.nzstart = (int)0;
    header.ispg = (int)0;
    header.nsymbt = (int)0;
    header.map[0]='M';
    header.map[1]='A';
    header.map[2]='P';
    header.map[3]=' ';
    header.machst[0]='M';
    header.machst[1]='A';
    header.machst[2]='U';
    header.machst[3]=' ';
    header.nlabl=0;
   readMrcImage(filenameMRC, I, header);

}



// **********************
// **********************
//  I/O
//
// READ A VECTOR on a image
// at a certain slice
//
// **********************
// **********************
template<typename T>
int readMrcSlice(const char * filenameMRC, T * I, MRCHeader & header, long int sliceNumber){
  //std::cerr<<"MRC read image\n";
  /*headerOut=header;
  headerOut.nz=1;
  headerOut.c=1;
  */
  if (sliceNumber>header.nz){
    sliceNumber=header.nz-1;
  }
  if (sliceNumber<0){
   sliceNumber=0;
  }
  
  //constant
  const int headerSizeMRC = 1024;
  const int _MAX_UNSIGNED_SHORT_ = 65535;

  //header
  FILE *m_fp = NULL;
  m_fp=fopen(filenameMRC,"r+");
  if(m_fp==NULL) return 1;

  //read file header
  rewind(m_fp);
  if(fread(&header,1,headerSizeMRC,m_fp)<headerSizeMRC)
	return 1;


  //READ THE FILE
  int         i;
  if ( ( abs( header.mode ) > _MAX_UNSIGNED_SHORT_ ) || ( abs(header.nx) > _MAX_UNSIGNED_SHORT_ ) ){

        std::cerr<<"Warning: Swapping header byte order for 4-byte types\n";
        int     extent = headerSizeMRC - 800; // exclude labels from swapping
        char _tmp;
	for ( i=0; i<extent; i+=4 ){
	    //swaps bytes
	    char* v = (char *) (&(header)+i);
	    for ( int iTmp=0; iTmp<2; iTmp++ )
	    {
		_tmp = v[iTmp];
		v[iTmp] = v[3-iTmp];
		v[3-iTmp] = _tmp;
	    }
	}
    }



    const unsigned long int nx =(int) header.nx;
    const unsigned long int ny =(int) header.ny;
    //const unsigned long int nz =(int) header.nz;
    const unsigned long int nxy = nx * ny;
    //const unsigned long int nxyz = nx * ny * nz;

    DataType datatype;
    size_t datatypesize;
    unsigned long int blockMemory=268435456; //256 Mb
    void * buf;


    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
	//std::cerr<<"File mode 0: UChar\n";
	datatype=UChar;
	datatypesize = 1;
        buf = new unsigned char [blockMemory/datatypesize];
    }else if (header.mode==1){ //image : 16-bit halfwords
	//std::cerr<<"File mode 1: Short\n";
	datatype=Short;
	datatypesize = 2;
	buf = new short [blockMemory/datatypesize];
    }else if (header.mode==2){ //image : 32-bit reals
	//std::cerr<<"File mode 2: Float\n";
	datatype=Float;
        datatypesize = 4;
	buf = new float [blockMemory/datatypesize];
    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
	//std::cerr<<"File mode 5: Char\n";
	datatype=SChar;
	datatypesize = 1;
	buf = new char [blockMemory/datatypesize];
    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
	//std::cerr<<"File mode 6: UShort\n";
	datatype=UShort;
	datatypesize = 2;
	buf = new unsigned short [blockMemory/datatypesize];
    }else{
	//std::cerr<<"Unsupported type\n";
	datatype=UChar;
	datatypesize = 1;
	buf = new unsigned char [blockMemory/datatypesize];
    }

    unsigned long int blockSize=(unsigned long int)blockMemory/(double)datatypesize; 

    //std::cerr<<"read header\n";
    int error_fseek = fseek( m_fp, headerSizeMRC+datatypesize*sliceNumber*nxy, SEEK_SET );
    if (error_fseek != 0)
           return -1;
    
    //std::cerr<<"read data\n";
    for ( unsigned long int tt=0; tt<nxy; tt+=blockSize ){
	unsigned long int endRead=tt+blockSize;
	if ( endRead >= nxy ){
		endRead=nxy;
	}
	unsigned long int toRead=endRead-tt;

	// IT MAY BE SLIGHTLY MORE EFFICIENT TO PUT IT IN DIRECTLY
	// BUT NEED TO MANAGE POLYMORPHISM in EXECUTION
	/*if (sizeof(T)==datatypesize ){
		std::cerr<<"same type t="<<tt<<"->"<<endRead<<"   (max="<<nxyz<<")\n";
		size_t result = fread( &I[tt], datatypesize, toRead, m_fp );
		if (result < 0 ){
			std::cerr<<"Error: wrong read\n";
        	        return -1;
		}
	}*/	
		size_t result = fread( buf, datatypesize, toRead, m_fp );
		if (result < 0 ){
			std::cerr<<"Error: wrong read\n";
        	        return -1;
		}
	
		//std::cerr<<" block ="<<tt<<"->"<<endRead<<"   (max="<<nxyz<<")\n";

		    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
			for (unsigned long int ii=0;ii<toRead;ii++)
				I[tt+ii]=(T)(((unsigned char*) buf)[ii]);
		    }else if (header.mode==1){ //image : 16-bit halfwords
			for (unsigned long int ii=0;ii<toRead;ii++)
				I[tt+ii]=(T)(((short*) buf)[ii]);
		    }else if (header.mode==2){ //image : 32-bit reals
			for (unsigned long int ii=0;ii<toRead;ii++)
				I[tt+ii]=(T)(((float*) buf)[ii]);
		    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
			for (unsigned long int ii=0;ii<toRead;ii++)
				I[tt+ii]=(T)(((char*) buf)[ii]);
		    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
			for (unsigned long int ii=0;ii<toRead;ii++)			
				I[tt+ii]=(T)(((unsigned short*) buf)[ii]);
		    }else{
			//std::cerr<<"Unsupported type\n";
			for (unsigned long int ii=0;ii<toRead;ii++)
				I[tt+ii]=(T)(((unsigned char*) buf)[ii]);
		    }
	}

    
    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
        delete [] (unsigned char *) buf;
    }else if (header.mode==1){ //image : 16-bit halfwords
	delete [] (short *) buf;
    }else if (header.mode==2){ //image : 32-bit reals
	delete [](float *) buf;
    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
	delete [](char *) buf;
    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
	delete [](unsigned short *) buf;
    }else{
	delete [] (unsigned char *) buf;
    }

  //close the file
  fclose(m_fp);
  return 0;
}



template<typename T>
int mmmHeader(T * I, MRCHeader & header){
  //check the header

    //detect if big endilan or little endian
    if ( htonl(47) == 47 ) { // Big endian
       header.machst[0] = header.machst[1] = 17;
    } else {// Little endian  
	header.machst[0] = 68;
	header.machst[1] = 65;
    }


    unsigned long int nx = header.nx;
    unsigned long int ny = header.ny;
    unsigned long int nz = header.nz;
    if (nx<1) nx=1;
    if (ny<1) ny=1;
    if (nz<1) nz=1;
    unsigned long int nxyz = nx*ny*nz;


    DataType datatype;
    size_t datatypesize;
    unsigned long int blockMemory=268435456; //256 Mb
    void * buf = NULL;

    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
	//std::cerr<<"File mode 0: UChar\n";
	datatype=UChar;
	datatypesize = 1;
        buf = new unsigned char [blockMemory/datatypesize];
    }else if (header.mode==1){ //image : 16-bit halfwords
	//std::cerr<<"File mode 1: Short\n";
	datatype=Short;
	datatypesize = 2;
	buf = new short [blockMemory/datatypesize];
    }else if (header.mode==2){ //image : 32-bit reals
	//std::cerr<<"File mode 2: Float\n";
	datatype=Float;
        datatypesize = 4;
	buf = new float [blockMemory/datatypesize];
    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
	//std::cerr<<"File mode 5: Char\n";
	datatype=SChar;
	datatypesize = 1;
	buf = new char [blockMemory/datatypesize];
    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
	//std::cerr<<"File mode 6: UShort\n";
	datatype=UShort;
	datatypesize = 2;
	buf = new unsigned short [blockMemory/datatypesize];
    }else{
	//std::cerr<<"Unsupported type\n";
	datatype=UChar;
	datatypesize = 1;
	buf = new unsigned char [blockMemory/datatypesize];
    }

   if ( !buf ){
	std::cerr<<"failed to allocate memory ("<< (int)blockMemory/1024 <<"Mb) for the buffer\n";
	return -1;
   }

    //compute and write file header with amended statistical values
    double dmean = 0;
    double dRMS = 0;
    double dmin = (double)I[0];
    double dmax = (double)I[0];
    for (unsigned long int ii=0;ii<nxyz;ii++){
	dmean+=(double)I[ii];
    }
    dmean/=nxyz;
    for (unsigned long int ii=0;ii<nxyz;ii++){
	double tmpI=(double)I[ii];
	dRMS+=pow(tmpI-dmean,2.0);
	if (dmin>tmpI) dmin=tmpI;
	if (dmax<tmpI) dmax=tmpI;
    }
    dRMS/=(nxyz);
	
    header.dmean = (float)dmean;
    header.drms   = (float)dRMS;
    header.dmin  = (float)dmin;
    header.dmax  = (float)dmax;
    return 0;
}



// **********************
// **********************
//  I/O
//
// Write A VECTOR on a image
//
// **********************
// **********************
template<typename T>
int writeMrcImage(const char * filenameMRC, T * I, MRCHeader & header){
  //check the header

    //detect if big endilan or little endian
    if ( htonl(47) == 47 ) { // Big endian
       header.machst[0] = header.machst[1] = 17;
    } else {// Little endian  
	header.machst[0] = 68;
	header.machst[1] = 65;
    }


    unsigned long int nx = header.nx;
    unsigned long int ny = header.ny;
    unsigned long int nz = header.nz;
    if (nx<1) nx=1;
    if (ny<1) ny=1;
    if (nz<1) nz=1;
    unsigned long int nxyz = nx*ny*nz;


    DataType datatype;
    size_t datatypesize;
    unsigned long int blockMemory=268435456; //256 Mb
    void * buf = NULL;

    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
	//std::cerr<<"File mode 0: UChar\n";
	datatype=UChar;
	datatypesize = 1;
        buf = new unsigned char [blockMemory/datatypesize];
    }else if (header.mode==1){ //image : 16-bit halfwords
	//std::cerr<<"File mode 1: Short\n";
	datatype=Short;
	datatypesize = 2;
	buf = new short [blockMemory/datatypesize];
    }else if (header.mode==2){ //image : 32-bit reals
	//std::cerr<<"File mode 2: Float\n";
	datatype=Float;
        datatypesize = 4;
	buf = new float [blockMemory/datatypesize];
    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
	//std::cerr<<"File mode 5: Char\n";
	datatype=SChar;
	datatypesize = 1;
	buf = new char [blockMemory/datatypesize];
    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
	//std::cerr<<"File mode 6: UShort\n";
	datatype=UShort;
	datatypesize = 2;
	buf = new unsigned short [blockMemory/datatypesize];
    }else{
	//std::cerr<<"Unsupported type\n";
	datatype=UChar;
	datatypesize = 1;
	buf = new unsigned char [blockMemory/datatypesize];
    }

   if ( !buf ){
	std::cerr<<"failed to allocate memory ("<< (int)blockMemory/1024 <<"Mb) for the buffer\n";
	return 1;
   }

    //compute and write file header with amended statistical values
    double dmean = 0;
    double dRMS = 0;
    double dmin = (double)I[0];
    double dmax = (double)I[0];
    for (unsigned long int ii=0;ii<nxyz;ii++){
	dmean+=(double)I[ii];
    }
    dmean/=nxyz;
    for (unsigned long int ii=0;ii<nxyz;ii++){
	double tmpI=(double)I[ii];
	dRMS+=pow(tmpI-dmean,2.0);
	if (dmin>tmpI) dmin=tmpI;
	if (dmax<tmpI) dmax=tmpI;
    }
    dRMS/=(nxyz);
	
    header.dmean = (float)dmean;
    header.drms   = (float)dRMS;
    header.dmin  = (float)dmin;
    header.dmax  = (float)dmax;

     FILE *m_fp = NULL;
     m_fp=fopen(filenameMRC,"w");
     if(m_fp==NULL) return 1;

     //write file header
     //rewind(m_fp);
     fseek(m_fp, 0, SEEK_SET);
     fwrite(&header,1,1024,m_fp);

     unsigned long int blockSize=(unsigned long int)blockMemory/(double)datatypesize;


    for ( unsigned long int tt=0; tt<nxyz; tt+=blockSize ){
	unsigned long int endWrite=tt+blockSize;
	if ( endWrite >= nxyz ){
		endWrite=nxyz;
	}
	unsigned long int toWrite=endWrite-tt;

		//std::cerr<<" block ="<<tt<<"->"<<endWrite<<"   (max="<<nxyz<<")\n";

		    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
			for (unsigned long int ii=0;ii<toWrite;ii++)
				((unsigned char*)buf)[ii]=(unsigned char)I[tt+ii];
		    }else if (header.mode==1){ //image : 16-bit halfwords
			for (unsigned long int ii=0;ii<toWrite;ii++)
				((short*)buf)[ii]=(unsigned long int)I[tt+ii];
		    }else if (header.mode==2){ //image : 32-bit reals
			for (unsigned long int ii=0;ii<toWrite;ii++)
				((float*)buf)[ii]=(float)I[tt+ii];
		    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
			for (unsigned long int ii=0;ii<toWrite;ii++)
				((char*)buf)[ii]=(char)I[tt+ii];
		    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
			for (unsigned long int ii=0;ii<toWrite;ii++)			
				((unsigned short*)buf)[ii]=(unsigned short)I[tt+ii];
		    }else{
			//std::cerr<<"Unsupported type\n";
			for (unsigned long int ii=0;ii<toWrite;ii++)
				((unsigned char*)buf)[ii]=(unsigned char)I[tt+ii];
		    }
		size_t result = fwrite( buf, datatypesize, toWrite, m_fp );
		if (result < 0 ){
			std::cerr<<"Error: wrong read\n";
        	        return 1;
		}
	}



    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
        delete [] (unsigned char *) buf;
    }else if (header.mode==1){ //image : 16-bit halfwords
	delete [] (short *) buf;
    }else if (header.mode==2){ //image : 32-bit reals
	delete [](float *) buf;
    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
	delete [](char *) buf;
    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
	delete [](unsigned short *) buf;
    }else{
	delete [] (unsigned char *) buf;
    }

  //close the file
  fclose(m_fp);
    return 0;
}


int writeEmptyMrcImage(const char * filenameMRC, MRCHeader & header){
  //check the header

    //detect if big endilan or little endian
    if ( htonl(47) == 47 ) { // Big endian
       header.machst[0] = header.machst[1] = 17;
    } else {// Little endian  
	header.machst[0] = 68;
	header.machst[1] = 65;
    }


    unsigned long int nx = header.nx;
    unsigned long int ny = header.ny;
    unsigned long int nz = header.nz;
    if (nx<1) nx=1;
    if (ny<1) ny=1;
    if (nz<1) nz=1;
    unsigned long int nxyz = nx*ny*nz;


    DataType datatype;
    size_t datatypesize;
    unsigned long int blockMemory=268435456; //256 Mb
    void * buf = NULL;

    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
	//std::cerr<<"File mode 0: UChar\n";
	datatype=UChar;
	datatypesize = 1;
        buf = new unsigned char [blockMemory/datatypesize];
    }else if (header.mode==1){ //image : 16-bit halfwords
	//std::cerr<<"File mode 1: Short\n";
	datatype=Short;
	datatypesize = 2;
	buf = new short [blockMemory/datatypesize];
    }else if (header.mode==2){ //image : 32-bit reals
	//std::cerr<<"File mode 2: Float\n";
	datatype=Float;
        datatypesize = 4;
	buf = new float [blockMemory/datatypesize];
    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
	//std::cerr<<"File mode 5: Char\n";
	datatype=SChar;
	datatypesize = 1;
	buf = new char [blockMemory/datatypesize];
    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
	//std::cerr<<"File mode 6: UShort\n";
	datatype=UShort;
	datatypesize = 2;
	buf = new unsigned short [blockMemory/datatypesize];
    }else{
	//std::cerr<<"Unsupported type\n";
	datatype=UChar;
	datatypesize = 1;
	buf = new unsigned char [blockMemory/datatypesize];
    }

   if ( !buf ){
	std::cerr<<"failed to allocate memory ("<< (int)blockMemory/1024 <<"Mb) for the buffer\n";
	return 1;
   }

    //compute and write file header with amended statistical values
    header.dmean = 0;
    header.drms   = 0;
    header.dmin  = 0;
    header.dmax  = 0;
    FILE *m_fp = NULL;
    m_fp=fopen(filenameMRC,"w");
    if(m_fp==NULL) return 1;

     //write file header
     //rewind(m_fp);
     fseek(m_fp, 0, SEEK_SET);
     fwrite(&header,1,1024,m_fp);

     unsigned long int blockSize=(unsigned long int)blockMemory/(double)datatypesize;


    for ( unsigned long int tt=0; tt<nxyz; tt+=blockSize ){
	unsigned long int endWrite=tt+blockSize;
	if ( endWrite >= nxyz ){
		endWrite=nxyz;
	}
	unsigned long int toWrite=endWrite-tt;

		    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
			for (unsigned long int ii=0;ii<toWrite;ii++)
				((unsigned char*)buf)[ii]=0;
		    }else if (header.mode==1){ //image : 16-bit halfwords
			for (unsigned long int ii=0;ii<toWrite;ii++)
				((short*)buf)[ii]=0;
		    }else if (header.mode==2){ //image : 32-bit reals
			for (unsigned long int ii=0;ii<toWrite;ii++)
				((float*)buf)[ii]=0;
		    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
			for (unsigned long int ii=0;ii<toWrite;ii++)
				((char*)buf)[ii]=0;
		    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
			for (unsigned long int ii=0;ii<toWrite;ii++)			
				((unsigned short*)buf)[ii]=0;
		    }else{
			//std::cerr<<"Unsupported type\n";
			for (unsigned long int ii=0;ii<toWrite;ii++)
				((unsigned char*)buf)[ii]=0;
		    }
		size_t result = fwrite( buf, datatypesize, toWrite, m_fp );
		if (result < 0 ){
			std::cerr<<"Error: wrong read\n";
        	        return 1;
		}
	}

    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
        delete [] (unsigned char *) buf;
    }else if (header.mode==1){ //image : 16-bit halfwords
	delete [] (short *) buf;
    }else if (header.mode==2){ //image : 32-bit reals
	delete [](float *) buf;
    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
	delete [](char *) buf;
    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
	delete [](unsigned short *) buf;
    }else{
	delete [] (unsigned char *) buf;
    }

  //close the file
  fclose(m_fp);
    return 0;
}



// *********************
template<typename T>
int replaceMrcSlice(const char * filenameMRC, T * I, MRCHeader & header, long int sliceNumber){
  //check the header

    //detect if big endilan or little endian
    if ( htonl(47) == 47 ) { // Big endian
       header.machst[0] = header.machst[1] = 17;
    } else {// Little endian  
	header.machst[0] = 68;
	header.machst[1] = 65;
    }
    const int headerSizeMRC = 1024;

    unsigned long int nx = header.nx;
    unsigned long int ny = header.ny;
    unsigned long int nz = header.nz;
    if (nx<1) nx=1;
    if (ny<1) ny=1;
    if (nz<1) nz=1;
    unsigned long int nxy = nx*ny;
    //unsigned long int nxyz = nx*ny*nz;


    DataType datatype;
    size_t datatypesize;
    unsigned long int blockMemory=268435456; //256 Mb
    void * buf = NULL;

    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
	//std::cerr<<"File mode 0: UChar\n";
	datatype=UChar;
	datatypesize = 1;
        buf = new unsigned char [blockMemory/datatypesize];
    }else if (header.mode==1){ //image : 16-bit halfwords
	//std::cerr<<"File mode 1: Short\n";
	datatype=Short;
	datatypesize = 2;
	buf = new short [blockMemory/datatypesize];
    }else if (header.mode==2){ //image : 32-bit reals
	//std::cerr<<"File mode 2: Float\n";
	datatype=Float;
        datatypesize = 4;
	buf = new float [blockMemory/datatypesize];
    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
	//std::cerr<<"File mode 5: Char\n";
	datatype=SChar;
	datatypesize = 1;
	buf = new char [blockMemory/datatypesize];
    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
	//std::cerr<<"File mode 6: UShort\n";
	datatype=UShort;
	datatypesize = 2;
	buf = new unsigned short [blockMemory/datatypesize];
    }else{
	//std::cerr<<"Unsupported type\n";
	datatype=UChar;
	datatypesize = 1;
	buf = new unsigned char [blockMemory/datatypesize];
    }

   if ( !buf ){
	std::cerr<<"failed to allocate memory ("<< (int)blockMemory/headerSizeMRC <<"Mb) for the buffer\n";
	return 1;
   }


    FILE *m_fp = NULL;
    m_fp=fopen(filenameMRC,"r+");
    if(m_fp==NULL) return 1;

    int error_fseek = fseek( m_fp, headerSizeMRC+datatypesize*sliceNumber*nxy, SEEK_SET );
    if (error_fseek != 0){
           std::cerr<<"ERROR: writing data on file\n";
           exit(1);
    }


     unsigned long int blockSize=(unsigned long int)blockMemory/(double)datatypesize;


    for ( unsigned long int tt=0; tt<nxy; tt+=blockSize ){
	unsigned long int endWrite=tt+blockSize;
	if ( endWrite >= nxy ){
		endWrite=nxy;
	}
	unsigned long int toWrite=endWrite-tt;

		//std::cerr<<" block ="<<tt<<"->"<<endWrite<<"   (max="<<nxyz<<")\n";

		    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
			for (unsigned long int ii=0;ii<toWrite;ii++)
				((unsigned char*)buf)[ii]=(unsigned char)I[tt+ii];
		    }else if (header.mode==1){ //image : 16-bit halfwords
			for (unsigned long int ii=0;ii<toWrite;ii++)
				((short*)buf)[ii]=(unsigned long int)I[tt+ii];
		    }else if (header.mode==2){ //image : 32-bit reals
			for (unsigned long int ii=0;ii<toWrite;ii++)
				((float*)buf)[ii]=(float)I[tt+ii];
		    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
			for (unsigned long int ii=0;ii<toWrite;ii++)
				((char*)buf)[ii]=(char)I[tt+ii];
		    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
			for (unsigned long int ii=0;ii<toWrite;ii++)			
				((unsigned short*)buf)[ii]=(unsigned short)I[tt+ii];
		    }else{
			//std::cerr<<"Unsupported type\n";
			for (unsigned long int ii=0;ii<toWrite;ii++)
				((unsigned char*)buf)[ii]=(unsigned char)I[tt+ii];
		    }
		size_t result = fwrite( buf, datatypesize, toWrite, m_fp );
		if (result < 0 ){
			std::cerr<<"Error: wrong read\n";
        	        return 1;
		}
	}






    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
        delete [] (unsigned char *) buf;
    }else if (header.mode==1){ //image : 16-bit halfwords
	delete [] (short *) buf;
    }else if (header.mode==2){ //image : 32-bit reals
	delete [](float *) buf;
    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
	delete [](char *) buf;
    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
	delete [](unsigned short *) buf;
    }else{
	delete [] (unsigned char *) buf;
    }



  //close the file
  fclose(m_fp);
  return 0;
}


// ***********
//
template<typename T>
int mmm( const char * filename){
    MRCHeader header;
    readHeaderMrc(filename, header);
    unsigned long int nx=header.nx;
    unsigned long int ny=header.ny;
    unsigned long int nz=header.nz;
    unsigned long int nxyz=nx*ny*nz;
    T * I = new T[nxyz];
    readMrcImage(filename, I, header);    
    mmmHeader ( I, header );
    writeMrcImage(filename, I, header);
    delete [] I;
    return 0;
}



// **********************
// **********************
//  I/O
//
// WRITE A VECTOR on a image
//
// **********************
// **********************
template<typename T>
void writeMrcImage(const char * filenameMRC, T * I, const unsigned long int nx, const unsigned long int ny, const unsigned long int nz, double APix=1.0, DataType datatype= Float){

   MRCHeader header;
    if (datatype==UChar){
	header.mode=0;
    }else if (datatype==Short){
	header.mode=1;
    }else if (datatype==Float){
	header.mode=2;
    }else if (datatype==SChar){
	header.mode=5;
    }else if (datatype==UShort){
	header.mode=6;
    }else{
	//std::cerr<<"Unsupported type\n";
	datatype=UChar;
    }

   header.nx=nx;
   header.ny=ny;
   header.nz=nz;
   header.mx=nx;
   header.my=ny;
   header.mz=nz;
   header.a=nx*APix;
   header.b=ny*APix;
   header.c=nz*APix;


    // fix some portions of the header
    header.mapc = 1;
    header.mapr = 2;
    header.maps = 3;
    header.alpha = (float)90.;
    header.beta = (float)90.;
    header.gamma = (float)90.;
    header.xOrigin = (float)0.;
    header.yOrigin = (float)0.;
    header.zOrigin = (float)0.;
    header.nxstart = (int)0;
    header.nystart = (int)0;
    header.nzstart = (int)0;
    header.ispg = (int)0;
    header.nsymbt = (int)0;
    header.map[0]='M';
    header.map[1]='A';
    header.map[2]='P';
    header.map[3]=' ';
    header.machst[0]='M';
    header.machst[1]='A';
    header.machst[2]='U';
    header.machst[3]=' ';
    header.nlabl=0;
   writeMrcImage(filenameMRC, I, header);

}



double calculateXPixelSpacing(const MRCHeader& header) {
    double spacing;
    // Calculate spacing for each axis
    spacing = header.mx > 0 ? header.a / header.mx : 0;
    return spacing;
}


#endif
