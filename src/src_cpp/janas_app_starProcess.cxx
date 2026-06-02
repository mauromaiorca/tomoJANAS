/*
 * File: janas_app_starProcess.h
 * (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
 */

/**
 * Assess Particles Program
 * \defgroup starProcess 
 */


#include <iostream>
#include <complex>
#include <fstream>
#include <list>
#include <iterator>     // std::back_inserter
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <sstream>
#include <string>
#include <regex>

		
//ITK

#define PI 3.14159265358979323846
#include "CsvStarReadWriteAnalyseLibs.h"
//#include "reconstructionLib.h"
//
//#include "fourierFunctions.h"
//#include "genericFilters.h"
//#include "msaEstimate.h"
//#include "msaLibs.h"
//#include "msaCsvReadWriteAnalyse.h"
//#include "plot2d.h"

#include "mrcIO.h"

#include "scores.h"


//#include "msaLibs.h"
//#include "scoringFunctionsLib.h"
//#include "preprocessingLib.h"
//#include "projectionsLib.h"
//#include "fiboLibs.h"
//#include "riaLibs.h"
//#include "analysisLib.h"
//#include "eulerAnalysisLibs.h"
//#include "visualizationLibs.h"




//




// #####################################
//AUSILIARY FUNCTION FOR SIMPLIFY NAMEFILES for COMPARISON
std::vector< std::string > simplifyListNamefiles (std::vector< std::string > inputList){
       std::vector< std::string > returnVector;
       for (unsigned long int ii=0;ii<inputList.size();ii++){
         //std::cerr<<inputList[ii]<<"  =>  ";
         std::string reducedMicrographName ( inputList[ii] );
         int lastDirChar=reducedMicrographName.find_last_of('/');
         if ( lastDirChar  < reducedMicrographName.size() ){
           reducedMicrographName=reducedMicrographName.substr(lastDirChar+1);
         }
         reducedMicrographName=reducedMicrographName.substr(0,reducedMicrographName.find_last_of('.'));
         std::string charsToRemove ("._-");
         for (char c: charsToRemove){
          reducedMicrographName.erase(std::remove(reducedMicrographName.begin(), reducedMicrographName.end(), c), reducedMicrographName.end());
         }
         //std::cerr<<reducedMicrographName<<"\n";
         returnVector.push_back(reducedMicrographName);
       }
       return returnVector;
   }


std::vector< std::string > reduceDirNamefiles (std::vector< std::string > inputList, int dept = 0){
  std::vector< std::string > returnVector;
  for (unsigned long int ii=0;ii<inputList.size();ii++){
         std::string reducedMicrographName ( inputList[ii] );
         int lastDirCharPos=-1;
         int tmpDept=dept+1;
         std::string tmpString (inputList[ii]);
         do{
           lastDirCharPos=tmpString.find_last_of('/');
           tmpString=tmpString.substr(0,lastDirCharPos-1);
           tmpDept--;
           //std::cerr<<tmpString<<"   ";
         } while (tmpDept>0 && lastDirCharPos -1 >=0 );
         
         if (lastDirCharPos < -1 || lastDirCharPos >= reducedMicrographName.size() -1 ){
          lastDirCharPos = -1;
         }
         
         //std::cerr<<"\n";
         //std::cerr<<"  || "<< inputList[ii] <<"  =>  " << reducedMicrographName.substr(lastDirCharPos+1)<<"\n";
         returnVector.push_back( reducedMicrographName.substr(lastDirCharPos+1) );
  }
  return returnVector;
}



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




template<typename T>
std::vector<long int> bubbleSortAscendingIndexes(const std::vector<T> valuesIn){
  long int size=valuesIn.size();
  T * values = new T [size];
  std::vector<long int> indexes;
  for (unsigned long int i=0; i<size; i++){
    values[i]=valuesIn[i];
    indexes.push_back(i);
  }
  for (long int i = (size - 1); i > 0; i--)
    {
      for (long int j = 1; j <= i; j++)
    {
      if (values[j - 1] > values[j])
        {
          T temp = values[j - 1];
          values[j - 1] = values[j];
          values[j] = temp;
          long int tmpIndex = indexes[j - 1];
          indexes[j - 1] = indexes[j];
          indexes[j] = tmpIndex;

        }
    }
    }
  delete [] values;
  return indexes;
}

double archDistanceSquaredApproximate(double phi, double theta, double phiV, double thetaV){
   double PI180=0.01745329251;
   double phi1=phi*PI180;
   double phi2=phiV*PI180;
   double theta1=theta*PI180;
   double theta2=thetaV*PI180;

   double point[]={sin(theta1)*cos(phi1), sin(theta1)*sin(phi1), sin(theta1)};
   double point2[]={sin(theta2)*cos(phi2), sin(theta2)*sin(phi2), sin(theta2)};
   double distanceSq=pow(point[0]-point2[0],2.0)+pow(point[1]-point2[1],2.0)+pow(point[2]-point2[2],2.0);
   return distanceSq;
}



//vector
//function that assignes for each particles index an index for micrograph line holding it. -1 if the index is not there
std::vector<long int> buildParticleMicrographMatchesIdx(char * sourceMicrographFileName, char * sourceParticlesStarFileName){
std::cerr<<sourceMicrographFileName<<"   "<<sourceParticlesStarFileName<<"\n";
       std::vector<long int> outputVector;
       std::string micrographNameTag("_rlnMicrographName");
       std::vector<std::string> originalMicrographItems;
       std::vector<std::string> referenceParticlesItems;
       readStar(originalMicrographItems, micrographNameTag, sourceMicrographFileName);
       readStar(referenceParticlesItems, micrographNameTag, sourceParticlesStarFileName);
       
// ###################################
// MAKE AN ANALYSIS OF THE HEADER
// ###################################
       std::vector<std::string> referenceStarParticleFields;
       std::vector<int> referenceStarParticleFieldsIdx;
       getStarHeaders(referenceStarParticleFields, referenceStarParticleFieldsIdx, sourceParticlesStarFileName);
       //for (unsigned long int ii=0; ii<referenceStarParticleFields.size(); ii++){
       //  std::cerr<< referenceStarParticleFields[ii]<< "\n";
       //}
       //std::cerr<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
       std::vector<std::string> referenceStarMicrographFields;
       std::vector<int> referenceStarMicrographFieldsIdx;
       getStarHeaders(referenceStarMicrographFields, referenceStarMicrographFieldsIdx, sourceMicrographFileName); 
       //for (unsigned long int ii=0; ii<referenceStarMicrographFields.size(); ii++){
       //  std::cerr<< referenceStarMicrographFields[ii]<< "\n";
       //} 
       
       std::vector< std::string > simplifiedReferenceParticlesItems= simplifyListNamefiles (referenceParticlesItems);
       //std::cerr<<"***********************************\n***********************************\n";
       std::vector< std::string > simplifiedOriginalMicrographItems= simplifyListNamefiles (originalMicrographItems);

       if (simplifiedOriginalMicrographItems.size() < 1){
         std::cerr<<"ERROR: original micrograph items is empty, exiting!\n";
         exit(0);
       }
       unsigned long int lastReferenceMicrographIdx = 0;
       //std::cerr<<"###########################\n###########################\n";
       for (unsigned long int jj=0; jj<simplifiedReferenceParticlesItems.size() ; jj++){
        
         bool foundReferenceMicrograph = false;
         if ( simplifiedReferenceParticlesItems[jj].find(simplifiedOriginalMicrographItems[lastReferenceMicrographIdx]) != std::string::npos ){
                 //std::cerr<< "||  ";
                 //std::cerr<<simplifiedOriginalMicrographItems[lastReferenceMicrographIdx]<<"    " <<simplifiedReferenceParticlesItems[jj]<<" ~~ ";
                 foundReferenceMicrograph=true;
                 //std::cerr<<"OK  \n";
         }else{
                 for (unsigned long int ii=0; ii<simplifiedOriginalMicrographItems.size() && !foundReferenceMicrograph; ii++){
                   //std::cerr<< "||  ";
                   //std::cerr<<simplifiedOriginalMicrographItems[ii]<<"    " <<simplifiedReferenceParticlesItems[jj]<<" ~~ ";
                   if ( simplifiedReferenceParticlesItems[jj].find(simplifiedOriginalMicrographItems[ii]) != std::string::npos ){
                     foundReferenceMicrograph=true;
                     lastReferenceMicrographIdx=ii;
                     //std::cerr<<"OK  \n";
                   }else{
                     //std::cerr<<"DIFFERENT  \n";
                   }
                 }
         }
         if (foundReferenceMicrograph){
           outputVector.push_back(lastReferenceMicrographIdx);
           //std::cerr << simplifiedReferenceParticlesItems[jj]<<"  =>  "<< simplifiedOriginalMicrographItems[referenceMicrographIdx]<<"\n";
         }else{
          outputVector.push_back(-1);
         }
       } 

   for (long int ii=0;ii<referenceParticlesItems.size(); ii++){
     std::cerr<<referenceParticlesItems[ii]<<"  --- >  ";
     if (outputVector[ii]>=0){
       std::cerr<<originalMicrographItems[outputVector[ii]]<<"\n";
     }else{
       std::cerr<<"NOT THERE \n";
     }
   }
   //readStar(originalMicrographItems, micrographNameTag, sourceMicrographFileName);
   //readStar(referenceParticlesItems, micrographNameTag, sourceParticlesStarFileName);

   return outputVector;
}

//
void updateMicrographStarFile1(char * sourceMicrographFileName, char * sourceParticlesStarFileName, char * outputFile, double bin = 1.0f){
       std::vector<std::string> referenceStarParticleFields;
       std::vector<int> referenceStarParticleFieldsIdx;
       getStarHeaders(referenceStarParticleFields, referenceStarParticleFieldsIdx, sourceParticlesStarFileName);

       std::vector<std::string> referenceStarMicrographFields;
       std::vector<int> referenceStarMicrographFieldsIdx;
       std::vector<std::string> selectedReferenceStarParticlesFields;
       std::vector<int> selectedReferenceStarParticlesFieldsIdx;
       getStarHeaders(referenceStarMicrographFields, referenceStarMicrographFieldsIdx, sourceMicrographFileName); 

       for (int jj=0 ; jj<referenceStarParticleFields.size(); jj++ ){
         //check if the element is in list micrographs
         bool found = false;
         for (int ii=0; ii<referenceStarMicrographFields.size() && !found; ii++){
           if(referenceStarMicrographFields[ii].find(referenceStarParticleFields[jj]) != std::string::npos ){
            //std::cerr<<std::string::npos<<"\n";
            found = true;
           }
         }
         if (!found){
            selectedReferenceStarParticlesFields.push_back(referenceStarParticleFields[jj]);
            selectedReferenceStarParticlesFieldsIdx.push_back(referenceStarParticleFieldsIdx[jj]);
         }
       }
       for (int jj=0 ; jj<selectedReferenceStarParticlesFields.size(); jj++ ){
         std::cerr<<selectedReferenceStarParticlesFields[jj]<<"\n";
       }
       

       std::ifstream fileMicrographs(sourceMicrographFileName);
       std::ifstream fileParticles(sourceParticlesStarFileName);
       std::ofstream fileOutput;
       fileOutput.open(outputFile);
       fileOutput.close();
       fileOutput.open (outputFile, std::ofstream::out | std::ofstream::app);   

       
       fileOutput<<"\ndata_\n\nloop_\n";
       for (int ii=0;ii<referenceStarMicrographFields.size();ii++){
          fileOutput<<referenceStarMicrographFields[ii]<<" #"<< ii + 1 <<"\n";
       }
       //fileOutput<<"######\n";
       for (int ii=0;ii<selectedReferenceStarParticlesFields.size();ii++){
          fileOutput<<selectedReferenceStarParticlesFields[ii]<<" #"<< referenceStarMicrographFields.size()+ii + 1 <<"\n";
       }

       std::vector<long int> indexMicrographs=buildParticleMicrographMatchesIdx(sourceMicrographFileName, sourceParticlesStarFileName);

       //get in line
       //fileMicrographs
       //std::vector<int> referenceStarMicrographFieldsIdx;
       //std::vector<std::string> selectedReferenceStarParticlesFields;
       long int startMicrograph=getStarStart(sourceMicrographFileName);
       long int startParticles=getStarStart(sourceParticlesStarFileName);
       unsigned long int mmCounter=0;
       unsigned long int ccCounter=0;
       unsigned long int counter=0;
       std::string strLine;
       
       //get micrograph lines
       std::vector<std::string> micrographLines;
       while ( std::getline(fileMicrographs, strLine) ){
         if (++mmCounter > startMicrograph){
          micrographLines.push_back(strLine);
         }
       }
       
       //TO USE: int findColumnItemPosition(std::string str, int idx)
       int imageNameIdx = -1;
       int CoordinateXIdx = -1;
       int CoordinateYIdx = -1;
       int OriginXIdx = -1;
       int OriginYIdx = -1;

       for (int ii=0; ii<selectedReferenceStarParticlesFields.size(); ii++ ){
         if ( selectedReferenceStarParticlesFields[ii].compare("_rlnImageName")==0){
           imageNameIdx=ii;
         }
         if ( selectedReferenceStarParticlesFields[ii].compare("_rlnCoordinateX")==0){
           CoordinateXIdx=ii;
         }
         if ( selectedReferenceStarParticlesFields[ii].compare("_rlnCoordinateY")==0){
           CoordinateYIdx=ii;
         }
         if ( selectedReferenceStarParticlesFields[ii].compare("_rlnOriginX")==0){
           OriginXIdx=ii;
         }
         if ( selectedReferenceStarParticlesFields[ii].compare("_rlnOriginY")==0){
           OriginYIdx=ii;
         }
       }
       //std::cerr<<imageNameIdx<< "   " << CoordinateXIdx << "   " << CoordinateYIdx << "\n\n";
       std::cerr<<"Binning="<<bin<<"\n";
       
       std::vector<std::string> originalMicrographItems;
       readStar(originalMicrographItems, "_rlnMicrographName", sourceMicrographFileName);
       
       
       while ( std::getline(fileParticles, strLine) ){
         if (++ccCounter > startParticles){
            if ( indexMicrographs[counter] >= 0 && indexMicrographs[counter] < micrographLines.size() ){
               fileOutput<< micrographLines[indexMicrographs[counter]] <<  "  ";
               //need to put the line for existing fields
               //selectedReferenceStarParticlesFieldsIdx.push_back(referenceStarMicrographFieldsIdx[jj]);
               std::vector<std::string> particleLine = stringLineToVector(strLine);
               for (int jj=0 ; jj<selectedReferenceStarParticlesFields.size(); jj++ ){
                 if (jj==imageNameIdx){//is imageName
                   std::vector<std::string> tmpVect = stringLineToVector(particleLine[selectedReferenceStarParticlesFieldsIdx[jj]],'@');
                   fileOutput << "   "<< tmpVect[0] << "@"<< originalMicrographItems[indexMicrographs[counter]];
                   //std::cerr<< micrographLines[indexMicrographs[counter]] << "\n";
                 }if ( (jj==CoordinateXIdx  || jj==CoordinateYIdx || jj==OriginXIdx || jj==OriginYIdx) && bin != 1.0f){
                   double binnedVal=stof(particleLine[selectedReferenceStarParticlesFieldsIdx[jj]]);
                   binnedVal/=(double)bin;
                   fileOutput << "   " << std::setfill(' ') << std::fixed << std::right << std::setw(8) << std::showpoint << std::setprecision(6)  <<  binnedVal;
                
                 }else{
                   fileOutput << "  "  << particleLine[selectedReferenceStarParticlesFieldsIdx[jj]] ;
                 }
               }
               fileOutput <<"\n";
            } else {
               std::cerr<<"WARNING: micrograph "<< sourceMicrographFileName << " not found \n";
            }
            counter++;
         }
         
       }
       
       
       fileOutput.close();

}




// #############################################
// 
//
void updateMicrographStarFile(char * sourceMicrographFileName, char * sourceParticlesStarFileName, char * outputFile, double bin = 1.0f){
       //std::cerr<<"eccolo\n";
       std::vector<std::string> referenceStarParticleFields;
       std::vector<int> referenceStarParticleFieldsIdx;
       getStarHeaders(referenceStarParticleFields, referenceStarParticleFieldsIdx, sourceParticlesStarFileName);

       std::vector<std::string> referenceStarMicrographFields;
       std::vector<int> referenceStarMicrographFieldsIdx;
       std::vector<std::string> selectedReferenceStarParticlesFields;
       std::vector<int> selectedReferenceStarParticlesFieldsIdx;
       getStarHeaders(referenceStarMicrographFields, referenceStarMicrographFieldsIdx, sourceMicrographFileName); 

       for (int jj=0 ; jj<referenceStarParticleFields.size(); jj++ ){
         //check if the element is in list micrographs
         bool found = false;
         for (int ii=0; ii<referenceStarMicrographFields.size() && !found; ii++){
           if(referenceStarMicrographFields[ii].find(referenceStarParticleFields[jj]) != std::string::npos ){
            //std::cerr<<std::string::npos<<"\n";
            found = true;
           }
         }
         if (!found){
            selectedReferenceStarParticlesFields.push_back(referenceStarParticleFields[jj]);
            selectedReferenceStarParticlesFieldsIdx.push_back(referenceStarParticleFieldsIdx[jj]);
         }
       }
       for (int jj=0 ; jj<selectedReferenceStarParticlesFields.size(); jj++ ){
         std::cerr<<selectedReferenceStarParticlesFields[jj]<<"\n";
       }
       

       std::ifstream fileMicrographs(sourceMicrographFileName);
       std::ifstream fileParticles(sourceParticlesStarFileName);
       std::ofstream fileOutput;
       fileOutput.open(outputFile);
       fileOutput.close();
       fileOutput.open (outputFile, std::ofstream::out | std::ofstream::app);   

       
       fileOutput<<"\ndata_\n\nloop_\n";
       for (int ii=0;ii<referenceStarMicrographFields.size();ii++){
          fileOutput<<referenceStarMicrographFields[ii]<<" #"<< ii + 1 <<"\n";
       }
       //fileOutput<<"######\n";
       for (int ii=0;ii<selectedReferenceStarParticlesFields.size();ii++){
          fileOutput<<selectedReferenceStarParticlesFields[ii]<<" #"<< referenceStarMicrographFields.size()+ii + 1 <<"\n";
       }

       std::vector<long int> indexMicrographs=buildParticleMicrographMatchesIdx(sourceMicrographFileName, sourceParticlesStarFileName);

       //get in line
       //fileMicrographs
       //std::vector<int> referenceStarMicrographFieldsIdx;
       //std::vector<std::string> selectedReferenceStarParticlesFields;
       long int startMicrograph=getStarStart(sourceMicrographFileName);
       long int startParticles=getStarStart(sourceParticlesStarFileName);
       unsigned long int mmCounter=0;
       unsigned long int ccCounter=0;
       unsigned long int counter=0;
       std::string strLine;
       
       //get micrograph lines
       std::vector<std::string> micrographLines;
       while ( std::getline(fileMicrographs, strLine) ){
         if (++mmCounter > startMicrograph){
          micrographLines.push_back(strLine);
         }
       }
       
       //TO USE: int findColumnItemPosition(std::string str, int idx)
       int imageNameIdx = -1;
       int CoordinateXIdx = -1;
       int CoordinateYIdx = -1;
       int OriginXIdx = -1;
       int OriginYIdx = -1;

       for (int ii=0; ii<selectedReferenceStarParticlesFields.size(); ii++ ){
         if ( selectedReferenceStarParticlesFields[ii].compare("_rlnImageName")==0){
           imageNameIdx=ii;
         }
         if ( selectedReferenceStarParticlesFields[ii].compare("_rlnCoordinateX")==0){
           CoordinateXIdx=ii;
         }
         if ( selectedReferenceStarParticlesFields[ii].compare("_rlnCoordinateY")==0){
           CoordinateYIdx=ii;
         }
         if ( selectedReferenceStarParticlesFields[ii].compare("_rlnOriginX")==0){
           OriginXIdx=ii;
         }
         if ( selectedReferenceStarParticlesFields[ii].compare("_rlnOriginY")==0){
           OriginYIdx=ii;
         }
       }
       //std::cerr<<imageNameIdx<< "   " << CoordinateXIdx << "   " << CoordinateYIdx << "\n\n";
       std::cerr<<"Binning="<<bin<<"\n";
       
       std::vector<std::string> originalMicrographItems;
       readStar(originalMicrographItems, "_rlnMicrographName", sourceMicrographFileName);
       
       
       while ( std::getline(fileParticles, strLine) ){
         if (++ccCounter > startParticles){
            if ( indexMicrographs[counter] >= 0 && indexMicrographs[counter] < micrographLines.size() ){
               fileOutput<< micrographLines[indexMicrographs[counter]] <<  "  ";
               //need to put the line for existing fields
               //selectedReferenceStarParticlesFieldsIdx.push_back(referenceStarMicrographFieldsIdx[jj]);
               std::vector<std::string> particleLine = stringLineToVector(strLine);
               for (int jj=0 ; jj<selectedReferenceStarParticlesFields.size(); jj++ ){
                 if (jj==imageNameIdx){//is imageName
                   std::vector<std::string> tmpVect = stringLineToVector(particleLine[selectedReferenceStarParticlesFieldsIdx[jj]],'@');
                   fileOutput << "   "<< tmpVect[0] << "@"<< originalMicrographItems[indexMicrographs[counter]];
                   //std::cerr<< micrographLines[indexMicrographs[counter]] << "\n";
                 }else if ( (jj==CoordinateXIdx  || jj==CoordinateYIdx || jj==OriginXIdx || jj==OriginYIdx) && bin != 1.0f){
                   double binnedVal=stof(particleLine[selectedReferenceStarParticlesFieldsIdx[jj]]);
                   binnedVal/=(double)bin;
                   fileOutput << "   " << std::setfill(' ') << std::fixed << std::right << std::setw(8) << std::showpoint << std::setprecision(6)  <<  binnedVal;
                 }else{
                   fileOutput << "  "  << particleLine[selectedReferenceStarParticlesFieldsIdx[jj]] ;
                 }
               }
               fileOutput <<"\n";
            } else {
               std::cerr<<"WARNING: micrograph "<< sourceMicrographFileName << " not found \n";
            }
            counter++;
         }
         
       }
       
       
       fileOutput.close();
       return;

}






  typedef float WorkingPixelType;
  typedef float OutputPixelType;
  const unsigned int Dimension = 3;




// *************************************
// inputParametersType
// **************************************
 /**
 * input parameters for starProcess
 *  \ingroup starProcess
 */
typedef struct inputParametersType {
    bool verboseOn;



    char * starFileIn;
    char * starFileOut;
    char * vemFileOut;
    char * csvFileOut;
    char * columnToDelete;
    char * RetainGroupMetadata;
    char * halfMapsTag;
    char * halfMap1Out;
    char * halfMap2Out;    
    char * referenceParticlesStarFile;
    bool micrographs;
    int depthDirMicrograph;
    float binningFactor;
    int RetainGroupNumber;
    char * coord_dir;
    char * coord_suffix;
    bool showInfo;
    bool showInfoDiff;
    bool stackNumbering;
    char * refinedFileName;
    char * refinedStackFileName;
    char * imodXfFile;
    char * imodXfFileOptions;

    char * micrographsWithCtf;
    bool sortImageName;
    char * commonImagesStarFile;   
    char * reducedImagesStarFile;   
    char * templateForUpdateStarFile;
    char * parametersToUpdate;
    char * starInfoDifferenceFile;
    char * tagToExport;
    char * fileToExportTag;
    char * stackFileNameRenamed;
    char * labelToMultiply;
    double valueLabelToMultiply;
    
    char * SourceFileForReplacingLabel;
    char * SourceLabelsNameCsvForReplacing;
    char * DestinationLabelNameCsvForReplacing;


    char * tag1_toInvert;
    char * tag2_toInvert;



    char * extractFromSimilarParametersTemplateFullFile;
    char * ParametersToCompare;
    char * outputImageTag;
    bool haveBackupImageNameTag;

    bool DoCheckForSimilarImages;


    char * newValueToSetLabel;
    char * labelNameToUpdate;

    char * SourceFilesForAveragingLabels;
    char * SourceLabelnameForAveraging;
    
    char * gtFileForAssessingLabels;
    char * gtLabelsnameForAssessing;
    bool updateRandomSubset;

    char * csvListFileToMerge;
    char * csvOutListFileToSplit;

    char * labelToSort;
    char * labelToNumberAsIndex;


    int numHomogeneousDistribution;
    char * simulatedStackName;



    char * starFileWithAutorefineUpdates;

    char * datAnglesEFFileOut;


    bool removeImageNamePrefix;

    char * classNameToSelect;

}inputParametersType;


/* ******************************************
 *  USAGE
 ***************************************** */
void usage(  char ** argv ){
    std::cerr<<"starProcess (c) Mauro Maiorca\n";
    std::cerr<<"\n";
    std::cerr<<"Usage: " << argv[0] << "\n";
    std::cerr<<"              --i starFileIn.star\n";
    std::cerr<<"              --o starFileOut.star\n";
    std::cerr<<"\n";    
    std::cerr<<"       Import from micrograph options:\n";    
    std::cerr<<"              --micrographs [depthDirMicrograph=0]\n";
    std::cerr<<"                     (produce the list of unique micrographs)\n";
    std::cerr<<"              --coordinates coord_dir coord_suffix [depthDirMicrograph=0]\n";
    std::cerr<<"              --rp referenceParticlesStarFile.star [binningFactor=1]\n";
    std::cerr<<"                   (to get the fields)\n";
    std::cerr<<"              --retrieveRefined refinedFileName.star refinedStackFileName.mrc\n";
    std::cerr<<"                   (retrieve refined information, except for the file name)\n";    
    std::cerr<<"              --importCtf micrographsWithCtf.star\n";
    std::cerr<<"                   (import the CTF from micrographsCTF file)\n";
    std::cerr<<"\n";
    std::cerr<<"       Image Processing based operations:\n";    
    std::cerr<<"              --checkForSimilarImages \n";
    std::cerr<<"                   (gets as output star file two extra tags with CC score and closest particle) \n";

    std::cerr<<"\n";

    
    std::cerr<<"       Export operations:\n";    
    std::cerr<<"              --hm halfMap1Out halfMap2Out halfMapsTag[=_rlnRandomSubset] (export to two different files)\n";
    std::cerr<<"              --g RetainGroupNumber [RetainGroupMetadata=_subset_group]\n";
    std::cerr<<"              --exportTag tagToExport fileToExportTag.csv\n";
    std::cerr<<"              --imodXf imodXfFile.xf =[shift,rotation]\n";
    std::cerr<<"                   (produce an imod xf file with translations from X,Y origin, default    )\n";
    std::cerr<<"                   (    --o produces the star file with X,Y origin at zero           )\n";
    std::cerr<<"\n";


    std::cerr<<"       Particle Stack Filename Operations:\n";    
    std::cerr<<"              --stackNumbering \n";
    std::cerr<<"              --stackRename stackFileNameRenamed.mrc\n";
    std::cerr<<"              --sortImageName \n";    
    std::cerr<<"              --simulateStarStack numHomogeneousDistribution simulatedStackName[=simulatedStack.mrcs] \n";    

    std::cerr<<"              --removeImageNamePrefix (useful for importing cryoSparc files)\n";    


    std::cerr<<"       Other operations:\n";
    std::cerr<<"              --multiplyLabelValues labelToMultiply valueLabelToMultiply\n";
    std::cerr<<"              --mergeStarFiles csvListFileToMerge\n";
    std::cerr<<"              --splitStarFiles csvOutListFileToSplit\n";
    std::cerr<<"              --setLabel newValueToSetLabel labelNameToUpdate[=_rlnRandomSubset]\n";
    std::cerr<<"              --extractMissingImages reducedImagesStarFile.star\n";
    std::cerr<<"                    (extract images in the starFileIn.star but missing in reducedImagesStarFile.star)\n";


    std::cerr<<"              --extractCommonImages commonImagesStarFile.star\n";
    std::cerr<<"                    (extract images in the starFileIn.star AND in commonImagesStarFile.star)\n";

    std::cerr<<"              --extractFromSimilarParameters extractFromSimilarParametersTemplateFullFile.star ParametersToCompare=[angles,origin,defocus] outputImageTag=[_janas_originalTag_rlnImageName]\n";
    std::cerr<<"                    (Useful for cryosparc's particle subtraction. It extract images in the starFileIn.star that has similar parameters from the template image, and put the filename of the template image.)\n";

    std::cerr<<"              --updateParameters templateForUpdateStarFile.star parametersToUpdate=[class,angles,euler,psi,origin,ctf,subset]\n";
    std::cerr<<"              --backupImageNameTag  outputImageTag=[_janas_originalTag_rlnImageName]\n";


    std::cerr<<"              --invertTagName tag1_toInvert tag2_toInvert\n";


    std::cerr<<"              --extractFromSimilarParameters extractFromSimilarParametersTemplateFullFile.star ParametersToCompare=[angles,origin,defocus] outputImageTag=[_janas_originalTag_rlnImageName]\n";

    std::cerr<<"              --classSelect classNameToSelect\n";


    std::cerr<<"              --autorefineUpdate starFileWithAutorefineUpdates.star\n";
    std::cerr<<"              --replaceLabel SourceFileForReplacingLabel.star SourceLabelsNameCsvForReplacing [DestinationLabelNameCsvForReplacing]\n";
    std::cerr<<"                      (assume label exists)";
    std::cerr<<"              --averageLabels SourceFilesForAveragingLabels.star SourceLabelnameForAveraging \n";
    std::cerr<<"              --assessEuler gtFileForAssessingLabels.star gtLabelsnameForAssessing \n";
    std::cerr<<"              --updateRandomSubset (random assign 1 or 2 to the label _rlnRandomSubset)\n";
    std::cerr<<"              --labelNumbering labelToNumberAsIndex\n";
    std::cerr<<"              --sortLabel labelToSort\n";
    std::cerr<<"              --eraseColumn columnToDelete\n";
    std::cerr<<"\n";
    std::cerr<<"       Conversion:\n";    
    std::cerr<<"              --vem vemFileOut.vem\n";
    std::cerr<<"              --csv csvFileOut.vem\n";
    std::cerr<<"              --dat datAnglesEFFileOut.dat\n";

    std::cerr<<"\n";    
    std::cerr<<"       Info options:\n";  
    std::cerr<<"              --info\n";    
    std::cerr<<"              --infoEuler\n";     
   std::cerr<<"              --infoDiff starInfoDifferenceFile.star\n";
    std::cerr<<"              --h (help)\n";
    std::cerr<<"              --silent\n";
    std::cerr<<"\n";
    exit(1);
}


// *************************************
//
// retrieveInputParameters
// *************************************
void retrieveInputParameters(inputParametersType * parameters, int argc, char** argv){
	//if ( argc < 2)
	//	usage(argv);


    parameters->starFileIn= NULL;
    parameters->starFileOut= NULL;
    parameters->vemFileOut= NULL;    
    parameters->csvFileOut= NULL;    
    parameters->RetainGroupMetadata= NULL;
    parameters->halfMap1Out= NULL;
    parameters->halfMap2Out= NULL;
    parameters->columnToDelete=NULL;
    parameters->halfMapsTag=NULL;
    parameters->RetainGroupNumber= -1;
    parameters->referenceParticlesStarFile= NULL;
    parameters->binningFactor=1.0f;
    parameters->showInfo=false;
    parameters->showInfoDiff=false;
    parameters->micrographs=false;
    parameters->depthDirMicrograph=0;
    parameters->coord_dir=NULL;
    parameters->coord_suffix=NULL;
    parameters->stackNumbering=false;
    parameters->refinedFileName=NULL;    
    parameters->refinedStackFileName=NULL;

    parameters->micrographsWithCtf=NULL;
    parameters->imodXfFile=NULL;
    parameters->imodXfFileOptions=(char *)"shift";
    parameters->sortImageName=false;

    parameters->templateForUpdateStarFile=NULL;
    parameters->parametersToUpdate=(char *)"angles,origin";


    parameters->extractFromSimilarParametersTemplateFullFile=NULL;
    parameters->ParametersToCompare=(char *)"angles,origin,defocus";
    parameters->outputImageTag=(char *)"_janas_originalTag_rlnImageName";



    parameters->numHomogeneousDistribution=-1;
    parameters->simulatedStackName=(char *)"simulatedStack.mrcs";


    parameters->DoCheckForSimilarImages=false;

    parameters->haveBackupImageNameTag=false;

    parameters->reducedImagesStarFile=NULL;
    parameters->commonImagesStarFile=NULL;
    parameters->starInfoDifferenceFile=NULL;

    parameters->tagToExport=NULL;
    parameters->fileToExportTag=NULL;

    parameters->labelToMultiply=NULL;
    parameters->valueLabelToMultiply=1.0;
    
    parameters->stackFileNameRenamed=NULL;
    parameters->SourceFileForReplacingLabel=NULL;
    parameters->SourceLabelsNameCsvForReplacing=NULL;
    parameters->DestinationLabelNameCsvForReplacing=NULL;

    parameters->newValueToSetLabel=NULL;
    parameters->labelNameToUpdate=NULL;

    parameters->SourceFilesForAveragingLabels=NULL;
    parameters->SourceLabelnameForAveraging=NULL;

    parameters->gtFileForAssessingLabels=NULL;
    parameters->gtLabelsnameForAssessing=NULL;

    parameters->csvListFileToMerge=NULL;
    parameters->csvOutListFileToSplit=NULL;


    parameters->labelToNumberAsIndex=NULL;
    parameters->labelToSort=NULL;


    parameters->updateRandomSubset=false;

    parameters->starFileWithAutorefineUpdates=NULL;

    parameters->datAnglesEFFileOut=NULL;


    parameters->tag1_toInvert=NULL;
    parameters->tag2_toInvert=NULL;

    parameters->verboseOn=true;




    parameters->classNameToSelect=NULL;

    parameters->removeImageNamePrefix=false;


    for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if (!optionStr.compare("--h") ){
            usage(argv);
        }
    }


  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--info") ){
          parameters->showInfo=true;
        }
  }    


  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--silent") ){
          parameters->verboseOn=false;
        }
  }    


    for (unsigned int ii=1;ii<argc;ii++){
          std::string optionStr(argv[ii]);
          if ( !optionStr.compare("--updateRandomSubset") ){
            parameters->updateRandomSubset=true;
          }
    }


  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--removeImageNamePrefix") ){
          parameters->removeImageNamePrefix=true;
        }
  }    

  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--stackNumbering") ){
          parameters->stackNumbering=true;
        }
  }  

  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--checkForSimilarImages") ){
          parameters->DoCheckForSimilarImages=true;
        }
  }  



  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--infoDiff") ){
          parameters->showInfoDiff=true;
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->starInfoDifferenceFile=argv[jj];
          }
          ii+=idx;
        }
  }


  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--simulateStarStack") ){
          parameters->showInfoDiff=true;
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->numHomogeneousDistribution=atoi(argv[jj]);
            if (idx==1) parameters->simulatedStackName=argv[jj];
          }
          ii+=idx;
        }
  }


  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( optionStr.compare("--sortImageName")==0 ){
          parameters->sortImageName=true;
        }
  }    


  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--labelNumbering") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->labelToNumberAsIndex=argv[jj];
          }
          ii+=idx;
        }
  }


  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--sortLabel") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->labelToSort=argv[jj];
          }
          ii+=idx;
        }
  }


  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--invertTagName") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->tag1_toInvert=argv[jj];
            if (idx==1) parameters->tag2_toInvert=argv[jj];
          }
          ii+=idx;
        }
  }

  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--dat") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->datAnglesEFFileOut=argv[jj];
          }
          ii+=idx;
        }
  }

  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--i") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->starFileIn=argv[jj];
          }
          ii+=idx;
        }
  }


  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--eraseColumn") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->columnToDelete=argv[jj];
          }
          ii+=idx;
        }
  }
    
  for (unsigned int ii=1;ii<argc;ii++){
      std::string optionStr(argv[ii]);
      if ( !optionStr.compare("--exportTag") ){
        int idx=0;
        for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
          std::string subparamStr(argv[jj]);
          if (!subparamStr.substr(0,2).compare("--")) break;
          if (idx==0) parameters->tagToExport=argv[jj];
          if (idx==1) parameters->fileToExportTag=argv[jj];
        }
        ii+=idx;
      }
  }


  for (unsigned int ii=1;ii<argc;ii++){
      std::string optionStr(argv[ii]);
      if ( !optionStr.compare("--classSelect") ){
        int idx=0;
        for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
          std::string subparamStr(argv[jj]);
          if (!subparamStr.substr(0,2).compare("--")) break;
          if (idx==0) parameters->classNameToSelect=argv[jj];
        }
        ii+=idx;
      }
  }



      for (unsigned int ii=1;ii<argc;ii++){
          std::string optionStr(argv[ii]);
          if ( !optionStr.compare("--multiplyLabelValues") ){
            int idx=0;
            for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
              std::string subparamStr(argv[jj]);
              if (!subparamStr.substr(0,2).compare("--")) break;
              if (idx==0) parameters->labelToMultiply=argv[jj];
              if (idx==1) parameters->valueLabelToMultiply=atof(argv[jj]);
            }
            ii+=idx;
          }
    }
    
      for (unsigned int ii=1;ii<argc;ii++){
          std::string optionStr(argv[ii]);
          if ( !optionStr.compare("--mergeStarFiles") ){
            int idx=0;
            for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
              std::string subparamStr(argv[jj]);
              if (!subparamStr.substr(0,2).compare("--")) break;
              if (idx==0) parameters->csvListFileToMerge=argv[jj];
            }
            ii+=idx;
          }
    }
      for (unsigned int ii=1;ii<argc;ii++){
          std::string optionStr(argv[ii]);
          if ( !optionStr.compare("--splitStarFiles") ){
            int idx=0;
            for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
              std::string subparamStr(argv[jj]);
              if (!subparamStr.substr(0,2).compare("--")) break;
              if (idx==0) parameters->csvOutListFileToSplit=argv[jj];
            }
            ii+=idx;
          }
    }


      for (unsigned int ii=1;ii<argc;ii++){
          std::string optionStr(argv[ii]);
          if ( !optionStr.compare("--setLabel") ){
            int idx=0;
            for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
              std::string subparamStr(argv[jj]);
              if (!subparamStr.substr(0,2).compare("--")) break;
              if (idx==0) parameters->newValueToSetLabel=argv[jj];
              if (idx==1) parameters->labelNameToUpdate=argv[jj];
            }
            ii+=idx;
          }
    }




  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--imodXf") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->imodXfFile=argv[jj];
            if (idx==1) parameters->imodXfFileOptions=argv[jj];
          }
          ii+=idx;
        }
  }

  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--extractMissingImages") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->reducedImagesStarFile=argv[jj];
          }
          ii+=idx;
        }
  }


  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--extractCommonImages") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->commonImagesStarFile=argv[jj];
          }
          ii+=idx;
        }
  }



  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--importCtf") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->micrographsWithCtf=argv[jj];
          }
          ii+=idx;
        }
  }






    for (unsigned int ii=1;ii<argc;ii++){
          std::string optionStr(argv[ii]);
          if ( !optionStr.compare("--stackRename") ){
            int idx=0;
            for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
              std::string subparamStr(argv[jj]);
              if (!subparamStr.substr(0,2).compare("--")) break;
              if (idx==0) parameters->stackFileNameRenamed=argv[jj];
            }
            ii+=idx;
          }
    }
    
    
    
  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--retrieveRefined") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->refinedFileName=argv[jj];
            if (idx==1) parameters->refinedStackFileName=argv[jj];
            
          }
          ii+=idx;
        }
  }



  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--updateParameters") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->templateForUpdateStarFile=argv[jj];
            if (idx==1) parameters->parametersToUpdate=argv[jj];
          }
          ii+=idx;
        }
  }


  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--extractFromSimilarParameters") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->extractFromSimilarParametersTemplateFullFile=argv[jj];
            if (idx==1) parameters->ParametersToCompare=argv[jj];
            if (idx==2) parameters->outputImageTag=argv[jj];
          }
          ii+=idx;
        }
  }



  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--backupImageNameTag") ){
          int idx=0;
          parameters->haveBackupImageNameTag=true;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->outputImageTag=argv[jj];
          }
          ii+=idx;
        }
  }


  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--autorefineUpdate") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->starFileWithAutorefineUpdates=argv[jj];
          }
          ii+=idx;
        }
  }


  for (unsigned int ii=1;ii<argc;ii++){
          std::string optionStr(argv[ii]);
          if ( !optionStr.compare("--replaceLabel") ){
            int idx=0;
            for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
              std::string subparamStr(argv[jj]);
              if (!subparamStr.substr(0,2).compare("--")) break;
              if (idx==0) parameters->SourceFileForReplacingLabel=argv[jj];
              if (idx==1) parameters->SourceLabelsNameCsvForReplacing=argv[jj];
              if (idx==2) parameters->DestinationLabelNameCsvForReplacing=argv[jj];
            }
            ii+=idx;
          }
    }

    for (unsigned int ii=1;ii<argc;ii++){
          std::string optionStr(argv[ii]);
          if ( !optionStr.compare("--averageLabels") ){
            int idx=0;
            for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
              std::string subparamStr(argv[jj]);
              if (!subparamStr.substr(0,2).compare("--")) break;
              if (idx==0) parameters->SourceFilesForAveragingLabels=argv[jj];
              if (idx==1) parameters->SourceLabelnameForAveraging=argv[jj];
            }
            ii+=idx;
          }
    }

    
    for (unsigned int ii=1;ii<argc;ii++){
          std::string optionStr(argv[ii]);
          if ( !optionStr.compare("--assessEuler") ){
            int idx=0;
            for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
              std::string subparamStr(argv[jj]);
              if (!subparamStr.substr(0,2).compare("--")) break;
              if (idx==0) parameters->gtFileForAssessingLabels=argv[jj];
              if (idx==1) parameters->gtLabelsnameForAssessing=argv[jj];
            }
            ii+=idx;
          }
    }


    
  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--rp") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->referenceParticlesStarFile=argv[jj];
            if (idx==1) parameters->binningFactor=atof(argv[jj]);
          }
          ii+=idx;
        }
  }
  
  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--coordinates") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->coord_dir=argv[jj];
            if (idx==1) parameters->coord_suffix=argv[jj];
            if (idx==2) parameters->depthDirMicrograph=atoi(argv[jj]);
          }
          ii+=idx;
        }
  }
      
  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--o") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->starFileOut=argv[jj];
          }
          ii+=idx;
        }
  }

    
    

  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--hm") ){
          parameters->halfMapsTag=(char *)"_rlnRandomSubset";
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->halfMap1Out=argv[jj];
            if (idx==1) parameters->halfMap2Out=argv[jj];            
            if (idx==2) parameters->halfMapsTag=argv[jj];
          }
          ii+=idx;
        }
  }

    
  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--micrographs") ){
          parameters->micrographs=true;
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->depthDirMicrograph=atoi(argv[jj]);
          }
          ii+=idx;

        }
  }


  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--g") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->RetainGroupNumber=atoi(argv[jj]);
            if (idx==1) parameters->RetainGroupMetadata=argv[jj];
          }
          ii+=idx;
        }
  }

  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--vem") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->vemFileOut=argv[jj];
          }
          ii+=idx;
        }
  }

  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--csv") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->csvFileOut=argv[jj];
          }
          ii+=idx;
        }
  }

    //if ( !parameters->starFileIn){
    //    usage(argv);
    //}
}


// **********************************
//
//    INT MAIN
//
// **********************************
int main( int argc, char **argv ){
//std::cerr<<"euler analysis (c) Mauro Maiorca\n";

	inputParametersType parameters;
	retrieveInputParameters(&parameters, argc, argv);






 
 if (parameters.starFileIn){
   
   char * starFileIn=parameters.starFileIn;
   
   
   if (parameters.micrographs && parameters.starFileOut){
     if (parameters.verboseOn) std::cerr<<"processing micrographs\n";
	   std::vector<std::string> rlnImageNameVector;
	   std::vector<std::string> rlnMicrographNameVector;
           readStar(rlnImageNameVector, "_rlnImageName", starFileIn);
           readStar(rlnMicrographNameVector, "_rlnMicrographName", starFileIn);
	   
	   //std::cerr<<"_rlnImageName ==>"<<rlnImageNameVector.size()<<"\n";
	   //std::cerr<<"_rlnMicrographName ==>"<<rlnMicrographNameVector.size()<<"\n";	   
	   std::vector<std::string> uniqueDataNames;
	   std::vector<unsigned long int> uniqueDataCounter;
	   
       unsigned long int size = rlnMicrographNameVector.size();
       std::vector<std::string> startData;
       if (size==0){
          size = rlnImageNameVector.size();
          startData=reduceDirNamefiles (rlnImageNameVector, parameters.depthDirMicrograph);
       }else{
         size = rlnMicrographNameVector.size();
         startData= reduceDirNamefiles(rlnMicrographNameVector, parameters.depthDirMicrograph);
       }
       
       for (unsigned long int ii=0; ii<size; ii++){
            //std::cerr<<startData[ii]<<"\n";
            std::string DataName;
            if(rlnMicrographNameVector.size()==0){
		  std::size_t start1 = rlnImageNameVector[ii].find("@");
		  std::string itemNo = rlnImageNameVector[ii].substr(0, start1);
		  //DataName = rlnImageNameVector[ii].substr(start1+1);
		  DataName = startData[ii].substr(start1+1);
	    }else{
	      //DataName=rlnMicrographNameVector[ii];
	      DataName = startData[ii];
	    } 
            bool found = false;
            //std::cerr<<DataName<<"\n";
	    for(unsigned long int kk=0; kk<uniqueDataNames.size() && !found; kk++){
			  if ( !DataName.compare(uniqueDataNames[kk]) ){
				found = true;
				uniqueDataCounter[kk]++;
			  }
		  }
		  if (!found){
			  uniqueDataNames.push_back(DataName);
			  uniqueDataCounter.push_back(1);
		  }
           }

           std::ofstream fileOut;
           fileOut.open(parameters.starFileOut);
           fileOut.close();
           fileOut.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);    
           fileOut<<"\ndata_\n\nloop_\n";
           fileOut<<"_rlnMicrographName #1\n";

           for(unsigned long int kk=0; kk<uniqueDataNames.size(); kk++){
	                   fileOut<<uniqueDataNames[kk]<<"\n";
           }
           fileOut.close();
   }

     // ###############################
     // ###############################
     // ##      export column           #
     if ( parameters.tagToExport && parameters.fileToExportTag ){
      if (parameters.verboseOn) std::cerr<<"export Column from a specific tag"<<"\n";
         std::vector<std::string> rlnTagToExportVector;
         readStar(rlnTagToExportVector, parameters.tagToExport, starFileIn);
         replaceAddCvsColumn(parameters.tagToExport, rlnTagToExportVector, parameters.fileToExportTag);
     }
     
     
   // ###############################
   // ###############################
   // ##      COORDINATES           #
   if ( parameters.coord_dir && parameters.coord_suffix ){
    if (parameters.verboseOn) std::cerr<<"processing coordinates\n";
	   std::vector<std::string> rlnImageNameVector;
	   readStar(rlnImageNameVector, "_rlnImageName", starFileIn);
	   std::vector<bool> doneParticle(rlnImageNameVector.size(), false);
           std::vector<std::string> rlnCoordinateX;
           std::vector<std::string> rlnCoordinateY;
           readStar(rlnCoordinateX, "_rlnCoordinateX", starFileIn);
           readStar(rlnCoordinateY, "_rlnCoordinateY", starFileIn);
//           for (unsigned long int ii=0; ii<rlnImageNameVector.size(); ii++){
//           std::cerr<<rlnCoordinateY[ii]<<"\n";
//           }
           
           std::vector<std::string> startData=reduceDirNamefiles (rlnImageNameVector, parameters.depthDirMicrograph);
           std::string command (std::string("mkdir -p ") + std::string(parameters.coord_dir));
           int resultCommand=system( command.c_str() );

       //std::vector<std::string> subdirs;
       std::string directoriesList("");
       for (unsigned long int ii=0; ii<rlnImageNameVector.size(); ii++){
            if (!doneParticle[ii]){
               std::string startDataTmp=startData[ii];
   	       
               std::string subdir=startDataTmp.substr(startDataTmp.find("@")+1);
               subdir=subdir.substr(0,subdir.find_last_of('/'));
               

               //std::cerr<<startData[ii]<< "  ->" << subdir << "\n";
               std::string dirToCreate (std::string(parameters.coord_dir) + std::string("/") + subdir);
               long int position=directoriesList.find(dirToCreate);
               if (position<0 || position >= directoriesList.size() ){
   	        std::string command1 (std::string(std::string("mkdir -p ") + dirToCreate));
                resultCommand=system( command1.c_str() );
                directoriesList += std::string(dirToCreate);
                std::cerr<<command1<<"\n";
               }
               
	       std::string DataName = startDataTmp.substr(startDataTmp.find_last_of('/')+1);
	       DataName=DataName.substr(0,DataName.find_last_of('.'));
	       std::string filenameOut=dirToCreate+std::string("/")+DataName+std::string(parameters.coord_suffix);
	       //std::cerr<<startData[ii]<< "  ->" << subdir << " => "<< filenameOut << "\n";
	       //

               std::ofstream fileOut;
               fileOut.open(filenameOut);
               fileOut.close();
               fileOut.open (filenameOut, std::ofstream::out | std::ofstream::app);    
               fileOut<<"\ndata_\n\nloop_\n";
               fileOut<<"_rlnCoordinateX #1\n_rlnCoordinateY #2\n";               
               for(unsigned long int kk=ii; kk<rlnImageNameVector.size(); kk++){
   	         std::string DataName1=startData[kk];
   	         //std::size_t start1 = DataName1.find_last_of('/');
	         DataName1 = DataName1.substr(DataName1.find_last_of('/')+1);
	         DataName1=DataName1.substr(0, DataName1.find_last_of('.') );
	         //std::cerr<< DataName << "   => " << DataName1 <<" \n ";// << rlnCoordinateX[kk]<<"   "<<rlnCoordinateX[kk]<<"\n";
                 if (!doneParticle[kk] && DataName.compare(DataName1) == 0 ){
                    //std::cerr<< DataName << "   => " << DataName1 << " " << rlnCoordinateX[kk]<<"   "<<rlnCoordinateX[kk]<<"\n";
                    fileOut<<rlnCoordinateX[kk]<<"   "<<rlnCoordinateY[kk]<<"\n";
                    doneParticle[kk]=true;
                 }
               }
               fileOut.close();
            }
        }

            /*bool found = false;
            //std::cerr<<DataName<<"\n";
	    for(unsigned long int kk=0; kk<uniqueDataNames.size() && !found; kk++){
		  if ( !DataName.compare(uniqueDataNames[kk]) ){
				found = true;
			  }
		  }
		  if (!found){
			  uniqueDataNames.push_back(DataName);
		  }
           }*/

   }

   // ########################################################
   // ########################################################
   // ##  sortImageName                                     #   
   if ( parameters.sortImageName && parameters.starFileOut){
     if (parameters.verboseOn) std::cerr<<"sort imageName\n";
       //acquire the filenames, and take the indexes
       std::cerr<<"sort imageName in " << starFileIn <<" ... ";
       std::vector<std::string> rlnImageNameVector;
       readStar(rlnImageNameVector, "_rlnImageName", starFileIn);
       std::vector<long int> sortedIdxVector;
       std::vector<long int> originalIdxnameVector;
       for (long int ii = 0; ii < rlnImageNameVector.size(); ii++){
        std::string targetNumber=rlnImageNameVector[ii].substr(0,rlnImageNameVector[ii].find("@"));
        //sortedIdxVector.push_back(ii);
        originalIdxnameVector.push_back(std::stol( targetNumber )-1);
        //std::cerr<<originalIdxVector[ii]<<"\n";
       }
       sortedIdxVector=bubbleSortAscendingIndexes(originalIdxnameVector);
       std::vector<std::string> linesByLines (sortedIdxVector.size(), "");
       //std::vector<std::string> linesByLines;


        std::ofstream fileOutput;
        fileOutput.open(parameters.starFileOut);
        fileOutput.close();
        fileOutput.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);   
        long int startMicrograph=getStarStart(starFileIn);
        std::ifstream fileParticles(starFileIn);
        std::string strLine;
        for (int counter=0;counter<startMicrograph;counter++){
           std::getline(fileParticles, strLine);
           fileOutput << strLine <<"\n";
        }
        for (unsigned long int ii=0, counter=0; ii<rlnImageNameVector.size(); ii++){
            std::getline(fileParticles, strLine);
            std::string str = strLine;
            std::regex rNoSpaces("\\s+");
            std::string strNoSpaces = std::regex_replace(str, rNoSpaces, "");    
            if (strNoSpaces.length()>0){
              if (sortedIdxVector[ii]>=0 && sortedIdxVector[ii]<linesByLines.size()){
                linesByLines[counter]=strLine;
                //linesByLines.push_back(strLine);
                counter++;
              }
            }
        }
        //std::cerr<<"\n\n*******************\n";
        for (unsigned long int ii=0; ii<linesByLines.size(); ii++){
             if (sortedIdxVector[ii]>=0 && sortedIdxVector[ii]<linesByLines.size()){
        
              //std::cerr<<ii<< " ===> " <<sortedIdxVector[ii]<<"\n";
              fileOutput << linesByLines[sortedIdxVector[ii]] <<"\n";
             }
        }
        fileOutput.close();


      std::cerr<<"   [DONE]\n";
   }


   // ########################################################
   // ########################################################
   // ##  labelNumbering                                          #   
   if ( parameters.labelToNumberAsIndex && parameters.starFileOut){
    //std::cerr<<"              --labelNumbering labelToNumberAsIndex\n";
         if (parameters.verboseOn) std::cerr<<"from file" << starFileIn <<"numbering this label: " << parameters.labelToNumberAsIndex << "\n";
      unsigned long int starFileSize=getNumItemsStar(starFileIn);
      std::vector<std::string> outIdx;
      for (unsigned long int gg=0; gg<starFileSize; gg++){
        outIdx.push_back(std::to_string(gg+1));
      }
      replaceAddValueStar(parameters.labelToNumberAsIndex, outIdx, starFileIn, parameters.starFileOut);

   }
   // ########################################################
   // ########################################################
   // ##  sortLabel                                          #   
   if ( parameters.labelToSort && parameters.starFileOut){
       //acquire the filenames, and take the indexes
       //std::string labelToSort(parameters.labelToSort);
       //if (parameters.labelToSort){
       // labelToSort=parameters.labelToSort;
       //}
//--sortLabel labelToSort
//      std::cerr<< "star file version="<<checkStarFileVersion (starFileIn)<<"\n";


       if (parameters.verboseOn) std::cerr<<"sort lines in file: " << starFileIn << "  according to label: " << parameters.labelToSort<< " ... ";
       //std::vector<std::string> rlnImageNameVector;

       std::vector<double> valuesToSortVector;
       readStar(valuesToSortVector, parameters.labelToSort, starFileIn);
       std::vector<long int> sortedIdxVector=bubbleSortAscendingIndexes(valuesToSortVector);
       std::vector<std::string> linesByLines (sortedIdxVector.size(), "");


  std::string startHeader = getStarHeader ( starFileIn);
  int countLines=std::count(startHeader.begin(), startHeader.end(), '\n');
  //std::cerr<<"\n\nstartMicrograph="<<countLines<<"\n";



        std::ofstream fileOutput;
        fileOutput.open(parameters.starFileOut);
        fileOutput.close();
        fileOutput.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);   
        long int startMicrograph=getStarStart(starFileIn);
        std::ifstream fileParticles(starFileIn);
        std::string strLine;
        for (int counter=0;counter<startMicrograph;counter++){
           std::getline(fileParticles, strLine);
           fileOutput << strLine <<"\n";
        }
        
//        std::cerr<<"\n";
        for (unsigned long int ii=0, counter=0; ii<valuesToSortVector.size(); ii++){
            //std::cerr<<" idx="<<valuesToSortVector[sortedIdxVector[ii]]<<"\n";
            std::getline(fileParticles, strLine);
            //std::cerr<<"ii="<< ii<<"  content="<<strLine<<"\n"; 
            std::string str = strLine;
            std::regex rNoSpaces("\\s+");
            std::string strNoSpaces = std::regex_replace(str, rNoSpaces, "");    
            if (strNoSpaces.length()>0){
                if (sortedIdxVector[ii]>=0 && sortedIdxVector[ii]<linesByLines.size()){
                  linesByLines[counter]=strLine;
                  //linesByLines.push_back(strLine);
                  counter++;
                }
                //std::cerr<<"ii="<< ii<<"  content="<<strLine<<"\n"; 
            }/*else{
               std::cerr<<"ECCHIMI: ii="<< ii<<"  content="<<strLine<<"\n"; 
            }*/
        }
        //std::cerr<<"\n\n*******************\n";
        for (unsigned long int ii=0; ii<linesByLines.size(); ii++){
             if (sortedIdxVector[ii]>=0 && sortedIdxVector[ii]<linesByLines.size()){
              //std::cerr<<ii<< " ===> " <<sortedIdxVector[ii]<<"\n";
              fileOutput << linesByLines[sortedIdxVector[ii]] <<"\n";
             }
        }
        fileOutput.close();


      std::cerr<<"   [DONE]\n";
   }


   // ########################################################
   // ########################################################
   // ##  converting dat for cryoEF                          #   
   if ( parameters.datAnglesEFFileOut ){
      if (parameters.verboseOn) std::cerr<<"convert start to dat (useful for cryoEF, for example)";
    	std::vector<double> rotAnglesVector;
      std::vector<double> tiltAnglesVector;
      readStar(rotAnglesVector, "_rlnAngleRot", starFileIn);
      readStar(tiltAnglesVector, "_rlnAngleTilt", starFileIn);
         std::ofstream fileOutput;
          fileOutput.open(parameters.datAnglesEFFileOut);
          fileOutput.close();
          fileOutput.open (parameters.datAnglesEFFileOut, std::ofstream::out | std::ofstream::app);
      for(unsigned long int uuu=0; uuu<rotAnglesVector.size(); uuu++){
        fileOutput << rotAnglesVector[uuu] <<" "<< tiltAnglesVector[uuu]<<"\n";
      }
      fileOutput.close();
   }

   // ########################################################
   // ########################################################
   // ##  stackNumbering                                     #   
   if ( parameters.stackNumbering && parameters.starFileOut){
    if (parameters.verboseOn) std::cerr<<"add numbering the filename for a specific stack\n";
	std::vector<std::string> rlnImageNameVector;
	readStar(rlnImageNameVector, "_rlnImageName", starFileIn);
	
	//get the file start
	//for each line get where the filename image start and ends
	//replace with the new characters
        std::ofstream fileOutput;
        fileOutput.open(parameters.starFileOut);
        fileOutput.close();
        fileOutput.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);   
        long int startMicrograph=getStarStart(starFileIn);
        std::ifstream fileParticles(starFileIn);
        std::string strLine;
        for (int counter=0;counter<startMicrograph;counter++){
           std::getline(fileParticles, strLine);
           fileOutput << strLine <<"\n";
        }
        std::string previousFileName(rlnImageNameVector[0].substr(rlnImageNameVector[0].find("@")+1));;
        for (unsigned long int ii=0,counter=1;ii<rlnImageNameVector.size();ii++){
            std::size_t start1 = rlnImageNameVector[ii].find("@");
            std::string currentFilename(rlnImageNameVector[ii].substr(start1+1));
            if ( previousFileName != currentFilename ){
              previousFileName = currentFilename;
              counter=1;
            }

            //std::cerr<<ii<<"  "<<counter<<"  previousFileName="<<previousFileName<< "   currentFileName="<<currentFilename<<"  \n";
            std::getline(fileParticles, strLine);
            int start=strLine.find(rlnImageNameVector[ii]);
            int end=start+rlnImageNameVector[ii].size();
            //std::cerr<<start<<"   ->"<<  rlnImageNameVector[ii]  <<" || "<< strLine<< "\n";
            
            if (start >=0  && start<strLine.size()){
               std::string tmpInitialString = strLine.substr(0,start);
               while(!tmpInitialString.empty() && std::isspace(*tmpInitialString.begin()))
                  tmpInitialString.erase(tmpInitialString.begin());  
            
              fileOutput << tmpInitialString;
              if (start>0){
                fileOutput << "   ";
              }
              fileOutput << std::setfill('0') << std::fixed << std::right << std::setw(8) << std::showpoint << std::setprecision(6)  <<  counter;
              fileOutput << "@" << rlnImageNameVector[ii].substr(start1+1);
              fileOutput << "   " << strLine.substr(end+1);
              fileOutput << "\n";
              

              counter++;

            }
        }
        fileOutput.close();
   }



     // ########################################################
     // ########################################################
     // ##  removeImageNamePrefix                                     #
     if ( parameters.removeImageNamePrefix && parameters.starFileOut){
      if (parameters.verboseOn) std::cerr<<"remove Image Name Prefix\n";

      std::vector<std::string> rlnImageNameVector;
      readStar(rlnImageNameVector, "_rlnImageName", starFileIn);
         //std::cerr<<rlnImageNameVector[0]<<"\n";
      //get the file start
      //for each line get where the filename image start and ends
      //replace with the new characters
         int labelIdx=getStarHeaderItemIdx("_rlnImageName", starFileIn);
         //std::cerr<< "labelIdx="<<labelIdx<<"\n";
         std::ofstream fileOutput;
          fileOutput.open(parameters.starFileOut);
          fileOutput.close();
          fileOutput.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);
          long int startMicrograph=getStarStart(starFileIn);
          std::ifstream fileParticles(starFileIn);
          std::string strLine;
          for (int counter=0;counter<startMicrograph;counter++){
             std::getline(fileParticles, strLine);
             fileOutput << strLine <<"\n";
          }
         for (unsigned long int ii=0,counter=1;ii<rlnImageNameVector.size();ii++){
             std::getline(fileParticles, strLine);
             std::string inputStr(strLine);

            //remove head spaces
            //strLine=strLine.erase(0, strLine.find_first_not_of(" \t\n\r\f\v"));
            //strLine = std::regex_replace(strLine, std::regex("^ +"), "");











             std::vector<long int> BeginAndLenght=getStringBeginAndLenghtAtStarPosition(strLine, labelIdx);
             std::string targetField(inputStr.substr(BeginAndLenght[0],BeginAndLenght[0]+BeginAndLenght[1]));
             std::string outputStr(inputStr.substr(0,BeginAndLenght[0]));
             if (labelIdx>0) outputStr+=std::string (" ");


            //replace
            size_t at = targetField.find('@');
            std::string second_part = targetField.substr(at + 1);
            size_t slash = second_part.rfind('/');
            size_t initial = at+slash + 2;
            std::string without_prefix = (slash != std::string::npos) ? second_part.substr(slash + 1) : second_part;
            size_t last_underscore = without_prefix.find_first_of("_");
            if (last_underscore != std::string::npos) {
                without_prefix = targetField.substr(0,initial)+without_prefix.substr(last_underscore + 1);
            }
            outputStr+=without_prefix;
             
             
             
             
             outputStr+=inputStr.substr(BeginAndLenght[0]+BeginAndLenght[1]);
             fileOutput << outputStr << "\n";
             counter++;
          }
         fileOutput.close();
     }




     // ########################################################
     // ########################################################
     // ##  rename stack                                     #
     // std::cerr<<"              --stackRename stackFileNameRenamed.mrc\n";
     if ( parameters.stackFileNameRenamed && parameters.starFileOut){
      if (parameters.verboseOn) std::cerr<<"rename Image Name to a certain stack\n";

      std::vector<std::string> rlnImageNameVector;
      readStar(rlnImageNameVector, "_rlnImageName", starFileIn);
         std::cerr<<rlnImageNameVector[0]<<"\n";
      //get the file start
      //for each line get where the filename image start and ends
      //replace with the new characters
         int labelIdx=getStarHeaderItemIdx("_rlnImageName", starFileIn);
         std::cerr<< "labelIdx="<<labelIdx<<"\n";
         std::ofstream fileOutput;
          fileOutput.open(parameters.starFileOut);
          fileOutput.close();
          fileOutput.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);
          long int startMicrograph=getStarStart(starFileIn);
          std::ifstream fileParticles(starFileIn);
          std::string strLine;
          for (int counter=0;counter<startMicrograph;counter++){
             std::getline(fileParticles, strLine);
             fileOutput << strLine <<"\n";
          }
         for (unsigned long int ii=0,counter=1;ii<rlnImageNameVector.size();ii++){
             std::getline(fileParticles, strLine);
             std::string inputStr(strLine);

            //remove head spaces
            //strLine=strLine.erase(0, strLine.find_first_not_of(" \t\n\r\f\v"));
            //strLine = std::regex_replace(strLine, std::regex("^ +"), "");




             std::vector<long int> BeginAndLenght=getStringBeginAndLenghtAtStarPosition(strLine, labelIdx);
             std::string targetField(inputStr.substr(BeginAndLenght[0],BeginAndLenght[0]+BeginAndLenght[1]));
             std::string outputStr(inputStr.substr(0,BeginAndLenght[0]));
             if (labelIdx>0) outputStr+=std::string (" ");
             outputStr+=targetField.substr(0,targetField.find("@")+1)+std::string (parameters.stackFileNameRenamed);
             outputStr+=inputStr.substr(BeginAndLenght[0]+BeginAndLenght[1]);
             fileOutput << outputStr << "\n";
             counter++;
          }
         fileOutput.close();
     }

     
     // ########################################################
     // ########################################################
     // ##  multiplyLabelValues                                     #
     // ##  std::cerr<<"              --multiplyLabelValues labelToMultiply valueLabelToMultiply\n";
     if ( parameters.labelToMultiply && parameters.starFileOut){
      if (parameters.verboseOn) std::cerr<<"multiply labels\n";
         std::vector<std::string> rlnToChangeLabelVector;
         readStar(rlnToChangeLabelVector, parameters.labelToMultiply, starFileIn);
         replaceAddValueStar(parameters.labelToMultiply, rlnToChangeLabelVector, starFileIn, parameters.starFileOut);
     }

     // ########################################################
     // ########################################################
     // ##  mergeStarFiles                                     #
     // ##  std::cerr<<"              --mergeStarFiles csvListFileToMerge\n";
     if ( parameters.csvListFileToMerge && parameters.starFileOut){
      if (parameters.verboseOn) std::cerr<<"merge star files\n";

         //std::vector<std::string> rlnToChangeLabelVector;
         //readStar(rlnToChangeLabelVector, parameters.labelToMultiply, starFileIn);
         //replaceAddValueStar(parameters.labelToMultiply, rlnToChangeLabelVector, starFileIn, parameters.starFileOut);
     
         mergeStarFiles(parameters.csvListFileToMerge, parameters.starFileOut);

     
     }
    
     // ########################################################
     // ########################################################
     // ##  splitStarFiles                                     #
     // ##  std::cerr<<"              --splitStarFiles csvOutListFileToSplit\n";
     if ( parameters.csvOutListFileToSplit){
          if (parameters.verboseOn) std::cerr<<"split star files\n";
          splitStarFiles(starFileIn, parameters.csvOutListFileToSplit);
     }

     // ########################################################
     // ########################################################
     // ##  setLabel                                     #
     // ##      std::cerr<<"              --setLabel newValueToSetLabel labelNameToUpdate[=_rlnRandomSubset]\n";
     if ( starFileIn && parameters.newValueToSetLabel){
        if (parameters.verboseOn) std::cerr<<"set new values to a certain label\n";

        unsigned long int numItems = countParticlesItems( starFileIn );
        std::vector<std::string> outputValuesVector;
        for (unsigned long int ii=0; ii< numItems; ii++){
            outputValuesVector.push_back( parameters.newValueToSetLabel );
        }
        if (parameters.starFileOut){
             replaceAddValueStar("_rlnRandomSubset", outputValuesVector, starFileIn, parameters.starFileOut);
        }else{
             replaceAddValueStar("_rlnRandomSubset", outputValuesVector, starFileIn, ".___TMP_FILE____.tmp");
             copyCvsFile(".___TMP_FILE____.tmp", starFileIn);
             removeCvsFile(".___TMP_FILE____.tmp");
        }

     }



   // ########################################################
   // ########################################################
   // ##  extractMissingImages                                     #  
   if ( parameters.reducedImagesStarFile && parameters.starFileOut){
        if (parameters.verboseOn) std::cerr<<"extract Missing Images from the template file\n";

	std::vector<std::string> rlnImageNameVector;
	readStar(rlnImageNameVector, "_rlnImageName", starFileIn);
	
	std::vector<std::string> rlnImageNameVectorReduced;
	readStar(rlnImageNameVectorReduced, "_rlnImageName", parameters.reducedImagesStarFile);
//	std::cerr<<"\nECCHIMI\n\n";
	
//    std::cerr<<"              --extractMissingImages reducedImagesStarFile.star\n";
	if (parameters.verboseOn) std::cerr<<"full size star="<<rlnImageNameVector.size()<<"\n";
	if (parameters.verboseOn) std::cerr<<"templete star file size="<<rlnImageNameVectorReduced.size()<<"\n";  
	
	
	
        std::vector<int> lookupExisting (rlnImageNameVector.size(), 1);        
        std::vector<unsigned long int> imageNumberReduced (rlnImageNameVectorReduced.size() );
        std::vector<std::string> stackNameReduced (rlnImageNameVectorReduced.size() );
	for (unsigned long int ii=0; ii<rlnImageNameVectorReduced.size(); ii++){
		  std::size_t start1 = rlnImageNameVectorReduced[ii].find("@");
		  imageNumberReduced[ii]=std::stol(rlnImageNameVectorReduced[ii].substr(0, start1));
		  stackNameReduced[ii]=rlnImageNameVectorReduced[ii].substr(start1+1);
	}
        

        

       
        
        /*
        for (unsigned long int ii=0; ii<rlnImageNameVectorReduced.size(); ii++){
          std::string referenceStr=rlnImageNameVectorReduced[ii].substr(0,rlnImageNameVectorReduced[ii].find("@"));
          unsigned long int referenceNo=std::stol(referenceStr)-1;
          //std::cerr<<referenceNo<<"\n";
          lookupExisting[referenceNo]=0;
        }*/
        for (unsigned long int ii=0; ii<rlnImageNameVector.size(); ii++){
	        std::size_t start1 = rlnImageNameVector[ii].find("@");
	        unsigned long int targetNumber = std::stol(rlnImageNameVector[ii].substr(0, start1));
	        std::string targetStackName ( rlnImageNameVector[ii].substr(start1+1) );
        	bool found=false;
	        for (unsigned long int kk=0; kk<rlnImageNameVectorReduced.size() && !found; kk++){
	        	if (targetNumber==imageNumberReduced[kk]){
	        		if(targetStackName==stackNameReduced[kk]){
	        			found=true;
	        			lookupExisting[ii]=0;
	        		}
	        	}
	        }
        }




	//get the file start
	//for each line get where the filename image start and ends
	//replace with the new characters
	
        std::ofstream fileOutput;
        fileOutput.open(parameters.starFileOut);
        fileOutput.close();
        fileOutput.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);   
        long int startMicrograph=getStarStart(starFileIn);
        std::ifstream fileParticles(starFileIn);
        std::string strLine;
        for (int counter=0;counter<startMicrograph;counter++){
           std::getline(fileParticles, strLine);
           fileOutput << strLine <<"\n";
        }
        bool foundFirst=false;
        
        for (unsigned long int ii=0;ii<rlnImageNameVector.size();ii++){
            std::getline(fileParticles, strLine);
            int start=strLine.find(rlnImageNameVector[ii]);
            int end=start+rlnImageNameVector[ii].size();
            //std::cerr<<start<<"   ->"<<  rlnImageNameVector[ii]  <<" || "<< strLine<< "\n";
            
            if (start >=0  && start<strLine.size()){
              if (lookupExisting[ii]>0){
                 if (!foundFirst){                 
                    foundFirst=true;
                 }else{
                    fileOutput <<"\n";
                 }
                 fileOutput << strLine;
              }
            }
        }
        
        fileOutput.close();
   }


// ########################################################
// ########################################################
// ##  extractCommonImages                                     #  
if ( parameters.commonImagesStarFile && parameters.starFileOut){
        if (parameters.verboseOn) std::cerr<<"extract common Images from two star files\n";

	std::vector<std::string> rlnImageNameVector;
	readStar(rlnImageNameVector, "_rlnImageName", starFileIn);
	
	std::vector<std::string> rlnImageNameVectorReduced;
	readStar(rlnImageNameVectorReduced, "_rlnImageName", parameters.commonImagesStarFile);
//	std::cerr<<"\nECCHIMI\n\n";
	
//    std::cerr<<"              --extractMissingImages commonImagesStarFile.star\n";
	std::cerr<<"rlnImageNameVector.size="<<rlnImageNameVector.size()<<"\n";
	std::cerr<<"rlnImageNameVectorReduced.size="<<rlnImageNameVectorReduced.size()<<"\n";  
	
	
	
  std::vector<int> lookupExisting (rlnImageNameVector.size(), 0);        
  std::vector<unsigned long int> imageNumberReduced (rlnImageNameVectorReduced.size() );
  std::vector<std::string> stackNameReduced (rlnImageNameVectorReduced.size() );
	for (unsigned long int ii=0; ii<rlnImageNameVectorReduced.size(); ii++){
		  std::size_t start1 = rlnImageNameVectorReduced[ii].find("@");
		  imageNumberReduced[ii]=std::stol(rlnImageNameVectorReduced[ii].substr(0, start1));
		  stackNameReduced[ii]=rlnImageNameVectorReduced[ii].substr(start1+1);
	}
        

        

        for (unsigned long int ii=0; ii<rlnImageNameVector.size(); ii++){
	        std::size_t start1 = rlnImageNameVector[ii].find("@");
	        unsigned long int targetNumber = std::stol(rlnImageNameVector[ii].substr(0, start1));
	        std::string targetStackName ( rlnImageNameVector[ii].substr(start1+1) );
        	bool found=false;
	        for (unsigned long int kk=0; kk<rlnImageNameVectorReduced.size() && !found; kk++){
	        	if (targetNumber==imageNumberReduced[kk]){
	        		if(targetStackName==stackNameReduced[kk]){
	        			found=true;
	        			lookupExisting[ii]=1;
	        		}
	        	}
	        }
        }




	//get the file start
	//for each line get where the filename image start and ends
	//replace with the new characters
	
        std::ofstream fileOutput;
        fileOutput.open(parameters.starFileOut);
        fileOutput.close();
        fileOutput.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);   
        long int startMicrograph=getStarStart(starFileIn);
        std::ifstream fileParticles(starFileIn);
        std::string strLine;
        for (int counter=0;counter<startMicrograph;counter++){
           std::getline(fileParticles, strLine);
           fileOutput << strLine <<"\n";
        }
        bool foundFirst=false;
        
        for (unsigned long int ii=0;ii<rlnImageNameVector.size();ii++){
            std::getline(fileParticles, strLine);
            int start=strLine.find(rlnImageNameVector[ii]);
            int end=start+rlnImageNameVector[ii].size();
            //std::cerr<<start<<"   ->"<<  rlnImageNameVector[ii]  <<" || "<< strLine<< "\n";
            
            if (start >=0  && start<strLine.size()){
              if (lookupExisting[ii]>0){
                 if (!foundFirst){                 
                    foundFirst=true;
                 }else{
                    fileOutput <<"\n";
                 }
                 fileOutput << strLine;
              }
            }
        }
        
        fileOutput.close();
        
   }



// ########################################################
// ########################################################
// ##  classSelect                                     #  
if ( parameters.classNameToSelect && parameters.starFileOut){
        if (parameters.verboseOn) std::cerr<<"select particles from a specific class\n";

	std::vector<std::string> classNameVector;
	readStar(classNameVector, "_rlnClassNumber", starFileIn);
	

	//get the file start
	//for each line get where the filename image start and ends
	//replace with the new characters
	
        std::ofstream fileOutput;
        fileOutput.open(parameters.starFileOut);
        fileOutput.close();
        fileOutput.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);   
        long int startMicrograph=getStarStart(starFileIn);
        std::ifstream fileParticles(starFileIn);
        std::string strLine;
        for (int counter=0;counter<startMicrograph;counter++){
           std::getline(fileParticles, strLine);
           fileOutput << strLine <<"\n";
        }
        bool foundFirst=false;
        
        for (unsigned long int ii=0;ii<classNameVector.size();ii++){
            std::getline(fileParticles, strLine);
            if ( (classNameVector[ii])==(parameters.classNameToSelect)){
              fileOutput << strLine << std::endl;
            }
        }
        fileOutput.close();
   }




   // ########################################################
   // ########################################################
   if ( parameters.imodXfFile){
        if (parameters.verboseOn) std::cerr<<"export parameters to imodXfFile\n";

        std::string inputParameters(parameters.imodXfFileOptions);

	std::vector<double> rlnCoordXVector;
	std::vector<double> rlnCoordYVector;
	std::vector<double> rlnPsi;
	bool doShift = false;
	bool doRotation = false;	


        if ( inputParameters.find("shift")!=std::string::npos  ){
	 readStar(rlnCoordXVector, "_rlnOriginX", starFileIn);
	 readStar(rlnCoordYVector, "_rlnOriginY", starFileIn);
	 doShift=true;
	}
	if (inputParameters.find("rotation")!=std::string::npos){
	 readStar(rlnPsi, "_rlnOriginY", starFileIn);
	 doRotation=true;
        }



        std::ofstream xfFileOutput;
        xfFileOutput.open(parameters.imodXfFile);
        xfFileOutput.close();
        xfFileOutput.open (parameters.imodXfFile, std::ofstream::out | std::ofstream::app);
        for (unsigned long int ii=0;ii<rlnCoordXVector.size();ii++){
            if (doRotation){
                    double angle=rlnPsi[ii]*PI/180.0;
                    double sinAngle=sin(angle), cosAngle=cos(angle);
                    double A=cosAngle, B=-sinAngle, C=sinAngle, D=cosAngle;
                    xfFileOutput << std::setfill('0') << std::fixed << std::right << std::setw(8) << std::showpoint << std::setprecision(6) <<  A << "   " << B  <<"  "<<  C << "   " << D <<"   ";
            }else{
                    xfFileOutput << "1.0000   0.00000  0.00000  1.0000   ";
            }
            if (doShift){
               xfFileOutput << std::setfill('0') << std::fixed << std::right << std::setw(8) << std::showpoint << std::setprecision(6)  <<  rlnCoordXVector[ii] << "  " << rlnCoordYVector[ii] <<"\n";
            }else{
               xfFileOutput << " 0.00000 0.00000\n";
            }
        }
        xfFileOutput.close();

        if (parameters.starFileOut){

	        int XIdx=getStarHeaderItemIdx("_rlnOriginX", starFileIn);
	        int YIdx=getStarHeaderItemIdx("_rlnOriginY", starFileIn);
	        std::cerr<<"XIdx="<<XIdx<<"    YIdx="<<YIdx<<"\n";
                int idx1 = XIdx;
                int idx2 = YIdx;
                if (idx2<idx1){
                  int tmp = idx2;
                  idx2 = idx1;
                  idx1 = tmp;
                }
                std::ofstream fileOutput;
                fileOutput.open(parameters.starFileOut);
                fileOutput.close();
                fileOutput.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);   
                long int startStar=getStarStart(starFileIn);
                std::ifstream fileParticles(starFileIn);
                std::string strLine;
                for (int counter=0;counter<startStar;counter++){
                   std::getline(fileParticles, strLine);
                   fileOutput << strLine <<"\n";
                }
                
                for (unsigned long int ii=0; ii<rlnCoordXVector.size();ii++){
                  std::getline(fileParticles, strLine);
                    //fileOutput << strLine <<"\n";
                    std::vector<long int> posIdx1=getStringBeginAndLenghtAtStarPosition(strLine, idx1);
                    std::vector<long int> posIdx2=getStringBeginAndLenghtAtStarPosition(strLine, idx2);
                    fileOutput << strLine.substr(0,posIdx1[0])<<"   ";
                    fileOutput << " 0.000000 0.000000  ";
                    fileOutput << strLine.substr(posIdx2[0]+posIdx2[1]);
                    
                    //fileOutput << strLine.substr(posIdx[0]+posIdx[1],posIdx2[0])<<"  0.02000  ";
                    //fileOutput << strLine.substr(posIdx2[0]+posIdx2[1]);
                    fileOutput << "\n";
                }
                fileOutput.close();
        
        }
   }
   // ###################################################################
   // ###################################################################
   // ##      create stack from star file                               #   
   if ( parameters.micrographsWithCtf && parameters.starFileOut){
        if (parameters.verboseOn) std::cerr<<"create stack from star file\n";


// FROM CTF FILE
/* 
   _rlnMicrographName #1 
_rlnCtfImage #2 
_rlnDefocusU #3 
_rlnDefocusV #4 
_rlnCtfAstigmatism #5 
_rlnDefocusAngle #6 
_rlnVoltage #7 
_rlnSphericalAberration #8 
_rlnAmplitudeContrast #9 
_rlnMagnification #10 
_rlnDetectorPixelSize #11 
_rlnCtfFigureOfMerit #12 
_rlnCtfMaxResolution #13 
*/
// FROM SHINY FILE
std::string cpFileMicrograph = std::string(std::string("relion_star_handler --remove_column rlnVoltage --i ")+ std::string(parameters.micrographsWithCtf) + std::string(" --o _tmp_sourceCTF.star"));
std::string cpFileParticles = std::string(std::string("relion_star_handler  --remove_column rlnDefocusU --i ")+ std::string(starFileIn) + std::string(" --o _tmp_particlesCTF.star"));
int outVar=system(cpFileMicrograph.c_str());
outVar=system("relion_star_handler --i  _tmp_sourceCTF.star --o  _tmp_sourceCTF.star --remove_column rlnCtfImage" );
outVar=system("relion_star_handler --i  _tmp_sourceCTF.star --o  _tmp_sourceCTF.star --remove_column rlnCtfMaxResolution" );

outVar=system(cpFileParticles.c_str());
outVar=system("relion_star_handler --i  _tmp_particlesCTF.star --o  _tmp_particlesCTF.star --remove_column rlnDefocusV" );
outVar=system("relion_star_handler --i  _tmp_particlesCTF.star --o  _tmp_particlesCTF.star --remove_column rlnCtfAstigmatism" );
outVar=system("relion_star_handler --i  _tmp_particlesCTF.star --o  _tmp_particlesCTF.star --remove_column rlnDefocusAngle" );
outVar=system("relion_star_handler --i  _tmp_particlesCTF.star --o  _tmp_particlesCTF.star --remove_column rlnSphericalAberration" );
outVar=system("relion_star_handler --i  _tmp_particlesCTF.star --o  _tmp_particlesCTF.star --remove_column rlnAmplitudeContrast" );
outVar=system("relion_star_handler --i  _tmp_particlesCTF.star --o  _tmp_particlesCTF.star --remove_column rlnMagnification" );
outVar=system("relion_star_handler --i  _tmp_particlesCTF.star --o  _tmp_particlesCTF.star --remove_column rlnDetectorPixelSize" );
outVar=system("relion_star_handler --i  _tmp_particlesCTF.star --o  _tmp_particlesCTF.star --remove_column rlnCtfFigureOfMerit" );
outVar=system("relion_star_handler --i  _tmp_particlesCTF.star --o  _tmp_particlesCTF.star --remove_column rlnOriginalParticleName" );


std::vector<std::string> MicrographNameSource;
readStar(MicrographNameSource, "_rlnMicrographName", "_tmp_sourceCTF.star");

std::vector<std::string> MicrographNameParticles;
readStar(MicrographNameParticles, "_rlnMicrographName", "_tmp_particlesCTF.star");

//simplifyListNamefiles (std::vector< std::string > inputList)
MicrographNameSource=reduceDirNamefiles (MicrographNameSource, 1);
MicrographNameParticles=reduceDirNamefiles (MicrographNameParticles, 1);
MicrographNameSource=simplifyListNamefiles (MicrographNameSource);
MicrographNameParticles=simplifyListNamefiles (MicrographNameParticles);

//std::vector<long int> getStringBeginAndLenghtAtStarPosition(std::string line,"_tmp_particlesCTF.star" int position){

//std::cerr<<"source    => "<<MicrographNameSource[0]<<"\n";
//std::cerr<<"particles => "<<MicrographNameParticles[0]<<"\n";

        std::vector<std::string> micrographCtfBlock=extractBlockLabels("_rlnDefocusU", "_rlnCtfFigureOfMerit", "_tmp_sourceCTF.star");

        std::ofstream fileOutput;
        fileOutput.open(parameters.starFileOut);
        fileOutput.close();
        fileOutput.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);   
        long int startMicrograph=getStarStart("_tmp_particlesCTF.star");
        std::ifstream fileParticles("_tmp_particlesCTF.star");
        std::string strLine;
        for (int counter=0;counter<startMicrograph;counter++){
           std::getline(fileParticles, strLine);
           fileOutput << strLine <<"\n";
        }
        int maxField=1+StarMaxFieldIdx("_tmp_particlesCTF.star");
        fileOutput << "_rlnDefocusU   #" << maxField++ <<"\n";
        fileOutput << "_rlnDefocusV   #" << maxField++ <<"\n";
        fileOutput << "_rlnCtfAstigmatism   #" << maxField++ <<"\n";
        fileOutput << "_rlnDefocusAngle   #" << maxField++ <<"\n";
        fileOutput << "_rlnSphericalAberration   #" << maxField++ <<"\n";
        fileOutput << "_rlnAmplitudeContrast   #" << maxField++ <<"\n";
        fileOutput << "_rlnMagnification   #" << maxField++ <<"\n";
        fileOutput << "_rlnDetectorPixelSize   #" << maxField++ <<"\n";
        fileOutput << "_rlnCtfFigureOfMerit   #" << maxField++ <<"\n";

        //std::string currentMicrograph=MicrographNameSource[0];
        unsigned long int currentMicrographIdx=0;
        for (unsigned long int ii=0;ii<MicrographNameParticles.size();ii++){
           std::getline(fileParticles, strLine);
           std::string stripStr = std::regex_replace(strLine, std::regex("\\s+"), "");
           if(stripStr.size()>0){
            fileOutput << strLine <<"  ";
            if (MicrographNameParticles[ii].compare(MicrographNameSource[currentMicrographIdx])==0){
               fileOutput << micrographCtfBlock[currentMicrographIdx];
            }else{
            bool found = false;
              for (unsigned long int jj=0;jj<MicrographNameSource.size() && ! found;jj++){
                 if (MicrographNameParticles[ii].compare(MicrographNameSource[jj])==0){
                   found=true;
                   currentMicrographIdx=jj;
                   fileOutput << micrographCtfBlock[currentMicrographIdx];
                 }
              }
            }
            fileOutput << "\n";
           }
        }



/* 
_rlnVoltage #1 
_rlnDefocusU #2 
_rlnDefocusV #3 
_rlnDefocusAngle #4 
_rlnSphericalAberration #5 
_rlnDetectorPixelSize #6 
_rlnCtfFigureOfMerit #7 
_rlnMagnification #8 
_rlnAmplitudeContrast #9 
_rlnImageName #10 

std::regex r("\\s+");
       std::string str0 = std::regex_replace(strLine, r, ",");
_rlnMicrographName #14 
*/
   
   } 
   // ###################################################################
   // ###################################################################
   // ##  update informations of input file from target file       
   if ( parameters.starFileWithAutorefineUpdates && parameters.starFileOut){
        //std::cerr<<"              --autorefineUpdate starFileWithAutorefineUpdates.star\n";

        const int original_relionStarFileVersion = checkStarFileVersion (starFileIn);
        const int update_relionStarFileVersion   = checkStarFileVersion (parameters.starFileWithAutorefineUpdates);
        if (original_relionStarFileVersion != update_relionStarFileVersion){
          std::cerr<<"ERROR: different star files versions... exiting\n ";
          exit(1);
        }


        std::vector<std::string> original_rlnImageNameVector;
        std::vector<std::string> update_rlnImageNameVector;
        readStar(original_rlnImageNameVector, "_rlnImageName", starFileIn);
        readStar(update_rlnImageNameVector, "_rlnImageName", parameters.starFileWithAutorefineUpdates);

        std::vector<double> update_ImageNumber;
        std::vector<std::string> update_ImageName;
        std::vector<bool> alreadyUpdatedItem;
        for (int ii=0;ii<update_rlnImageNameVector.size(); ii++){
            std::string targetFile=update_rlnImageNameVector[ii].substr(update_rlnImageNameVector[ii].find("@")+1);
            long int numItem=std::stol (update_rlnImageNameVector[ii].substr(0,update_rlnImageNameVector[ii].find("@"))) - 1;
            update_ImageNumber.push_back(numItem+1);
            update_ImageName.push_back(targetFile);
            alreadyUpdatedItem.push_back(false);
            //std::cerr<< update_ImageNumber [ii]<<"  -> "<<update_ImageName[ii] <<"\n";
        } 


//std::cerr<<"\nORIGINAL=\n";
        std::vector<double> original_ImageNumber;
        std::vector<std::string> original_ImageName;
        for (int ii=0;ii<original_rlnImageNameVector.size(); ii++){
            std::string targetFile=original_rlnImageNameVector[ii].substr(original_rlnImageNameVector[ii].find("@")+1);
            long int numItem=std::stol (original_rlnImageNameVector[ii].substr(0,original_rlnImageNameVector[ii].find("@"))) - 1;
            original_ImageNumber.push_back(numItem+1);
            original_ImageName.push_back(targetFile);
            //std::cerr<< original_ImageNumber [ii]<<"  -> "<<original_ImageName[ii] <<"\n";
        } 



        std::ofstream fileOutput;
        fileOutput.open(parameters.starFileOut);
        fileOutput.close();
        fileOutput.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);   
        long int startFileIn=getStarStart(starFileIn);
        long int startFileTemplate=getStarStart(parameters.starFileWithAutorefineUpdates);
        std::ifstream inputFileStream(starFileIn);
        std::ifstream templateFileStream(parameters.starFileWithAutorefineUpdates);
        std::string strLineIn;
        for (int counter=0;counter<startFileIn;counter++){
           std::getline(inputFileStream, strLineIn);
           fileOutput << strLineIn <<"\n";
        }
        std::string strLineTemplate;
        for (int counter=0;counter<startFileTemplate;counter++){
           std::getline(templateFileStream, strLineTemplate);
        }


        
        std::vector<std::string> templateFileLineByLine;
        while ( std::getline(templateFileStream, strLineTemplate)){
              std::string strCheck=strLineTemplate;
              if (strCheck.erase(0, strCheck.find_first_not_of(" \t\n\r\f\v")).size() != std::string::npos) { //check the line is not empty
                templateFileLineByLine.push_back(strLineTemplate);
              }
        }

        //std::cerr<<" filling output:\n";
        //const long int newItemInputStar = getNumItemsStar(starFileIn);

        long int templateIdx=-1;
        long int currentLine = 0;
        while (  std::getline(inputFileStream, strLineIn) ){

           std::string strCheck=strLineIn;
           strCheck=strCheck.erase(0, strCheck.find_first_not_of(" \t\n\r\f\v"));
           long int targetLine=0;
          //std::cerr<<"    ===  "<< strLineIn <<" \n";

          
           if (strCheck.size()>1) { //if okay
               //if(strLineIn.find_first_not_of(' ') != std::string::npos){//check the line is not empty
                bool found = false;
                for (unsigned long kk=0; kk<update_ImageNumber.size() && !found; kk++){
                    if (!alreadyUpdatedItem[kk]){
                      if (original_ImageNumber[currentLine]==update_ImageNumber[kk]){
                        if( original_ImageName[currentLine].compare(update_ImageName[kk])==0 ){ //if those strings are equal
                            found=true;
                            alreadyUpdatedItem[kk]=true;
                            targetLine=kk;
                        }

                      }
                    }
                }
                if(found){
                        fileOutput << templateFileLineByLine[targetLine] <<"\n";
                }else{
                        //just place the original image
                        fileOutput << strLineIn <<"\n";
                }
                currentLine++;
            }

               
           }
           fileOutput.close();





   }

   // ###################################################################
   // ###################################################################
   // ##  update informations of input file from target file       
   if ( parameters.templateForUpdateStarFile && parameters.starFileOut){
        if (parameters.verboseOn) std::cerr<<"update informations of input file from target file\n";


     //std::cerr<<"              --updateParameters templateForUpdateStarFile.star parametersToUpade=[angles,origin,ctf]\n";
        std::string inputParameters(parameters.parametersToUpdate);
       std::vector <int> listTagIdxInput;
       std::vector <int> listTagIdxTemplate;
       const int relionStarFileVersion=checkStarFileVersion (parameters.templateForUpdateStarFile);


       if ( inputParameters.find("euler")!=std::string::npos  ){
         listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnAngleRot",starFileIn));
         listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnAngleTilt",starFileIn));
         listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnAngleRot",parameters.templateForUpdateStarFile));
         listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnAngleTilt",parameters.templateForUpdateStarFile));
       }
       if ( inputParameters.find("class")!=std::string::npos ){
         listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnClassNumber",starFileIn));
         listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnClassNumber",parameters.templateForUpdateStarFile));
       }
       if ( inputParameters.find("subset")!=std::string::npos ){
         listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnRandomSubset",starFileIn));
         listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnRandomSubset",parameters.templateForUpdateStarFile));
       }       

       if ( inputParameters.find("psi")!=std::string::npos  ){
         listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnAnglePsi",starFileIn));
         listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnAnglePsi",parameters.templateForUpdateStarFile));
       }

       if ( inputParameters.find("angles")!=std::string::npos  ){
         listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnAngleRot",starFileIn));
         listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnAngleTilt",starFileIn));
         listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnAnglePsi",starFileIn));
         listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnAngleRot",parameters.templateForUpdateStarFile));
         listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnAngleTilt",parameters.templateForUpdateStarFile));
         listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnAnglePsi",parameters.templateForUpdateStarFile));
       }
       if ( inputParameters.find("origin")!=std::string::npos ){
           if (relionStarFileVersion==3100){
               listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnOriginXAngst",starFileIn));
               listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnOriginYAngst",starFileIn));
               listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnOriginXAngst",parameters.templateForUpdateStarFile));
               listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnOriginYAngst",parameters.templateForUpdateStarFile));
           }else{
             listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnOriginX",starFileIn));
             listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnOriginY",starFileIn));
             listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnOriginX",parameters.templateForUpdateStarFile));
             listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnOriginY",parameters.templateForUpdateStarFile));
           }
       }
       if ( inputParameters.find("ctf")!=std::string::npos ){
         listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnSphericalAberration",starFileIn));
         listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnAmplitudeContrast",starFileIn));
         listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnCtfFigureOfMerit",starFileIn));
         listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnDefocusU",starFileIn));
         listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnDefocusV",starFileIn));
         listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnDefocusAngle",starFileIn));
         listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnCtfMaxResolution",starFileIn));
         listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnSphericalAberration",parameters.templateForUpdateStarFile));
         listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnAmplitudeContrast",parameters.templateForUpdateStarFile));
         listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnCtfFigureOfMerit",parameters.templateForUpdateStarFile));
         listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnDefocusU",parameters.templateForUpdateStarFile));
         listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnDefocusV",parameters.templateForUpdateStarFile));
         listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnDefocusAngle",parameters.templateForUpdateStarFile));
         listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnCtfMaxResolution",parameters.templateForUpdateStarFile));
       }


       std::vector<std::string> intputImageTagVector;
       readStar(intputImageTagVector, "_rlnImageName", starFileIn);

       
       
       std::vector<long int> intputImageNumberVector;
       for (long int ii=0; ii<intputImageTagVector.size(); ii++){
            int startIdx=intputImageTagVector[ii].find("@");
            intputImageNumberVector.push_back(std::stol(intputImageTagVector[ii].substr(0,startIdx)));
            intputImageTagVector[ii]=intputImageTagVector[ii].substr(startIdx+1);
       }

       std::vector<std::string> templateImageTagVector;
       readStar(templateImageTagVector, "_rlnImageName", parameters.templateForUpdateStarFile);
       std::vector<long int> templateImageNumberVector;
       for (long int ii=0; ii<templateImageTagVector.size(); ii++){
            int startIdx=templateImageTagVector[ii].find("@");
            templateImageNumberVector.push_back(std::stol(templateImageTagVector[ii].substr(0,startIdx)));
            templateImageTagVector[ii]=templateImageTagVector[ii].substr(startIdx+1);
       }

       //load all the lines of the inputFileIn memory


        std::ofstream fileOutput;
        fileOutput.open(parameters.starFileOut);
        fileOutput.close();
        fileOutput.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);   
        long int startFileIn=getStarStart(starFileIn);
        long int startFileTemplate=getStarStart(parameters.templateForUpdateStarFile);
        std::ifstream inputFileStream(starFileIn);
        std::ifstream templateFileStream(parameters.templateForUpdateStarFile);
        std::string strLineIn;
        for (int counter=0;counter<startFileIn;counter++){
           std::getline(inputFileStream, strLineIn);
           fileOutput << strLineIn <<"\n";
        }
        std::string strLineTemplate;
        for (int counter=0;counter<startFileTemplate;counter++){
           std::getline(templateFileStream, strLineTemplate);
        }

        std::vector<std::string> templateFileLineByLine;
        for (long int ii=0; ii<templateImageTagVector.size();ii++){
          if (std::getline(templateFileStream, strLineTemplate)){
              std::string strCheck=strLineTemplate;
              if (strCheck.erase(0, strCheck.find_first_not_of(" \t\n\r\f\v")).size() != std::string::npos) { //check the line is not empty
                templateFileLineByLine.push_back(strLineTemplate);
              }
          }
        }

        //std::cerr<<"qui boh?\n";

        long int templateIdx=-1;
        for (unsigned long int ii=0; ii<intputImageTagVector.size();ii++){
           std::getline(inputFileStream, strLineIn);
           std::string strCheck=strLineIn;
           strCheck=strCheck.erase(0, strCheck.find_first_not_of(" \t\n\r\f\v"));
           if (strCheck.size()>1) {      
               //if(strLineIn.find_first_not_of(' ') != std::string::npos){//check the line is not empty
                bool found = false;
                for (unsigned long kk=templateIdx+1; kk<templateImageTagVector.size() && !found; kk++){
                  //std::cerr<<kk<<", ";
                  if ( templateImageNumberVector[kk]==intputImageNumberVector[ii] ){
                        if ( templateImageTagVector[kk].compare(intputImageTagVector[ii])==0  ){ //if those strings are equal
                          templateIdx=kk;
                          found=true;
                        }
                  }
                }
                for (unsigned long kk=0; kk<templateIdx && !found && kk<templateImageTagVector.size(); kk++){
                  //std::cerr<<kk<<", ";
                  if ( templateImageNumberVector[kk]==intputImageNumberVector[ii] ){
                        if ( templateImageTagVector[kk].compare(intputImageTagVector[ii])==0  ){ //if those strings are equal
                          templateIdx=kk;
                          found=true;
                        }
                  }
                }
                //std::cerr<<templateIdx <<"  ";
                if(found){
                        std::string inputStr=strLineIn;
                        for (int hhh=0;hhh<listTagIdxInput.size(); hhh++){
                                std::vector<long int> templatePos=getStringBeginAndLenghtAtStarPosition(templateFileLineByLine[templateIdx], listTagIdxTemplate[hhh]);
                                std::string valueToReplace = templateFileLineByLine[templateIdx].substr(templatePos[0],templatePos[1]);
                                //std::cerr<<listTagIdxTemplate[hhh]<<"->("<<valueToReplace<<") ";
                                inputStr=replaceValueStrlineStarFile(inputStr, listTagIdxInput[hhh], valueToReplace);
                        }
                        fileOutput << inputStr <<"\n";
                }else{
                        //just place the original image
                        fileOutput << strLineIn <<"\n";
                }
           }

        }
   }



   // ###################################################################
   // ###################################################################
   // ##  It extract images in the starFileIn.star that has similar parameters from the template       
   // ##  Useful for cryosparc's particle subtraction. 
if ( parameters.extractFromSimilarParametersTemplateFullFile && parameters.starFileOut){
        if (parameters.verboseOn) std::cerr<<"extract images in the starFileIn.star that has similar parameters from target file\n";

      int maxHeaderLabel=getMaxHeaderLabelNumber (parameters.extractFromSimilarParametersTemplateFullFile);
      std::vector<std::string> rlnImageNameVector;
      readStar(rlnImageNameVector, "_rlnImageName", starFileIn);
       std::string ParametersToCompare(parameters.ParametersToCompare);
       std::vector <int> listTagIdxInput;
       std::vector <int> listTagIdxTemplate;
       const int relionStarFileVersion=checkStarFileVersion (parameters.extractFromSimilarParametersTemplateFullFile);


       if ( ParametersToCompare.find("angles")!=std::string::npos  ){
         listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnAngleRot",starFileIn));
         listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnAngleTilt",starFileIn));
         listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnAnglePsi",starFileIn));
         listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnAngleRot",parameters.extractFromSimilarParametersTemplateFullFile));
         listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnAngleTilt",parameters.extractFromSimilarParametersTemplateFullFile));
         listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnAnglePsi",parameters.extractFromSimilarParametersTemplateFullFile));
       }
       if ( ParametersToCompare.find("origin")!=std::string::npos ){
           if (relionStarFileVersion==3100){
               listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnOriginXAngst",starFileIn));
               listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnOriginYAngst",starFileIn));
               listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnOriginXAngst",parameters.extractFromSimilarParametersTemplateFullFile));
               listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnOriginYAngst",parameters.extractFromSimilarParametersTemplateFullFile));
           }else{
             listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnOriginX",starFileIn));
             listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnOriginY",starFileIn));
             listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnOriginX",parameters.extractFromSimilarParametersTemplateFullFile));
             listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnOriginY",parameters.extractFromSimilarParametersTemplateFullFile));
           }
       }
       if ( ParametersToCompare.find("defocus")!=std::string::npos ){
         listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnDefocusU",starFileIn));
         listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnDefocusV",starFileIn));
         listTagIdxInput.push_back(getStarHeaderItemIdx("_rlnDefocusAngle",starFileIn));
         listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnDefocusU",parameters.extractFromSimilarParametersTemplateFullFile));
         listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnDefocusV",parameters.extractFromSimilarParametersTemplateFullFile));
         listTagIdxTemplate.push_back(getStarHeaderItemIdx("_rlnDefocusAngle",parameters.extractFromSimilarParametersTemplateFullFile));
       }
       std::vector<std::string> intputImageTagVector;
       readStar(intputImageTagVector, "_rlnImageName", starFileIn);

       std::vector<long int> intputImageNumberVector;
       for (long int ii=0; ii<intputImageTagVector.size(); ii++){
            int startIdx=intputImageTagVector[ii].find("@");
            intputImageNumberVector.push_back(std::stol(intputImageTagVector[ii].substr(0,startIdx)));
            intputImageTagVector[ii]=intputImageTagVector[ii].substr(startIdx+1);
       }

       std::vector<std::string> templateImageTagVector;
       readStar(templateImageTagVector, "_rlnImageName", parameters.extractFromSimilarParametersTemplateFullFile);
       std::vector<long int> templateImageNumberVector;
       for (long int ii=0; ii<templateImageTagVector.size(); ii++){
            int startIdx=templateImageTagVector[ii].find("@");
            templateImageNumberVector.push_back(std::stol(templateImageTagVector[ii].substr(0,startIdx)));
            templateImageTagVector[ii]=templateImageTagVector[ii].substr(startIdx+1);
       }
       //load all the lines of the inputFileIn memory


        std::ofstream fileOutput;
        fileOutput.open(parameters.starFileOut);
        fileOutput.close();
        fileOutput.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);   
        long int startFileIn=getStarStart(starFileIn);
        long int startFileTemplate=getStarStart(parameters.extractFromSimilarParametersTemplateFullFile);
        std::ifstream inputFileStream(starFileIn);
        std::ifstream templateFileStream(parameters.extractFromSimilarParametersTemplateFullFile);
        std::string strLineIn;
        for (int counter=0;counter<startFileIn;counter++){
           std::getline(inputFileStream, strLineIn);
           fileOutput << strLineIn <<"\n";
        }
        std::string strLineTemplate;
        for (int counter=0;counter<startFileTemplate;counter++){
           std::getline(templateFileStream, strLineTemplate);
        }

        std::string newLabel("_janas_backup_rlnImageName #");
        newLabel+=std::to_string(maxHeaderLabel+1);
        fileOutput << newLabel <<"\n";

        std::vector<std::string> InputFileLineByLine;
        for (long int ii=0; ii<intputImageTagVector.size();ii++){
          if (std::getline(inputFileStream, strLineIn)){
              std::string strCheck=strLineIn;
              if (strCheck.erase(0, strCheck.find_first_not_of(" \t\n\r\f\v")).size() != std::string::npos) { //check the line is not empty
                InputFileLineByLine.push_back(strLineIn);
              }
          }
        }


        std::vector<std::string> rlnResultNameVector(templateImageTagVector.size());
        long int inputImageNameIdx=getStarHeaderItemIdx("_rlnImageName",starFileIn);
        long int templateImageNameIdx=getStarHeaderItemIdx("_rlnImageName",parameters.extractFromSimilarParametersTemplateFullFile);

        std::vector<bool> assignedLine(rlnImageNameVector.size(), false);
        long int inputFileIdx=-1;
        for (unsigned long int ii=0; ii<templateImageTagVector.size();ii++){
          //std::cerr<< ii <<"\n";
           std::cout << "\rProgress: " << (ii*100)/templateImageTagVector.size() << "%" << std::flush;
           std::getline(templateFileStream, strLineIn);
           std::string strCheck=strLineIn;
           strCheck=strCheck.erase(0, strCheck.find_first_not_of(" \t\n\r\f\v"));
           if (strCheck.size()>1) {      
               //if(strLineIn.find_first_not_of(' ') != std::string::npos){//check the line is not empty
                bool found = false;
                //list of values
                std::vector<int> listTemplateValues( listTagIdxTemplate.size() );
                for (int hhh=0;hhh<listTagIdxInput.size(); hhh++){
                    std::vector<long int> templatePos=getStringBeginAndLenghtAtStarPosition(strLineIn, listTagIdxTemplate[hhh]);
                    std::string valueRetrieved = strLineIn.substr(templatePos[0],templatePos[1]);
                    double num1 = std::stod(valueRetrieved);
                    int roundedNum = static_cast<int>(std::round(num1));
                    listTemplateValues[hhh]=roundedNum;
                }

                for (unsigned long int jj=0; jj<rlnImageNameVector.size() && !found ;jj++){
                  if(!assignedLine[jj]){
                      //here the comparison happens:
                      int numHits = 0;
                      for (int hhh=0;hhh<listTagIdxInput.size(); hhh++){
                          std::vector<long int> templatePos=getStringBeginAndLenghtAtStarPosition(InputFileLineByLine[jj], listTagIdxTemplate[hhh]);
                          std::string valueRetrieved = InputFileLineByLine[jj].substr(templatePos[0],templatePos[1]);
                          double num1 = std::stod(valueRetrieved);
                          int roundedNum = static_cast<int>(std::round(num1));
                          if ( listTemplateValues[hhh] == roundedNum) numHits++;
                      }
                      if (numHits==listTagIdxInput.size()){
                        found=true;
                        inputFileIdx=jj;
                      }
                  }
                }
                //std::cerr<<"input str ="<<strLineIn<<"\n";
                //std::cerr<<"target str="<<InputFileLineByLine[inputFileIdx]<<"\n";
                //std::cerr<<"\n\n";
           
                if(found){
                        rlnResultNameVector[ii]=rlnImageNameVector[inputFileIdx];
//                        std::string inputStr=strLineIn;
//                        InputFileLineByLine[inputFileIdx];
//                        std::vector<long int>inputPos=getStringBeginAndLenghtAtStarPosition(InputFileLineByLine[inputFileIdx], templateImageNameIdx);
//                        std::string valueToReplace = templateImageTagVector[ii];
//                        std::string tmp_outputStr=replaceValueStrlineStarFile(inputStr, inputPos, valueToReplace);
                        fileOutput << strLineIn << "   " << rlnImageNameVector[inputFileIdx] <<"\n";
                        assignedLine[inputFileIdx]=true;

                    }else{
                              //don't write anythin
                              //fileOutput << strLineIn <<"\n";
                        rlnResultNameVector[ii]=templateImageTagVector[ii];
                }
           }

          //std::cerr<<"==>"<< templateImageTagVector[ii]<<"  => "  << rlnResultNameVector[ii] <<"\n";
        }
        std::cerr<< "\n";

   }










   // ###################################################################
   // ###################################################################
   // ##  It extract images in the starFileIn.star that has similar parameters from the template       
   // ##  Useful for cryosparc's particle subtraction. 
if ( parameters.DoCheckForSimilarImages && parameters.starFileOut){
      if (parameters.verboseOn) std::cerr<<"Do Check For Similar Images (based on parameters)\n";

      int maxHeaderLabel=getMaxHeaderLabelNumber (parameters.extractFromSimilarParametersTemplateFullFile);
      std::vector<std::string> rlnImageNameVector;
      readStar(rlnImageNameVector, "_rlnImageName", starFileIn);

       std::vector<long> intputImageNumberVector(rlnImageNameVector.size());
       std::vector<std::string> intputImageTagVector(rlnImageNameVector.size());
       for (long int ii=0; ii<rlnImageNameVector.size(); ii++){
            int startIdx=rlnImageNameVector[ii].find("@");
            intputImageNumberVector[ii]=std::stol(rlnImageNameVector[ii].substr(0,startIdx));
            intputImageTagVector[ii]=rlnImageNameVector[ii].substr(startIdx+1);
       }


        std::ofstream fileOutput;
        fileOutput.open(parameters.starFileOut);
        fileOutput.close();
        fileOutput.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);   
        long int startFileIn=getStarStart(starFileIn);
        long int startFileTemplate=getStarStart(parameters.extractFromSimilarParametersTemplateFullFile);
        std::ifstream inputFileStream(starFileIn);
        std::string strLineIn;
        for (int counter=0;counter<startFileIn;counter++){
           std::getline(inputFileStream, strLineIn);
           fileOutput << strLineIn <<"\n";
        }

        std::string newLabel("_janas_similarImage_CC #");
        newLabel+=std::to_string(maxHeaderLabel+1);
        fileOutput << newLabel <<"\n";
        std::string newLabel2("_janas_similarImage_targetImage #");
        newLabel2+=std::to_string(maxHeaderLabel+2);
        fileOutput << newLabel2 <<"\n";

        MRCHeader header;
        std::string firstStackFilename(rlnImageNameVector[0].substr(rlnImageNameVector[0].find("@")+1));        
        readHeaderMrc(firstStackFilename.c_str(), header);

        std::cerr<<" "<<header.nx<<"   "<<header.ny<<"\n";
        unsigned long int nxy=header.nx*header.ny;
        float * image1 = new float [nxy];
        float * image2 = new float [nxy];

        long int inputFileIdx=-1;
        for (unsigned long int ii=0; ii<rlnImageNameVector.size();ii++){
          //std::cerr<< ii <<"\n";
           std::cout << "\rProgress: " << (ii*100)/rlnImageNameVector.size() << "%" << std::flush;
           std::getline(inputFileStream, strLineIn);
           //std::cerr << "image=  " << intputImageNumberVector[ii] << "    " << rlnImageNameVector[ii] <<"\n";

          readMrcSlice(intputImageTagVector[ii].c_str(), image1, header, intputImageNumberVector[ii]-1);

           double maxScore=0;
           unsigned long int targetIdx=ii;
           for (unsigned long int jj=0; jj<rlnImageNameVector.size() ;jj++){
                  if( ii != jj ){
                      readMrcSlice(intputImageTagVector[jj].c_str(), image2, header, intputImageNumberVector[jj]-1);
                      double score=crossCorrelationDistance(image1, image2, nxy);
                      if (score > maxScore){
                        maxScore=score;
                        targetIdx=jj;
                      }
                  }
            }
            fileOutput << strLineIn<< "   " << std::to_string(maxScore) << "   "<< rlnImageNameVector[ii] <<"\n";
            //std::cerr<<"==>"<< templateImageTagVector[ii]<<"  => "  << rlnResultNameVector[ii] <<"\n";
        }
        std::cerr<< "\n";
        delete [] image1;
        delete [] image2;
        return 0;

   }



   // ###################################################################
   // ###################################################################
   // ##  backupImageNameTag. 
if ( parameters.haveBackupImageNameTag && parameters.starFileOut){
      if (parameters.verboseOn) std::cerr<<"backup Image Name Tag\n";

      std::vector<std::string> rlnImageNameVector;
      readStar(rlnImageNameVector, "_rlnImageName", starFileIn);
      int maxHeaderLabel=getMaxHeaderLabelNumber (starFileIn);
        std::ofstream fileOutput;
        fileOutput.open(parameters.starFileOut);
        fileOutput.close();
        fileOutput.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);   
        long int startFileIn=getStarStart(starFileIn);
        long int startFileTemplate=getStarStart(parameters.extractFromSimilarParametersTemplateFullFile);
        std::ifstream inputFileStream(starFileIn);
        std::ifstream templateFileStream(parameters.extractFromSimilarParametersTemplateFullFile);
        std::string strLineIn;
        for (int counter=0;counter<startFileIn;counter++){
           std::getline(inputFileStream, strLineIn);
           fileOutput << strLineIn <<"\n";
        }

        //std::string newLabel("_janas_backup_rlnImageName #");
        std::string newLabel(parameters.outputImageTag);
        newLabel+=std::string(" #");
        newLabel+=std::to_string(maxHeaderLabel+1);
        fileOutput << newLabel <<"\n";
        for (unsigned long int ii=0; ii<rlnImageNameVector.size();ii++){
           std::getline(inputFileStream, strLineIn);
           std::string strCheck=strLineIn;
           fileOutput << strLineIn << "   " << rlnImageNameVector[ii] <<"\n";
        }
           
   }










   // ###################################################################
   // ###################################################################
   // ##  invertTagName. 
if ( parameters.tag1_toInvert && parameters.tag2_toInvert && parameters.starFileOut){
      if (parameters.verboseOn) std::cerr<<"invert Tag Names\n";

    std::string header = getStarHeader ( starFileIn );
    //std::cerr<<"header="<<header;
    std::string tag1_toInvert("\n");
    tag1_toInvert+=std::string(parameters.tag1_toInvert);
    tag1_toInvert+=std::string(" ");
    std::string tag2_toInvert("\n");
    tag2_toInvert+=std::string(parameters.tag2_toInvert);
    tag2_toInvert+=std::string(" ");

    size_t pos1 = header.find(tag1_toInvert);
    size_t pos2 = header.find(tag2_toInvert);

    if (pos1 == std::string::npos ) {
      std::cerr << "ERROR: invertTagName One of the tags was not found in the header." << std::endl;
      return 1;
    }


    std::string headerReverted("");
    std::cerr<<"pos1="<<pos1<<"\n";
    std::cerr<<"pos2="<<pos2<<"\n";
    std::cerr<<"tag1_toInvert=$$$"<<tag1_toInvert<<"$$$\n";
    std::cerr<<"tag2_toInvert=$$$"<<tag2_toInvert<<"$$$\n";

    if (pos1 < pos2) {
        headerReverted = header.substr(0, pos1) + tag2_toInvert + header.substr(pos1 + tag1_toInvert.size(), pos2 - pos1 - tag1_toInvert.size()) + tag1_toInvert + header.substr(pos2 + tag2_toInvert.size());
    } else if (pos1 > pos2) {
        headerReverted = header.substr(0, pos2) + tag1_toInvert + header.substr(pos2 + tag2_toInvert.size(), pos1 - pos2 - tag2_toInvert.size()) + tag2_toInvert + header.substr(pos1 + tag1_toInvert.size());
    } else{
    	std::cerr<<"WARNING: one of the two headers is missing, doing nothing;\n";
	headerReverted=header;
    }

    //std::cerr<<"headerReverted="<<headerReverted;
    long int starStart=getStarStart(starFileIn);

    std::istringstream stream(headerReverted);
    std::string line;
    std::vector<std::string> lines;
    // Store each line in the vector and count them
    while (std::getline(stream, line)) {
        lines.push_back(line);
        //std::cerr<<line<<"\n";
    }
    int numOfLines = lines.size();
    std::ifstream inputFileStream(starFileIn);
    std::string str;
    std::ofstream fileOutput;
   
    fileOutput.open(parameters.starFileOut, std::ios::out | std::ios::trunc);
    fileOutput.close();
    fileOutput.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);
    unsigned long int counter = 0;
    
    for (int ii =0 ; ii<starStart; ii++ ){
      std::getline(inputFileStream, str);
    }
    for (int ii =0 ; ii<lines.size(); ii++){
      fileOutput << lines[ii]  <<"\n";
    }
    while (std::getline(inputFileStream, str) ){
           fileOutput << str  <<"\n";
    }
    inputFileStream.close();

   }








     // ###################################################################
     // ###################################################################
     // ##  updateRandomSubset  #
     if ( parameters.updateRandomSubset){
      if (parameters.verboseOn) std::cerr<<"update Random Subset\n";


         unsigned long int numItems = countParticlesItems( starFileIn );
         std::vector<int> randomValuesVector;
         for (unsigned long int ii=0; ii<numItems; ii++){
             randomValuesVector.push_back(rnd_01()+1);
         }
         if (parameters.starFileOut){
             replaceAddValueStar("_rlnRandomSubset", randomValuesVector, starFileIn, parameters.starFileOut);
         }else{
             replaceAddValueStar("_rlnRandomSubset", randomValuesVector, starFileIn, ".___TMP_FILE____.tmp");
             copyCvsFile(".___TMP_FILE____.tmp", starFileIn);
             removeCvsFile(".___TMP_FILE____.tmp");
         }
      }
     
   // ###################################################################
   // ###################################################################
   // ##  retrieve all information from refined file, except file name  #   
   if ( parameters.refinedFileName && parameters.starFileOut){
      if (parameters.verboseOn) std::cerr<<"retrieve all information from refined file, except file name\n";

//     std::cerr<<"              --retrieveRefined refinedFileName.star\n";

	std::vector<std::string> rlnMicrographNameVector;
	readStar(rlnMicrographNameVector, "_rlnMicrographName", starFileIn);

	std::vector<double> rlnMicrographCoordXVector;
	readStar(rlnMicrographCoordXVector, "_rlnCoordinateX", starFileIn);

	std::vector<double> rlnMicrographCoordYVector;
	readStar(rlnMicrographCoordYVector, "_rlnCoordinateY", starFileIn);

	std::vector<std::string> rlnImageNameNameVector;
	readStar(rlnImageNameNameVector, "_rlnImageName", starFileIn);



	std::vector<std::string> rlnRefinedMicrographNameVector;
	readStar(rlnRefinedMicrographNameVector, "_rlnMicrographName", parameters.refinedFileName);

	std::vector<double> rlnRefinedMicrographCoordXVector;
	readStar(rlnRefinedMicrographCoordXVector, "_rlnCoordinateX", parameters.refinedFileName);

	std::vector<double> rlnRefinedMicrographCoordYVector;
	readStar(rlnRefinedMicrographCoordYVector, "_rlnCoordinateY", parameters.refinedFileName);

	int refinedImageNameidx=getStarHeaderItemIdx("_rlnImageName", parameters.refinedFileName);
	//std::cerr<<idx<<"\n";
	
	std::vector<bool> processed (rlnRefinedMicrographNameVector.size(), false);
	std::vector<long int> foundLine (rlnMicrographNameVector.size(), -1);
	
	for (unsigned long int ii=0; ii<rlnMicrographNameVector.size();ii++){
   	   for (unsigned long int jj=0, found=false; jj<rlnRefinedMicrographNameVector.size() && !found ;jj++){
	     if (!processed[jj]){
	       if(rlnMicrographCoordXVector[ii]==rlnRefinedMicrographCoordXVector[jj]){
	         if ( rlnMicrographCoordYVector[ii]==rlnRefinedMicrographCoordYVector[jj]){
	          if(rlnMicrographNameVector[ii].compare(rlnRefinedMicrographNameVector[jj]) == 0 ){
                       processed[jj]=true;
                       foundLine[jj]=ii+1;
                       found=true;
	          }
	         }
	       }
	     }
	   }
	   
	}


        std::ofstream fileOutput;
        fileOutput.open(parameters.starFileOut);
        fileOutput.close();
        fileOutput.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);   
        long int startMicrograph=getStarStart(parameters.refinedFileName);
        std::ifstream fileRefinedParticles(parameters.refinedFileName);
        std::string strLine;
        for (int counter=0;counter<startMicrograph;counter++){
           std::getline(fileRefinedParticles, strLine);
           fileOutput << strLine <<"\n";
        }
        for (unsigned long int ii=0,counter=1;ii<rlnMicrographNameVector.size();ii++){
          std::getline(fileRefinedParticles, strLine);
          if(foundLine[ii]>=0){
            //fileOutput << strLine <<"\n";
            std::vector<long int> posIdx=getStringBeginAndLenghtAtStarPosition(strLine, refinedImageNameidx);
            fileOutput << strLine.substr(0,posIdx[0])<<"   ";
            fileOutput << std::setfill('0') << std::fixed << std::right << std::setw(8) << std::showpoint << std::setprecision(6)  <<  foundLine[ii];
            fileOutput << "@" << parameters.refinedStackFileName << " ";
            fileOutput << strLine.substr(posIdx[0]+posIdx[1]);
            fileOutput << "\n";
          }
        }
        fileOutput.close();
   }   






     // ###################################################################
     // ###################################################################
     // ##  updateRandomSubset  #
     if ( parameters.updateRandomSubset){
      if (parameters.verboseOn) std::cerr<<"update Random Subset\n";

         unsigned long int numItems = countParticlesItems( starFileIn );
         std::vector<int> randomValuesVector;
         for (unsigned long int ii=0; ii<numItems; ii++){
             randomValuesVector.push_back(rnd_01()+1);
         }
         if (parameters.starFileOut){
             replaceAddValueStar("_rlnRandomSubset", randomValuesVector, starFileIn, parameters.starFileOut);
         }else{
             replaceAddValueStar("_rlnRandomSubset", randomValuesVector, starFileIn, ".___TMP_FILE____.tmp");
             copyCvsFile(".___TMP_FILE____.tmp", starFileIn);
             removeCvsFile(".___TMP_FILE____.tmp");
         }
      }
     
   // ###################################################################
   // ###################################################################
   // ##  retrieve all information from refined file, except file name  #   
   if ( parameters.refinedFileName && parameters.starFileOut){
      if (parameters.verboseOn) std::cerr<<"retrieve refined information, except for the file name\n";

//     std::cerr<<"              --retrieveRefined refinedFileName.star\n";

	std::vector<std::string> rlnMicrographNameVector;
	readStar(rlnMicrographNameVector, "_rlnMicrographName", starFileIn);

	std::vector<double> rlnMicrographCoordXVector;
	readStar(rlnMicrographCoordXVector, "_rlnCoordinateX", starFileIn);

	std::vector<double> rlnMicrographCoordYVector;
	readStar(rlnMicrographCoordYVector, "_rlnCoordinateY", starFileIn);

	std::vector<std::string> rlnImageNameNameVector;
	readStar(rlnImageNameNameVector, "_rlnImageName", starFileIn);



	std::vector<std::string> rlnRefinedMicrographNameVector;
	readStar(rlnRefinedMicrographNameVector, "_rlnMicrographName", parameters.refinedFileName);

	std::vector<double> rlnRefinedMicrographCoordXVector;
	readStar(rlnRefinedMicrographCoordXVector, "_rlnCoordinateX", parameters.refinedFileName);

	std::vector<double> rlnRefinedMicrographCoordYVector;
	readStar(rlnRefinedMicrographCoordYVector, "_rlnCoordinateY", parameters.refinedFileName);

	int refinedImageNameidx=getStarHeaderItemIdx("_rlnImageName", parameters.refinedFileName);
	//std::cerr<<idx<<"\n";
	
	std::vector<bool> processed (rlnRefinedMicrographNameVector.size(), false);
	std::vector<long int> foundLine (rlnMicrographNameVector.size(), -1);
	
	for (unsigned long int ii=0; ii<rlnMicrographNameVector.size();ii++){
   	   for (unsigned long int jj=0, found=false; jj<rlnRefinedMicrographNameVector.size() && !found ;jj++){
	     if (!processed[jj]){
	       if(rlnMicrographCoordXVector[ii]==rlnRefinedMicrographCoordXVector[jj]){
	         if ( rlnMicrographCoordYVector[ii]==rlnRefinedMicrographCoordYVector[jj]){
	          if(rlnMicrographNameVector[ii].compare(rlnRefinedMicrographNameVector[jj]) == 0 ){
                       processed[jj]=true;
                       foundLine[jj]=ii+1;
                       found=true;
	          }
	         }
	       }
	     }
	   }
	   
	}


        std::ofstream fileOutput;
        fileOutput.open(parameters.starFileOut);
        fileOutput.close();
        fileOutput.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);   
        long int startMicrograph=getStarStart(parameters.refinedFileName);
        std::ifstream fileRefinedParticles(parameters.refinedFileName);
        std::string strLine;
        for (int counter=0;counter<startMicrograph;counter++){
           std::getline(fileRefinedParticles, strLine);
           fileOutput << strLine <<"\n";
        }
        for (unsigned long int ii=0,counter=1;ii<rlnMicrographNameVector.size();ii++){
          std::getline(fileRefinedParticles, strLine);
          if(foundLine[ii]>=0){
            //fileOutput << strLine <<"\n";
            std::vector<long int> posIdx=getStringBeginAndLenghtAtStarPosition(strLine, refinedImageNameidx);
            fileOutput << strLine.substr(0,posIdx[0])<<"   ";
            fileOutput << std::setfill('0') << std::fixed << std::right << std::setw(8) << std::showpoint << std::setprecision(6)  <<  foundLine[ii];
            fileOutput << "@" << parameters.refinedStackFileName << " ";
            fileOutput << strLine.substr(posIdx[0]+posIdx[1]);
            fileOutput << "\n";
          }
        }
        fileOutput.close();
   }   



   // ########################################################
   // ########################################################
   // ##  replaceLabel
   //          std::cerr<<"              --replaceLabel SourceFileForReplacingLabel.star SourceLabelsNameCsvForReplacing DestinationLabelNameCsvForReplacing\n";
   if (  parameters.starFileOut && parameters.SourceFileForReplacingLabel && parameters.SourceLabelsNameCsvForReplacing ){
     if (parameters.DestinationLabelNameCsvForReplacing == NULL){
       parameters.DestinationLabelNameCsvForReplacing = parameters.SourceLabelsNameCsvForReplacing;
     }
           if (parameters.verboseOn) std::cerr<<"replace Label\n";

     replaceLabel(starFileIn, parameters.starFileOut, parameters.SourceFileForReplacingLabel, parameters.SourceLabelsNameCsvForReplacing, parameters.DestinationLabelNameCsvForReplacing);
   }
     
     // ########################################################
     // ########################################################
     // ##  averageLabels
     //          std::cerr<<"              --averageLabels SourceFilesForAveragingLabels.star SourceLabelnameForAveraging DestinationLabelnameForAveraging\n";
     if (  parameters.starFileOut && parameters.SourceFilesForAveragingLabels && parameters.SourceLabelnameForAveraging ){
           if (parameters.verboseOn) std::cerr<<"average label values\n";
           averageLabels(starFileIn, parameters.starFileOut, parameters.SourceFilesForAveragingLabels, parameters.SourceLabelnameForAveraging);
     }

     // ########################################################
     // ########################################################
     // ##  assessEuler
     // std::cerr<<"              --assessEuler gtFileForAssessingLabels.star gtLabelsnameForAssessing \n";
     if (  parameters.gtFileForAssessingLabels  ){
           if (parameters.verboseOn) std::cerr<<"assessEuler\n";

         std::vector<double> AngleRot;
         std::vector<double> AngleTilt;
         std::vector<double> gtValuesRot;
         std::vector<double> gtValuesTilt;
         
         readStar(AngleRot, "_rlnAngleRot", starFileIn);
         readStar(gtValuesRot, "_rlnAngleRot", parameters.gtFileForAssessingLabels);
         readStar(AngleTilt, "_rlnAngleTilt", starFileIn);
         readStar(gtValuesTilt, "_rlnAngleTilt", parameters.gtFileForAssessingLabels);

         double mean=0.0;
         for (unsigned long int ii=0; ii<AngleRot.size(); ii++){
          mean+=archDistanceSquaredApproximate(AngleRot[ii], AngleTilt[ii], gtValuesRot[ii], gtValuesTilt[ii]);
         }
         mean /= AngleRot.size();
         
         std::cerr<< mean << "\n";
     }
  
     
   // ########################################################
   // ########################################################
   // ##  importing micrograph names from reference          #   
   if ( parameters.referenceParticlesStarFile && parameters.starFileOut){
        if (parameters.verboseOn) std::cerr<<"processing referenceParticles\n";
       //std::cerr<<"              --rp referenceParticlesStarFile.star\n";

       
       /*for (unsigned long int ii=0; ii<simplifiedOriginalMicrographItems.size(); ii++){
         if (str.find(lastMicrographItem) != string::npos){
             
         }else{
           std::cerr<<simplifiedOriginalMicrographItems[ii]<<"  ->  "
 
         }
         if (foundReferenceMicrograph)
       }*/
       updateMicrographStarFile(starFileIn, parameters.referenceParticlesStarFile, parameters.starFileOut, parameters.binningFactor);
   }
   
   if (parameters.RetainGroupNumber>=0){
            if (parameters.verboseOn) std::cerr<<"Retain Specific Group\n";

           std::string itemType ="_subset_group";
           if (parameters.RetainGroupMetadata){
             itemType=std::string(parameters.RetainGroupMetadata);
           }
           std::string valueToRetain = std::to_string(parameters.RetainGroupNumber);
           std::cerr<<"valueToRetain = "<<valueToRetain << "  ("<< itemType<<")\n";
           removeLinesStar(itemType, valueToRetain, starFileIn, parameters.starFileOut);
           if (parameters.starFileOut){
            starFileIn=parameters.starFileOut;
           }else{
             starFileIn=parameters.starFileIn;
           }
           //int idx=readStar(phiList, "_rlnCoordinateY", parameters.metadataStar);
           //std::cerr<<"index of _rlnCoordinateY is " << idx << "\n";
   }
   if(parameters.halfMapsTag){
       removeLinesStar("_rlnRandomSubset", std::string("1"), parameters.starFileIn, parameters.halfMap1Out);
       removeLinesStar("_rlnRandomSubset", std::string("2"), parameters.starFileIn, parameters.halfMap2Out);
   }

   if( parameters.columnToDelete ){
     if (parameters.starFileOut == NULL){
      parameters.starFileOut=parameters.starFileIn;
     }
     std::string columnsToDelete ( parameters.columnToDelete );
     int positionWidlchard = columnsToDelete.find("*");
     if (positionWidlchard < 0 ){
       //std::cerr<<"no wildcard character\n";
       std::cerr<<"deleting tag:"<< parameters.columnToDelete <<"\n";
       removeColumnStar( parameters.columnToDelete, starFileIn, parameters.starFileOut);
     }else{
       //std::cerr << "wildcard character at " << positionWidlchard << "\n";
       std::vector<std::string> columnsToDeleteLists;
       std::vector<int> fieldsIdx;
       getStarHeaders(columnsToDeleteLists, fieldsIdx, parameters.starFileIn);
       std::string basenameToDelete = columnsToDelete.substr (0, positionWidlchard);
       std::string TmpFilename = generateTmpFilename(parameters.starFileIn).c_str();
       copyCvsFile(parameters.starFileIn, TmpFilename.c_str());

      //  std::cerr<<"TmpFilename2="<<TmpFilename<<"\n";
       for (int kk=columnsToDeleteLists.size()-1; kk >= 0; kk--){
        //std::cerr << std::string(columnsToDeleteLists[kk]) << "\n";
        std::string tmpString = std::string(columnsToDeleteLists[kk]).substr (0, positionWidlchard);
        if ( basenameToDelete == tmpString){
          std::cerr<<"deleting tag:"<< columnsToDeleteLists[kk] <<"\n";
//          std::string TmpFilename1 = generateTmpFilename(parameters.starFileIn).c_str();
          removeColumnStar( columnsToDeleteLists[kk].c_str(), TmpFilename.c_str(), TmpFilename.c_str());
//          copyCvsFile(TmpFilename1.c_str(), TmpFilename.c_str());
        }
       }
       //std::cerr<<"TmpFilename3="<<TmpFilename<<"\n";
       //std::cerr<<"parameters.starFileOut="<<parameters.starFileOut<<"\n";
       copyCvsFile(TmpFilename.c_str(), parameters.starFileOut);
       removeCvsFile(TmpFilename.c_str());
     }
   }
   if (parameters.vemFileOut){
    StarToCsv(starFileIn, parameters.vemFileOut, true);
   }
   if (parameters.csvFileOut){
    StarToCsv(starFileIn, parameters.csvFileOut, false);
   }
// ***********************
// ***********************
// SHOW SIMPLE INFO
// ***********************
   if(parameters.showInfo){
    if (parameters.verboseOn) std::cerr<<"show Info\n";

    unsigned long int numItems=countParticlesItems(parameters.starFileIn);
    std::cerr<<"\n**************\n****  num particles= "<< numItems <<"\n**************\n";
    std::vector<double> randomSubsetList;
    readStar(randomSubsetList, "_rlnRandomSubset", parameters.starFileIn);
    unsigned long int count=0;
    std::cerr<<"subsets:";   
    if (randomSubsetList.size()>0){
      unsigned long int count1=0;
      unsigned long int count2=0;
      double mean1 =0;
      double mean2 =0;
      for (unsigned long int ii=0; ii<randomSubsetList.size(); ii++){
        if(randomSubsetList[ii]==1.0){
           count1++;
        }else if(randomSubsetList[ii]==2.0){
           count2++;
        }
      }
      std::cerr<< "\n";
      std::cerr<< "    "<< count1 << "  (subset1)\n";
      std::cerr<< "    "<< count2 << "  (subset2)\n";
    }else{
       std::cerr<<" NONE\n"; 
    }
   }



// *************************
// *************************
//  INFO DIFFERENCE
// *************************
   if (parameters.starInfoDifferenceFile){
    if (parameters.verboseOn) std::cerr<<"star Info Difference File\n";

            std::vector<double> randomSubsetList;
            readStar(randomSubsetList, "_rlnRandomSubset", parameters.starFileIn);

            //euler
	    std::vector<double> eulerDiff;
            std::vector<double> eulerDiff_h1;
	    std::vector<double> eulerDiff_h2;
	    std::vector<double> phiListParticle1;
            std::vector<double> thetaListParticle1;
            std::vector<double> phiListParticle2;
            std::vector<double> thetaListParticle2;
            readMetadataValues(phiListParticle1, "_rlnAngleRot", parameters.starFileIn);
            readMetadataValues(thetaListParticle1, "_rlnAngleTilt", parameters.starFileIn);
            readMetadataValues(phiListParticle2, "_rlnAngleRot", parameters.starInfoDifferenceFile);
            readMetadataValues(thetaListParticle2, "_rlnAngleTilt", parameters.starInfoDifferenceFile);
	    if (phiListParticle1.size()!=phiListParticle2.size() || thetaListParticle1.size()!=thetaListParticle2.size()){
		std::cerr<<"ERROR: different size euler.. EXIT\n";
                exit(0);
            }
	    double mean=0;
            double meanH1=0;
            double meanH2=0;
	    double sd=0;
            double sdH1=0;
            double sdH2=0;

            for (unsigned long int ii=0; ii<randomSubsetList.size(); ii++){
             double distance=archDistanceSquaredApproximate(phiListParticle1[ii], thetaListParticle1[ii], phiListParticle2[ii], thetaListParticle2[ii]);
             distance=pow(distance,0.5)*(180.0/PI);
             eulerDiff.push_back(distance);
             mean+=distance;
             if ( ii<randomSubsetList.size() ){
		     if(randomSubsetList[ii]==1.0){
		       eulerDiff_h1.push_back(distance);
		       meanH1+=distance;
		     }else if(randomSubsetList[ii]==2.0){
		       eulerDiff_h2.push_back(distance);
                       meanH2+=distance;
		     }
	     }
            }
            if (randomSubsetList.size()-1>0){
               mean/=randomSubsetList.size();
	       for (unsigned long int ii=0; ii<randomSubsetList.size(); ii++){
                  sd+=pow(eulerDiff[ii]-mean,2.0);
               }
               sd=pow(sd/(randomSubsetList.size()-1.0),0.5);
	    }
	    if ( eulerDiff_h1.size()>0 ){
               meanH1/=eulerDiff_h1.size();
	       for (unsigned long int ii=0; ii<eulerDiff_h1.size(); ii++){
                  sdH1+=pow(eulerDiff_h1[ii]-mean,2.0);
               }
               sdH1=pow(sdH1/(eulerDiff_h1.size()-1.0),0.5);
            }
	    if ( eulerDiff_h2.size()>0 ){
               meanH2/=eulerDiff_h2.size();
	       for (unsigned long int ii=0; ii<eulerDiff_h2.size(); ii++){
                  sdH2+=pow(eulerDiff_h2[ii]-mean,2.0);
               }
               sdH2=pow(sdH2/(eulerDiff_h2.size()-1.0),0.5);
            }
	    std::cerr<<"euler stats:\n";
	    std::cerr<<"     mean(SD): "<< mean<<" ("<< sd<<")\n";
	    std::cerr<<"       h1(SD): "<< meanH1<<" ("<< sdH1 <<")\n";
	    std::cerr<<"       h2(SD): "<< meanH2<<" ("<< sdH2 <<")\n";

//origin
	    std::vector<double> originDiff;
            std::vector<double> originDiff_h1;
	    std::vector<double> originDiff_h2;
	    std::vector<double> XListParticle1;
            std::vector<double> YListParticle1;
            std::vector<double> XListParticle2;
            std::vector<double> YListParticle2;
            readMetadataValues(XListParticle1, "_rlnOriginX", parameters.starFileIn);
            readMetadataValues(YListParticle1, "_rlnOriginY", parameters.starFileIn);
            readMetadataValues(XListParticle2, "_rlnOriginX", parameters.starInfoDifferenceFile);
            readMetadataValues(YListParticle2, "_rlnOriginY", parameters.starInfoDifferenceFile);
	    if (XListParticle1.size()!=YListParticle1.size() || XListParticle2.size()!=YListParticle2.size()){
		std::cerr<<"ERROR: different size origin.. EXIT\n";
                exit(0);
            }
	    mean=0;
            meanH1=0;
            meanH2=0;
	    sd=0;
            sdH1=0;
            sdH2=0;

            for (unsigned long int ii=0; ii<randomSubsetList.size(); ii++){
             double distance=pow(pow(XListParticle1[ii]-XListParticle2[ii],2.0)+pow(YListParticle1[ii]-YListParticle2[ii],2.0),0.5);
             originDiff.push_back(distance);
             mean+=distance;
             if ( ii<randomSubsetList.size() ){
		     if(randomSubsetList[ii]==1.0){
		       originDiff_h1.push_back(distance);
		       meanH1+=distance;
		     }else if(randomSubsetList[ii]==2.0){
		       originDiff_h2.push_back(distance);
                       meanH2+=distance;
		     }
	     }
            }
            if (randomSubsetList.size()-1>0){
               mean/=randomSubsetList.size();
	       for (unsigned long int ii=0; ii<randomSubsetList.size(); ii++){
                  sd+=pow(originDiff[ii]-mean,2.0);
               }
               sd=pow(sd/(randomSubsetList.size()-1.0),0.5);
	    }
	    if ( originDiff_h1.size()>0 ){
               meanH1/=originDiff_h1.size();
	       for (unsigned long int ii=0; ii<originDiff_h1.size(); ii++){
                  sdH1+=pow(originDiff_h1[ii]-mean,2.0);
               }
               sdH1=pow(sdH1/(originDiff_h1.size()-1.0),0.5);
            }
	    if ( originDiff_h2.size()>0 ){
               meanH2/=originDiff_h2.size();
	       for (unsigned long int ii=0; ii<originDiff_h2.size(); ii++){
                  sdH2+=pow(originDiff_h2[ii]-mean,2.0);
               }
               sdH2=pow(sdH2/(originDiff_h2.size()-1.0),0.5);
            }
	    std::cerr<<"origin stats:\n";
	    std::cerr<<"     mean(SD): "<< mean<<" ("<< sd<<")\n";
	    std::cerr<<"       h1(SD): "<< meanH1<<" ("<< sdH1 <<")\n";
	    std::cerr<<"       h2(SD): "<< meanH2<<" ("<< sdH2 <<")\n";

           //psi
	    std::vector<double> psiDiff;
            std::vector<double> psiDiff_h1;
	    std::vector<double> psiDiff_h2;
	    std::vector<double> psiListParticle1;
            std::vector<double> psiListParticle2;
            readMetadataValues(psiListParticle1, "_rlnAnglePsi", parameters.starFileIn);
            readMetadataValues(psiListParticle2, "_rlnAnglePsi", parameters.starInfoDifferenceFile);
	    if (psiListParticle1.size()!=psiListParticle2.size() ){
		std::cerr<<"ERROR: different size psi list.. EXIT\n";
                exit(0);
            }
	    mean=0;
            meanH1=0;
            meanH2=0;
	    sd=0;
            sdH1=0;
            sdH2=0;

            for (unsigned long int ii=0; ii<randomSubsetList.size(); ii++){
             double distance=abs(psiListParticle1[ii]-psiListParticle2[ii]);
             psiDiff.push_back(distance);
             mean+=distance;
             if ( ii<randomSubsetList.size() ){
		     if(randomSubsetList[ii]==1.0){
		       psiDiff_h1.push_back(distance);
		       meanH1+=distance;
		     }else if(randomSubsetList[ii]==2.0){
		       psiDiff_h2.push_back(distance);
                       meanH2+=distance;
		     }
	     }
            }
            if (randomSubsetList.size()-1>0){
               mean/=randomSubsetList.size();
	       for (unsigned long int ii=0; ii<randomSubsetList.size(); ii++){
                  sd+=pow(psiDiff[ii]-mean,2.0);
               }
               sd=pow(sd/(randomSubsetList.size()-1.0),0.5);
	    }
	    if ( psiDiff_h1.size()>0 ){
               meanH1/=psiDiff_h1.size();
	       for (unsigned long int ii=0; ii<psiDiff_h1.size(); ii++){
                  sdH1+=pow(psiDiff_h1[ii]-mean,2.0);
               }
               sdH1=pow(sdH1/(psiDiff_h1.size()-1.0),0.5);
            }
	    if ( psiDiff_h2.size()>0 ){
               meanH2/=psiDiff_h2.size();
	       for (unsigned long int ii=0; ii<psiDiff_h2.size(); ii++){
                  sdH2+=pow(psiDiff_h2[ii]-mean,2.0);
               }
               sdH2=pow(sdH2/(psiDiff_h2.size()-1.0),0.5);
            }
	    std::cerr<<"psi stats:\n";
	    std::cerr<<"     mean(SD): "<< mean<<" ("<< sd<<")\n";
	    std::cerr<<"       h1(SD): "<< meanH1<<" ("<< sdH1 <<")\n";
	    std::cerr<<"       h2(SD): "<< meanH2<<" ("<< sdH2 <<")\n";
  }
 } //parameters star file in

// #############################################################
// #############################################################
//        HERE it is not necessary to have a star file in


   // ###################################################################
   // ###################################################################
   // Simulate a starFile with homogeneous particle distribution
  else if ( parameters.numHomogeneousDistribution > 0 && parameters.starFileOut){
            if (parameters.verboseOn) std::cerr<<"Do Simulate a starFile with homogeneous particle distribution\n";
            std::vector<double> phiV;
            std::vector<double> thetaV;
            std::vector<double> subsetV;


            const unsigned long int numActualViews=parameters.numHomogeneousDistribution+1;
            const double ga = (3.0 - pow(5.0,0.5)) * PI; // golden angle  
            double K=-1.0;
            long int views=0;
            for (unsigned long int ii=0; ii<numActualViews; ii++, K+=2.0/((double)numActualViews-1.0)){
              if ( K <= 1.0 ){
                      //wrapping angle
                      double angle=ii*ga;
                      double range=2.0*PI;
                      if (angle<0){
                        double tmp=range+range-angle;
                        long int tmp1=ceil(tmp/range);
                        angle=angle+range*tmp1;
                      }
                      double wrappedPhi=fmod(angle,range);
                      double anglePhi = wrappedPhi*180.0/PI;
                      double angleTheta = acos(K)*180.0/PI;
                      thetaV.push_back( angleTheta );
                      phiV.push_back( anglePhi );
                      views++;
              }
            }
            std::ofstream fileOutput;
            fileOutput.open(parameters.starFileOut);
            fileOutput.close();
            fileOutput.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);  
                        std::string strLine;


              fileOutput << "# version 30001\n";
              fileOutput << "\n";
              fileOutput << "data_optics\n";
              fileOutput << "\n";
              fileOutput << "loop_\n"; 
              fileOutput << "_rlnVoltage #1\n"; 
              fileOutput << "_rlnImagePixelSize #2\n"; 
              fileOutput << "_rlnSphericalAberration #3\n"; 
              fileOutput << "_rlnAmplitudeContrast #4\n"; 
              fileOutput << "_rlnOpticsGroup #5\n"; 
              fileOutput << "_rlnImageSize #6\n";
              fileOutput << "_rlnImageDimensionality #7\n"; 
              fileOutput << "_rlnOpticsGroupName #8\n";
              fileOutput << "  300.000000     1.000000     2.700000     0.100000            1          128            2 opticsGroup5 \n";
              fileOutput << "\n";
              fileOutput << "\n";
              fileOutput << "# version 30001\n";
              fileOutput << "\n";
              fileOutput << "data_particles\n";
              fileOutput << "\n";
              fileOutput << "loop_\n";
              fileOutput << "_rlnImageName #1\n";
              fileOutput << "_rlnAngleRot #2\n";
              fileOutput << "_rlnAngleTilt #3\n";
              fileOutput << "_rlnAnglePsi #4\n";
              fileOutput << "_rlnOriginXAngst #5\n";
              fileOutput << "_rlnOriginYAngst #6\n"; 
              fileOutput << "_rlnDefocusU #7\n";
              fileOutput << "_rlnDefocusV #8\n";
              fileOutput << "_rlnDefocusAngle #9\n";
              fileOutput << "_rlnPhaseShift #10\n";
              fileOutput << "_rlnCtfBfactor #11\n";
              fileOutput << "_rlnOpticsGroup #12\n";
              fileOutput << "_rlnRandomSubset #13\n";
              fileOutput << "_rlnClassNumber #14\n";
              for (unsigned long int ii=0, counterOutput=1; ii<phiV.size(); ii++, counterOutput++ ){
                double phi=phiV[ii];
                double psi=rand()%360;
                double theta=thetaV[ii];
                double originX = 0;
                double originY = 0;
                double defocusU=1000.000;
                double defocusV=1000.000;
                double defocusAngle=10.00;
                double PhaseShift=0;
                double CtfBfactor=0;
                int OpticsGroup=1;
                int RandomSubset=ii%2+1;
                int classNumber=1;

                fileOutput << counterOutput<<"@"<< parameters.simulatedStackName <<" "
                  << phi<<" "
                  << psi <<" "
                  << theta << " "
                  << originX << " "
                  << originY << " "
                  << defocusU << " "
                  << defocusV << " "
                  << defocusAngle << " "
                  << PhaseShift << " "
                  << CtfBfactor << " "
                  << OpticsGroup << " "
                  << RandomSubset << " "
                  << classNumber << "\n";
              }
  }

// #############################################################
// #############################################################
//        if nothing is covered, then send the usage message
else{
  usage(argv);
}



 return 0;

}
