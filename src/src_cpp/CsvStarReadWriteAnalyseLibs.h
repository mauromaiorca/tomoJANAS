/*
 * File: CsvStarReadWriteAnalyseLibs.h
 * (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
 *
 */

#ifndef __CVS_STAR_READ_WRITE_ANALISE__H___
#define __CVS_STAR_READ_WRITE_ANALISE__H___

#define __MAX_STAR_HEADER_SIZE__ 130

#include <fstream>
#include <string>
#include <cstdio>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <iomanip>      // std::setprecision
#include <regex>
#include "mrcIO.h"
#include "randomLibs.h"

//some prototypes
void replaceAddCvsColumn(std::string columnName, std::vector<std::string> & itemVector, const char * filename);
void replaceAddCvsColumn(std::string columnName, std::vector<double> & itemVector, const char * filename);

//###############################
// ACCESSORY PARAMETERS PARSING CSV OPTIONS
// GETS PARAMETERS FROM OPTIONS
std::vector<std::vector<std::string> > getCsvOptionParameters(char * startingString, int numItemsPerRow){
                std::stringstream SigmaWidthListStream(startingString);
                std::vector<std::vector<std::string> > outputList;
                while( SigmaWidthListStream.good() ){
                  std::string substr;
                  getline( SigmaWidthListStream, substr, ',' );
                  std::stringstream SigmaWidthListStream2(substr);
                  int ccounter = 0;
                  std::vector<std::string> values;                   
                  while( SigmaWidthListStream2.good() ){
                    std::string substr2;
                    getline( SigmaWidthListStream2, substr2, ':' );
                    //values.push_back(std::stod(substr2));
                    values.push_back(substr2);
                    ccounter ++;
                  }
                  if (values.size()==numItemsPerRow){
                    outputList.push_back(values);
                  }
                }
                return outputList;
}
std::vector<std::vector<double> > getCsvOptionDoubleParameters(char * startingString, int numItemsPerRow){
    std::vector<std::vector<std::string> > strParameters=getCsvOptionParameters(startingString, numItemsPerRow);
    std::vector<std::vector<double> > outputList;
    for (int ii=0; ii<strParameters.size(); ii++){
      std::vector<double> tmpVal;
      for (int jj=0; jj<(strParameters[ii]).size(); jj++){
        tmpVal.push_back(std::stod(strParameters[ii][jj]));
      }
      outputList.push_back(tmpVal);
    }
    return outputList;
}
std::vector<double> getColonSeparatedValues(char * startingString){
                std::stringstream SigmaWidthListStream(startingString);
                std::vector<double> outputList;
                while( SigmaWidthListStream.good() ){
                  std::string substr;
                  getline( SigmaWidthListStream, substr, ':' );
                  outputList.push_back(atof(substr.c_str()));
                }
                return outputList;
}
// ******************************

//###############################
// ACCESSORY STRING FUNCTIONS
// (useful for star file)
std::vector<long int> getStringBeginAndLenghtAtStarPosition(std::string line, int position){

  long int startIdx = 0;
  long int length = 0;
  const long int lineSize = line.size();
  long int currentPosition = 0;

  //search for the space at the beginning
  for (bool keepGoing = true ;startIdx < lineSize && keepGoing; ){
   if ( !( line[startIdx]==' ' || line[startIdx]=='\t' ) ){
     keepGoing=false;
   }else{
     startIdx++;
   }
  }

  for (int ii=startIdx, foundStart=false; ii<lineSize - 1 && !foundStart;ii++){
   if (currentPosition==position){
     foundStart=true;
   }else{
    if ( (line[ii]==' ' || line[ii]=='\t' ) &&
         !(line[ii+1]==' ' || line[ii+1]=='\t' || line[ii+1]=='\n')  ){
       currentPosition++;
       startIdx=ii+1;
    }
   }
  }

  for (int jj=startIdx,foundEnd=false; jj<lineSize && !foundEnd; jj++){
   if ( line[jj]==' ' || line[jj]=='\t' || line[jj]=='\n'  ){
       foundEnd=true;
   }else{
     length++;
   }
  }
   
 std::vector<long int> result;
 result.push_back(startIdx);
 result.push_back(length);
 return result;
}

// ****************
// Replace a value of a starfile at a certain position idx
std::string replaceValueStrlineStarFile(std::string inputStr, int StarPositionIdx, std::string valueToReplace){
   std::vector<long int> BeginAndLenght=getStringBeginAndLenghtAtStarPosition(inputStr, StarPositionIdx);
   std::string outputStr(inputStr.substr(0,BeginAndLenght[0]));
   if (StarPositionIdx>0){
     outputStr+=std::string (" ");
   }
   outputStr+=valueToReplace+std::string (" ");
   outputStr+=inputStr.substr(BeginAndLenght[0]+BeginAndLenght[1]);
   //outputStr+=std::string (" --> ")+valueToReplace;
   return outputStr;
}

//add leading zeros to a string
std::string zeros_lead_to_string(const long int value, const int numFields){
   std::ostringstream oss;
   oss << std::setw(numFields)<<std::setfill('0')<<value;
   return oss.str();
}

// ACCESSORY FILESYSTEM FUNCTIONS

std::string generateTmpFilename(const char * extension){
  std::string dir ("");
  std::string extensionString(extension);
  int slashPosition = extensionString.find_last_of("/");
  if (slashPosition>=0){
    dir  = extensionString.substr (0, slashPosition) + std::string("/");
    extensionString = extensionString.substr (slashPosition+1);
    //std::cerr<<dir<<"\n";
    //std::cerr<<extensionString<<"\n";
  }
  std::string outFile ("___TMP___");
  outFile=dir+outFile;
  outFile=outFile+std::to_string(rnd_Int()) + std::string("___")+ extensionString;
  //std::cerr<<outFile<<"\n";
  return outFile;
}


int renameFile(const char * filenameSrc, const char * filenameDst){
   std::ifstream infile(filenameSrc); //check fileExists
   if ( infile.good() ){
    std::string mvStr=std::string("mv ")+std::string(filenameSrc)+std::string(" ")+std::string(filenameDst);
    int ff=system(mvStr.c_str());
    return ff;
   }else{
     std::cerr<<"ERROR mv: reading file "<<filenameSrc<<"  \n";
   }
   return 1;
}


// ***************************
// copyCvsFile
// NEED TO BE IMPROVED (no system call)
// ***************************
int copyCvsFile(const char * filenameSrc, const char * filenameDst){

//   std::ifstream  srcFile(parameters.inputList);
//   std::ofstream  dstFile(parameters.outputList);
//   dstFile << srcFile.rdbuf();
   std::ifstream infile(filenameSrc); //check fileExists
   if ( infile.good() ){
    std::string copyStr=std::string("cp ")+std::string(filenameSrc)+std::string("  ")+std::string(filenameDst);
    int ff=system(copyStr.c_str());
    return ff;
   }else{
     std::cerr<<"ERROR cp: reading file "<<filenameSrc<<"  \n";
   }
   return 1;
}

int removeCvsFile(const char * filename){

//   std::ifstream  srcFile(parameters.inputList);
//   std::ofstream  dstFile(parameters.outputList);
//   dstFile << srcFile.rdbuf();
    std::ifstream infile(filename); //check fileExists
    if ( infile.good() ){
       std::string rmStr=std::string("rm ")+std::string(filename);
       int ff=system(rmStr.c_str());
       return ff;
    }else{
     std::cerr<<"ERROR rm: reading file "<<filename<<"  \n";
   }
    return 1;
}

std::string filesystem_basename(const std::string& s) {
   char sep = '/';
#ifdef _WIN32
   sep = '\\';
#endif
   size_t i = s.rfind(sep, s.length());
   if (i != std::string::npos) {
      return(s.substr(i+1, s.length() - i));
   }
   return("");
}

int filesystem_create_directory(const std::string& s) {
   std::string mkdirStr("mkdir ");
   mkdirStr=mkdirStr+s;
   int rr=system(mkdirStr.c_str());
   return rr;
}



// ////////////////
//
// accessory functions
//
int stringPositionVector(const std::vector<std::string> vecStr, const char * strItem){
         size_t index=std::distance(vecStr.begin(),std::find(vecStr.begin(),vecStr.end(), strItem));
         if (index<vecStr.size()){
          return index;
         }else{
          return -1;
         }
}




std::string getColumnAtPosition(std::string valueStr, int index, const char * separatorChar = ":"){

   std::string tmpStartStr=valueStr;
   std::string tmpStr=valueStr;
   int startIdx=0;
   int endIdx=valueStr.length();
   for (int ii = 0; ii<=index && startIdx<tmpStartStr.length(); ii++){
    tmpStartStr=tmpStartStr.substr(startIdx);
    int separatorIdx=tmpStartStr.find(separatorChar);
    tmpStr=tmpStartStr.substr(0,separatorIdx);
    startIdx=separatorIdx+1;
   }
   return tmpStr;
}

std::string filelineToString(const char * filename){
        std::string item_name;
        std::ifstream nameFileout;
        nameFileout.open(filename);
        std::getline(nameFileout, item_name);
        nameFileout.close();
        return item_name;
}


// ////////////////
//
// readCvs
//
int readCvsColumns(std::vector<std::string> & itemVector, const char * filename){
    itemVector.clear();
    std::ifstream file(filename);
    std::string str;
    bool csvTitleFound = false;
    //int position = -1;
    while (std::getline(file, str) && !csvTitleFound){
	//str.erase(std::remove(str.begin(),str.end(),' '),str.end()); //remove spaces
	//str.erase(std::remove(str.begin(),str.end(),'\t'),str.end()); //remove tabs
	if (str[0]=='#'){
			str.erase(0, 1); //remove first character (#)
			std::stringstream ss;
			ss.str(str);
			char delim=',';
			std::string item;
			bool found = false;
			while (std::getline(ss, item, delim)) {
        			itemVector.push_back(item);
			}
			csvTitleFound=true;
        }
    }
    return itemVector.size();
}

int indexCvsColumn(std::string itemType, const char * filename){
  std::vector<std::string> columns;
  readCvsColumns(columns, filename);
  for (unsigned long int ii=0; ii<columns.size();ii++){
	if ((columns[ii]).compare(itemType) == 0){
		return ii;
	}
  }
  return -1;
}


// ////////////////
//
// readCvs
//
int readCvs(std::vector<double> & itemVector, std::string itemType, const char * filename)
{

    std::transform(itemType.begin(), itemType.end(), itemType.begin(), ::tolower);
    std::ifstream file(filename);
    std::string str;
    bool csvTitleFound = false;
    int position = -1;
//    if (!itemType.compare("")){ //empty string
//      csvTitleFound = true;
//      position = 0;
//    }
    

    while (std::getline(file, str))
    {
	str.erase(std::remove(str.begin(),str.end(),' '),str.end()); //remove spaces
	str.erase(std::remove(str.begin(),str.end(),'\t'),str.end()); //remove tabs
	if (str[0]=='#'){
		if(!csvTitleFound){
			std::transform(str.begin(), str.end(), str.begin(), ::tolower); //to lower
			str.erase(0, 1); //remove first character (#)
			std::stringstream ss;
			ss.str(str);
			char delim=',';
			std::string item;
			bool found = false;
			int positionTmp=position;
			while (std::getline(ss, item, delim) && !found) {
				positionTmp++;
				if (item.compare(itemType) == 0){
					found = true;
				}
			}
			if (found)
				position=positionTmp;
			csvTitleFound=true;
		}
	}else if (position>=0){
		std::stringstream ss;
		ss.str(str);
		char delim=',';
		std::string item;
		bool found = false;
		int tmpPosition = 0;
		while (std::getline(ss, item, delim) && tmpPosition<=position) {
			if (tmpPosition==position){
				itemVector.push_back(atof(item.c_str()));
			}
			tmpPosition++;
		}

	}
    }
    return position;
}


int readCvs(std::vector<std::string> & itemVector, std::string itemType, const char * filename)
{

    std::transform(itemType.begin(), itemType.end(), itemType.begin(), ::tolower);
    std::ifstream file(filename);
    std::string str;
    bool csvTitleFound = false;
    int position = -1;
//    if (!itemType.compare("")){ //empty string
//      csvTitleFound = true;
//      position = 0;
//    }
    

    while (std::getline(file, str))
    {
	str.erase(std::remove(str.begin(),str.end(),' '),str.end()); //remove spaces
	str.erase(std::remove(str.begin(),str.end(),'\t'),str.end()); //remove tabs
	if (str[0]=='#'){
		if(!csvTitleFound){
			std::transform(str.begin(), str.end(), str.begin(), ::tolower); //to lower
			str.erase(0, 1); //remove first character (#)
			std::stringstream ss;
			ss.str(str);
			char delim=',';
			std::string item;
			bool found = false;
			int positionTmp=position;
			while (std::getline(ss, item, delim) && !found) {
				positionTmp++;
				if (item.compare(itemType) == 0){
					found = true;
				}
			}
			if (found)
				position=positionTmp;
			csvTitleFound=true;
		}
	}else if (position>=0){
		std::stringstream ss;
		ss.str(str);
		char delim=',';
		std::string item;
		bool found = false;
		int tmpPosition = 0;
		while (std::getline(ss, item, delim) && tmpPosition<=position) {
			if (tmpPosition==position){
				itemVector.push_back(item.c_str());
			}
			tmpPosition++;
		}

	}
    }
    return position;
}


int readCvs(std::vector<long int> & itemVector, std::string itemType, const char * filename)
{

    std::transform(itemType.begin(), itemType.end(), itemType.begin(), ::tolower);
    std::ifstream file(filename);
    std::string str;
    bool csvTitleFound = false;
    int position = -1;
    

    while (std::getline(file, str))
    {
	str.erase(std::remove(str.begin(),str.end(),' '),str.end()); //remove spaces
	str.erase(std::remove(str.begin(),str.end(),'\t'),str.end()); //remove tabs
	if (str[0]=='#'){
		if(!csvTitleFound){
			std::transform(str.begin(), str.end(), str.begin(), ::tolower); //to lower
			str.erase(0, 1); //remove first character (#)
			std::stringstream ss;
			ss.str(str);
			char delim=',';
			std::string item;
			bool found = false;
			int positionTmp=position;
			while (std::getline(ss, item, delim) && !found) {
				positionTmp++;
				if (item.compare(itemType) == 0){
					found = true;
				}
			}
			if (found)
				position=positionTmp;
			csvTitleFound=true;
		}
	}else if (position>=0){
		std::stringstream ss;
		ss.str(str);
		char delim=',';
		std::string item;
		bool found = false;
		int tmpPosition = 0;
		while (std::getline(ss, item, delim) && tmpPosition<=position) {
			if (tmpPosition==position){
				itemVector.push_back(std::stoi(item));
			}
			tmpPosition++;
		}

	}
    }
    return position;
}



// *****************
//
//  star file versione
//
// *****************
//
int checkStarFileVersion (const char * filename){

    //std::vector<double> itemVector;
    std::ifstream file(filename);
    std::string strLine;
    unsigned long int cc = 0;
    bool got_data_images = false;
    bool got_loop = false;
    int targetID = -1;
    unsigned long int headerLines = 0;

    int dataBlockCount=0;
    int loopBlockCount=0;
    //int relion3000_hints=0;
    int relion3100_hints=0;
    while (std::getline(file, strLine) && ++cc < __MAX_STAR_HEADER_SIZE__){
        std::regex rNoSpaces("\\s+");
        std::string strNoSpaces = std::regex_replace(strLine, rNoSpaces, ",");
        std::regex r("\\s+");
        std::string str0 = std::regex_replace(strLine, r, ",");
       if ( strNoSpaces.length()<1 ){
         //do nothing
       }else if ( strNoSpaces.substr(0,1).compare("#")==0 ){
         //it is a comment, do nothing
       }else if( strstr(strNoSpaces.substr(0,5).c_str(),"data_") ) {
           dataBlockCount++;
           if( strstr(strNoSpaces.substr(0,14).c_str(),"data_particles") ) relion3100_hints++;
       }else if( strstr(str0.substr(0,5).c_str(),"loop_") ){
            loopBlockCount++;
       }
       /*else if( strstr(str0.substr(0,21).c_str(),"_rlnDetectorPixelSize") ){
           relion3000_hints++;
       }else if( strstr(str0.substr(0,11).c_str(),"_rlnOriginX") ){
           relion3000_hints++;
       }else if( strstr(str0.substr(0,11).c_str(),"_rlnOriginY") ){
           relion3000_hints++;
       }
       else if( strstr(str0.substr(0,18).c_str(),"_rlnImagePixelSize") ){
            relion3100_hints++;
       } else if( strstr(str0.substr(0,16).c_str(),"_rlnOriginXAngst") ){
            relion3100_hints++;
       } else if( strstr(str0.substr(0,16).c_str(),"_rlnOriginYAngst") ){
                  relion3100_hints++;
       }*/
    }
    if (dataBlockCount>1 && loopBlockCount>1 && relion3100_hints > 0){
        return 3100;
    }else if (dataBlockCount==1 && loopBlockCount==1) {
        return 3000;
    }else{
        return 0;
    }
    
}



std::string relionDataStartLabel(const char * filename) {
    int starFileVersion=checkStarFileVersion (filename); //return 3100; 3000
    std::string dataLabelStr("data_");
    if (starFileVersion==3100){
        dataLabelStr="data_particles";
    }
    return dataLabelStr;
}

// *****************
//
//  Read star file
//
// *****************
std::string findColumnItem(std::string str, int idx){
  long int posStart=-1;
  long int  posEnd=str.length()-1;
  
  for (int ii=0, found=false; ii<idx && ! found; ii++){
	  if (str[ii]==' ' || str[ii]=='\t'){
		  posStart++;
	  }else{
		  found=true;
	  }
  }
  for (int ii=0; ii<idx; ii++){
   posStart+=1+str.substr(posStart+1,posEnd).find(",");
  }
  posEnd=str.substr(posStart+1,posEnd).find(",");
  return str.substr(posStart+1,posEnd);
}


// EXTRACT A VECTOR FROM A LINE
std::vector<std::string> stringLineToVector(const std::string str0, const char delimiter=','){
  std::string str (str0);
  str=str.erase(0, str.find_first_not_of(" \t\n\r\f\v"));
  std::regex r("\\s+");
  str = std::regex_replace(str, r, ",");
  std::stringstream ss(str);
  std::vector<std::string> result;
  std::string tmp;
  while(getline(ss, tmp, delimiter)){
   result.push_back(tmp);
  }
  return result;
}


//_rlnOpticsGroupName
int readDataOptics(std::vector<std::string> & itemVector, std::string itemLabel, const char * filename){
//find startData
//find endData
    itemVector.clear();
    std::string data_opticsStr("data_optics");
    std::ifstream file(filename);
    std::string str;
    unsigned long int cc = 0;
    bool got_data_optics = false;
    bool got_loop = false;
    bool got_data_optics_Header = false;
    std::vector<std::string> fields;
    std::vector<int> fieldsIdx;
    int targetID = -1;
    int opticsGroupID=-1;
    unsigned long int headerLines = 0;
    unsigned long int dataOpticsLines = 0;
    std::string resultStr="";
    while (std::getline(file, str) && ++cc < __MAX_STAR_HEADER_SIZE__){
       std::regex r("\\s+");
       str = std::regex_replace(str, r, ",");
       std::regex rNoSpaces("\\s+");
       std::string strNoSpaces = std::regex_replace(str, rNoSpaces, ",");
       if (strNoSpaces.length()<1 && !got_data_optics_Header){
         //do nothing
       }else if (str.substr(0,1).compare("#")==0 && !got_data_optics_Header){
         //it is a comment, do nothing
       }else if (!got_data_optics){
           if( strstr(strNoSpaces.substr(0,data_opticsStr.size()).c_str(),data_opticsStr.c_str()) ) {
            got_data_optics=true;
          }
       }else if (!got_loop){
          if( strstr(str.substr(0,5).c_str(),"loop_") ){
            got_loop=true;
            //std::cerr<<"got_loop\n";
          }
       }else if( strstr(strNoSpaces.substr(0,1).c_str(),"_") ){
        got_data_optics_Header = true;
        headerLines=cc+1;
        std::size_t found = str.find(",");
        std::string headerItem=str.substr(0,found);
        std::size_t start1=str.find("#");
        //std::size_t start1 = headerItem.find(",");
        std::string columnItem=str.substr(start1+1, str.length()-1);
        std::size_t start2=columnItem.find(",");
        columnItem=str.substr(start1+1, start2);
        int columnItemValue=atoi(columnItem.c_str())-1;
        fields.push_back(headerItem);
        fieldsIdx.push_back(columnItemValue);
        //std::cerr<<"--->"<<headerItem<<" column="<<columnItemValue<<"\n";
        //        fields.push_back(std::string("_rln"));
        //substr (size_t pos = 0, size_t len = npos)
        if (headerItem.compare(itemLabel)==0){
          targetID=columnItemValue;
        }
        if (headerItem.compare("_rlnOpticsGroupName")==0){
             opticsGroupID=columnItemValue;
        }

       }else if ( got_data_optics_Header == true && strNoSpaces.length()>1 && str.substr(0,1).compare("#")!=0){
           dataOpticsLines++;
       }else{
           cc = __MAX_STAR_HEADER_SIZE__ + 1;
       }
    }

    //std::cerr<<"dataOpticsLines="<<dataOpticsLines<<"\n";
    
    unsigned long int ccCounter=0;
    unsigned long int hhCounter=0;
    file.clear();
    file.seekg(0,std::ios::beg);
    //std::cerr<<"headerLines="<<headerLines<<"\n";
    while (std::getline(file, str) && ccCounter < __MAX_STAR_HEADER_SIZE__){
     if (++ccCounter >= headerLines){
        //remove head spaces
        str=str.erase(0, str.find_first_not_of(" \t\n\r\f\v"));
        if (str.size()>0 && ++hhCounter<=dataOpticsLines) { //check the line is not empty
                //std::cerr<<str;//<<"HUH\n";
                std::regex r("\\s+");
                str = std::regex_replace(str, r, ",");
                std::string result=findColumnItem(str, targetID);
                itemVector.push_back(result);
/*
                std::string opticsGroupTmp=findColumnItem(str, opticsGroupID);
                if (opticsGroupTmp.compare(OpticsGroupName)==0){
                    resultStr=result;
                }
 */
        }
      }
    }
    return dataOpticsLines;
}
    
int readStar(std::vector<double> & itemVector, std::string itemType, const char * filename){
    std::string dataLabelStr(relionDataStartLabel(filename));
    int dataLabelStrSize=dataLabelStr.size();
    
    std::ifstream file(filename);
    std::string str;
    unsigned long int cc = 0;
    bool got_data_images = false;
    bool got_loop = false;
    std::vector<std::string> fields;
    std::vector<int> fieldsIdx;
    int targetID = -1;
    unsigned long int headerLines = 0;
    while (std::getline(file, str) && ++cc < __MAX_STAR_HEADER_SIZE__){
       std::regex r("\\s+");
       str = std::regex_replace(str, r, ",");
       std::regex rNoSpaces("\\s+");
       std::string strNoSpaces = std::regex_replace(str, rNoSpaces, ",");    
       if (strNoSpaces.length()<1){
         //do nothing
       }else if (str.substr(0,1).compare("#")==0){
         //it is a comment, do nothing
       }else if (!got_data_images){
           if( strstr(strNoSpaces.substr(0,dataLabelStrSize).c_str(),dataLabelStr.c_str()) ) {
            got_data_images=true;
            //std::cerr<<"star file: data_images\n";
          }
       }else if (!got_loop){
          if( strstr(str.substr(0,5).c_str(),"loop_") ){
            got_loop=true;
            //std::cerr<<"got_loop\n";
          }
       }else if( strstr(strNoSpaces.substr(0,1).c_str(),"_") ){
        headerLines=cc+1;
        std::size_t found = str.find(",");
        std::string headerItem=str.substr(0,found);
        std::size_t start1=str.find("#");
        //std::size_t start1 = headerItem.find(",");
        std::string columnItem=str.substr(start1+1, str.length()-1);
        std::size_t start2=columnItem.find(",");
        columnItem=str.substr(start1+1, start2);
        int columnItemValue=atoi(columnItem.c_str())-1;
        fields.push_back(headerItem);
        fieldsIdx.push_back(columnItemValue);
        //std::cerr<<"--->"<<headerItem<<" column="<<columnItemValue<<"\n";
        //        fields.push_back(std::string("_rln"));
        //substr (size_t pos = 0, size_t len = npos)
        if (headerItem.compare(itemType)==0){
          targetID=columnItemValue;
        }
       }
    }

    unsigned long int ccCounter=0;
    file.clear();
    file.seekg(0,std::ios::beg);
    //std::cerr<<"headerLines="<<headerLines<<"\n";
    while (std::getline(file, str) ){
     if (++ccCounter >= headerLines){
		//remove head spaces
        str=str.erase(0, str.find_first_not_of(" \t\n\r\f\v"));
        if (str.size()>0) { //check the line is not empty        
                //std::cerr<<str;//<<"HUH\n";
                std::regex r("\\s+");
	            str = std::regex_replace(str, r, ",");
                std::string result=findColumnItem(str, targetID);
                itemVector.push_back(atof(result.c_str()));
        }
      }
    }
    return targetID;
}



int readStar(std::vector<std::string> & itemVector, std::string itemType, const char * filename){
    
    std::string dataLabelStr(relionDataStartLabel(filename));
    int dataLabelStrSize=dataLabelStr.size();
    
    std::ifstream file(filename);
    std::string str;
    unsigned long int cc = 0;
    bool got_data_images = false;
    bool got_loop = false;
    std::vector<std::string> fields;
    std::vector<int> fieldsIdx;
    int targetID = -1;
    unsigned long int headerLines = 0;
    while (std::getline(file, str) && ++cc < __MAX_STAR_HEADER_SIZE__){
       std::regex r("\\s+");
       str = std::regex_replace(str, r, ",");
       std::regex rNoSpaces("\\s+");
       std::string strNoSpaces = std::regex_replace(str, rNoSpaces, ",");    
       if (strNoSpaces.length()<1){
         //do nothing
       }else if (str.substr(0,1).compare("#")==0){
         //it is a comment, do nothing
       }else if (!got_data_images){
          if( strstr(strNoSpaces.substr(0,dataLabelStrSize).c_str(),dataLabelStr.c_str()) ) {
            got_data_images=true;
            //std::cerr<<"star file: data_images\n";
          }
       }else if (!got_loop){
          if( strstr(str.substr(0,5).c_str(),"loop_") ){
            got_loop=true;
            //std::cerr<<"got_loop\n";
          }
       }else if( strstr(strNoSpaces.substr(0,1).c_str(),"_") ){
        headerLines=cc+1;
        std::size_t found = str.find(",");
        std::string headerItem=str.substr(0,found);
        std::size_t start1=str.find("#");
        //std::size_t start1 = headerItem.find(",");
        std::string columnItem=str.substr(start1+1, str.length()-1);
        std::size_t start2=columnItem.find(",");
        columnItem=str.substr(start1+1, start2);
        int columnItemValue=atoi(columnItem.c_str())-1;
        fields.push_back(headerItem);
        fieldsIdx.push_back(columnItemValue);
        //std::cerr<<"--->"<<headerItem<<" column="<<columnItemValue<<"\n";
        //        fields.push_back(std::string("_rln"));
        //substr (size_t pos = 0, size_t len = npos)
        if (headerItem.compare(itemType)==0){
          targetID=columnItemValue;
        }
       }
    }

    unsigned long int ccCounter=0;
    file.clear();
    file.seekg(0,std::ios::beg);
    //std::cerr<<"headerLines="<<headerLines<<"\n";
    //        	std::cerr<<str<<"\n";

    while (std::getline(file, str) ){
     if (++ccCounter >= headerLines){
		//remove head spaces
		str=str.erase(0, str.find_first_not_of(" \t\n\r\f\v"));
        if (str.size()>0) { //check the line is not empty        
          std::regex r("\\s+");
	  str = std::regex_replace(str, r, ",");
	  //std::cerr<<"\n"<<str;
          std::string result=findColumnItem(str, targetID);
          itemVector.push_back(result);
        }
      }
    }
    return targetID;
}




// *****************
//
//  get max header label Number
//
// *****************
//
int getMaxHeaderLabelNumber (const char * filename){
  int maxHeaderLabel=0;
  std::string dataLabelStr(relionDataStartLabel(filename));
    int dataLabelStrSize=dataLabelStr.size();
    
    std::ifstream file(filename);
    std::string str;
    unsigned long int cc = 0;
    bool got_data_images = false;
    bool got_loop = false;
    int targetID = -1;
    unsigned long int headerLines = 0;
    while (std::getline(file, str) && ++cc < __MAX_STAR_HEADER_SIZE__){
       std::regex r("\\s+");
       str = std::regex_replace(str, r, ",");
       std::regex rNoSpaces("\\s+");
       std::string strNoSpaces = std::regex_replace(str, rNoSpaces, ",");    
       if (strNoSpaces.length()<1){
         //do nothing
       }else if (str.substr(0,1).compare("#")==0){
         //it is a comment, do nothing
       }else if (!got_data_images){
           if( strstr(strNoSpaces.substr(0,dataLabelStrSize).c_str(),dataLabelStr.c_str()) ) {
            got_data_images=true;
            //std::cerr<<"star file: data_images\n";
          }
       }else if (!got_loop){
          if( strstr(str.substr(0,5).c_str(),"loop_") ){
            got_loop=true;
            //std::cerr<<"got_loop\n";
          }
       }else if( strstr(strNoSpaces.substr(0,1).c_str(),"_") ){
        headerLines=cc+1;
        std::size_t found = str.find(",");
        std::string headerItem=str.substr(0,found);
        std::size_t start1=str.find("#");
        //std::size_t start1 = headerItem.find(",");
        std::string columnItem=str.substr(start1+1, str.length()-1);
        std::size_t start2=columnItem.find(",");
        columnItem=str.substr(start1+1, start2);
        int columnItemValue=atoi(columnItem.c_str());
        if (columnItemValue>maxHeaderLabel)
          maxHeaderLabel=columnItemValue;
        //maxHeaderLabel
        //std::cerr<<"--->"<<headerItem<<" column="<<columnItemValue<<"\n";
        //        fields.push_back(std::string("_rln"));
        //substr (size_t pos = 0, size_t len = npos)
       }
    }

    return maxHeaderLabel;


}





// *************************************
// Returns header from star file
//
std::string getStarHeader (char * starFile){

    std::string dataLabelStr(relionDataStartLabel(starFile));
    int dataLabelStrSize=dataLabelStr.size();
    std::string outHeader;
    std::ifstream file(starFile);
    std::string strLine;
    unsigned long int cc = 0;
    bool got_data_images = false;
    bool got_loop = false;
    std::vector<std::string> fields;
    std::vector<int> fieldsIdx;
    int columnItemValue = -1;
    unsigned long int headerLines = 0;
    int lastLabelPosition=0;

    while (std::getline(file, strLine) && ++cc < __MAX_STAR_HEADER_SIZE__){
        std::string tmpOutLine=strLine+std::string("\n");
        
        std::regex r("\\s+");
	std::string str0 = std::regex_replace(strLine, r, ",");
       std::regex rNoSpaces("\\s+");
	std::string strNoSpaces = std::regex_replace(strLine, rNoSpaces, ",");    
       if (strNoSpaces.length()<1){
         //do nothing
         if( !got_loop || !got_data_images ){
          headerLines++;
          outHeader+=tmpOutLine;
         }
       }else if (str0.substr(0,1).compare("#")==0){
         //it is a comment, do nothing
         headerLines++;
         outHeader+=tmpOutLine;
       }else if (!got_data_images){
          headerLines++;
          outHeader+=tmpOutLine;
          if( strstr(strNoSpaces.substr(0,dataLabelStrSize).c_str(),dataLabelStr.c_str()) ) {
            got_data_images=true;
          }
       }else if (!got_loop){
          headerLines++;
          outHeader+=tmpOutLine;
          if( strstr(str0.substr(0,5).c_str(),"loop_") ){
            got_loop=true;
            //std::cerr<<"got_loop\n";
          }
       }else if( strstr(str0.substr(0,1).c_str(),"_") ){
        headerLines++;
        outHeader+=tmpOutLine;
        lastLabelPosition=headerLines;
       }
    }
    return outHeader;
}




void getStarHeaders(std::vector<std::string> & fields, std::vector<int> & fieldsIdx, const char * filename){
    std::string dataLabelStr(relionDataStartLabel(filename));
    int dataLabelStrSize=dataLabelStr.size();
    
    //std::vector<double> itemVector;
    std::ifstream file(filename);
    std::string strLine;
    unsigned long int cc = 0;
    bool got_data_images = false;
    bool got_loop = false;
    int targetID = -1;
    unsigned long int headerLines = 0;

    while (std::getline(file, strLine) && ++cc < __MAX_STAR_HEADER_SIZE__){
        std::regex rNoSpaces("\\s+");
	std::string strNoSpaces = std::regex_replace(strLine, rNoSpaces, ",");
        std::regex r("\\s+");
	std::string str0 = std::regex_replace(strLine, r, ",");
       if ( strNoSpaces.length()<1 ){
         //do nothing
       }else if ( strNoSpaces.substr(0,1).compare("#")==0 ){
         //it is a comment, do nothing
       }else if (!got_data_images){
          if( strstr(strNoSpaces.substr(0,dataLabelStrSize).c_str(),dataLabelStr.c_str()) ) {
            got_data_images=true;
          }
       }else if (!got_loop){
          if( strstr(str0.substr(0,5).c_str(),"loop_") ){
            got_loop=true;
          }
       }else if( strstr(str0.substr(0,1).c_str(),"_") ){
        headerLines=cc+1;
        std::size_t found = str0.find(",");
        std::string headerItem=str0.substr(0,found);
        std::size_t start1=str0.find("#");
        //std::size_t start1 = headerItem.find(",");
        std::string columnItem=str0.substr(start1+1, str0.length()-1);
        std::size_t start2=columnItem.find(",");
        columnItem=str0.substr(start1+1, start2);
        int columnItemValue=atoi(columnItem.c_str())-1;
        fields.push_back(headerItem);
        fieldsIdx.push_back(columnItemValue);
       }else{
         cc=__MAX_STAR_HEADER_SIZE__+1;
       }
    }
}

int getStarHeaderItemIdx(std::string targetStr, const char * filename){
        int targetIndex = -1;
        std::vector<std::string> fields;
        std::vector<int> fieldsIdx;
        getStarHeaders( fields, fieldsIdx, filename);
        for (unsigned long int ii=0;ii<fields.size();ii++){
          if ( targetStr.compare(fields[ii])==0 ){
            targetIndex=fieldsIdx[ii];
          }
        }
        return targetIndex;
}

bool checkLabelExists(std::string label, const char * filename){
          std::vector<std::string> headerFields;
          std::vector<int> fieldsIdx;
          getStarHeaders(headerFields, fieldsIdx, filename);
          if ( std::find(headerFields.begin(), headerFields.end(), label) == headerFields.end() ){
            //std::cerr<<"ERROR: no occurrence of label "<< labelNormalizedScoresForSubmaps<<" in file  "<<inputStarFileForSubmaps<<" ... EXIT! \n";
            //exit(1);
            return false;
          }
          return true;
}


// *****************************
// getStarStart
long int getStarStart(const char * filenameIn){
  std::string startHeader = getStarHeader ( (char *)filenameIn);
  int countLines=std::count(startHeader.begin(), startHeader.end(), '\n');
  return countLines;

}

/*
// *****************************
// getStarStart
long int getStarStart(const char * filenameIn){
    std::string dataLabelStr(relionDataStartLabel(filenameIn));
    int dataLabelStrSize=dataLabelStr.size();
    
    std::ifstream file(filenameIn);
    std::string strLine;
    unsigned long int cc = 0;
    bool got_data_images = false;
    bool got_loop = false;
    std::vector<std::string> fields;
    std::vector<int> fieldsIdx;
    int columnItemValue = -1;
    unsigned long int headerLines = 0;

    while (std::getline(file, strLine) && ++cc < __MAX_STAR_HEADER_SIZE__){
        std::regex r("\\s+");
	std::string str0 = std::regex_replace(strLine, r, ",");
       std::regex rNoSpaces("\\s+");
	std::string strNoSpaces = std::regex_replace(strLine, rNoSpaces, ",");    
       if (strNoSpaces.length()<1){
         //do nothing
         headerLines++;
       }else if (str0.substr(0,1).compare("#")==0){
         //it is a comment, do nothing
         headerLines++;
       }else if (!got_data_images){
          headerLines++;
          if( strstr(strNoSpaces.substr(0,dataLabelStrSize).c_str(),dataLabelStr.c_str()) ) {
            got_data_images=true;
          }
       }else if (!got_loop){
          headerLines++;
          if( strstr(str0.substr(0,5).c_str(),"loop_") ){
            got_loop=true;
            //std::cerr<<"got_loop\n";
          }
       }else if( strstr(str0.substr(0,1).c_str(),"_") ){
        headerLines++;
       }
    }

    return headerLines;
}
*/

// *****************************
// getNumItemsStar
int getNumItemsStar(const char * filename){

    long int headerLines = getStarStart(filename)+1;
    std::ifstream file(filename);
    unsigned long int ccCounter=0;
    unsigned long int lineCounter=0;
    file.seekg(0,std::ios::beg);
    std::string str;
    while (std::getline(file, str) ){
     if (++ccCounter >= headerLines){
        str=str.erase(0, str.find_first_not_of(" \t\n\r\f\v"));
        if (str.size()>0) { //check the line is not empty        
                lineCounter++;
        }
      }
    }
    return lineCounter;

}


unsigned long int addHeaderFieldsNewStar(std::vector<std::string> newFields, const char * filenameIn, const char * filenameOut){

    //check the strings
    std::string fileInStr(filenameIn);
    std::string fileOutStr(filenameOut);
    std::string tmpFile="";
    if (fileInStr.compare(filenameOut)==0){
      tmpFile=generateTmpFilename(".star");
      copyCvsFile(fileInStr.c_str(), tmpFile.c_str() );
      fileInStr=tmpFile;
    }

    std::string dataLabelStr(relionDataStartLabel(fileInStr.c_str()));
    int dataLabelStrSize=dataLabelStr.size();
    
    std::ifstream file(fileInStr.c_str());
    std::ofstream fileOut;
    fileOut.open(filenameOut);
    fileOut.close();
    fileOut.open (filenameOut, std::ofstream::out | std::ofstream::app);    
    std::string strLine;
    unsigned long int cc = 0;
    bool got_data_images = false;
    bool got_loop = false;
    std::vector<std::string> fields;
    std::vector<int> fieldsIdx;
    int columnItemValue = -1;
    unsigned long int headerLines = 0;

    while (std::getline(file, strLine) && ++cc < __MAX_STAR_HEADER_SIZE__){
        std::regex r("\\s+");
	std::string str0 = std::regex_replace(strLine, r, ",");
       std::regex rNoSpaces("\\s+");
	std::string strNoSpaces = std::regex_replace(strLine, rNoSpaces, ",");    
       if (strNoSpaces.length()<1){
         //do nothing
         fileOut<<strLine<<"\n";
       }else if (str0.substr(0,1).compare("#")==0){
         //it is a comment, do nothing
         fileOut<<strLine<<"\n";
       }else if (!got_data_images){
          fileOut<<strLine<<"\n";
          if( strstr(strNoSpaces.substr(0,dataLabelStrSize).c_str(),dataLabelStr.c_str()) ) {
            got_data_images=true;
            //std::cerr<<"star file: data_images\n";
          }
       }else if (!got_loop){
          fileOut<<strLine<<"\n";
          if( strstr(str0.substr(0,5).c_str(),"loop_") ){
            got_loop=true;
            //std::cerr<<"got_loop\n";
          }
       }else if( strstr(str0.substr(0,1).c_str(),"_") ){
        fileOut<<strLine<<"\n";
        headerLines=cc+1;
        std::size_t found = str0.find(",");
        std::string headerItem=str0.substr(0,found);
        std::size_t start1=str0.find("#");
        //std::size_t start1 = headerItem.find(",");
        std::string columnItem=str0.substr(start1+1, str0.length()-1);
        std::size_t start2=columnItem.find(",");
        columnItem=str0.substr(start1+1, start2);
        columnItemValue=atoi(columnItem.c_str())-1;
       }
    }
    columnItemValue++;
    for (unsigned long int ii=0; ii<newFields.size(); ii++ ){      
      if ( getStarHeaderItemIdx(newFields[ii], fileInStr.c_str()) < 0 ){
       columnItemValue++;
       headerLines++;
       fileOut<<newFields[ii]<< " #"<< std::to_string(columnItemValue) <<" \n";
      }
    }
    fileOut.close();
    if (tmpFile.length()>0){
      removeCvsFile(tmpFile.c_str());
    }

    return headerLines;
}




unsigned long int removeLinesStar(std::string itemType, std::string valueToRetain, const char * filenameIn, const char * filenameOut){

            std::vector<std::string> subset;
            readStar(subset, itemType, filenameIn);
    
            std::ofstream fileOutput;
            fileOutput.open(filenameOut);
            fileOutput.close();
            fileOutput.open (filenameOut, std::ofstream::out | std::ofstream::app);
            long int startMicrograph=getStarStart(filenameIn);
            std::ifstream fileParticles(filenameIn);
            std::string strLine;
            for (int ccc=0;ccc<startMicrograph;ccc++){
               std::getline(fileParticles, strLine);
               fileOutput << strLine <<"\n";
            }
          unsigned long int counter = 0;
             for (unsigned long int ii=0; ii<subset.size(); ii++){
               std::getline(fileParticles, strLine);
                 
               if ( subset[ii] ==  valueToRetain){
                   std::string inputStr(strLine);
                   fileOutput << inputStr << "\n";
               }
              }
        
        fileOutput.close();
    return subset.size();
    
}
    








// ******************
//extract blockFields
std::vector<std::string> extractBlockLabels(std::string startLabel, std::string endLabel, const char * starFilename){
        std::vector<std::string> returnVector;
        std::ofstream fileOutput;
        int startIdx=getStarHeaderItemIdx(startLabel, starFilename);
	int endIdx=getStarHeaderItemIdx(endLabel, starFilename);


        long int startMicrograph=getStarStart(starFilename);
        std::ifstream fileRefinedParticles(starFilename);
        std::string strLine;
        for (int counter=0;counter<startMicrograph;counter++){
           std::getline(fileRefinedParticles, strLine);
        }
        bool done=false;
        for (unsigned long int ii=0,counter=1;!done; ii++){
           std::getline(fileRefinedParticles, strLine);
           if(!fileRefinedParticles.eof()){
             std::string stripStr = std::regex_replace(strLine, std::regex("\\s+"), "");
             if(stripStr.size()>0){
              
              std::vector<long int> beginVec=getStringBeginAndLenghtAtStarPosition(strLine, startIdx);
              std::vector<long int> endVec=getStringBeginAndLenghtAtStarPosition(strLine, endIdx);
              int startPoint = beginVec[0];
              int totalLength=endVec[0]-beginVec[0]+endVec[1];
              //std::cerr<< strLine.substr(startPoint,totalLength) << "\n";
              returnVector.push_back(strLine.substr(startPoint,totalLength));
             }
           }else{
             done=true;
           }
        }
        return returnVector;
}


int findColumnItemPosition(std::string str, int idx){
          int beforePos;
          long int posStart=-1;
          long int  posEnd=str.length()-1;
          for (int ii=0; ii<idx; ii++){
           posStart+=1+str.substr(posStart+1,posEnd).find(",");
          }
          posEnd=str.substr(posStart+1,posEnd).find(",");
          beforePos = posStart+1;
          return beforePos;
}



unsigned long int StarMaxFieldIdx(const char * filenameIn){
    std::string dataLabelStr(relionDataStartLabel(filenameIn));
    int dataLabelStrSize=dataLabelStr.size();
    
    std::ifstream file(filenameIn);
    std::string strLine;
    unsigned long int cc = 0;
    bool got_data_images = false;
    bool got_loop = false;
    std::vector<std::string> fields;
    std::vector<int> fieldsIdx;
    int columnItemValue = -1;
    unsigned long int headerLines = 0;

    long int maxItemValue = -1;
    while (std::getline(file, strLine) && ++cc < __MAX_STAR_HEADER_SIZE__){
       std::regex r("\\s+");
       std::string str0 = std::regex_replace(strLine, r, ",");
       std::regex rNoSpaces("\\s+");
       std::string strNoSpaces = std::regex_replace(strLine, rNoSpaces, ",");    
       if (strNoSpaces.length()<1){
         //do nothing
       }else if (str0.substr(0,1).compare("#")==0){
         //it is a comment, do nothing
       }else if (!got_data_images){
          if( strstr(strNoSpaces.substr(0,dataLabelStrSize).c_str(),dataLabelStr.c_str()) ) {
            got_data_images=true;
          }
       }else if (!got_loop){
          if( strstr(str0.substr(0,5).c_str(),"loop_") ){
            got_loop=true;
          }
       }else if( strstr(str0.substr(0,1).c_str(),"_") ){
        headerLines=cc+1;
        std::size_t found = str0.find(",");
        std::string headerItem=str0.substr(0,found);
        std::size_t start1=str0.find("#");
        std::string columnItem=str0.substr(start1+1, str0.length()-1);
        std::size_t start2=columnItem.find(",");
        columnItem=str0.substr(start1+1, start2);
        columnItemValue=atoi(columnItem.c_str())-1;
        if (columnItemValue > maxItemValue){
         maxItemValue=columnItemValue;
        }
       }
    }
    return maxItemValue+1;
}



unsigned long int replaceAddValueStar(std::string itemType, std::vector<std::string> replacingValuesString, const char * filenameIn, const char * filenameOut){
    std::string dataLabelStr(relionDataStartLabel(filenameIn));
    int dataLabelStrSize=dataLabelStr.size();
    //int numExistingLabels = StarMaxFieldIdx( filenameIn );
    //removeCvsFile(filenameOut);



    //check the strings
    std::string fileInStr(filenameIn);
    std::string fileOutStr(filenameOut);
    std::string tmpFile="";
    if (fileInStr.compare(filenameOut)==0 ){
      tmpFile=generateTmpFilename(".star");
      copyCvsFile(filenameIn, tmpFile.c_str() );
      fileInStr=tmpFile;
    }
    std::ifstream file(fileInStr);
    std::ofstream fileOut;
    fileOut.open(filenameOut);
    fileOut.close();
    fileOut.open (filenameOut, std::ofstream::out | std::ofstream::app);    
    std::string strLine;
    unsigned long int cc = 0;
    bool got_data_images = false;
    bool got_loop = false;
    std::vector<std::string> fields;
    std::vector<int> fieldsIdx;
    int targetID = -1;
    unsigned long int headerLines = 0;

    //input file is empty
/*    
    if (numExistingLabels==0){
      std::cerr<<"no labels!"<<"\n";
      fileOut<<strLine<<"data_\n\nloop_\n";
      fileOut<<itemType<<"\n";
      for (unsigned long int ii=0;ii<replacingValuesString.size(); ii++){
        fileOut<<replacingValuesString[ii]<<"\n";
      }
      fileOut.close();
      file.close();
      if (tmpFile.length()>0){
          removeCvsFile(tmpFile.c_str());
      }
      return itemType.size();
    }
    
    else if (numExistingLabels==1){
      //std::string tmpHeader = getStarHeader (filenameIn);
      if (checkLabelExists(itemType, fileInStr.c_str()) ){
        std::cerr<<"exists!"<<"\n";
        fileOut<<strLine<<"data_\n\nloop_\n";
        fileOut<<itemType<<" #1\n";
        for (unsigned long int ii=0;ii<replacingValuesString.size(); ii++){
          fileOut<<replacingValuesString[ii]<<"\n";
        }
        fileOut.close();
        file.close();
        if (tmpFile.length()>0){
          removeCvsFile(tmpFile.c_str());
        }
        return itemType.size();
      }else{
        std::cerr<<"does not exists!"<<"\n";
        std::vector<std::string> fields;
        std::vector<int> fieldsIdx;
        getStarHeaders(fields, fieldsIdx, filenameIn);
        fileOut<<"data_\n\nloop_\n";
        fileOut<<fields[0]<<"\n";
        std::vector<std::string> existingItem;
        readStar(existingItem, fields[0], filenameIn);
        for (unsigned long int ii=0;ii<existingItem.size() && ii<replacingValuesString.size(); ii++){
          fileOut<<existingItem[ii]<<"  "<<replacingValuesString[ii]<<"\n";
        }
        fileOut.close();
        file.close();
        if (tmpFile.length()>0){
          removeCvsFile(tmpFile.c_str());
        }
        return itemType.size();
      }

    }
*/


    bool newLabel = false;
    if ( getStarHeaderItemIdx(itemType, fileInStr.c_str()) < 0 ){
      newLabel = true;
    }
    
    while (std::getline(file, strLine) && ++cc < __MAX_STAR_HEADER_SIZE__){
        std::regex r("\\s+");
	std::string str0 = std::regex_replace(strLine, r, ",");
        std::regex rNoSpaces("\\s+");
	std::string strNoSpaces = std::regex_replace(strLine, rNoSpaces, ",");    
       if (strNoSpaces.length()<1){
         //do nothing
         fileOut<<strLine<<"\n";
       }else if (str0.substr(0,1).compare("#")==0){
         //it is a comment, do nothing
         fileOut<<strLine<<"\n";
       }else if (!got_data_images){
          fileOut<<strLine<<"\n";
          if( strstr(strNoSpaces.substr(0,dataLabelStrSize).c_str(),dataLabelStr.c_str()) ) {
            got_data_images=true;
            //std::cerr<<"star file: data_images\n";
          }
       }else if (!got_loop){
          fileOut<<strLine<<"\n";
          if( strstr(str0.substr(0,5).c_str(),"loop_") ){
            got_loop=true;
            //std::cerr<<"got_loop\n";
          }
       }else if( strstr(str0.substr(0,1).c_str(),"_") ){ //was _rln
        fileOut<<strLine<<"\n";
        headerLines=cc+1;
        std::size_t found = str0.find(",");
        std::string headerItem=str0.substr(0,found);
        std::size_t start1=str0.find("#");
        //std::size_t start1 = headerItem.find(",");
        std::string columnItem=str0.substr(start1+1, str0.length()-1);
        std::size_t start2=columnItem.find(",");
        columnItem=str0.substr(start1+1, start2);
        int columnItemValue=atoi(columnItem.c_str())-1;
        fields.push_back(headerItem);
        fieldsIdx.push_back(columnItemValue);
        if (headerItem.compare(itemType)==0){
          targetID=columnItemValue;
        }
       }
    }

    if ( newLabel ){
      int nextIdx = StarMaxFieldIdx(fileInStr.c_str())+1;
      fileOut<<itemType<< " #"<< std::to_string(nextIdx)  <<" \n";
    }  

    unsigned long int ccCounter=0;
    file.clear();
    file.seekg(0,std::ios::beg);
    unsigned long int counter = 0;
    while (std::getline(file, strLine) ){
     if (++ccCounter >= headerLines){
        std::regex r("\\s+");
	std::string str0 = std::regex_replace(strLine, r, ",");
        std::string result=findColumnItem(str0, targetID);

        std::regex r1("[,]+");
	std::string str1 = std::regex_replace(str0,r1," ");

        std::regex r2("\\s+");
	std::string str2 = std::regex_replace(str1, r2, "");        

        if ( newLabel ){
           fileOut<<str1<<" ";
                if (counter<replacingValuesString.size()){
                  fileOut<<replacingValuesString[counter];
                }else{
                  if (str2.size()>0){
                   fileOut<<"0";
                  }
                }
            fileOut<<" ";
            counter++;
        }else{
                int beforePos=findColumnItemPosition(str0, targetID);
	        if (beforePos>=0){
          	  std::string beforeStr = str1.substr(0,beforePos);
	          fileOut<<beforeStr;
	        }
	        if (beforePos>0){
	          fileOut<<" ";
	        }
                //OK fileOut<<replacingString;
                if (counter<replacingValuesString.size()){
                  fileOut<<replacingValuesString[counter];
                }else{
                  if (str2.size()>0){
                   fileOut<<"0";
                  }
                }
	        fileOut<<" ";
	        counter++;
	        std::string afterStr = str1.substr( beforePos+result.length() );
	        fileOut<<afterStr;
        }
        fileOut<<"\n";
      }
    }
    fileOut.close();

    if (tmpFile.length()>0){
      removeCvsFile(tmpFile.c_str());
    }

    return counter;
}

unsigned long int replaceAddValueStar(std::string itemType, std::vector<double> replacingValuesDouble, const char * filenameIn, const char * filenameOut){
    std::vector<std::string> replacingValuesString;
    for (unsigned long int ii=0; ii<replacingValuesDouble.size();ii++){
        replacingValuesString.push_back(std::to_string(replacingValuesDouble[ii]));
    }
    return replaceAddValueStar(itemType,replacingValuesString,filenameIn,filenameOut);
}

unsigned long int replaceAddValueStar(std::string itemType, std::vector<int> replacingValuesInt, const char * filenameIn, const char * filenameOut){
    std::vector<std::string> replacingValuesString;
    for (unsigned long int ii=0; ii<replacingValuesInt.size();ii++){
        replacingValuesString.push_back(std::to_string(replacingValuesInt[ii]));
    }
    return replaceAddValueStar(itemType,replacingValuesString,filenameIn,filenameOut);
}



unsigned long int removeColumnStar(std::string itemType, const char * filenameIn, const char * filenameOut){
    std::string dataLabelStr(relionDataStartLabel(filenameIn));
    int dataLabelStrSize=dataLabelStr.size();

    //check the strings
    std::string fileInStr(filenameIn);
    std::string fileOutStr(filenameOut);
    bool createdTmpFile = false;
    std::string tmpFile="";
    if (fileInStr.compare(filenameOut)==0){
      tmpFile=generateTmpFilename(".star");
      copyCvsFile(filenameIn, tmpFile.c_str() );
      fileInStr=tmpFile;
      createdTmpFile=true;
    }
    std::ifstream file(fileInStr);
    std::ofstream fileOut;
    fileOut.open(filenameOut);
    fileOut.close();
    fileOut.open (filenameOut, std::ofstream::out | std::ofstream::app);    
    std::string strLine;
    unsigned long int cc = 0;
    bool got_data_images = false;
    bool got_loop = false;
    std::vector<std::string> fields;
    std::vector<int> fieldsIdx;
    int targetID = -1;
    unsigned long int headerLines = 0;

    bool newLabel = false;
    if ( getStarHeaderItemIdx(itemType, fileInStr.c_str()) < 0 ){
      newLabel = true;
      copyCvsFile(fileInStr.c_str(), filenameOut);
      if (createdTmpFile){
        removeCvsFile(tmpFile.c_str());
      }
      return 0;
    }
    

    while (std::getline(file, strLine) && ++cc < __MAX_STAR_HEADER_SIZE__){
        std::regex r("\\s+");
	std::string str0 = std::regex_replace(strLine, r, ",");
        std::regex rNoSpaces("\\s+");
	std::string strNoSpaces = std::regex_replace(strLine, rNoSpaces, ",");    
       if (strNoSpaces.length()<1){
         //do nothing
         fileOut<<strLine<<"\n";
       }else if (str0.substr(0,1).compare("#")==0){
         //it is a comment, do nothing
         fileOut<<strLine<<"\n";
       }else if (!got_data_images){
          fileOut<<strLine<<"\n";
          if( strstr(strNoSpaces.substr(0,dataLabelStrSize).c_str(),dataLabelStr.c_str()) ) {
            got_data_images=true;
            //std::cerr<<"star file: data_images\n";
          }
       }else if (!got_loop){
          fileOut<<strLine<<"\n";
          if( strstr(str0.substr(0,5).c_str(),"loop_") ){
            got_loop=true;
            //std::cerr<<"got_loop\n";
          }
       }else if( strstr(str0.substr(0,1).c_str(),"_") ){ //was _rln
        headerLines=cc+1;
        std::size_t found = str0.find(",");
        std::string headerItem=str0.substr(0,found);
        std::size_t start1=str0.find("#");
        //std::size_t start1 = headerItem.find(",");
        std::string columnItem=str0.substr(start1+1, str0.length()-1);
        std::size_t start2=columnItem.find(",");
        columnItem=str0.substr(start1+1, start2);
        int columnItemValue=atoi(columnItem.c_str())-1;
        fields.push_back(headerItem);
        fieldsIdx.push_back(columnItemValue);
        if (headerItem.compare(itemType)==0){
          targetID=columnItemValue;
        }else{
          fileOut<<strLine<<"\n";
        }
       }
    }

    /*
    if ( newLabel ){
      int nextIdx = StarMaxFieldIdx(fileInStr.c_str())+1;
      fileOut<<itemType<< " #"<< std::to_string(nextIdx)  <<" \n";
    }   
    */
    unsigned long int ccCounter=0;
    file.clear();
    file.seekg(0,std::ios::beg);
    unsigned long int counter = 0;
    while (std::getline(file, strLine) ){
     if (++ccCounter >= headerLines){
        std::regex r("\\s+");
	std::string str0 = std::regex_replace(strLine, r, ",");
        std::string result=findColumnItem(str0, targetID);

        std::regex r1("[,]+");
	std::string str1 = std::regex_replace(str0,r1," ");

        //std::regex r2("\\s+");
	//std::string str2 = std::regex_replace(str1, r2, "");        

       
                int beforePos=findColumnItemPosition(str0, targetID);
	        if (beforePos>=0){
          	  std::string beforeStr = str1.substr(0,beforePos);
	          fileOut<<beforeStr;
	        }
	        //if (beforePos>0){
	        //  fileOut<<" ";
	        //}
                //OK fileOut<<replacingString;
                /*
                //not add current value
                if (counter<replacingValuesString.size()){
                  fileOut<<replacingValuesString[counter];
                }else{
                  if (str2.size()>0){
                   fileOut<<"0";
                  }
                }
                */
	        fileOut<<" ";
	        counter++;
	        std::string afterStr = str1.substr( beforePos+result.length() );
	        fileOut<<afterStr;
        
        fileOut<<"\n";
      }
    }
    fileOut.close();

   if (createdTmpFile){
     removeCvsFile(tmpFile.c_str());
   }

    return counter;
}





bool isStarFile (const char * filename){
  std::string filenameStr (filename);
  std::string starExt (".star");
  bool isStarFileBool = false;
  if (filenameStr.length()>=starExt.length()){
    if ( starExt.compare(filenameStr.substr(filenameStr.length()-5))==0 ){
      isStarFileBool = true;
    }
  }
  return isStarFileBool;
}


//
//
//
/*
template<typename T>
long int refineStarFileByIndexes(std::vector<T> indexesLabelled, const char * filenameIn, const char * filenameOut, int selection=-1){
    std::string dataLabelStr(relionDataStartLabel(filenameIn));
    int dataLabelStrSize=dataLabelStr.size();
    
    //selection =-1 all
    //selection = 0 even
    //selection = 1 odd
    std::string itemType;
    std::ifstream file(filenameIn);
    std::ofstream fileOut;
    fileOut.open(filenameOut);
    fileOut.close();
    fileOut.open (filenameOut, std::ofstream::out | std::ofstream::app);    
    std::string strLine;
    unsigned long int cc = 0;
    bool got_data_images = false;
    bool got_loop = false;
    std::vector<std::string> fields;
    std::vector<int> fieldsIdx;
    int targetID = -1;
    unsigned long int headerLines = 0;

    while (std::getline(file, strLine) && ++cc < __MAX_STAR_HEADER_SIZE__){
       std::regex r("\\s+");
       std::string str0 = std::regex_replace(strLine, r, ",");
       std::regex rNoSpaces("\\s+");
       std::string strNoSpaces = std::regex_replace(strLine, rNoSpaces, ",");    
       if (strNoSpaces.length()<1){
         //do nothing
         fileOut<<strLine<<"\n";
       }else if (str0.substr(0,1).compare("#")==0){
         //it is a comment, do nothing
         fileOut<<strLine<<"\n";
       }else if (!got_data_images){
          fileOut<<strLine<<"\n";
          if( strstr(strNoSpaces.substr(0,dataLabelStrSize).c_str(),dataLabelStr.c_str()) ) {
            got_data_images=true;
            //std::cerr<<"star file: data_images\n";
          }
       }else if (!got_loop){
          fileOut<<strLine<<"\n";
          if( strstr(str0.substr(0,5).c_str(),"loop_") ){
            got_loop=true;
            //std::cerr<<"got_loop\n";
          }
       }else if( strstr(str0.substr(0,1).c_str(),"_") ){ //was _rln
        fileOut<<strLine<<"\n";
        headerLines=cc+1;
        std::size_t found = str0.find(",");
        std::string headerItem=str0.substr(0,found);
        std::size_t start1=str0.find("#");
        //std::size_t start1 = headerItem.find(",");
        std::string columnItem=str0.substr(start1+1, str0.length()-1);
        std::size_t start2=columnItem.find(",");
        columnItem=str0.substr(start1+1, start2);
        int columnItemValue=atoi(columnItem.c_str())-1;
        fields.push_back(headerItem);
        fieldsIdx.push_back(columnItemValue);
        if (headerItem.compare("_rlnImageName")==0){
          targetID=columnItemValue;
        }        
       }
    }

    unsigned long int ccCounter=0;
    file.clear();
    file.seekg(0,std::ios::beg);
    //std::cerr<<"headerLines="<<headerLines<<"\n";
    unsigned long int counter = 0;
  
    unsigned long int index=0;
    while (std::getline(file, strLine) ){
    if (++ccCounter >= headerLines){
         if ( indexesLabelled[index] > 0){
           counter++;
           if (selection == -1){
             fileOut<<strLine<<"\n";
           }else if (selection == 0){
             if (counter%2==0){
               fileOut<<strLine<<"\n";
             }
           }else if (selection == 1){
             if (counter%2==1){
               fileOut<<strLine<<"\n";
             }
           }
         }
         index++;
        }
    }

    fileOut.close();
    return counter;
}
*/



template<typename T>
long int refineStarFileByIndexes(std::vector<T> indexesLabelled, const char * filenameIn, const char * filenameOut, int selection=-1, bool useStarFileRandomSubset = false){
    
    std::vector<double> subset;
    if (useStarFileRandomSubset){
        readStar(subset, "_rlnRandomSubset", filenameIn);
        if ( subset.size() < 0 ){
            useStarFileRandomSubset=false;
        }
    }

    
    //if ( subset.size() !=  indexesLabelled.size() ){
    //}
    std::string dataLabelStr(relionDataStartLabel(filenameIn));
    int dataLabelStrSize=dataLabelStr.size();

    //get the file start
    //for each line get where the filename image start and ends
    //replace with the new characters
       //int labelIdx=getStarHeaderItemIdx("_rlnImageName", starFileIn);
       //std::cerr<< "labelIdx="<<labelIdx<<"\n";
       std::ofstream fileOutput;
        fileOutput.open(filenameOut);
        fileOutput.close();
        fileOutput.open (filenameOut, std::ofstream::out | std::ofstream::app);
        long int startMicrograph=getStarStart(filenameIn);
        std::ifstream fileParticles(filenameIn);
        std::string strLine;
        for (int ccc=0;ccc<startMicrograph;ccc++){
           std::getline(fileParticles, strLine);
           fileOutput << strLine <<"\n";
        }
      unsigned long int counter = 0;
       if ( selection < 0 ){
         for (unsigned long int ii=0; ii<indexesLabelled.size(); ii++){
           std::getline(fileParticles, strLine);
           if ( (int)indexesLabelled[ii] > 0 ){
               std::string inputStr(strLine);
               fileOutput << inputStr << "\n";
           }
          }
       }else if ( useStarFileRandomSubset ){
//           std::cerr<<"predefined halfmaps\n";
           for (unsigned long int ii=0; ii < indexesLabelled.size() ; ii++){
             std::getline(fileParticles, strLine);
             if ( (int)indexesLabelled[ii] > 0 && (int)subset[ii]==selection+1){
                 counter++;
                 std::string inputStr(strLine);
                 fileOutput << inputStr << "\n";
             }
            }
       } else {
//           std::cerr<<"induce halfmaps\n";
           for (unsigned long int ii=0, cc=0; ii < indexesLabelled.size() ; ii++){
             std::getline(fileParticles, strLine);
             if ( (int)indexesLabelled[ii] > 0){
                 std::string inputStr(strLine);
                 cc++;
                 if ( cc%2 == selection ){
                     counter++;
                     fileOutput<<inputStr<<"\n";
                 }
             }
            }
        }
    
    fileOutput.close();
    return counter;
}



template<typename T>
long int refineVemFileByIndexes(std::vector<T> indexesLabelled, const char * filenameIn, const char * filenameOut){
    std::string itemType;
    std::ifstream file(filenameIn);
    std::ofstream fileOut;
    fileOut.open(filenameOut);
    fileOut.close();
    fileOut.open (filenameOut, std::ofstream::out | std::ofstream::app);    
    std::string strLine;
    std::getline(file, strLine);
    fileOut<<strLine<<"\n";
    unsigned long int counter = 0;
    while (std::getline(file, strLine) ){
         if ( indexesLabelled[counter] > 0){
           fileOut<<strLine<<"\n";
         }
         counter++;
    }
    fileOut.close();
    return counter;
}





void replaceAddMetadataColumn( std::string itemType,  std::vector<std::string> & itemVector, const char * filenameIn, const char * filenameOut){
  if (isStarFile(filenameIn)){
    replaceAddValueStar(itemType, itemVector, filenameIn, filenameOut);
  }else{
    replaceAddCvsColumn(itemType, itemVector, filenameIn);
  }
}



// ********************************
//  Count items
// ********************************
unsigned long int countParticlesItems(const char * filenameIn){
    std::string dataLabelStr(relionDataStartLabel(filenameIn));
    int dataLabelStrSize=dataLabelStr.size();
    
    std::ifstream file(filenameIn);
    std::string strLine;
    unsigned long int cc = 0;
    bool got_data_images = false;
    bool got_loop = false;
    std::vector<std::string> fields;
    std::vector<int> fieldsIdx;
    int targetID = -1;
    unsigned long int headerLines = 0;
    std::string header="#\n";
    while (std::getline(file, strLine) && ++cc < __MAX_STAR_HEADER_SIZE__){
        std::regex r("\\s+");
	std::string str0 = std::regex_replace(strLine, r, ",");
        std::regex rNoSpaces("\\s+");
	std::string strNoSpaces = std::regex_replace(strLine, rNoSpaces, ",");    
       if (strNoSpaces.length()<1){
         //do nothing
       }else if (str0.substr(0,1).compare("#")==0){
         //it is a comment, do nothing
       }else if (!got_data_images){
          //fileOut<<strLine<<"\n";
          //std::cerr<<strstr(strNoSpaces.substr(0,5).c_str(),"data_") << "<<-----\n";
          if( strstr(strNoSpaces.substr(0,dataLabelStrSize).c_str(),dataLabelStr.c_str()) ) {
            got_data_images=true;
            //std::cerr<<"star file: data_images\n";
          }
       }else if (!got_loop){
          //fileOut<<strLine<<"\n";
          if( strstr(str0.substr(0,5).c_str(),"loop_") ){
            got_loop=true;
            //std::cerr<<"got_loop\n";
          }
       }else if( strstr(str0.substr(0,1).c_str(),"_") ){ //was _rln
        //fileOut<<strLine<<"\n";
        headerLines=cc+1;
        std::size_t found = str0.find(",");
        std::string headerItem=str0.substr(0,found);
        std::size_t start1=str0.find("#");
        //std::size_t start1 = headerItem.find(",");
        std::string columnItem=str0.substr(start1+1, str0.length()-1);
        std::size_t start2=columnItem.find(",");
        columnItem=str0.substr(start1+1, start2);
        int columnItemValue=atoi(columnItem.c_str())-1;
        if (headerItem.compare("_rlnAngleRot")==0){
           headerItem=std::string("Phi");
        }else if (headerItem.compare("_rlnAngleTilt")==0){
           headerItem=std::string("Theta");
        } else if (headerItem.compare("_rlnAnglePsi")==0){
           headerItem=std::string("Psi");
        }
        
        fields.push_back(headerItem);
        fieldsIdx.push_back(columnItemValue);
        header=header.substr(0, header.length()-1);
        std::regex rA("\\s+");
	header = std::regex_replace(header, rA, "");
        if (header.length()>2){
         header=header+std::string(",");
        }
        header=header+headerItem+std::string("\n");
      }
    }
    unsigned long int ccCounter=0;
    unsigned long int objectsCounter=0;
    file.clear();
    file.seekg(0,std::ios::beg);
    while (std::getline(file, strLine) ){
     if (++ccCounter >= headerLines){
        objectsCounter++;
      }
    }
    return objectsCounter;
}



// ********************************
//
// ********************************
void readMetadataValues(std::vector<double> & itemVector, std::string itemType,const char * filename){
  if ( isStarFile(filename) ){
    readStar(itemVector, itemType, filename);
  }else{
   //std::cerr<< filename <<" is NOT a star file\n";
   readCvs(itemVector, itemType, filename);
  }
}
void readMetadataValues(std::vector<std::string> & itemVector, std::string itemType,const char * filename){
  if ( isStarFile(filename) ){
    readStar(itemVector, itemType, filename);
  }else{
   //std::cerr<< filename <<" is NOT a star file\n";
   readCvs(itemVector, itemType, filename);
  }
}

unsigned long int StarToCsv(const char * filenameIn, const char * filenameOut, bool changeHeaderToVem = false){
    std::string dataLabelStr(relionDataStartLabel(filenameIn));
    int dataLabelStrSize=dataLabelStr.size();

    std::ifstream file(filenameIn);
    std::ofstream fileOut;
    fileOut.open(filenameOut);
    fileOut.close();
    fileOut.open (filenameOut, std::ofstream::out | std::ofstream::app);    
    std::string strLine;
    unsigned long int cc = 0;
    bool got_data_images = false;
    bool got_loop = false;
    std::vector<std::string> fields;
    std::vector<int> fieldsIdx;
    int targetID = -1;
    unsigned long int headerLines = 0;
    std::string header="#\n";
    while (std::getline(file, strLine) && ++cc < __MAX_STAR_HEADER_SIZE__){
        std::regex r("\\s+");
	std::string str0 = std::regex_replace(strLine, r, ",");
        std::regex rNoSpaces("\\s+");
	std::string strNoSpaces = std::regex_replace(strLine, rNoSpaces, ",");    
       if (strNoSpaces.length()<1){
         //do nothing
       }else if (str0.substr(0,1).compare("#")==0){
         //it is a comment, do nothing
       }else if (!got_data_images){
          //fileOut<<strLine<<"\n";
          //std::cerr<<strstr(strNoSpaces.substr(0,5).c_str(),"data_") << "<<-----\n";
          if( strstr(strNoSpaces.substr(0,dataLabelStrSize).c_str(),dataLabelStr.c_str()) ) {
            got_data_images=true;
            //std::cerr<<"star file: data_images\n";
          }
       }else if (!got_loop){
          //fileOut<<strLine<<"\n";
          if( strstr(str0.substr(0,5).c_str(),"loop_") ){
            got_loop=true;
            //std::cerr<<"got_loop\n";
          }
       }else if( strstr(str0.substr(0,1).c_str(),"_") ){ //was _rln
        //fileOut<<strLine<<"\n";
        headerLines=cc+1;
        std::size_t found = str0.find(",");
        std::string headerItem=str0.substr(0,found);
        std::size_t start1=str0.find("#");
        //std::size_t start1 = headerItem.find(",");
        std::string columnItem=str0.substr(start1+1, str0.length()-1);
        std::size_t start2=columnItem.find(",");
        columnItem=str0.substr(start1+1, start2);
        int columnItemValue=atoi(columnItem.c_str())-1;
        if (changeHeaderToVem){
         if (headerItem.compare("_rlnAngleRot")==0){
           headerItem=std::string("Phi");
         }else if (headerItem.compare("_rlnAngleTilt")==0){
           headerItem=std::string("Theta");
         } else if (headerItem.compare("_rlnAnglePsi")==0){
           headerItem=std::string("Psi");
         }
        }
        
        fields.push_back(headerItem);
        fieldsIdx.push_back(columnItemValue);
        header=header.substr(0, header.length()-1);
        std::regex rA("\\s+");
	header = std::regex_replace(header, rA, "");
        if (header.length()>2){
         header=header+std::string(",");
        }
        header=header+headerItem+std::string("\n");
      }
    }
    fileOut<<header;
    unsigned long int ccCounter=0;
    file.clear();
    file.seekg(0,std::ios::beg);
    unsigned long int counter = 0;
    while (std::getline(file, strLine) ){
     if (++ccCounter >= headerLines){
        std::string str00 = strLine.erase(0, strLine.find_first_not_of(" \t\n\r\f\v"));
        std::regex r("\\s+");
	std::string str0 = std::regex_replace(str00, r, ",");
	str0=str0.substr(0, str0.length()-1);
        fileOut<<str0<<"\n";
      }
    }
    fileOut.close();
    return counter;
}



// ////////////////
//
// writeCvs
//
int writeCvs(std::vector<double> & itemVector, std::string itemType, const char * filename, bool append = true) {
    std::ifstream file(filename);
    std::string str;
    std::string outstr("");
    int counterLine = 0;
    int counterItems=0;

    file.clear();
    file.seekg (0, std::ios::beg);

    std::getline(file, str);
    str.erase(std::remove(str.begin(),str.end(),' '),str.end()); //remove spaces
    str.erase(std::remove(str.begin(),str.end(),'\t'),str.end()); //remove tabs

    //insert header
    if (str[0] != '#'){
          outstr.append("#");
          outstr.append(itemType);
          outstr.append("\n");
          file.clear();
          file.seekg (0, std::ios::beg);
    } else {
      	if(append){
      			outstr.append(str);
      			outstr.append(",");
        }else{
            outstr = "#";
        }
        outstr.append(itemType);
        outstr.append("\n");
    }


    for (int counterItems = 0; counterItems < itemVector.size(); ){
      std::string tmpOutStr ("");
      bool endOfFileReached = true;
      if (std::getline(file, str)){
          endOfFileReached = false;
      }
      bool commentLine = ( str[0]=='#' ) ? true : false;
      if (commentLine){
         	tmpOutStr=str;
        	tmpOutStr.append("\n");
      }else {
        if( !endOfFileReached && append && str.size() > 0 ){
         			tmpOutStr=str;
        			tmpOutStr.append(",");
        }

  			std::ostringstream ss;
  			ss << std::setprecision(11) <<std::fixed<< itemVector[counterItems];
        tmpOutStr.append(std::string (ss.str()));
  			outstr.append( tmpOutStr );
  			outstr.append("\n");
        counterItems++;
      }

    }


    std::ofstream out(filename);
    out << outstr;
    out.close();
    return counterItems;
}



int writeCvs(std::vector<std::string> & itemVector, std::string itemType, const char * filename, bool append = true) {
    std::ifstream file(filename);
    std::string str;
    std::string outstr("");
    int counterLine = 0;
    int counterItems=0;

    file.clear();
    file.seekg (0, std::ios::beg);

    std::getline(file, str);
    str.erase(std::remove(str.begin(),str.end(),' '),str.end()); //remove spaces
    str.erase(std::remove(str.begin(),str.end(),'\t'),str.end()); //remove tabs

    //insert header
    if (str[0] != '#'){
          outstr.append("#");
          outstr.append(itemType);
          outstr.append("\n");
          file.clear();
          file.seekg (0, std::ios::beg);
    } else {
      		if(append){
      			outstr.append(str);
      			outstr.append(",");
      		}else{
            outstr = "#";
          }
          outstr.append(itemType);
          outstr.append("\n");
    }


    for (int counterItems = 0; counterItems < itemVector.size(); ){
      std::string tmpOutStr ("");
      bool endOfFileReached = true;
      if (std::getline(file, str)){
          endOfFileReached = false;
      }
      bool commentLine = ( str[0]=='#' ) ? true : false;
      if (commentLine){
         	tmpOutStr=str;
        	tmpOutStr.append("\n");
      }else {
        if( !endOfFileReached && append && str.size() > 0 ){
         			tmpOutStr=str;
        			tmpOutStr.append(",");
        }

                        tmpOutStr.append(std::string (  itemVector[counterItems]  ));
  			outstr.append( tmpOutStr );
  			outstr.append("\n");
        counterItems++;
      }

    }


    std::ofstream out(filename);
    out << outstr;
    out.close();
    return counterItems;
}


// //////
//
// test
//
int test(){
  std::vector<double> lisResolutions;
  lisResolutions.push_back(0.1);
  lisResolutions.push_back(0.2);
  lisResolutions.push_back(0.3);
  lisResolutions.push_back(0.4);
  lisResolutions.push_back(0.5);
  lisResolutions.push_back(0.6);
  writeCvs(lisResolutions, "resolutions", "foo.txt", false);

  std::vector<double> listCSV;
  listCSV.push_back(110);
  listCSV.push_back(30);
  listCSV.push_back(3);
  listCSV.push_back(10);
  listCSV.push_back(12);
  listCSV.push_back(6);
  writeCvs(listCSV, "msa", "foo.txt", true);

  std::vector<double> resultList;
  readCvs(resultList, std::string("msa"), "foo.txt");
	for (int i=0;i<resultList.size();i++){
		std::cerr<<resultList[i]<<"\n";
  }
  std::cerr<<"\n";
  //writeCvs(listMSA, "eccoli", "foo.txt");
    return 0;
}



//ADD and REMOVE columns
int removeCvsColumn(std::string columnName, const char * filename){
 std::vector<std::string> columns;
 readCvsColumns(columns, filename);
 int indexToDelete=indexCvsColumn(columnName, filename);
 //std::cerr<<"indexToDelete="<<indexToDelete<<"   filename="<<filename<<"\n";
 if (indexToDelete<0){
   return -1;
 }

 //string with header 
 std::string header="#";
 for (int ii=0; ii<columns.size(); ii++){
  if ( indexToDelete!=ii ){
    header=header+columns[ii];
    if (ii<columns.size()-1 && !( indexToDelete==columns.size()-1 && ii==columns.size()-2) ){
      header=header+std::string(",");
    }
  }
 }


    std::ifstream file(filename);
    std::string str;
    std::string outstr(header);
    int counterLine = 0;
    int counterItems=0;

    file.clear();
    file.seekg (0, std::ios::beg);

    std::getline(file, str);
    str.erase(std::remove(str.begin(),str.end(),' '),str.end()); //remove spaces
    str.erase(std::remove(str.begin(),str.end(),'\t'),str.end()); //remove tabs


    bool endOfFileReached = false;
    unsigned int counter;
    for (; !endOfFileReached; ){
      if (std::getline(file, str)){
        bool commentLine = ( str[0]=='#' ) ? true : false;
        if( !commentLine && !endOfFileReached && str.size() > 0 ){
           std::string tmpOutStr ("");
           for (unsigned long int tmpStart=0, nexItemPos = 0, counter=0; nexItemPos<str.length(); counter++){
             std::string stringTmp = str.substr (tmpStart);
             nexItemPos=stringTmp.find(",");
             stringTmp=str.substr(tmpStart,nexItemPos);
             tmpStart+=nexItemPos+1;
             if (counter != indexToDelete && counter > 0 && !(indexToDelete==0&&counter==1)){
               tmpOutStr.append(",");
             }
             if(counter != indexToDelete){
              tmpOutStr.append(stringTmp);
             }
  	   }
  	   outstr.append("\n");
  	   outstr.append( tmpOutStr );
        }
      }else{
       endOfFileReached = true;
      }
    }

    std::ofstream out(filename);
    out << outstr;
    out.close();
    return counterItems;
}





// ***************************
// replaceAddCvsColumn
// ***************************
void replaceAddCvsColumn(std::string columnName, std::vector<std::string> & itemVector, const char * filename){
 //std::cerr<<"filename="<<filename<<"   columnName="<<columnName<< "\n";
 int index=indexCvsColumn(columnName, filename);
 if (index >= 0){
   //column needs to be deleted
   removeCvsColumn(columnName, filename);
 }
 writeCvs(itemVector, columnName, filename, true);
}


void replaceAddCvsColumn(std::string columnName, std::vector<double> & itemVector, const char * filename){
 //std::cerr<<"filename="<<filename<<"   columnName="<<columnName<< "\n";
 int index=indexCvsColumn(columnName, filename);
 if (index >= 0){
   //column needs to be deleted
   removeCvsColumn(columnName, filename);
 }
 writeCvs(itemVector, columnName, filename, true);
}



/*
data_
loop_
_rlnImageName
_rlnMicrographName
_rlnAngleRot
_rlnAngleTilt
_rlnAnglePsi
_rlnOriginX
_rlnOriginY
_rlnCoordinateX
_rlnCoordinateY
_rlnDetectorPixelSize
*/
void cvsToStarFile(const char * filenameSrc, const char * filenameDst, const char * projectionfileName="projectionfileName"){
  std::vector<std::string> existingColumns;
  readCvsColumns(existingColumns, filenameSrc);

  int activeIdx=stringPositionVector(existingColumns, "active");
  int projectionfileNameIdx=stringPositionVector(existingColumns, projectionfileName);
  int originXIdx=stringPositionVector(existingColumns, "_rlnOriginX");
  int originYIdx=stringPositionVector(existingColumns, "_rlnOriginY");
  int coordinateXIdx=stringPositionVector(existingColumns, "_rlnCoordinateX");
  int coordinateYIdx=stringPositionVector(existingColumns, "_rlnCoordinateY");
  int DetectorPixelSizeIdx=stringPositionVector(existingColumns, "_rlnDetectorPixelSize");

  int rotIdx=stringPositionVector(existingColumns, "Phi");
  int tiltIdx=stringPositionVector(existingColumns, "Theta");
  int psiIdx=stringPositionVector(existingColumns, "Psi");

  
  if (projectionfileNameIdx<0){
   std::cerr<<"WARNING: bad file ("<< filenameSrc<<")\n";
   return;
  }
  std::vector<std::string> filenames;
  readCvs(filenames, projectionfileName, filenameSrc);
  if (filenames.size()<1){
     return;
  }

  std::vector<double> rotVector;
  std::vector<double> tltVector;  
  std::vector<double> psiVector;  
  if (rotIdx>=0){
   readCvs(rotVector, "Phi", filenameSrc);
  }
  if (tiltIdx>=0){
   readCvs(tltVector, "Theta", filenameSrc);
  }  
  if (psiIdx>=0){
   readCvs(psiVector, "Psi", filenameSrc);
  }  
  
  MRCHeader stackHeader;
  readHeaderMrc(filenames[0].c_str(), stackHeader);
  unsigned long int nx=stackHeader.nx;
  unsigned long int ny=stackHeader.ny;
  unsigned long int nz=1;
  unsigned long int nxyz=nx*ny*nz;
  
  //check the header of cvs file
  std::string header ("data_\nloop_\n_rlnImageName\n_rlnMicrographName\n_rlnOriginX\n_rlnOriginY\n_rlnCoordinateX\n_rlnCoordinateY\n_rlnDetectorPixelSize\n_rlnAngleRot\n_rlnAngleTilt\n_rlnAnglePsi\n");
  std::ofstream outfile;
  outfile.open(filenameDst);
  outfile << header;
  outfile.close();
  std::string space(" ");
  for (unsigned long int ii=0;ii<filenames.size(); ii++){
   std::string lineStr=std::string("0001@")+filenames[ii]+space;
   lineStr=lineStr+filenames[ii]+space;
   lineStr=lineStr+std::to_string(nx)+space+std::to_string(ny)+space+std::to_string(nx/2)+space+std::to_string(ny/2)+space;
   //lineStr=lineStr+std::to_string(0)+space+std::to_string(0)+space+std::to_string(nx/2)+space+std::to_string(ny/2)+space;
   lineStr=lineStr+std::to_string(1)+space; //pixel size (1 for the moment)
   if (ii<rotVector.size()){
    lineStr=lineStr+std::to_string(rotVector[ii])+space;
   }else{
    lineStr=lineStr+std::string("0")+space;
   }
   if (ii<tltVector.size()){
    lineStr=lineStr+std::to_string(tltVector[ii])+space;
   }else{
    lineStr=lineStr+std::string("0")+space;
   }   
   if (ii<psiVector.size()){
    lineStr=lineStr+std::to_string(psiVector[ii])+space;
   }else{
    lineStr=lineStr+std::string("0")+space;
   }
   lineStr=lineStr+std::string("\n");
   outfile.open(filenameDst, std::ios_base::app);
   outfile << lineStr;
   outfile.close();
  }
}




// *********************************
//
int replaceLabel(char * inputFile, char * outputFile, char * SourceFileForReplacingLabel, char * SourceLabelnameForReplacing, char * DestinationLabelnameForReplacing){
    
    std::stringstream sourceLabelsString(SourceLabelnameForReplacing);
    std::stringstream dstLabelsString(DestinationLabelnameForReplacing);
    std::vector<std::string> srcLabels;
    std::vector<std::string> dstLabels;
    while( sourceLabelsString.good() )
    {
        std::string substr;
        getline( sourceLabelsString, substr, ',' );
        srcLabels.push_back( substr );
    }
    while( dstLabelsString.good() )
    {
        std::string substr;
        getline( dstLabelsString, substr, ',' );
        dstLabels.push_back( substr );
    }
    if ( srcLabels.size() != dstLabels.size() ){
        std::cerr<<"ERROR: icoherent label sizes... exit\n";
        exit(1);
    }
    if ( srcLabels.size() < 1 ){
        std::cerr<<"ERROR: no label to replace... exit\n";
        return 2;
    }
    std::vector< std::vector<std::string> > sourceData;
    
    
    std::vector<int> labelDstIdxVector;
    for (int jjj=0;jjj<srcLabels.size();jjj++){
        std::vector<std::string> dataColumn;
        readMetadataValues(dataColumn, srcLabels[jjj],SourceFileForReplacingLabel);
        sourceData.push_back(dataColumn);
        int labelDstIdx=getStarHeaderItemIdx(dstLabels[jjj], inputFile);
        labelDstIdxVector.push_back(labelDstIdx);
    }

    
        std::ofstream fileOutput;
        fileOutput.open(outputFile);
        fileOutput.close();
        fileOutput.open (outputFile, std::ofstream::out | std::ofstream::app);
        long int startMicrograph=getStarStart(inputFile);
        std::ifstream fileParticles(inputFile);
        std::string strLine;
        for (int counter=0;counter<startMicrograph;counter++){
           std::getline(fileParticles, strLine);
           fileOutput << strLine <<"\n";
        }
       for (unsigned long int ii=0,counter=1;ii<(sourceData[0]).size();ii++){
           std::getline(fileParticles, strLine);
           std::string inputStr(strLine);
           std::string ouptputString (inputStr);
           for (int jjj=0;jjj<srcLabels.size();jjj++){
               ouptputString = replaceValueStrlineStarFile(inputStr, labelDstIdxVector[jjj], sourceData[jjj][ii]);
               inputStr=ouptputString;
           }
           fileOutput << ouptputString << "\n";
           counter++;
        }
       fileOutput.close();
    return 0;
}











// *********************************
//
 int splitStarFiles(char * inputFile, char * ListOutFilesForSplitting){
        int numLinesInput=getNumItemsStar(inputFile);

        std::stringstream sourceSplittingFiles(ListOutFilesForSplitting);
        std::vector<std::string> outFiles;
        std::cerr<<"from: " << inputFile << " output splitted files: ";
        while( sourceSplittingFiles.good() )
        {
            std::string substr;
            getline( sourceSplittingFiles, substr, ',' );
            outFiles.push_back( substr );
            std::cerr<<substr<<"  ";
        }
        std::cerr<<"\n";
        if ( outFiles.size() < 1 ){
            std::cerr<<"WARNING: no files to split...\n";
            return 2;
        } else {
          //  std::cerr<<"processing..n";
          std::vector<std::ofstream> fileOutputVector ( outFiles.size() );
          for (int fff=0;fff<outFiles.size(); fff++){
              fileOutputVector[fff].open( outFiles[fff] );
              fileOutputVector[fff].close();
              fileOutputVector[fff].open (outFiles[fff] , std::ofstream::out | std::ofstream::app);
          }
          long int startMicrograph=getStarStart(inputFile);
          std::ifstream fileInputStream(inputFile);
          std::string strLine;
          for (int counter=0;counter<startMicrograph;counter++){
            std::getline(fileInputStream, strLine);
            for (int fff=0;fff<outFiles.size(); fff++){
                fileOutputVector[fff] << strLine <<"\n";
            }          
          }

          const long int maxBlockSize= numLinesInput/outFiles.size();
          //std::cerr<<"maxBlockSize="<<maxBlockSize<<"\n";
          long int currentBlock = 0;
          long int counter = 0;
          std::string strLineTmp;
          while (std::getline(fileInputStream, strLineTmp) ){

          	std::string str(strLineTmp);
            str=str.erase(0, str.find_first_not_of(" \t\n\r\f\v"));
            if (str.size()>0) { //check the line is not empty
                if (counter>=maxBlockSize){
                    currentBlock++;
                    counter = 0;
                }
                //std::cerr<<"currentBlock="<<currentBlock<<"   "<<"counter="<<counter<<"   ";
                //std::cerr<<"\n";
                if ( currentBlock<fileOutputVector.size() ){
                   fileOutputVector[currentBlock] << strLineTmp <<"\n";
                }else{
                  fileOutputVector[0] << strLineTmp <<"\n";
                }
                counter++;
            }
          }

          for (int fff=0;fff<outFiles.size(); fff++){
              fileOutputVector[fff].close();
          }          

        }

    return numLinesInput;
 }




 int mergeStarFiles(char * ListFilesForMerging, char * outputFile ){
 
        std::stringstream sourceLabelFiles(ListFilesForMerging);
        std::vector<std::string> srcFiles;
        std::cerr<<"files to merge: ";
        while( sourceLabelFiles.good() )
        {
            std::string substr;
            getline( sourceLabelFiles, substr, ',' );
            srcFiles.push_back( substr );
            std::cerr<<substr<<"  ";
        }
        std::cerr<<"\n";
        if ( srcFiles.size() < 1 ){
            std::cerr<<"WARNING: no files with label to replace...\n";
            return 2;
        } else {
          //  std::cerr<<"processing..n";
          const char * fileWithHeader=(srcFiles[0]).c_str();
          std::ofstream fileOutput;
          fileOutput.open(outputFile);
          fileOutput.close();
          fileOutput.open (outputFile, std::ofstream::out | std::ofstream::app);
          long int startMicrograph=getStarStart(fileWithHeader);
          std::ifstream fileTemplate(fileWithHeader);
          std::string strLine;
          for (int counter=0;counter<startMicrograph;counter++){
            std::getline(fileTemplate, strLine);
            fileOutput << strLine <<"\n";
          }
          for (int fff=0;fff<srcFiles.size(); fff++){
            //std::cerr<<" "<<srcFiles[fff]<<"  \n ";
            std::ifstream fileTmpToMerge(srcFiles[fff]);
            std::string strLineFlush;
            long int startMicrographTmp=getStarStart(srcFiles[fff].c_str());
            for (int counter=0;counter<startMicrographTmp;counter++){
              std::getline(fileTmpToMerge, strLineFlush);
            }
            std::string strLineTmp;
            while (std::getline(fileTmpToMerge, strLineTmp) ){
          		  std::string str(strLineTmp);
                str=str.erase(0, str.find_first_not_of(" \t\n\r\f\v"));
                if (str.size()>0) { //check the line is not empty
                  //std::cerr<<"   ***"<<strLineTmp<<"***\n";
                  fileOutput << strLineTmp <<"\n";
                }
            }
          }
        }
     return 0;
 }




// *********************************
//
int averageLabels(char * inputFile, char * outputFile, char * SourceFilesForAveragingLabels, char * SourceLabelnameForAveraging){
    
    
    std::stringstream sourceLabelFiles(SourceFilesForAveragingLabels);
    std::vector<std::string> srcFiles;
    std::cerr<<"files to average labels from: ";
    while( sourceLabelFiles.good() )
    {
        std::string substr;
        getline( sourceLabelFiles, substr, ',' );
        srcFiles.push_back( substr );
        std::cerr<<substr<<"  ";
    }
    if ( srcFiles.size() < 1 ){
        std::cerr<<"WARNING: no files with label to replace...\n";
        return 2;
    }
    std::cerr<<"\n";

    std::stringstream sourceLabels(SourceLabelnameForAveraging);
    std::vector<std::string> srcLabels;
    std::cerr<<"labels to average: ";
    while( sourceLabels.good() )
    {
        std::string substr;
        getline( sourceLabels, substr, ',' );
        srcLabels.push_back( substr );
        std::cerr<<substr<<"  ";
    }
    if ( srcLabels.size() < 1 ){
        std::cerr<<"WARNING: no label to replace...\n";
        return 2;
    }
    std::cerr<<"\n";
    
    //OK
    //std::vector< double > sourceData;
    //readMetadataValues(sourceData, SourceLabelnameForAveraging, srcFiles[0].c_str() );
    std::vector< double > * sourceDataVector=new std::vector< double > [srcLabels.size()];
    //std::vector< double > * sourceDataVectorTmp=new std::vector< double > [srcLabels.size()];
    
    for (unsigned int ll=0;ll<srcLabels.size();ll++){
        readMetadataValues(sourceDataVector[ll], srcLabels[ll].c_str(), srcFiles[0].c_str() );
        //std::cerr<<sourceDataVector[ll][0]<<"  ";
    }
    //std::cerr<<"\n";
    
    for (unsigned int ii=1;ii<srcFiles.size();ii++){
        for (unsigned int ll=0;ll<srcLabels.size();ll++){
            std::vector< double > sourceDataTmp;
            readMetadataValues(sourceDataTmp, srcLabels[ll].c_str(), srcFiles[ii].c_str() );
            //std::cerr<<sourceDataTmp[0]<<"  ";
            if (sourceDataTmp.size() != sourceDataVector[0].size()){
                std::cerr<<"ERROR: incompatible size of data in "<<srcFiles[ii]<<"  and " << srcFiles[0]<<"  for label:"<< srcLabels[ll] <<"  ..EXIT";
                exit(1);
            }
            for (unsigned int kk=0;kk<sourceDataTmp.size();kk++){
                (sourceDataVector[ll])[kk]+=sourceDataTmp[kk];
            }
        }
    }
    for (unsigned int ll=0;ll<srcLabels.size();ll++){
        for (unsigned int kk=0;kk<sourceDataVector[ll].size();kk++){
            (sourceDataVector[ll])[kk]/=srcFiles.size();
        }
        //std::cerr<<sourceDataVector[ll][0]<<"  ";
    }
    //std::cerr<<"\n";
    
        //int labelDstIdx=getStarHeaderItemIdx(DestinationLabelnameForAveraging, inputFile);
        std::vector< int > labelDstIdxVector;
        for (unsigned int ll=0;ll<srcLabels.size();ll++){
            labelDstIdxVector.push_back( getStarHeaderItemIdx(srcLabels[ll].c_str(), inputFile) );
        }
    
        std::ofstream fileOutput;
        fileOutput.open(outputFile);
        fileOutput.close();
        fileOutput.open (outputFile, std::ofstream::out | std::ofstream::app);
        long int startMicrograph=getStarStart(inputFile);
        std::ifstream fileParticles(inputFile);
        std::string strLine;
        for (int counter=0;counter<startMicrograph;counter++){
           std::getline(fileParticles, strLine);
           fileOutput << strLine <<"\n";
        }

       for (unsigned long int ii=0,counter=1;ii<(sourceDataVector[0]).size();ii++){
           std::getline(fileParticles, strLine);
           std::string inputStr(strLine);
           std::string ouptputString (inputStr);
           for (int ll=0;ll<srcLabels.size();ll++){
               ouptputString = replaceValueStrlineStarFile(inputStr, labelDstIdxVector[ll], std::to_string((sourceDataVector[ll])[ii]));
               inputStr=ouptputString;
           }
           fileOutput << ouptputString << "\n";
           counter++;
        }

    fileOutput.close();
    
    delete [] sourceDataVector;
    //delete [] sourceDataVectorTmp;
    return 0;
}


#endif


