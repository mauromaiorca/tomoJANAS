/*
 * File: euler_libs.h
 * (C) 2022 Mauro Maiorca 
 */


#ifndef ___EULER_LIBS___H__
#define ___EULER_LIBS___H__

#include <iostream>
//#include <complex>
//#include <fstream>
//#include <list>
//#include <iterator>     // std::back_inserter
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <sstream>
#include <string>
//#include <regex>


double wrapAngles(double angle, double range=360.0){
    double positiveAngles=fmod(range+fmod(angle,range),range);
    return positiveAngles;
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


template<typename T>
std::vector<long int> bubbleSortDecendingIndexes(const std::vector<T> valuesIn, std::vector<long int> & indexes){
  unsigned long int size=valuesIn.size();
  T * values = new T [size];
  for (unsigned long int i=0; i<size; i++){
    values[i]=valuesIn[i];
  }
  for (long int i = (size - 1); i > 0; i--)
    {
      for (long int j = 1; j <= i; j++)
    {
      if (values[j - 1] < values[j])
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


// ******************************
//
//   score normalization
// ******************************
std::vector<double> scoreNormalizationHomogeneousViews ( std::vector<double> & phiListParticle,
                                           std::vector<double> & thetaListParticle,
                                           const std::vector<double> scores,
//                                           const std::vector<double> _randomSubset,
                                           const unsigned long int numViews, bool doNormalization=true) {


  const unsigned long int numActualViews=numViews;
  //get the number of test projections by
  //long int numActualViews=20;
  const double ga = (3 - pow(5.0,0.5)) * M_PI; // golden angle
  
  std::vector<double> phiV;
  std::vector<double> thetaV;
  double K=-1.0;
  double numEffectiveViews=0;
  for (unsigned long int ii=0; ii<numActualViews; ii++, K+=2.0/((double)numActualViews-1.0)){
    if ( K <= 1.0 ){
            double anglePhi = wrapAngles(ii*ga, 2*M_PI)*180.0/M_PI;
            double angleTheta = acos(K)*180.0/M_PI;
            thetaV.push_back( angleTheta );
            phiV.push_back( anglePhi );
            numEffectiveViews++;
    }
  }

/*
    removeCvsFile("angoli.csv");
    std::ifstream inCsvFileTmp("angoli.csv");
    inCsvFileTmp.close();
    replaceAddCvsColumn("Phi", phiV, "angoli.csv");
    replaceAddCvsColumn("Theta", thetaV, "angoli.csv");
*/
    
  //fill lookup Matrix
  std::vector<long int> * lookupMatrixIdx = new std::vector<long int> [ phiV.size() ];
  std::vector<long int> numParticlesPerView = std::vector<long int> ( phiV.size(), 0.0 );
  long int maxParticlesPerViews=0;
  long int minParticlesPerViews=0;
  //FOR EACH PARTICLE, FIND THE CLOSEST IN THE HOMOGENEOUS VIEWS
  for (unsigned long int ii=0;ii<phiListParticle.size();ii++){
    long int target = -1;
    //EXTREMELY INEFFICIENT PROCEDURE
    // BEGIN
    double minDistance = 999999;
    for (unsigned long int jj=0; jj<thetaV.size(); jj++){
        double tmpDistance=archDistanceSquaredApproximate(phiListParticle[ii], thetaListParticle[ii], phiV[jj], thetaV[jj]);
        if ( tmpDistance  < minDistance ){
            minDistance=tmpDistance;
            target=jj;
        }
    }
    //END EXTREMELY INEFFICIENT PROCEDURE
    if (target>=0 && target<(long int)phiV.size()){
     lookupMatrixIdx[target].push_back(ii);
     //numParticlesPerView[target]++;
    }
  }
  for (unsigned long int jj=0; jj<thetaV.size(); jj++){
      numParticlesPerView[jj]=lookupMatrixIdx[jj].size();
      if (numParticlesPerView[jj]>maxParticlesPerViews){
          maxParticlesPerViews = numParticlesPerView[jj];
      }
      if (numParticlesPerView[jj] < minParticlesPerViews){
          minParticlesPerViews = numParticlesPerView[jj];
      }
  }
  //std::cerr<<"minParticlesPerViews="<< minParticlesPerViews <<"\n";
  //std::cerr<<"maxParticlesPerViews="<< maxParticlesPerViews <<"\n";
  //const double maxRangeScore = 1000.0;
  const unsigned long int eulerLayers=maxParticlesPerViews;
  //std::cerr<<"got the reference target\n";
  std::vector<long int> * lookupLayersIdx = new std::vector<long int> [ eulerLayers ];

  //sort indexes for each section/view;
  //Zanetti's comment

  for (unsigned long int jj=0; jj<thetaV.size(); jj++ ){
      if(lookupMatrixIdx[jj].size()>0){
          std::vector<double> scoresSection;
          for (unsigned long int kk=0; kk<lookupMatrixIdx[jj].size();kk++){
              long int idx=lookupMatrixIdx[jj][kk];
              scoresSection.push_back(scores[idx]);
          }
          bubbleSortDecendingIndexes(scoresSection, lookupMatrixIdx[jj]);
          //OK std::cerr<<"\n********\nview="<<jj<<"  (elements="<<lookupMatrixIdx[jj].size()<<")\n";
          //for (unsigned long int kk=0; kk<lookupMatrixIdx[jj].size();kk++){
          //    long int idx=lookupMatrixIdx[jj][kk];
              //OK std::cerr<<lookupMatrixIdx[jj][kk]<<"("<<scores[idx]<<")   ";
          //}
          //OK std::cerr<<"\n";
        }
  }
    
//OK std::cerr<<"###########################\n";
//OK std::cerr<<"###########################\n";
    
  std::vector<double> scoresNormalized = scores;
  //pick values for each section/view, normalize and sort
  for ( unsigned long int kk=0; kk<eulerLayers; kk++ ){
      for (long int jj=0; jj<(long int)thetaV.size(); jj++ ){
          if ( lookupMatrixIdx[jj].size() > kk ){
              long int idx=lookupMatrixIdx[jj][kk];
              lookupLayersIdx[kk].push_back( idx );
              double maxScoreView = scores[ lookupMatrixIdx[jj][0] ];
              //double minScoreView = score(lookupMatrixIdx[jj][lookupMatrixIdx[jj].size()-1]);
              if ( maxScoreView > 0 && kk > 0 && doNormalization){
                  scoresNormalized[idx]/=maxScoreView;
              }else { //do not make sense normalize negative scores
                  scoresNormalized[idx]=scores[idx];
              }
        }
      }
  }


    for ( unsigned long int kk=0; kk<eulerLayers; kk++ ){
        std::vector<double> scoresLayer;
        for (long int jj=0; jj<(long int)lookupLayersIdx[kk].size(); jj++ ){
            long int idx=lookupLayersIdx[kk][jj];
            scoresLayer.push_back(scoresNormalized[idx]);
        }
        bubbleSortDecendingIndexes(scoresLayer, lookupLayersIdx[kk]);
    }
     
    std::vector<double> outRank (phiListParticle.size(), 0.0);

    for ( unsigned long int kk=0, rankNum=0; kk<eulerLayers; kk++ ){
        //OK std::cerr<<"\n********layer="<<kk<<"  (elements="<<lookupLayersIdx[kk].size()<<")\n";
        for (long int jj=0; jj<(long int)lookupLayersIdx[kk].size(); jj++, rankNum++ ){
            long int idx=lookupLayersIdx[kk][jj];
            //outRank[idx]=double(1.0)/double(rankNum+1.0);
            outRank[idx]= double(rankNum)/double(phiListParticle.size());
            //std::cerr <<idx<<" => " << rankNum+1 << " => " << outRank[idx] <<"\n";
            //outRank[idx]=rankNum+1.0;
        }
    }

  delete [] lookupMatrixIdx;
  delete [] lookupLayersIdx;


  return outRank;
}




// ******************************
//
//   half maps scoreNormalizationHomogeneousViews
// ******************************
std::vector<double> scoreNormalizationHomogeneousViews_halfMaps ( std::vector<double> & phiListParticle,
                                           std::vector<double> & thetaListParticle,
                                           const std::vector<double> scores,
                                           const std::vector<double> randomSubset,
                                           const unsigned long int numViews, bool doNormalization=true) {

    std::vector<long int> idx_h1;
    std::vector<long int> idx_h2;
    std::vector<double> phiListParticle_h1;
    std::vector<double> phiListParticle_h2;
    std::vector<double> thetaListParticle_h1;
    std::vector<double> thetaListParticle_h2;
    std::vector<double> scores_h1;
    std::vector<double> scores_h2;
    std::vector<double> outrank;
    
    
    for (unsigned long int ii=0;ii<randomSubset.size(); ii++){
        outrank.push_back(0.0);
        if(randomSubset[ii]==1){
            idx_h1.push_back(ii);
            phiListParticle_h1.push_back(phiListParticle[ii]);
            thetaListParticle_h1.push_back(thetaListParticle[ii]);
            scores_h1.push_back(scores[ii]);
        }else{
            idx_h2.push_back(ii);
            phiListParticle_h2.push_back(phiListParticle[ii]);
            thetaListParticle_h2.push_back(thetaListParticle[ii]);
            scores_h2.push_back(scores[ii]);
        }
    }

    std::vector<double> outrank_h1 = scoreNormalizationHomogeneousViews ( phiListParticle_h1,thetaListParticle_h1,scores_h1,numViews,doNormalization);
    std::vector<double> outrank_h2 = scoreNormalizationHomogeneousViews ( phiListParticle_h2,thetaListParticle_h2,scores_h2,numViews,doNormalization);
    for (unsigned long int ii=0;ii<idx_h1.size(); ii++){
        long int idx=idx_h1[ii];
        phiListParticle[idx]=phiListParticle_h1[ii];
        thetaListParticle[idx]=thetaListParticle_h1[ii];
        outrank[idx]=outrank_h1[ii];
    }
    for (unsigned long int ii=0;ii<idx_h2.size(); ii++){
        long int idx=idx_h2[ii];
        phiListParticle[idx]=phiListParticle_h2[ii];
        thetaListParticle[idx]=thetaListParticle_h2[ii];
        outrank[idx]=outrank_h2[ii];
    }
/*
    for (int ii=0;ii<10;ii++){
        std::cerr<<"phi="<<phiListParticle[ii]
            <<"  theta="<<thetaListParticle[ii]
            <<"  scores="<<scores[ii]
            <<"  randomSubset="<<randomSubset[ii]
            <<"  outrank="<<outrank[ii]
            <<"\n";
    }
    */
    return outrank;
 }


// ******************************
//
//   Particles Group
// ******************************
std::vector<std::string> getEulerClassGroup ( std::vector<double> & phiListParticle,
                                           std::vector<double> & thetaListParticle,
                                           const unsigned long int numViews) {


  const unsigned long int numActualViews=numViews;
  //get the number of test projections by
  //long int numActualViews=20;
  const double ga = (3 - pow(5.0,0.5)) * M_PI; // golden angle
  
  std::vector<double> phiV;
  std::vector<double> thetaV;
  double K=-1.0;
  double numEffectiveViews=0;
  for (unsigned long int ii=0; ii<numActualViews; ii++, K+=2.0/((double)numActualViews-1.0)){
    if ( K <= 1.0 ){
            double anglePhi = wrapAngles(ii*ga, 2*M_PI)*180.0/M_PI;
            double angleTheta = acos(K)*180.0/M_PI;
            thetaV.push_back( angleTheta );
            phiV.push_back( anglePhi );
            numEffectiveViews++;
    }
  }


  //fill lookup Matrix
  std::vector<std::string> targetIdx;
  std::vector<long int> numParticlesPerView = std::vector<long int> ( phiV.size(), 0.0 );
  //long int maxParticlesPerViews=0;
  //long int minParticlesPerViews=0;
  //FOR EACH PARTICLE, FIND THE CLOSEST IN THE HOMOGENEOUS VIEWS
  for (unsigned long int ii=0;ii<phiListParticle.size();ii++){
    long int target = -1;
    //EXTREMELY INEFFICIENT PROCEDURE
    // BEGIN
    double minDistance = 999999;
    for (unsigned long int jj=0; jj<thetaV.size(); jj++){
        double tmpDistance=archDistanceSquaredApproximate(phiListParticle[ii], thetaListParticle[ii], phiV[jj], thetaV[jj]);
        if ( tmpDistance  < minDistance ){
            minDistance=tmpDistance;
            target=jj;
        }
    }
    //END EXTREMELY INEFFICIENT PROCEDURE
    targetIdx.push_back(std::to_string(target));
  }

  return targetIdx;
}



#endif
