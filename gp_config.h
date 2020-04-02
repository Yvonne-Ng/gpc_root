#ifndef GPCONFIG_H
#define GPCONFIG_H
#include <iostream>
#include <string> 
#include <vector>
#include <iostream>
#include <stdlib.h>

#include "ndlexceptions.h"

using namespace std;
//Yvonne made this on 3-10-2020
class gp_config{
  public:
    //This function is only used for TH1 root reading
  // GetCurrentArgument returns char*
  // Initialization, copied from CCLgp::learn()
  // constructorr
    gp_config(){
      ;

    }

    enum {
      KERNEL_USAGE_BACK,
      KERNEL_USAGE_FWD,
      KERNEL_USAGE_DYN
    };
    double tol=1e-6;
    string optimiser="scg";

    double gettol(){
      return tol;
    }

    vector<string> kernelTypes;
    vector<unsigned int> kernelUsageFlag;
    vector<double> ratQuadAlphas;
    vector<double> rbfInvWidths;
    vector<double> weightVariances;
    vector<double> biasVariances;
    vector<double> variances;
    vector<double> degrees;
    vector<bool> selectInputs;
    bool centreData=true;
    bool scaleData=false;
    bool outputScaleLearnt=false;
    int activeSetSize = -1;
    int approxType = -1;
    string approxTypeStr = "ftc";
    bool labelsProvided = true;
    double signalVariance = 0.0;
    vector<int> labels;
    int iters=1000;
    string modelFileName="gp_model";
    // Change these variables on your own, when CClgp::learnTH1 reads in the loop, initializes its variables value to the ones specified in this configuration
    void setkernel(char * kernel){

        kernelTypes.push_back(kernel); 
        kernelUsageFlag.push_back(KERNEL_USAGE_FWD);
        ratQuadAlphas.push_back(-1.0);
		//rbfInvWidths.push_back(-1.0);
		rbfInvWidths.push_back(5);
        weightVariances.push_back(-1.0);
        biasVariances.push_back(-1.0);
        variances.push_back(-1.0);
        degrees.push_back(-1.0);
        selectInputs.push_back(false);

    }
    void setGamma (double gamma ){
        if(kernelTypes.size()==0)
          exitError("Inverse width specification must come after covariance function type is specified.");
        if(kernelTypes[kernelTypes.size()-1]!="rbf" && kernelTypes[kernelTypes.size()-1]!="exp" && kernelTypes[kernelTypes.size()-1]!="ratquad")
          exitError("Inverse width parameter only valid for RBF, exponential and rational quadratic covariance function.");
        rbfInvWidths[rbfInvWidths.size()-1]=2*gamma;

    }
      void setDegree(double degree){

        if(kernelTypes.size()==0)
          exitError("Polynomial degree specification must come after covariance function type is specified.");
        if(kernelTypes[kernelTypes.size()-1]!="poly")
          exitError("Polynomial degree parameter only valid for poly covariance function.");
        degrees[degrees.size()-1]=degree;
        }
      void setWeight(double weight){

        if(kernelTypes.size()==0)
          exitError("`Weight variance' parameter specification must come after covariance function type is specified.");
        if(kernelTypes[kernelTypes.size()-1]!="poly" 
           && kernelTypes[kernelTypes.size()-1]!="mlp")
          exitError("`Weight variance' parameter only valid for polynomial and MLP covariance function.");
        weightVariances[weightVariances.size()-1]=weight;
        }

      void setBias(double bias){

        if(kernelTypes.size()==0)
          exitError("`Bias variance' parameter specification must come after covariance function type is specified.");
        if(kernelTypes[kernelTypes.size()-1]!="poly" 
           && kernelTypes[kernelTypes.size()-1]!="mlp")
          exitError("`Bias variance' parameter only valid for polynomial and MLP covariance function.");
        biasVariances[biasVariances.size()-1]=bias;
      }

      void setVariance(double variance){

        if(kernelTypes.size()==0) exitError("Variance parameter specification must come after covariance function type is specified.");
        variances[variances.size()-1]=variance;
      }
      void setInput(bool input){

        if(kernelTypes.size()==0)
          exitError("Input selection flag must come after covariance function type is specified.");
        selectInputs[selectInputs.size()-1]=input;
      }
    private:

      void exitError(const string error)
      {
        cerr << error << endl << endl;
        exit(1);
      }
};

#else /* GPCONFIG_H */
#endif /* GPCONFIG_H */

