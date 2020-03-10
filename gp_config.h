#ifndef GPCONFIG_H
#define GPCONFIG_H
#include <iostream>

//Yvonne made this on 3-10-2020
class CClgp : public gp_config{
  public:
    //This function is only used for TH1 root reading
  // GetCurrentArgument returns char*
  // Initialization, copied from CCLgp::learn()

    enum {
      KERNEL_USAGE_BACK,
      KERNEL_USAGE_FWD,
      KERNEL_USAGE_DYN
    };
    double tol=1e-6;
    string optimiser="scg";

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
}

