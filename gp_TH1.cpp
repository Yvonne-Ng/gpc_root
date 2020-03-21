#include "gp_TH1.h"
#include "gp_config.h"

using namespace std;

gp_TH1::gp_TH1(gp_config* config)
{
    // initializing the config
     Config = config;
}
  
void gp_TH1::readFromTH1(TH1D *&hist){

  int maxFeat=1;  //only need 1d gaussian process right now
  int numData = hist->GetNbinsX();
  cout<<"got here"<<endl;

  X.resize(numData, maxFeat);
  X.zeros();
  y.resize(numData, 1);
  y.zeros();
  cout<<"number of bins:"<<hist->GetNbinsX()<<endl;

  for (int i =0; i< hist->GetNbinsX(); i++){
    cout<<"got here"<<endl;
    X.setVal(hist->GetBinCenter(i), i, 0);
	y.setVal(hist->GetBinContent(i),i);
    cout<<"pt :"<<i <<endl;
    cout<<"X : "<< hist->GetBinCenter(i)<<endl;
    cout<<"Y : "<< hist->GetBinContent(i)<<endl;

  }


}



void gp_TH1::exitError(const string error){

  cerr << error << endl << endl;
  exit(1);
}

void gp_TH1::writeGpToFile(const CGp& model, const string modelFileName, const string comment)
{
  model.toFile(modelFileName, comment);
}

void gp_TH1::learn(){

  //Setting all variables to equal to the ones from gp_config
  //
  enum {
    KERNEL_USAGE_BACK,
    KERNEL_USAGE_FWD,
    KERNEL_USAGE_DYN
  };

  double tol=Config->tol;
  string optimiser=Config->optimiser;

  
  vector<string> kernelTypes=Config->kernelTypes;
  vector<unsigned int> kernelUsageFlag= Config->kernelUsageFlag;
  vector<double> ratQuadAlphas= Config->ratQuadAlphas;
  vector<double> rbfInvWidths= Config->rbfInvWidths;
  vector<double> weightVariances=Config->weightVariances;
  vector<double> biasVariances=Config->biasVariances;
  vector<double> variances=Config->variances;
  vector<double> degrees=Config->degrees;
  vector<bool> selectInputs=Config->selectInputs;
  bool centreData=Config->centreData;
  bool scaleData=Config->scaleData;
  bool outputScaleLearnt=Config->outputScaleLearnt;
  int activeSetSize = Config->activeSetSize;
  int approxType = Config->approxType;
  string approxTypeStr = Config->approxTypeStr;
  bool labelsProvided = Config->labelsProvided;
  double signalVariance = Config->signalVariance;
  vector<int> labels=Config->labels;
  int iters=Config->iters;
  string modelFileName=Config->modelFileName;

  int inputDim = X.getCols();    

  CCmpndKern kern(X);
  vector<CKern*> kernels;
  for(int i=0; i<kernelTypes.size(); i++) {
    CMatrix *M = 0;
    M = &X;
    if(kernelTypes[i]=="lin") {
      if(selectInputs[i])
        kernels.push_back(new CLinardKern(*M));
      else
        kernels.push_back(new CLinKern(*M));
      if(variances[i]!=-1.0)
        kernels[i]->setParam(variances[i], 0); // set variance parameter as specified.
    }
    else if(kernelTypes[i]=="poly") {
      if(selectInputs[i]) {
        kernels.push_back(new CPolyardKern(*M));
        if(degrees[i]!=-1.0)
          ((CPolyardKern*)kernels[i])->setDegree(degrees[i]);
      }
      else {
        kernels.push_back(new CPolyKern(*M));	
        if(degrees[i]!=-1.0)
          ((CPolyKern*)kernels[i])->setDegree(degrees[i]);
      }
      if(weightVariances[i]!=-1.0)
        kernels[i]->setParam(weightVariances[i], 0); // set `weight variance' parameter as specified.
      if(biasVariances[i]!=-1.0)
        kernels[i]->setParam(biasVariances[i], 1); // set `bias variance' parameter as specified.
      if(variances[i]!=-1.0)
        kernels[i]->setParam(variances[i], 2); // set variance parameter as specified.
    }
    else if(kernelTypes[i]=="rbf") {
      if(selectInputs[i])
        kernels.push_back(new CRbfardKern(*M));
      else
        kernels.push_back(new CRbfKern(*M));
      if(rbfInvWidths[i]!=-1.0)
        kernels[i]->setParam(rbfInvWidths[i], 0); /// set rbf inverse width as specified.
      if(variances[i]!=-1.0)
        kernels[i]->setParam(variances[i], 1); /// set variance parameter as specified.
	}
    else if(kernelTypes[i]=="exp") {
      if(selectInputs[i])
        exitError("Exponential covariance function not available with input selection yet.");
      else
        kernels.push_back(new CExpKern(*M));
      if(rbfInvWidths[i]!=-1.0)
        kernels[i]->setParam(rbfInvWidths[i], 0); /// set exp inverse width as specified.
      if(variances[i]!=-1.0)
        kernels[i]->setParam(variances[i], 1); /// set variance parameter as specified.
	}
    else if(kernelTypes[i]=="ratquad") {
      if(selectInputs[i])
        exitError("Rational quadratic covariance function not available with input selection yet.");
      else
        kernels.push_back(new CRatQuadKern(*M));
      if(ratQuadAlphas[i]!=-1.0)
        kernels[i]->setParam(ratQuadAlphas[i], 0); /// set rat quad length scale as specified.
	  if(rbfInvWidths[i]!=-1.0)
        kernels[i]->setParam(1/sqrt(rbfInvWidths[i]), 1); /// set rat quad length scale as specified.
      if(variances[i]!=-1.0)
        kernels[i]->setParam(variances[i], 2); /// set variance parameter as specified.
    }
    else if(kernelTypes[i] == "mlp") {
      if(selectInputs[i])
        kernels.push_back(new CMlpardKern(*M));
      else
        kernels.push_back(new CMlpKern(*M));
      if(weightVariances[i]!=-1.0)
        kernels[i]->setParam(weightVariances[i], 0); // set `weight variance' parameter as specified.
      if(biasVariances[i]!=-1.0)
        kernels[i]->setParam(biasVariances[i], 1); // set `bias variance' parameter as specified.
      if(variances[i]!=-1.0)
        kernels[i]->setParam(variances[i], 2); // set variance parameter as specified.
    }
    else if(kernelTypes[i] == "bias" && kernelUsageFlag[i]!=KERNEL_USAGE_FWD) {
      // fwd covariance function always has bias component
      kernels.push_back(new CBiasKern(*M));
      if(variances[i]!=-1.0)
        kernels[i]->setParam(variances[i], 0); // set variance parameter as specified.
    }
    else if(kernelTypes[i] == "white" && kernelUsageFlag[i]!=KERNEL_USAGE_FWD) {
      // fwd covariance function always includes a white noise component
      kernels.push_back(new CWhiteKern(*M));
      if(variances[i]!=-1.0)
        kernels[i]->setParam(variances[i], 0); // set variance parameter as specified.
    }
    else {
      exitError("Unknown covariance function type: " + kernelTypes[i]);
    }
    switch (kernelUsageFlag[i]) {
    case KERNEL_USAGE_FWD:
      kern.addKern(kernels[i]);
      break;
    case KERNEL_USAGE_BACK:
      break;
    case KERNEL_USAGE_DYN:
      break;
    }
  }
  // if no covariance function was specified, add an RBF.
  if(kern.getNumKerns()==0)
  {
    CKern* defaultKern = new CRbfKern(X);
    kern.addKern(defaultKern);
  }
  CKern* biasKern = new CBiasKern(X);
  CKern* whiteKern = new CWhiteKern(X);
  kern.addKern(biasKern);
  kern.addKern(whiteKern);

  if(approxTypeStr=="ftc") {
    approxType = CGp::FTC;
    activeSetSize = -1;
  }
  else {
    if(approxTypeStr=="dtc") {
      approxType = CGp::DTC;
    }
    else if(approxTypeStr=="dtcvar") {
      approxType = CGp::DTCVAR;
      cout << "Warning: numerical stabilities exist in DTCVAR approximation." << endl;
    }
    else if(approxTypeStr=="fitc") {
      approxType = CGp::FITC;
      exitError("FITC Approximation currently not working.");
    }
    else if(approxTypeStr=="pitc") {
      approxType = CGp::PITC;
    }
    else {
      exitError("Unknown sparse approximation type: " + approxTypeStr + ".");
    }
    if(activeSetSize==-1)
      exitError("You must choose an active set size (option -a) for the command learn.");
  }


  
  CGaussianNoise noise(&y);
  noise.setBias(0.0);
  // Remove scales and center if necessary.
  CMatrix scale(1, y.getCols(), 1.0);
  CMatrix bias(1, y.getCols(), 0.0);
  if(centreData)
    bias.deepCopy(meanCol(y));
  if(scaleData)
    scale.deepCopy(stdCol(y));
  

  CGp* pmodel;
  CMatrix bK(1,1,0.0);
  pmodel = new CGp(&kern, &noise, &X, approxType, activeSetSize, getVerbosity());
  if(optimiser=="scg")
    pmodel->setDefaultOptimiser(CGp::SCG);
  else if(optimiser=="conjgrad")
    pmodel->setDefaultOptimiser(CGp::CG);
  else if(optimiser=="graddesc")
    pmodel->setDefaultOptimiser(CGp::GD);
  else if(optimiser=="quasinew")
    pmodel->setDefaultOptimiser(CGp::BFGS);
  else
    exitError("Unrecognised optimiser type: " + optimiser);
  pmodel->setBetaVal(1); //
  pmodel->setScale(scale);
  pmodel->setBias(bias);
  pmodel->updateM();
  pmodel->setOutputScaleLearnt(outputScaleLearnt);
  //writeGpToFile(*pmodel, "c:\\gp_model", "Write for testing of model");  
  pmodel->optimise(iters);
  string comment="";
//comment += + ".";
  writeGpToFile(*pmodel, modelFileName, comment);
}

void gp_TH1::printXY(){
  cout<<"X: "<<X<<endl;
  cout<<"Y: "<<y<<endl;
}
  



