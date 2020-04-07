#include "gp_TH1.h"
#include "gp_config.h"
#include <stdlib.h>


using namespace std;

gp_TH1::gp_TH1(gp_config* config)
{
    // initializing the config
     Config = config;

}
  
void gp_TH1::readFromTH1(TH1D *&hist, double minX, double maxX){

  int maxFeat=1;  //only need 1d gaussian process right now
  int numData = hist->GetNbinsX();

  X.resize(numData, maxFeat);
  X.zeros();
  y.resize(numData, 1);
  y.zeros();

  int minBin;
  int maxBin;

  if (minX==0 and maxX==0){
    minBin=0; 
    maxBin=hist->GetNbinsX();
  }
  else{
    minBin=hist->FindBin(minX);
    maxBin=hist->FindBin(maxX);
  }
  input_minX=minBin;
  input_maxX=maxBin;

  for (int i =minBin; i<maxBin; i++){
    X.setVal(hist->GetBinCenter(i), i, 0);
	y.setVal(hist->GetBinContent(i),i);

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

void gp_TH1::learn(string comment){

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
  //for (auto kernel : kernelTypes){
  //  cout<<"kernel"<<kernel;

  //}
  cout<<"got here "<<endl;

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
  activeSetSize = Config->activeSetSize;
  approxType = Config->approxType;
  string approxTypeStr = Config->approxTypeStr;
  bool labelsProvided = Config->labelsProvided;
  double signalVariance = Config->signalVariance;
  vector<int> labels=Config->labels;
  int iters=Config->iters;
  string modelFileName=Config->modelFileName;


  cout<<"got here 0.5"<<endl;
  int inputDim = X.getCols();    


  cout<<"got here 0.7"<<endl;
  //CCmpndKern kern(X);
  //cannot call constructor again
  //kern(X);

  //update X
  kern.updateX(X);

  cout<<"got here 1"<<endl;
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
        kernels[i]->setParam(rbfInvWidths[i], 1); /// set rbf inverse width as specified.
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

    cout<<"got here 2"<<endl;
    }
    else {
      exitError("Unknown covariance function type: " + kernelTypes[i]);
    }

    //for (auto k :kernels){
    //  cout<<"kernels from kernels vector"<<endl;
    //  cout<<k<<endl;
    //}

    cout<<"got here 3"<<endl;
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

    cout<<"got here 4"<<endl;
    if(activeSetSize==-1)
      exitError("You must choose an active set size (option -a) for the command learn.");
  }
  
  noise.setTarget(&y);
  noise.initVals();
  noise.initParams(); 

  noise.setBias(0.0);
  CMatrix scale(1, y.getCols(), 1.0);
  CMatrix bias(1, y.getCols(), 0.0);
  if(centreData)
    bias.deepCopy(meanCol(y));
  if(scaleData)
    scale.deepCopy(stdCol(y));
  
    cout<<"got here 5"<<endl;
  //replacing the function pmodel by the class pmodel
  //CGp* pmodel;
  CMatrix bK(1,1,0.0);
  pmodel = new CGp(&kern, &noise, &X, approxType, activeSetSize, getVerbosity());

  //cout<<"at initialization: "<<typeid(pmodel->getKernel()).name()<<endl;
  
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

    cout<<"got here 6"<<endl;
    pmodel->setBetaVal(1); //
    pmodel->setScale(scale);
    pmodel->setBias(bias);
    pmodel->updateM();
    pmodel->setOutputScaleLearnt(outputScaleLearnt);
    //writeGpToFile(*pmodel, "c:\\gp_model", "Write for testing of model");  
    //
    cout<<"got here 7"<<endl;
    pmodel->setVerbosity(3);

    //pmodel->optimise(iters);
    //
    auto nudgeLearningScale = [&](){

      //Can't multiply matrix directly
        //pmodel->setBetaVal(1*rand()); //
        //pmodel->setScale(scale*rand());
        //pmodel->setBias(bias*rand());
        //pmodel->updateM();
        //pmodel->setOutputScaleLearnt(outputScaleLearnt*rand());
        ;

    };
    auto nudgeInitParams = [&] () {

        for(int i=0; i<kernelTypes.size(); i++) {
          //TODO for all kerneltypes other than rbf too
          if(kernelTypes[i]=="rbf") {
            if(rbfInvWidths[i]!=-1.0)
              kernels[i]->setParam(rbfInvWidths[i], 1*rand()); 
            if(variances[i]!=-1.0)
              kernels[i]->setParam(variances[i], 1*rand()); 
          } 
        }
    };
    
    int n=0; 

    cout<<"before optimise"<<endl;
    while (n< 100){
      try {
      pmodel->optimise(1);
      }
      catch (int e){
        cout<<"an exception occured in optimization"<<endl;
        cout<<e<<endl;

        nudgeInitParams();
        nudgeLearningScale();

      //This changes the gradiant of the kernel training

    }
    }
    // writing the file to be saved//
    // TODO uncomment this 
    //
    cout<<"after optimise"<<endl;
    writeGpToFile(*pmodel, modelFileName, comment);

    cout<<"got here 8"<<endl;

  //cout<<"at end of learn kernel type: "<<typeid(pmodel->getKernel()).name()<<endl;

  //cout<<"type of pkern: "<<typeid(*pkern).name()<<endl;

  //cout<<"type of pkern: "<<typeid(*pmodel->getKernel()).name()<<endl;
  //cout<<"pkern address: "<<pmodel->getKernel()<<endl;
  //
}
  


TH1D* gp_TH1::fitAndOutput(bool output_dat){

  double pointSize = 2;
  double lineWidth = 2;
  int resolution = 80;

  string name = "gp";

  string modelFileName="gp_model";
  string labelFileName="";


  pmodel->py=&y;
  pmodel->updateM();
  pmodel->pX=&X;

  //checking for invalid combination of input
  if(pmodel->getNoiseType()!="gaussian" && pmodel->getInputDim()!=2) {
    exitError("Incorrect number of model inputs.");
  }
  if(pmodel->getNoiseType()=="gaussian" && pmodel->getInputDim()>2) {
    exitError("Incorrect number of model inputs.");
  }

  if(X.getCols()!=pmodel->getInputDim()) {
    exitError("Incorrect dimension of input data.");
  }


  //if pmodel noise type is gaussian 
  if(pmodel->getNoiseType()=="gaussian") {
    switch(pmodel->getApproximationType())
    {
    case CGp::DTC:
    case CGp::DTCVAR:
    case CGp::FITC:
    case CGp::PITC:
      CMatrix scatterActive(pmodel->X_u.getRows(), pmodel->X_u.getCols()+1);
      CMatrix scatterOut(pmodel->X_u.getRows(), pmodel->getOutputDim());     
      pmodel->out(scatterOut, pmodel->X_u);

      for(int i=0; i<pmodel->X_u.getRows(); i++) 
      {
        for(int j=0; j<pmodel->X_u.getCols(); j++)

        scatterActive.setVal(pmodel->X_u.getVal(i, j), i, j);
        scatterActive.setVal(scatterOut.getVal(i, 0), i, pmodel->X_u.getCols());
      }
      
      scatterActive.toUnheadedFile(name+"_active_set.dat");
    }
    CMatrix scatterData(X.getRows(), X.getCols()+1);
    for(int i=0; i<X.getRows(); i++) 
    {
      for(int j=0; j<X.getCols(); j++)
      {
      scatterData.setVal(X.getVal(i, j), i, j);
      }
      scatterData.setVal(y.getVal(i, 0), i, X.getCols());
    }

    scatterData.toUnheadedFile(name+"_scatter_data.dat");

	 
    if(pmodel->X_u.getCols()==1) // one dimensional input.
    {
      double outLap=0.25;
      //Change this resolution to = to maxX-minX /totalx
      //int numx=resolution;
      //int numx=X.getRows();
      int numx=input_maxX-input_minX+1;


      CMatrix minVals(1, X.getCols());
      CMatrix maxVals(1, X.getCols());
      X.maxRow(maxVals);
      X.minRow(minVals);
      

      double xspan=maxVals.getVal(0, 0)-minVals.getVal(0, 0);

      maxVals.setVal(maxVals.getVal(0, 0)+outLap*xspan, 0, 0);
      minVals.setVal(minVals.getVal(0, 0)-outLap*xspan, 0, 0);
      xspan=maxVals.getVal(0, 0)-minVals.getVal(0, 0);

      double xdiff=xspan/(numx-1);
      CMatrix Xinvals(numx, 1);
      CMatrix regressOut(numx, 2);
      CMatrix errorBarPlus(numx, 2);
      CMatrix errorBarMinus(numx, 2);

      double x;
      int j;
      double xbins[numx];
      //cout<<"input_minX: "<<input_minX<<endl;
      //cout<<"input_maxX: "<<input_maxX<<endl;
      //cout<<"X: "<<X<<endl;
      int i;
      for (j=input_minX; j< input_maxX-1;j++){
          i=j-input_minX;
          
          xbins[i]=X.getVal(j, 0);

          Xinvals.setVal(X.getVal(j, 0), i, 0);
          regressOut.setVal(X.getVal(j, 0), i, 0);
          errorBarPlus.setVal(X.getVal(j, 0), i, 0);
          errorBarMinus.setVal(X.getVal(j, 0), i, 0);
      }

      TH1D* output_hist=new TH1D("gp_prediction", "gp_prediction", numx, xbins);
      
      CMatrix outVals(Xinvals.getRows(), pmodel->getOutputDim());
      CMatrix stdVals(Xinvals.getRows(), pmodel->getOutputDim());
//      cout<<"xinvals row: "<<Xinvals.getRows()<<endl;
//
//      cout<<"output dimension: "<< pmodel->getOutputDim()<<endl;
//      cout<<"input X dimension row: "<<X.getRows()<<endl;
//      cout<<"input X dimension column: "<<X.getCols()<<endl;
//      cout<<"input y dimension rows: "<<y.getRows()<<endl;
//      cout<<"input y dimension column: "<<y.getCols()<<endl;
//


      pmodel->out(outVals, stdVals, Xinvals);

      for(int j=0; j<numx; j++) {
        double val = outVals.getVal(j);
        regressOut.setVal(val, j, 1);
        errorBarPlus.setVal(val + 2*stdVals.getVal(j), j, 1);
        errorBarMinus.setVal(val - 2*stdVals.getVal(j), j, 1);
        output_hist->SetBinContent(j, val);
        output_hist->SetBinError(j, stdVals.getVal(j));
      }

      if (output_dat){
      string lineFile = name + "_line_data.dat";
      regressOut.toUnheadedFile(lineFile);
      string errorFile = name + "_error_bar_data.dat";
      ofstream out(errorFile.c_str());

      if(!out) throw ndlexceptions::FileWriteError(errorFile);
      out << setiosflags(ios::scientific);
      out << setprecision(17);
      out << "# Prepared plot of model file " << endl;
      errorBarPlus.toUnheadedStream(out);
      out << endl;
      errorBarMinus.toUnheadedStream(out);

      string plotFileName = name + "_plot.gp";
      ofstream outGnuplot(plotFileName.c_str());
      if(!outGnuplot) throw ndlexceptions::FileWriteError(plotFileName);
      outGnuplot << "plot \"" << name << "_line_data.dat\" with lines lw " << lineWidth;
      outGnuplot << ", \"" << name << "_scatter_data.dat\" with points ps " << pointSize;

      if(pmodel->isSparseApproximation())
      {
        outGnuplot << ", \"" << name << "_active_set.dat\" with points ps " << pointSize;
      }

      outGnuplot << ", \"" << name << "_error_bar_data.dat\" with lines lw " << lineWidth << endl;
      outGnuplot << "pause -1";
      outGnuplot.close();

    }// if output_dat

    return output_hist;
    }// if 1d data
    else{

    exitError("not 1d data .");
    TH1D * h1;
    return h1;
    }
    
  }//if noise type==gaussian
  else 
  {
    exitError("Unknown noise model for gnuplot output.");
    TH1D * h1;
    return h1;
  }
}

void gp_TH1::printXY(){
  cout<<"X: "<<X<<endl;
  cout<<"Y: "<<y<<endl;
}
