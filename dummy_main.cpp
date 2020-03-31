#include <iostream>
#include "gp_config.h"
#include "gp_TH1.h"
#include "CMatrix.h"
#include "TFile.h"
#include "gp_TH1.h"

using namespace std;



int main(int argc, char* argv[])
{

  gp_config * config = new gp_config();
//TH1D* h1 = new TH1D("h1", "h1", 10, 0, 10);
  gp_TH1 * GP_TH1 = new gp_TH1(config);
  config->modelFileName="test.model";

  // -- first test config works, hello world
  //cout<<config->gettol()<<endl;
  //cout<<h1<<endl;
  //cout<<"hello world"<<endl;
  //
  // -- second test:  open root file and check the X and y matrix

  TFile * f1 = TFile::Open("/afs/cern.ch/work/y/ywng/workspace/dimuon/gpc/python/landau_histogram.root");

  f1->ls();

  TH1D * hist = (TH1D*)f1->Get("hist");

  cout<<"hisitogram: "<<hist<<endl;
  GP_TH1->readFromTH1(hist);
  GP_TH1->printXY();
  GP_TH1->learn("testlalala");
  GP_TH1->fitAndOutput();

}
