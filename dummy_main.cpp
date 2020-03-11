#include <iostream>
#include "gp_config.h"

using namespace std;



int main(int argc, char* argv[])
{

  gp_config * config = new gp_config();
  cout<<config->gettol()<<endl;
  cout<<"hello world"<<endl;

}
