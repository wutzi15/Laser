#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
//#include "TROOT"
#include "TMath.h"
#include "TGraph.h"
//#include <ifstream>
#include <fstream>
#include <iostream>
#include "TApplication.h"
#include "TImage.h"
#include <string>
#include <sstream>
#include "TStyle.h"
#include "TROOT.h"
#include "TGraph2D.h"
#include "boost/filesystem.hpp"

void graph(std::string name,bool big);
void init();
void AtlasStyle();
