#include "pref.h"
#include "graph.h"
TH2F *big;

inline void init(){
	big = new TH2F();
}


Int_t npeaks = 10;



Double_t fpeaks(Double_t *x, Double_t *par) {
	Double_t result = par[0] + par[1]*x[0];
	for (Int_t p=0;p<npeaks;p++) {
		Double_t norm  = par[3*p+2];
		Double_t mean  = par[3*p+3];
		Double_t sigma = par[3*p+4];
		result += norm*TMath::Gaus(x[0],mean,sigma);
	}
	return result;
}


void graph(std::string name, bool big){
	Float_t x[5001],y[5001];
	std::ifstream in;
	TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,500);
	if(!big){
		in.open(name.c_str());
		
		for(Int_t i = 0 ; i< 5001 ;i++){
			in >> x[i] >> y[i];
		} 
		
		
		c1->cd(1);
		TGraph *gr = new TGraph(5001,x,y);
		gr->SetMinimum(-60.);
		gr->Draw("A*");
		TH1F *h = gr->GetHistogram();
		for(Int_t i = 0 ; i< 5001 ;i++){
			if(y[i]>= -60.)
				h->SetBinContent(i,y[i]);
		}
		
		h->SetYTitle("Intensity in dB");
		h->SetXTitle("#lambda in nm");
		h->SetTitle("Sectrum");
		h->Draw();
		
		TSpectrum *s = new TSpectrum(1000);
		Int_t nfound = s->Search(h,1,"new");
		std::cout <<"Found " << nfound << " candiate peaks to fit\n";
		c1->Update();
		
			//estimate linear background
		Double_t par[3000];
		TF1 *fline = new TF1("fline","pol1",842,852);
		h->Fit("fline","qn");
		
		par[0] = fline->GetParameter(0);
		par[1] = fline->GetParameter(1);
			//loop on all found peaks. Eliminate peaks at the background level
		Float_t *xpeaks = s->GetPositionX();
		for (Int_t p=0;p<nfound;p++) {
			Float_t xp = xpeaks[p];
			Int_t bin = h->GetXaxis()->FindBin(xp);
			Float_t yp = h->GetBinContent(bin);
			if (yp-TMath::Sqrt(yp) < fline->Eval(xp)) continue;
			par[3*npeaks+2] = yp;
			par[3*npeaks+3] = xp;
			par[3*npeaks+4] = 3;
			npeaks++;
		}
		c1->Update();
		
		TImage *img = TImage::Create();
		
		img->FromPad(c1);
		std::stringstream _name;
		_name << name << ".png";
		_name >> name;
		boost::filesystem::path path(name);
		std::string stringtmp = path.parent_path().string() +"/"+path.stem().string()+".png" ;
		std::cout <<"\n \n stem \n \n"<<stringtmp<< '\t' << path.stem().string() << std::endl;
		img->WriteImage(stringtmp.c_str());
		return;
	} 
	TGraph2D *g = new TGraph2D("tmp.dat");
	g->Draw();
	TImage *img = TImage::Create();
	img->FromPad(c1);
	img->WriteImage("big1.png");
	c1->Destructor();
}


void AtlasStyle() {
	
	std::cout << "\nApplying ATLAS style settings...\n" << std::endl ;
	
	TStyle *atlasStyle = new TStyle("ATLAS","Atlas style");
	
		// use plain black on white colors
	Int_t icol=0; // WHITE
	atlasStyle->SetFrameBorderMode(icol);
	atlasStyle->SetFrameFillColor(icol);
	atlasStyle->SetCanvasBorderMode(icol);
	atlasStyle->SetCanvasColor(icol);
	atlasStyle->SetPadBorderMode(icol);
	atlasStyle->SetPadColor(icol);
	atlasStyle->SetStatColor(icol);
		// don't use: white fill color floa *all* objects
	
		// set the paper & margin sizes
	atlasStyle->SetPaperSize(20,26);
	atlasStyle->SetPadTopMargin(0.05);
	atlasStyle->SetPadRightMargin(0.05);
	atlasStyle->SetPadBottomMargin(0.16);
	atlasStyle->SetPadLeftMargin(0.16);
	
		// use large fonts
		// Helvetica italics
	Int_t font=42; // Helvetica
	Double_t tsize=0.05;
	atlasStyle->SetTextFont(font);
	
	atlasStyle->SetTextSize(tsize);
	atlasStyle->SetLabelFont(font,"x");
	atlasStyle->SetTitleFont(font,"x");
	atlasStyle->SetLabelFont(font,"y");
	atlasStyle->SetTitleFont(font,"y");
	atlasStyle->SetLabelFont(font,"z");
	atlasStyle->SetTitleFont(font,"z");
	
	atlasStyle->SetLabelSize(tsize,"x");
	atlasStyle->SetTitleSize(tsize,"x");
	atlasStyle->SetLabelSize(tsize,"y");
	atlasStyle->SetTitleSize(tsize,"y");
	atlasStyle->SetLabelSize(tsize,"z");
	atlasStyle->SetTitleSize(tsize,"z");
	
		// use bold lines and markers
	atlasStyle->SetMarkerStyle(20);
	atlasStyle->SetMarkerSize(0.25);
	atlasStyle->SetHistLineWidth(2.);
	atlasStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
	
		// get rid of X error bars and y error bar caps
	
		// do not display any of the standard histogram decorations
	atlasStyle->SetTitleBorderSize(1);
	atlasStyle->SetOptTitle(1);
	atlasStyle->SetTitleFillColor(icol);
	atlasStyle->SetTitleX(0.5);
	atlasStyle->SetTitleY(1);
	atlasStyle->SetOptStat(0);
	atlasStyle->SetOptFit(0);
	
		// put tick marks on top and RHS of plots
	atlasStyle->SetPadTickX(1);
	atlasStyle->SetPadTickY(1);
	
	
	
	gROOT->SetStyle("ATLAS");
	gROOT->ForceStyle();
	
}

