#include "pref.h"
#include "graph.h"
#include "defines.h"
#include "helper.h"





void gradient(double *asy, double *center,double *integral, int n,TCanvas *canv){
	bool warning = false;
	for(int i = 0; i < n -2; i++){
		
		double gradient_asy;
		double gradient_asy2;
		double gradient_centr ;
		double gradient_centr2;
		double gradient_int;
		double gradient_int2;
		gradient_asy = (asy[i+1]-asy[i]);
		if(i <= ((sizeof(asy)/sizeof(asy[0]))) -2){
			gradient_asy2 =  (asy[i+1]-asy[i+2]);
			gradient_int2 = (integral[i+1]-integral[i+2]);
		}else{
			gradient_asy2 = 0;
		}
		gradient_centr = (center[i+1]-center[i]);
		gradient_int = (integral[i+1]-integral[i]);
		gradient_centr2 = (center[i+1]-center[i+2]);
		if (gradient_asy == INFINITY || gradient_asy2== INFINITY || gradient_centr == INFINITY|| gradient_centr2 == INFINITY || gradient_int == INFINITY || gradient_int2 == INFINITY) {
			continue;
		}
			//Warnings
		if(fabs(gradient_asy) >= WARNING_PERC_ASY || fabs(gradient_centr)>= WARNING_PERC_CEN || fabs(gradient_int)>= WARNING_PERC_INT ){
			warning=true;
			std::cout << "\n\n WARNING!!! \n\n";
			std::cout << "Gradient detected @:"<<i << std::endl;
			std::cout << "Asym. Grad: " << gradient_asy << '\t' << "Center Grad: " << gradient_centr << '\t'<<"Inten. Grad: " << gradient_int<<std::endl;
			label(canv,i,WARNING);
			continue;
		}
		if(fabs(gradient_asy) >= NOTICE_PER_ASY || fabs(gradient_centr)>= NOTICE_PER_CEN || fabs(gradient_int)>= NOTICE_PERC_INT ){
			std::cout << "Gradient detected @:"<<i << std::endl;
			std::cout << "Asym. Grad: " << gradient_asy << '\t' << "Center Grad: " << gradient_centr << '\t'<<"Inten. Grad: " << gradient_int<<std::endl;
			label(canv,i,NOTICE);
			continue;
		}
	}
	if (warning) {
		TCanvas *warn = new TCanvas("Warning", "Warning", 200, 100);
		TText *warn_text = new TText(0.1, 0.4, "Warning!");
		warn_text->SetTextSize(0.5);
		warn_text->SetTextColor(kRed);
		warn->cd(1);
		warn_text->Draw();
		warn->Update();
	}
	
}

int main(int argc , char* argv[]){
	Double_t startwl, stopwl;
	startwl = boost::lexical_cast<double>(argv[1]);
	stopwl = boost::lexical_cast<double>(argv[2]);
	std::cout << startwl << '\t' << stopwl << std::endl;

	Double_t max = -210;
	Double_t maxwl = 0;
	int _argc = argc;
	TApplication *t = new TApplication("big",&_argc,argv);
	std::cout << "Running with boost" <<std::endl;
	std::vector<double> _x,_y;
	Double_t x[LINES], y[LINES], _inta[LINES], _intb[LINES]; 
	
	Double_t *cmp_int = new Double_t[argc];
	Double_t *argc_ary = new Double_t[argc];
	Double_t *cmp_int_root = new Double_t[argc];
	Double_t *asymmety_ary = new Double_t[argc];
	Double_t *width_ary = new Double_t [argc];
	
	std::ofstream of;
	std::ofstream integral_of;
	integral_of.open("integral.txt");
	TGraph2D *gr = new TGraph2D(LINES*(argc-1));
		//Setting up canvas for plot of all sectrums (is it called spectrums? ;) )
	TCanvas *c1 = new TCanvas("All Plots","All Plots",10,10,3000,1500);
	TH1F *integral_hist = new TH1F("Asymmerty", "Asymmetry", 100,0, 100);
	
	std::cout << argc%ROWS<< std::endl;
	if(!(argc % ROWS)){
		c1->Divide(argc/ROWS,ROWS);
		
	}else{
		c1->Divide(argc/ROWS+(argc %ROWS -1),ROWS);
	}
	
	of.open("tmp.dat");
	for (Int_t i = NUM_ARGS +1; i < argc ; i++){
		try{ 
			max = -210;
			maxwl = 0;
			argc_ary[i] = i-NUM_ARGS;
			
			std::ifstream in;
			in.seekg(0, std::ios::beg);
				// voodoo keep this;
			char **arg1 = t->Argv() ;
			std::string tmp = arg1[i];
			in.open(tmp.c_str());
			std::cout<< "file: " << tmp << std::endl;
			std::string line;
			int cline = 0;
			std::vector<double> a,b, inta, intb;
			
				//reading file
			while(getline(in,line)){
				std::string first,sec;
				
				std::string TR = line.c_str();
				std::replace(TR.begin(),TR.end(),',','\t');
				std::stringstream sstr(TR);
				sstr >> first >> sec;
				if(check(first) &&check(sec)){
					double _first = boost::lexical_cast<double>(first);
					double _sec = boost::lexical_cast<double>(sec);
					a.push_back(_first);
					b.push_back(_sec);
					cline++;
					inta.push_back(_first);
					_sec < -70.0 ? intb.push_back(0) : intb.push_back(_sec+70.0);
					if ( *(intb.end() - 1) < 0) {
						std::cout << "??? " << *(intb.end() - 1) << std::endl;
					}
				}
			}
			if (cline < LINES){
				for(int i = cline ; i < LINES ; i++){
					a.push_back(100);
					b.push_back(-210);
				}
			}
			std::cout<< "\n\n cline " << cline<< std::endl; 
			cline =(cline > LINES) ? LINES :cline;
			
			for(Int_t j = 0; j <LINES ;j++){
				x[j] = a[j];
				y[j] = b[j];
				_inta[j] = inta[j];
				_intb[j]= (intb[j] < 0)? 0:intb[j];
				of << x[j]<<'\t' << y[j]<<'\t' <<  '\n';
			}
			double s_integral = 0;
			std::cout <<"size of int " <<  intb.size()<< std::endl;
			for (size_t it = 0; it < intb.size() - 1; it++){
				double y_val = (intb[it]+intb[it+1])/2;
				assert (y_val >= 0);
				double area = 0.002*y_val;
				if(area > 0 )
					s_integral += area;
			}
			
			
			std::cout << "Simpson integral: " <<s_integral <<std::endl;
			integral_of << i << '\t' << s_integral << std::endl;
			integral_hist->Fill(s_integral);
			cmp_int[i] = s_integral;
			Int_t lines = (Int_t)intb.size();
			TGraph *r_integral = new TGraph(lines, _inta, _intb);
			
			std::cout << "ROOT integral: " << r_integral->Integral() << std::endl;
			cmp_int_root[i] = r_integral->Integral();
				
				//expanding
			expand(y, THRS_EXPAND, RATIO_EXPAND, LINES);
				
				
				//Filling TGraph2D
			
			for(Int_t j = 0; j <LINES ; j++){
				if (y[j] > max){
					max = y[j];
					maxwl = x[j];
				}
				gr->SetPoint(j+i*LINES, x[j],i,y[j]);
			}
			in.seekg(0, std::ios::beg);
			in.close();
			
				//Plotting each spectrum
			TGraph *_gr = new TGraph(LINES,x,y);
			_gr->GetHistogram()->GetXaxis()->SetTitle("#lambda in nm");
			_gr->GetHistogram()->GetYaxis()->SetTitle("Intensity in dB");
			c1->cd(i-NUM_ARGS);
			Double_t integral = 0;
			TH1F* hist = _gr->GetHistogram();
			for (Int_t k =0 ; k < hist->GetNbinsX();k++){
				integral += hist->GetBinContent(k);
			}
			_gr->Draw("AP");
			_gr->GetYaxis()->SetRangeUser(-80.,-10.);
			_gr->GetXaxis()->SetRangeUser(startwl,stopwl);
			_gr->SetTitle(tmp.c_str());
			c1->Update();
			TImage *img = TImage::Create();
			img->FromPad(c1);
			img->WriteImage("all.png");
			
			
				//Calculating asymmetry
			std::cout << "maximum: " << max << std::endl;
			double leftlimit, rightlimit = 1;
			leftlimit = findlower(x,y, max);
			rightlimit = findupper(x,y, max);
			if (leftlimit != 1 && rightlimit != 1){
				width_ary[i] = (leftlimit +rightlimit)/2;
			}else{
				width_ary[i] = maxwl;
			}
			double calced_asy = (maxwl-leftlimit)/(rightlimit-maxwl);
			asymmety_ary[i-3] = calced_asy;
			std::string asy_text = boost::lexical_cast<std::string>(calced_asy);
			std::cout << "Asymmetry: " << calced_asy << std::endl;
			
		}catch(std::exception e){
			std::cout << e.what()<< std::endl;
		}
	}
	of.close();
		//Setting style for 3D Plot
	TCanvas *d = new TCanvas("big","big",10,10,1500,800);
	d->Divide(2,2);
	d->cd(1);
	TGraph *the_ints = new TGraph(argc-1,argc_ary,cmp_int);
	the_ints->Draw("A*");
	the_ints->SetTitle("My Ints");
	d->Update();
	d->cd(2);
	std::cout << "Fitting\n\n";
	integral_hist->SetFillColor(kBlue);
	integral_hist->Draw();
	integral_hist->Fit("gaus");
	d->Update();
	d->cd(3);
	
	TGraph *roots_int = new TGraph(argc-1, argc_ary, cmp_int_root);
	roots_int->SetTitle("ROOTS Int");
	roots_int->Draw("A*");
	d->Update();
	d->cd(4);
	d->Update();
	gROOT->SetStyle("modern");
	gr->SetTitle("big");
	gr->GetHistogram("empty")->GetXaxis()->SetTitle("#lambda in nm");
	gr->GetHistogram("empty")->GetXaxis()->SetLimits(startwl,stopwl);
	gr->GetHistogram("empty")->GetYaxis()->SetTitle("Messurement"); 
	gr->GetHistogram("empty")->GetZaxis()->SetTitle("Intensity in dB");
	gr->GetHistogram("empty")->GetXaxis()->SetTitleOffset(1.5);
	gr->GetHistogram("empty")->GetYaxis()->SetTitleOffset(1.5);
	gr->GetHistogram("empty")->GetZaxis()->SetTitleOffset(1.5);
	gr->GetHistogram("empty")->GetZaxis()->SetRangeUser(-70.,max);
	gr->GetHistogram("empty")->GetXaxis()->CenterTitle();
	gr->GetHistogram("empty")->GetYaxis()->CenterTitle();
	gr->GetHistogram("empty")->GetZaxis()->CenterTitle();
	gr->Draw("PCOL");
	d->SetFillColor(16);
	
	
#ifdef RENDER
		//Render 3D animation
	const Int_t kUPDATE = 1;
	TSlider *slider = 0;
	
	for (Int_t i = 1; i <= 125; i++){
		TView3D *v = new TView3D();
		v->RotateView(5+i,45+i,d);
			//d->Update();
		
		if(i && (i%kUPDATE)== 0){
			if (i == kUPDATE){
				gr->Draw("PCOL");
				d->Update();
				slider = new TSlider("slider","test",850,-70,856,max);
			}
			if (slider) slider->SetRange(0,Float_t(i)/10000.);
			d->Modified();
			d->Update();
			d->Print("3d.gif+");
		} 
	}
	d->Update();
	d->Print("3d.gif++");
#endif
	

		//Saving image
	TImage *img = TImage::Create();
	boost::filesystem::path p(t->Argv(3));
	std::string file = p.parent_path().string();
	file += "big.png";
	img->FromPad(d);
	img->WriteImage(file.c_str());
		//cleaning
	boost::filesystem::remove(boost::filesystem::path("tmp.dat"));
	TCanvas *e = new TCanvas("Asymmetry","Asymmetry",10,10,1500,800);
	e->Divide(2,1);
	TGraph *asy_plot = new TGraph(argc-1, argc_ary, asymmety_ary);
	e->cd(1);
	asy_plot->SetTitle("Asymmetry");
	asy_plot->GetHistogram()->GetXaxis()->SetTitle("# Meassurement");
	asy_plot->GetHistogram()->GetYaxis()->SetTitle("Asymmetry");
	asy_plot->GetHistogram()->GetXaxis()->SetRange(1, argc);
	asy_plot->Draw("A*");
	e->Update();
	e->cd(2);
	
	
	TGraph *center_plot = new TGraph(argc-1 , argc_ary, width_ary);
	center_plot->GetHistogram()->GetXaxis()->SetTitle("# Meassurement");
	center_plot->GetHistogram()->GetYaxis()->SetTitle("Center in nm");
	center_plot->GetHistogram()->GetYaxis()->SetRangeUser(startwl, stopwl);
	center_plot->SetTitle("Center");
	center_plot->Draw("A*");
	e->Update();
		//Saving Images
	TImage *secimg = TImage::Create();
	boost::filesystem::path p2(t->Argv(3));
	file = p2.parent_path().string();
	file += "asy_cent.png";
	secimg->FromPad(e);
	secimg->WriteImage(file.c_str());
	
	TImage *thrdimg = TImage::Create();
	boost::filesystem::path p3(t->Argv(3));
	file = p3.parent_path().string();
	file += "allplots.png";
	thrdimg->FromPad(c1);
	thrdimg->WriteImage(file.c_str());
	
		//detecting Gradients
	gradient(asymmety_ary, width_ary,cmp_int, argc-1,c1);
	std::cout << "\n\n\nDone !!\nYou can quit now using CTRL+C \n" ;
	
	t->Run();
	
	
	delete[] cmp_int;
	delete[] argc_ary; 
	delete[] cmp_int_root;
	delete[] asymmety_ary;
	delete[] width_ary;
	return 0;
}
