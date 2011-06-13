#include "pref.h"
#include "graph.h"
#include "defines.h"
#include "helper.h"
int NUM_ARGS = 0;




void gradient(double *asy, double *center,double *integral, int n,TCanvas *canv){

	bool warning = false;
	for(int i = 0; i < n -2; i++){
			//foo
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

void read_file(std::string &line, std::vector<double> &a, std::vector<double> &b, std::vector<double> &inta, std::vector<double> &intb  ){
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
			//cline++;
		inta.push_back(_first);
		_sec < -70.0 ? intb.push_back(0) : intb.push_back(_sec+70.0);
		if ( *(intb.end() - 1) < 0) {
			std::cout << "??? " << *(intb.end() - 1) << std::endl;
		}
	}
	
}

bool check_extensions(int count , char* files[]){
	for (int i = NUM_ARGS +1; i < count ; i++){
		std::string tmp = files[i];
		boost::filesystem::path error_path(tmp);
		std::string extension = error_path.extension().c_str();
		int txt = extension.compare(".txt");
		int csv = extension.compare(".csv");
		int CSV = extension.compare(".CSV");
		bool filetype_ok = false;
		if (!txt) {
			filetype_ok = true;
		}
		if (!csv) {
			filetype_ok = true;
		}
		if (!CSV) {
			filetype_ok = true;
		}
		if (!filetype_ok) {
			std::cerr << "Filetype? " << files[i]<<'\t'<<i <<  '\t' << NUM_ARGS <<  std::endl;
			return false;
		}
	}
	return true;
	
}

int main(int argc , char* argv[]){
		
		//Program Options

	po::options_description desc("Allowed Options");
	desc.add_options()
		 ("help,h", "Produce this help message")
		 ("startwl,s",po::value<double>(),"Set the start Wavelength for the Analysis")
		 ("stopwl,p",po::value<double>(),"Set the stop Wavelength for the Analysis")
		 ("non-interactive,n","Runs the program in Noninteractive mode. It quits when it's finished")
		 ("version,v","Prints Version")
	;
		 
	po::variables_map vm;
		 po::store(po::parse_command_line(argc,argv,desc),vm);
		 po::notify(vm);
	if (vm.count("help")) {
		std::cout << desc<< std::endl;
		return 3;
	}
	if (vm.count("version")) {
		std::cout << "VCSEL Laser Analysis Version " << _VERSION << std::endl;
		std::cout << "Using ROOT version " << _ROOT_VERSION << " and Boost version " << _BOOST_VERSION << std::endl;
		return 0;
	}
	
	if (argc < 4) {
		std::cout << desc;
		return 2;
	}	
	double startwl, stopwl;
	startwl = 842.;
	stopwl = 860.;
	bool run = true;
	if (vm.count("startwl")) {
		startwl = vm["startwl"].as<double>();
		NUM_ARGS +=2;
	}
	if (vm.count("stopwl")) {
		double tmp =  vm["stopwl"].as<double>();
		stopwl =tmp;
		NUM_ARGS +=2;
	}
	if (vm.count("non-interactive")) {
		run = false;
		NUM_ARGS++;
	}
	
	
	//checking filetypes must be txt, csv or CSV
	if (!check_extensions(argc, argv)) {
		return 1;
	}
	std::cout <<"startwl: "<< startwl << '\t' << "stopwl: " << stopwl << std::endl;
	
	Double_t max = -210;
	Double_t maxwl = 0;
	int _argc = argc;
	TApplication *t = new TApplication("big",&_argc,argv);
	std::cout << "Running with boost and ROOT" <<std::endl;
	std::vector<double> _x,_y;
	Double_t x[LINES], y[LINES], _inta[LINES], _intb[LINES]; 
	
	Double_t *cmp_int = new Double_t[argc];
	Double_t *argc_ary = new Double_t[argc];
	Double_t *cmp_int_root = new Double_t[argc];
	Double_t *asymmety_ary = new Double_t[argc];
	Double_t *width_ary = new Double_t [argc];
	

	
	TGraph2D *gr = new TGraph2D(LINES*(argc-1));
		//Setting up canvas for plot of all sectrums (is it called spectrums? ;) )
	TCanvas *c1 = new TCanvas("All Plots","All Plots",10,10,3000,1500);
	TH1F *integral_hist = new TH1F("Asymmerty", "Asymmetry", 100,0, 100);
	

	if(!(argc % ROWS)){
		c1->Divide(argc/ROWS,ROWS);
		
	}else{
		c1->Divide(argc/ROWS+(argc %ROWS -1),ROWS);
	}
	
	for (Int_t i = NUM_ARGS +1; i < argc ; i++){
		try{ 
			
			max = -211;
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
				read_file(line, a, b, inta, intb);
				cline++;
			}
			
			if (cline < LINES){
				for(int i = cline ; i < LINES ; i++){
					a.push_back(100);
					b.push_back(-70);
				}
			}
			std::cout<< "\n\ncline: " << cline<< std::endl; 
			cline =(cline > LINES) ? LINES :cline;
			
			for(Int_t j = 0; j <LINES ;j++){
				x[j] = a[j];
				y[j] = b[j];
				_inta[j] = inta[j];
				_intb[j]= (intb[j] < 0)? 0:intb[j];
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
			integral_hist->Fill(s_integral);
			cmp_int[i] = s_integral;
			Int_t lines = (Int_t)intb.size();
			TGraph *r_integral = new TGraph(lines, _inta, _intb);
			
			std::cout << "ROOT integral: " << r_integral->Integral() << std::endl;
			cmp_int_root[i] = r_integral->Integral();
			
				//expanding
				//expand(y, THRS_EXPAND, RATIO_EXPAND, LINES);
			
			
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
			_gr->Draw("AP");
			_gr->GetYaxis()->SetRangeUser(-80.,-10.);
			_gr->GetXaxis()->SetRangeUser(startwl,stopwl);
			_gr->SetTitle(tmp.c_str());
			c1->Update();
			
			
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
			asymmety_ary[i-NUM_ARGS] = calced_asy;
			
			std::cout << "Asymmetry: " << calced_asy << std::endl;
			
		}catch(std::exception e){
			std::cout << e.what()<< std::endl;
		}
	}
	
	
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
		//settig everything to print fitresuts
	gStyle->SetOptStat(1211);
	gStyle->SetOptFit(1111);
	integral_hist->Draw();
	integral_hist->Fit("gaus","W","" ,10,100);
		//integral_hist->Draw("SAME");
	d->Update();
	d->cd(3);
	
	TGraph *roots_int = new TGraph(argc-1, argc_ary, cmp_int_root);
	roots_int->SetTitle("ROOTS Int");
	roots_int->Draw("A*");
	d->Update();
	d->cd(4);
	d->Update();
		//gROOT->SetStyle("modern");
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
	file += "_big.png";
	img->FromPad(d);
	img->WriteImage(file.c_str());
		//cleaning
	
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
	file += "_asy_cent.png";
	secimg->FromPad(e);
	secimg->WriteImage(file.c_str());
	
	TImage *thrdimg = TImage::Create();
	boost::filesystem::path p3(t->Argv(3));
	file = p3.parent_path().string();
	file += "_allplots.png";
	thrdimg->FromPad(c1);
	thrdimg->WriteImage(file.c_str());
	
		//detecting Gradients
	gradient(asymmety_ary, width_ary,cmp_int, argc-1,c1);
	std::cout << "\n\n\nDone !!\nYou can quit now using CTRL+C \n" ;
	
	if (run == true){
		t->Run();
	}
	std::cout << "With \n" ;
	
	delete[] cmp_int;
	delete[] argc_ary; 
	delete[] cmp_int_root;
	delete[] asymmety_ary;
	delete[] width_ary;
	return 0;
}
