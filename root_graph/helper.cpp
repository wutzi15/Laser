//
//  helper.cpp
//  root_graph
//
//  Created by Wutzi on 26.05.11.
//  Copyright 2011 me. All rights reserved.
//
#include "pref.h"
#include "helper.h"
#include "defines.h"

bool check(std::string a){
		//Check if string is a good double
	try{
		double d = boost::lexical_cast<double>(a);
		(void)d;
		return true;
	}catch(boost::bad_lexical_cast){
		return false;
	}
}

double findlower(double *x,double *y, double max ){
	for(int k = 0;k< LINES; k++){
		if(y[k] >= max - 30){
			return  x[k];
			
		}
	}
	return 1;
}


double findupper(double *x,double *y , double max){
	for(int l = LINES-1 ;l >= 0 ; l--){
		if(y[l]>= max -30){
			return x[l];
		}
	}
	return 1;
}

inline bool ispm(double num, double pm){
	if (num >= pm+PM_OFFSET || num <=PM_OFFSET-pm) {
		return true;
	}
	return false;
}

void label(TCanvas *canv, int pos, Level level){
	if (level == WARNING) {
		TPad *pad =(TPad*)canv->cd(pos);
		std::string title = pad->GetTitle();
		title += " WARNING!";
		pad->SetTitle(title.c_str());
		
	}
	if (level == NOTICE) {
		TPad *pad =(TPad*)canv->cd(pos);
		std::string title = pad->GetTitle();
		title += " NOTICE!";
		pad->SetTitle(title.c_str());
	}
}

void expand(double *y,double threshold,double ratio, int count){
	for (int i = 0; i< count; i++) {
		if (y[i] < threshold) {
			y[i] -= fabs(y[i] -threshold)*ratio;
		}
	}
}







