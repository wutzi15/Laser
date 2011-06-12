//
//  helper.h
//  root_graph
//
//  Created by Wutzi on 26.05.11.
//  Copyright 2011 me. All rights reserved.
//
#ifndef _HELPER
#define _HELPER

#include "defines.h"
bool check(std::string a);
double findlower(double *x,double *y, double max );
double findupper(double *x,double *y , double max);
inline bool ispm(double num, double pm);
void label(TCanvas *canv, int pos, Level level);
void expand(double *y,double threshold,double ratio,int count);

#endif