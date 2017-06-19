//
//  birthdeath_pt2.h
//  try1
//double *mean_concentration(Node_type *list,Common_parms *parms,int N,Node_type *node_out,double baseline)

//  Created by prince agyeman on 7/1/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef try1_birthdeath_pt2_h
#define try1_birthdeath_pt2_h
#include "deconvolution_main.h"
#include "randgen.c"
#include "hash.c"

void birth_death_l(Subject_type *sublist,double **ts,Common_parms *parms,int N,
                 unsigned long *seed,int iter);
void birth_death_f(Subject_type *sublist,double **ts,Common_parms *parms,int N,
                 unsigned long *seed,int iter);
void mean_contribution(Node_type *node,double **ts,Common_parms *parms,int N,double halflife);
double *mean_concentration(Node_type *list,Common_parms *parms,int N,Node_type *node_out,double baseline);
double likelihood(Subject_type *sublist,double **ts,Common_parms *parms,int N,
                  Node_type *node_out);
double likelihood2(Node_type *list,double **ts,Common_parms *parms,int N,
                  Node_type *node_out,double baseline);
double *calc_death_rate(Node_type *list,int num_node,double *partial_likelihood,double full_likelihood,double Birth_rate,double r);




#endif
