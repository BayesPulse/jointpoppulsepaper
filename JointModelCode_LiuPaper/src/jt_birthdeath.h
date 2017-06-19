//
//  jt_birthdeath.h
//

#include "jt_deconvolution_main.h"
#include "randgen.h"
#include "jt_hash.h"

void birth_death_l(Subject_type *sublist,double **ts,Common_parms *parms,int N,
                 unsigned long *seed,int iter);
void birth_death_f(Subject_type *sublist, double **ts, Common_parms *parms, Common_parms *, int N,
                 unsigned long *seed,int iter);
void mean_contribution(Node_type *node,double **ts,Common_parms *parms,int N,double halflife);
double *mean_concentration(Node_type *list,Common_parms *parms,int N,Node_type *node_out,double baseline);
double likelihood(Subject_type *sublist,double **ts,Common_parms *parms,int N,
                  Node_type *node_out);
double likelihood2(Node_type *list,double **ts,Common_parms *parms,int N,
                  Node_type *node_out,double baseline);
double *calc_death_rate_f(Node_type *parent, Node_type *list, int num_node, double *partial_likelihood, double full_likelihood, double Birth_rate, Common_parms *parms_cox);
double *calc_death_rate_l(Node_type *parent, Node_type *list, int num_node, double *partial_likelihood, double full_likelihood, double Birth_rate, Common_parms *parms_cox);
double *calc_death_rate(Node_type *list, int num_node, double *partial_likelihood, double full_likelihood, double Birth_rate, double r);

double *skernel(Node_type *list, Node_type *parents, Common_parms *parms_cox, int npulse);
double skernel_1(Node_type *fnode, Node_type *parents, Common_parms *parms_cox);


