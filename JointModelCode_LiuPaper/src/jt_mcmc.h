
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"

void mcmc(Subject_type *sublist,Common_parms *parms_l,Common_parms *parms_f,double **ts_l,double **ts_f,long iter,int N,
          Priors *priors,unsigned long *seed,char *filel,char *file2,double propvar[],Hyper_priors *hyper);


void draw_fe_prior_a_mean(Subject_type *sublist,Priors *priors,Common_parms *parms_f,unsigned long *seed,Hyper_priors *hyper);
void draw_fe_prior_w_mean(Subject_type *sublist,Priors *priors,Common_parms *parms_f,unsigned long *seed,Hyper_priors *hyper);
void draw_fe_prior_a_var(Subject_type *sublist,Priors *priors,Common_parms *parms_f,unsigned long *seed,Hyper_priors *hyper);
void draw_bh_mean_f(Subject_type *sublist,Priors *priors,Common_parms *parms,unsigned long *seed,Hyper_priors *hyper);
void draw_bh_mean_l(Subject_type *sublist,Priors *priors,Common_parms *parms,unsigned long *seed,Hyper_priors *hyper);
void draw_fe_priors_w_var(Subject_type *sublist,Priors *priors,double v1,double v2,unsigned long *seed,Hyper_priors *hyper);

void draw_bh_var_l(Subject_type *sublist,Priors *priors,double v1,double v2,unsigned long *seed,Hyper_priors *hyper);
void draw_bh_var_f(Subject_type *sublist,Priors *priors,double v1,double v2,unsigned long *seed,Hyper_priors *hyper);

void draw_fixed_mass(Subject_type *sublist,Priors *priors,Common_parms *parms_l,Common_parms *parms_f,unsigned long *seed,double **pmean_var);
void draw_fixed_width_l(Subject_type *sublist,Priors *priors,Common_parms *parms,double v,unsigned long *seed);
void draw_fixed_width_f(Subject_type *sublist,Priors *priors,Common_parms *parms,double v,unsigned long *seed);

void draw_fe_precision_l(Subject_type *sublist,Priors *priors,Common_parms *parms,double v1,double v2,unsigned long *seed);
void draw_fe_precision_f(Subject_type *sublist,Priors *priors,Common_parms *parms,double v1,double v2,unsigned long *seed);

void draw_bh_l(Subject_type *sublist,Common_parms *parms,Priors *priors,double **ts,
           int N,unsigned long *seed,double **var);
void draw_bh_f(Subject_type *sublist,Common_parms *parms,Priors *priors,double **ts,
           int N,unsigned long *seed,double **var);
void draw_times_f(Subject_type *sublist,Common_parms *parms,Common_parms *parms_cox,double **ts,
             int N,unsigned long *seed,double v) ;
void draw_times_l(Subject_type *sublist,Common_parms *parms,double **ts,
             int N,unsigned long *seed,double v) ;
void draw_random_effects_l(double **ts,Subject_type *sublist,Common_parms *parms,int N,double v1, double v2,unsigned long *seed);
void draw_random_effects_f(double **ts,Subject_type *sublist,Common_parms *parms,int N,double v1, double v2,unsigned long *seed);

double error_squared_l(double **ts,Subject_type *sublist,Common_parms *parms,int N); 
double error_squared_f(double **ts,Subject_type *sublist,Common_parms *parms,int N) ;

void draw_eta_l(Subject_type *sublist,Common_parms *parms,unsigned long *seed,double *v);
void draw_eta_f(Subject_type *sublist,Common_parms *parms,unsigned long *seed,double *v);

double Phi(double y);
void adjust_acceptance(double x,double *X);
void adjust2_acceptance(double x,double **X,double corr);
double gamm(double x) ;

double gammaPdf(double x, double a, double b);
double phi(double, double, double);



