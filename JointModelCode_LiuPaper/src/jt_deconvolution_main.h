#pragma once
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
/*****THIS INCLUDES THE information of each pulse********/
typedef struct node_tag{
	struct node_tag *succ;
	struct node_tag *pred;
	double time;  /*tau_ki*/
	double theta[2]; /*a_ki, w_ki*/
	/*  double width;*/
	double *mean_contrib;
	double eta[2]; /*kappa for t distribution*/
	double lambda; /*for fsh pulse only, lambda, denomsum  
	                */
} Node_type;


/*mean of the mean and variance random effects of each pulse*/
/*inverse of the variance of the mean and variance pulse random effects*/
/*this includes the information of each subject */
typedef struct{
	double sigma; /*error*/
	double lsigma;
    double Rinter1; /*(hard core)*/
	double Rinter2; /*struass*/
	double *re_precision;  /*nu_a,nu_w,draw_fe_precision same for each subject*/
	int numsub;
	int subindex;
	int iter;
/*	double nprior;*/
	double rho; /*cluster size, used for lh only,on the log scale*/
	double lrho;
	double lnu;
	double nu; /*cluster width, used for lh only,it is the square not standard deviation*/
	double beta; /*intensity function for strauss process*/
	double gamma; /*repulsion function for struass */
} Common_parms;
/*******the prior distribution is collasping all information of fsh and lh*******/

typedef struct{
	double meanbh_f[2]; /*population mean on baseline and halflife, fsh * u_b, u_h,draw_bh_mean*/
	double meanbh_l[2]; /*population mean on baseline and halflife, lh * u_b, u_h,draw_bh_mean*/

	double varbh_f[2]; /*population variance on baseline and halflife, fsh,draw_bh_var */
	double varbh_l[2]; /*population variance  on  baseline and halflife, lh,draw_bh_var */

	double *fe_mean_f;/*ua,uw,,population mean of pulse mass and width draw_fe_priors1*/
	double *fe_mean_l;
	double **fe_precision; /*the variance matrix of correlated pulse mass on the population level */

	double **re_var; /*inverse of fe_precision*/

	double fe_precision_wf;/*sigma_w，drawe_fe_priors2*/
	double fe_precision_wl;/*sigma_w，drawe_fe_priors2*/

	double *re_var_f; /*specified, the half-cauchy or uniform prior for subject-level variance */
	double *re_var_l;
	double alpha_f;  /* prior parameters for model error,fsh */
	double beta_f;
	double alpha_l;  /* prior parameters for model error,lh */
	double beta_l;
    double rho_prior[2]; /*log scale*/
	double nu_prior[2];  /*log scale*/

} Priors;

/*the hyper is also defined together for both the driver and response */
typedef struct {
	double hmean_f[2]; /*ma,mw*/
	double hmean_l[2]; /*ma,mw*/

	double hvar_f; /*vw*/
	double hvar_l; /*vw*/
	
	double precision_wf;
	double precision_wl;
	double **sig_a_inv;

	double **prec;  /* prior variance-covariance of pulse mass,,F-1*/
	double meanmeanbh_l[2];  /*mhn,prior mean fo baseline and half life*/
	double meanvarbh_l[2]; /*prior variance for baseline and halflife */
	double varbh_l[2]; /*prior for variance for baseline and halflife*/
	double meanmeanbh_f[2];  /*mhn,prior mean fo baseline and half life*/
	double meanvarbh_f[2]; /*prior variance for baseline and halflife */
	double varbh_f[2]; /*prior for variance for baseline and halflife*/

} Hyper_priors;

typedef struct subject_tag {
	struct subject_tag *succ;
	struct subject_tag *pred;
	Node_type *driver;
	Node_type *response;
	int numnode_f;
	int numnode_l;
	double theta_f[2]; /*u_ai,u_wi,draw_fixed_effect*/
	double theta_l[2];
	double basehalf_f[2]; /*u_bi, u_hi,draw_bhdd*/
	double basehalf_l[2];
	double decay_f;
	double decay_l;
	char *common_f;
	char *pulse_f;
	FILE *csub_f;
	FILE *psub_f;
	char *common_l;
	char *pulse_l;
	FILE *csub_l;
	FILE *psub_l;
} Subject_type;
