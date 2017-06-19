/*******************************************************************/
/*************************linklistv2.c *****************************/
/*******************************************************************/

#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"


/*******************************************************************
*******************GLOBAL VARIABLE DEFINITIONS**********************

mmm: Order statistic used for distribution of pulse locations.
This is inputted by the user and is typically 3
adelta_f: counter that increments when a halflife estimate is accetped in the MH draw
ndelta_f: counter that increments for each loop of the halflife MH draw
adelta_l: counter that increments when a halflife estimate is accetped in the MH draw
ndelta_l: counter that increments for each loop of the halflife MH draw

atime_l: counter that increments when a pulse location estimate is accetped in the MH draw
ntime_l: counter that increments for each loop of the pulse location MH draw

atime_f: counter that increments when a pulse location estimate is accetped in the MH draw
ntime_f: counter that increments for each loop of the pulse location MH draw

arem_f: counter that increments when an individual pulse mass estimate is accetped in the MH draw
nrem_f: counter that increments for each loop of the individual pulse mass MH draw
arem_l: counter that increments when an individual pulse mass estimate is accetped in the MH draw
nrem_l: counter that increments for each loop of the individual pulse mass MH draw


arew_l: counter that increments when an individual pulse width estimate is accetped in the MH draw
nrew_l: counter that increments for each loop of the individual pulse width MH draw
arew_f: counter that increments when an individual pulse width estimate is accetped in the MH draw
nrew_f: counter that increments for each loop of the individual pulse width MH draw

afepmv_l: counter that increments when an overall mean pulse mass standard deviation (sma) is accetped in the MH draw
nfepmv_l: counter that increments for each loop of the overall mean pulse mass standard deviation (sma) MH draw
afepwv_l: counter that increments when an overall mean pulse width standard deviation (smw) is accetped in the MH draw
nfepwv_l: counter that increments for each loop of the overall mean pulse width standard deviation (smw) MH draw

afepmv_f: counter that increments when an overall mean pulse mass standard deviation (sma) is accetped in the MH draw
nfepmv_f: counter that increments for each loop of the overall mean pulse mass standard deviation (sma) MH draw
afepwv_f: counter that increments when an overall mean pulse width standard deviation (smw) is accetped in the MH draw
nfepwv_f: counter that increments for each loop of the overall mean pulse width standard deviation (smw) MH draw

afepm: counter that increments when an overall mean pulse mass standard deviation (sma) is accetped in the MH draw
nfepm: counter that increments for each loop of the overall mean pulse mass standard deviation (sma) MH draw
afepwv_l: counter that increments when an overall mean pulse width standard deviation (smw) is accetped in the MH draw
nfepwv_l: counter that increments for each loop of the overall mean pulse width standard deviation (smw) MH draw

afepwv_f: counter that increments when an overall mean pulse width standard deviation (smw) is accetped in the MH draw
nfepwv_f: counter that increments for each loop of the overall mean pulse width standard deviation (smw) MH draw

afemv: counter that increments when an overall pulse mass standard deviation (sa) is accetped in the MH draw
nfemv: counter that increments for each loop of the overall pulse mass standard deviation (sa) MH draw


afewv_l: counter that increments when an overall pulse width standard deviation (sw) is accetped in the MH draw
nfewv_l: counter that increments for each loop of the overall pulse width standard deviation (sw) MH draw
afewv_f: counter that increments when an overall pulse width standard deviation (sw) is accetped in the MH draw
nfewv_f: counter that increments for each loop of the overall pulse width standard deviation (sw) MH draw

afebv_l: counter that increments when a baseline standard deviation (sb) is accetped in the MH draw
nfebv_l: counter that increments for each loop of the baseline standard deviation (sb) MH draw
afehv_l: counter that increments when a half-life standard deviation (sh) is accetped in the MH draw
nfehv_l: counter that increments for each loop of the half-life standard deviation (sh) MH draw
afebv_f: counter that increments when a baseline standard deviation (sb) is accetped in the MH draw
nfebv_f: counter that increments for each loop of the baseline standard deviation (sb) MH draw
afehv_f: counter that increments when a half-life standard deviation (sh) is accetped in the MH draw
nfehv_f: counter that increments for each loop of the half-life standard deviation (sh) MH draw
ae_f:    for drawing the eta
ne_f:
ae_l:
ne_l:
fitstart: The first time in hours that a pulse may occur
fitend: The last time in hours that a pulse may occur

*********************************************************************/

extern int mmm;
long adelta_l = 0;
long ndelta_l = 0;
long atime_l = 0;
long ntime_l = 0;
long arem_l = 0;
long nrem_l = 0;
long arew_l = 0;
long nrew_l = 0;
long afepw_l = 0;
long nfepw_l = 0;
long ae_ml = 0;
long ne_ml = 0;
long ae_wl = 0;
long ne_wl = 0;

long afepmv_l = 0;
long nfepmv_l = 0;
long afepwv_l = 0;
long nfepwv_l = 0;
long afewv_l = 0;
long nfewv_l = 0;
long afebv_l = 0;
long nfebv_l = 0;
long afehv_l = 0;
long nfehv_l = 0;
long afepm = 0;
long nfepm = 0;
long arho = 0;
long nrho = 0;
long anu = 0;
long nnu = 0;


long adelta_f = 0;
long ndelta_f = 0;
long atime_f = 0;
long ntime_f = 0;
long arem_f = 0;
long nrem_f = 0;
long arew_f = 0;
long nrew_f = 0;
long afepw_f = 0;
long nfepw_f = 0;
long afepmv_f = 0;
long nfepmv_f = 0;
long afepwv_f = 0;
long nfepwv_f = 0;
long afewv_f = 0;
long nfewv_f = 0;
long afebv_f = 0;
long nfebv_f = 0;
long afehv_f = 0;
long nfehv_f = 0;
extern double fitstart;
extern double fitend;
long afemv = 0;
long nfemv = 0;
long ae_mf = 0;
long ne_mf = 0;
long ae_wf = 0;
long ne_wf = 0;

/********************************************************************/
/*SUBROUTINES THAT EXIST IN THIS PROGRAM

 mcmc
 draw_fe_priors
 draw_fe_priors2
 draw_bh_mean
 draw_bh_var
 draw_fixed_effects
 draw_fe_precision
 draw_bh
 draw_times
 draw_random_effects
 error_squared
 adjust_acceptance
 adjust2_acceptance

 **********************************************************************/


/*********************************************************************/
/*START OF MCMC SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*mcmc: this runs the BDMCMC process
	ARGUMENTS: Subject_type *sublist; this is the current list of subjects that exist;
	Common_parms *parms; the current values of the common parameters;
	double **ts; this is the matrix of observed data (a column of
	times and S columns of log(concentration);
	long iter; the number of iterations to run;
	int N; the number of observations in **ts;
	Priors *priors; the parameters of the prior distributions;
	unsigned long *seed; seed values needed for the randon number generator;
	char *file1; the output file name for common parameters
	double propvar[]; vector of proposal variances input by user
	Hyper_priors *hyper; hyperpriors input into algorithm by user
	RETURNS: None; all updates are made internally
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i,j,k,l: generic counters
 num_node: counter of num_node for loops
 num_node2: number of pulses
 NN=50: Output every NNth iteration to files
 NNN=5000: Output every NNNth interation to screen or output file
 vtime: proposal variance for individual pulse locations
 vrem: proposal variance for individual masses
 vrew: proposal variance for individual widths
 vfepmv: proposal variance for standard deviation of overall mean mass, sma
 vfepwv: proposal variance for standard deviation of overall mean width, smw
 vfemv: proposal variance for standard deviation of pulse mass, sa
 vfewv: proposal variance for standard deviation of pulse width, sw
 vfebv: proposal variance for standard deviation of baseline, sb
 vfehv: proposal variance for standard deviation of half-life, sh
 ssq: sum of squared differences between log(concentration) and expected value
 **pmd_var: proposal variance-covariance matrix for baseline and halflife
 **pmd_vch: cholesky decomposed proposal variance-covariance matrix for
 baseline and halflife
 *new_node: index used to go through nodes in file output
 *subject: index used to go through subjects in file output
 *common: name given to inputted file name for common parameters
 *parm: name given to inputted file name for pulse parameters

 SUBROUTINES USED
 inverse_gamma: found in randgen.h; draws from the inverse gamma distribution
 error_squared: found in this file; calculates sum of squared error
 draw_bh: found in this file; draws baseline and half-life
 draw_times: found in this file; draws individual pulse locations
 birth_death: found in birthdeath.c; runs birthdeath process
 draw_fixed_effects: found in this file; draws subject specific mean mass and width
 draw_fe_precision: found in this file; draws st dev of pulse mass and width
 draw_random_effects: found in this file; draws individual masses and widths
 draw_bh_mean: found in this file; draws mean of baseline and half-life
 draw_bh_var: found in this file; draws st dev of baseline and half-life
 draw_fe_priors: found in this file; draws mean of mean pulse mass and width
 draw_fe_priors2: found in this file; draws st dev of mean pulse mass and width
 adjust_acceptance: found in this file; adjusts proposal variance
 adjust_acceptance2: found in this file; adjusts 2-d proposal variance*/
/***********************************************************************/
void mcmc(Subject_type *sublist, Common_parms *parms_l, Common_parms *parms_f, double **ts_l, double **ts_f, long iter, int N,
	Priors *priors, unsigned long *seed, char *filel, char *file2, double propvar[], Hyper_priors *hyper)
{
	int i, j, k, l, num_node, num_node2, NN = 50, NNN = 5000;
	double vfepmv_l, **sigma_p, vfepwv_l, vrem_l, vrew_l, vfemv_l, vfepw_l, vnu, vrho, vfepw_f, vfewv_l, vfebv_l, vfehv_l, veta_l[2], veta_f[2], vfepmv_f, vfepwv_f, vrem_f, vrew_f, vfemv_f, vfewv_f, vfebv_f, vfehv_f, vtime_l, ssq_l, vtime_f, ssq_f, **pmd_var_l, **pmd_vch_l, **pmd_var_f, **pmd_vch_f, **pmean_var;
	Node_type *new_node, *tdriver, *tresponse;
	Subject_type *subject;
	FILE *common_l, *common_f;
	int cholesky_decomp(double **A, int num_col);
	double *skernel(Node_type *, Node_type *, Common_parms *, int);
	double phi(double, double, double);
	double **cholesky_invert(int, double **);

	Common_parms *tparms_l, *tparms_f;
	Priors *tpriors;
	Subject_type *tsubject, *tsublist, *ttsublist, *ttsubject;

	double inverse_gamma(double alpha, double beta, unsigned long *seed);
	double error_squared_l(double **, Subject_type *, Common_parms *, int);
	double error_squared_f(double **, Subject_type *, Common_parms *, int);
	void draw_fe_prior_a_mean(Subject_type *sublist, Priors *priors, Common_parms *parms_f, unsigned long *seed, Hyper_priors *hyper);
	void draw_fe_prior_w_mean(Subject_type *sublist, Priors *priors, Common_parms *parms_f, unsigned long *seed, Hyper_priors *hyper);
	void draw_fe_prior_a_var(Subject_type *sublist, Priors *priors, Common_parms *parms_f, unsigned long *seed, Hyper_priors *hyper);
	void draw_bh_mean_f(Subject_type *sublist, Priors *priors, Common_parms *parms, unsigned long *seed, Hyper_priors *hyper);
	void draw_bh_mean_l(Subject_type *sublist, Priors *priors, Common_parms *parms, unsigned long *seed, Hyper_priors *hyper);
	void draw_fe_priors_w_var(Subject_type *sublist, Priors *priors, double v1, double v2, unsigned long *seed, Hyper_priors *hyper);

	void draw_bh_l(Subject_type *, Common_parms *, Priors *, double **,
		int, unsigned long *, double **);
	void draw_bh_f(Subject_type *, Common_parms *, Priors *, double **,
		int, unsigned long *, double **);


	void draw_times_f(Subject_type *, Common_parms *, Common_parms *, double **,
		int, unsigned  long *, double);
	void draw_times_l(Subject_type *, Common_parms *, double **,
		int, unsigned long *, double);
	void draw_times_old(Subject_type *, Common_parms *, double **,
		int, unsigned long *, double);
	void draw_times_old_f(Subject_type *, Common_parms *, double **,
		int, unsigned long *, double);

	void draw_bh_var_l(Subject_type *, Priors *, double, double, unsigned long *, Hyper_priors *);
	void draw_bh_var_f(Subject_type *, Priors *, double, double, unsigned long *, Hyper_priors *);


	void birth_death_l(Subject_type *, double **, Common_parms *, int,
		unsigned long *, int);
	void birth_death_f(Subject_type *, double **, Common_parms *, Common_parms *, int,
		unsigned long *, int);

	void draw_fixed_mass(Subject_type *, Priors *, Common_parms *, Common_parms *, unsigned long *, double **);
	void draw_fixed_width_l(Subject_type *, Priors *, Common_parms *, double, unsigned long *);
	void draw_fixed_width_f(Subject_type *, Priors *, Common_parms *, double, unsigned long *);
	void draw_fe_precision_l(Subject_type *, Priors *, Common_parms *, double, double, unsigned long *);
	void draw_fe_precision_f(Subject_type *, Priors *, Common_parms *, double, double, unsigned long *);

	void draw_bh_l(Subject_type *, Common_parms *, Priors *, double **,
		int, unsigned long *, double **);
	void draw_bh_f(Subject_type *, Common_parms *, Priors *, double **,
		int, unsigned long *, double **);
	void draw_random_effects_l(double **, Subject_type *, Common_parms *, int, double, double, unsigned long *);
	void draw_random_effects_f(double **, Subject_type *, Common_parms *, int, double, double, unsigned long *);


	void draw_eta_l(Subject_type *, Common_parms *, unsigned long *, double *);
	void draw_eta_f(Subject_type *, Common_parms *, unsigned long *, double *);
	void destroy_sublist(Subject_type *);
	Subject_type *initialize_subject(void);
	void insert_subject(Subject_type *, Subject_type *);

	double Phi(double y);
	void adjust_acceptance(double x, double *X);
	void adjust2_acceptance(double x, double **X, double corr);
	double gamm(double x);
	double gammalog(double x, double a, double b);
	double gammaPdf(double x, double a, double b);

	void adjust_acceptance(double, double *);
	void adjust2_acceptance(double, double **, double);
	void mh_logalphamean(Subject_type *, Common_parms *, Priors *, unsigned long *, double);
	void mh_logsigmean(Subject_type *, Common_parms *, Priors *, unsigned long *, double);



	
 
strcat(filel, ".out");
	common_l = fopen(filel, "w");

	strcat(file2, ".out");
	common_f = fopen(file2, "w");

	subject = sublist->succ;
	while (subject != NULL){
		subject->csub_l = fopen(subject->common_l, "w");
		subject->psub_l = fopen(subject->pulse_l, "w");
		subject->csub_f = fopen(subject->common_f, "w");
		subject->psub_f = fopen(subject->pulse_f, "w");
		subject = subject->succ;
	}
	sigma_p = (double **)calloc(2, sizeof(double *));
	for (i = 0; i < 2; i++)
		sigma_p[i] = (double *)calloc(2, sizeof(double));
	/*proposal variances for poulation variance of pulse width,sigma_wl,sigma_wf*/
	vfewv_l = propvar[0];
	vfewv_f = propvar[1];

	/*proposal variances for population variance of baseline and halflife,sigma_b, sigma_h)*/
	vfebv_l = propvar[2];
	vfehv_l = propvar[3];
	vfebv_f = propvar[4];
	vfehv_f = propvar[5];
	/*proosal variancee for subject pulse mass pair*/
	/*read in three parameters, lh, fsh, and correlation*/
	/*the output is and 2 by 2 matrix,big_sigma_m*/
	pmean_var = (double **)calloc(2, sizeof(double *));
	for (i = 0; i < 2; i++)
		pmean_var[i] = (double *)calloc(2, sizeof(double));
	pmean_var[0][0] = propvar[6];
	pmean_var[1][1] = propvar[7];
	pmean_var[0][1] = pmean_var[1][0] = sqrt(pmean_var[1][1] * pmean_var[0][0])* propvar[8];
	/*proposal variances for individual width mean,mu_wi*/
	vfepw_l = propvar[9];
	vfepw_f = propvar[10];

	/*proposal variances for common variances of subject level v_pm */
	vfepmv_l = propvar[11];
	vfepwv_l = propvar[12];
	vfepmv_f = propvar[13];
	vfepwv_f = propvar[14];

	/*proposal variances for subject baselines and halflifes*/
	pmd_var_l = (double **)calloc(2, sizeof(double *));
	for (i = 0; i < 2; i++)
		pmd_var_l[i] = (double *)calloc(2, sizeof(double));
	pmd_var_f = (double **)calloc(2, sizeof(double *));
	for (i = 0; i < 2; i++)
		pmd_var_f[i] = (double *)calloc(2, sizeof(double));


	pmd_var_l[0][0] = propvar[15];
	pmd_var_l[1][1] = propvar[16];
	pmd_var_l[0][1] = pmd_var_l[1][0] = propvar[17] * sqrt(pmd_var_l[0][0])*sqrt(pmd_var_l[1][1]);
	pmd_var_f[0][0] = propvar[18];
	pmd_var_f[1][1] = propvar[19];
	pmd_var_f[0][1] = pmd_var_f[1][0] = propvar[20] * sqrt(pmd_var_f[0][0])*sqrt(pmd_var_f[1][1]);

	/*proposa variance for pulse mass and width*/
	vrem_l = propvar[21];
	vrew_l = propvar[22];
	vrem_f = propvar[23];
	vrew_f = propvar[24];


	/*proposal variance for pulse locations*/
	vtime_l = propvar[25];
	vtime_f = propvar[26];

	/*eta*/
	veta_l[0] = propvar[27];
	veta_l[1] = propvar[28];

	veta_f[0] = propvar[29];
	veta_f[1] = propvar[30];
	vrho = propvar[31];
	vnu = propvar[32];

	/***cholesky decompose the proposal var-covar matrix for b and hl***/
	pmd_vch_l = (double **)calloc(2, sizeof(double *));
	for (i = 0; i < 2; i++)
		pmd_vch_l[i] = (double *)calloc(2, sizeof(double));

	for (i = 0; i < 2; i++)
	for (j = 0; j < 2; j++)
		pmd_vch_l[i][j] = pmd_var_l[i][j];

	if (!cholesky_decomp(pmd_vch_l, 2)){
		printf("not PSD matrix A\n");
		exit(0);
	}


	pmd_vch_f = (double **)calloc(2, sizeof(double *));
	for (i = 0; i < 2; i++)
		pmd_vch_f[i] = (double *)calloc(2, sizeof(double));

	for (i = 0; i < 2; i++)
	for (j = 0; j < 2; j++)
		pmd_vch_f[i][j] = pmd_var_f[i][j];

	if (!cholesky_decomp(pmd_vch_f, 2)){
		printf("not PSD matrix A\n");
		exit(0);
	}

	/*********type in the true value  ***********/

	for (i = 0; i < iter; i++) {

		/* run the birth-death algorithm to select pulses */

		parms_f->iter = i;
		parms_l->iter = i;



		/* draw fixed effects priors; mua and muw,the population mean*/

		draw_fe_prior_a_mean(sublist, priors, parms_f, seed, hyper);
		draw_fe_prior_w_mean(sublist, priors, parms_f, seed, hyper);



		/* draw fixed effects prior variances; sma and smw*/

		draw_fe_prior_a_var(sublist, priors, parms_f, seed, hyper);
		draw_fe_priors_w_var(sublist, priors, vfewv_l, vfewv_f, seed, hyper);  /*proposal variance */

		/* draw mean baseline and halflife; mub and muh*/

		draw_bh_mean_f(sublist, priors, parms_f, seed, hyper);
		draw_bh_mean_l(sublist, priors, parms_l, seed, hyper);

		///* draw variance for baseline and halflife; sb and sh*/

		draw_bh_var_l(sublist, priors, vfebv_l, vfehv_l, seed, hyper);
		draw_bh_var_f(sublist, priors, vfebv_f, vfehv_f, seed, hyper);

		/////* draw subject specific means; muak and muwk*/

		draw_fixed_mass(sublist, priors, parms_l, parms_f, seed, pmean_var);
		draw_fixed_width_l(sublist, priors, parms_l, vfepw_l, seed);
		draw_fixed_width_f(sublist, priors, parms_f, vfepw_f, seed);

		/* draw mass and width variances; sa and sw*/

		draw_fe_precision_l(sublist, priors, parms_l, vfepmv_l, vfepwv_l, seed);
		draw_fe_precision_f(sublist, priors, parms_f, vfepmv_f, vfepwv_f, seed);

		birth_death_l(sublist, ts_l, parms_l, N, seed, i);
		//
		
			draw_times_l(sublist, parms_l, ts_l,
				N, seed, vtime_l);
		
			////////
			//draw_times_old(tsublist, tparms_l, ts_l,
			//	N, seed, vtime_l);
			draw_random_effects_l(ts_l, sublist, parms_l, N, vrem_l, vrew_l, seed);
			draw_eta_l(sublist, parms_l, seed, veta_l);
			

		draw_bh_l(sublist, parms_l, priors, ts_l,
			N, seed, pmd_var_l);
		//////

		/////* draw pulse locations; tauki*/
		birth_death_f(sublist, ts_f, parms_f, parms_l, N, seed, i);
		////draw_times_old_f(tsublist, tparms_f, ts_f,
		////	N, seed, vtime_f);
		
		draw_times_f(sublist, parms_f, parms_l, ts_f,
			N, seed, vtime_f);
		draw_eta_f(sublist, parms_f, seed, veta_f);
	
		//
		draw_random_effects_f(ts_f, sublist, parms_f, N, vrem_f, vrew_f, seed);

		draw_bh_f(sublist, parms_f, priors, ts_f,
			N, seed, pmd_var_f);


		/////* draw pulse mass and width; Aki and s2pki*/

		
		mh_logalphamean(sublist, parms_l, priors, seed, vrho);
		mh_logsigmean(sublist, parms_l, priors, seed, vnu);


		/* draw model error; s2e*/
		ssq_l = error_squared_l(ts_l, sublist, parms_l, N);
		parms_l->sigma = inverse_gamma(priors->alpha_l + (parms_l->numsub*N) / 2, priors->beta_l + 0.5*ssq_l, seed);
		parms_l->lsigma = log(parms_l->sigma);
		ssq_f = error_squared_l(ts_f, sublist, parms_f, N);
		parms_f->sigma = inverse_gamma(priors->alpha_f + (parms_f->numsub*N) / 2, priors->beta_f + 0.5*ssq_f, seed);
		parms_f->lsigma = log(parms_f->sigma);

		fflush(stdout);

		/***************************************************************/
	
		/* Every 50th iteration, print estimates to files*/
		if (!(i%NN)) {
			subject = sublist->succ;
			while (subject != NULL){
				num_node = 0;
				new_node = subject->driver->succ;
				while (new_node != NULL){
					num_node++;
					new_node = new_node->succ;
				}
				/*Print subject specific parameters to c1sk.out files*/

				fprintf(subject->csub_l, "%d %d %lf %lf %lf %lf %lf %lf \n", i, num_node, subject->theta_l[0], subject->theta_l[1], subject->basehalf_l[0], subject->basehalf_l[1], parms_l->rho, parms_l->nu);
				fflush(subject->csub_l);
				new_node = subject->driver->succ;
				num_node2 = 0;
				/*Print pulse specific parameters to p1sk.out files*/
				while (new_node != NULL) {
					fprintf(subject->psub_l, "%d %d %d %lf %lf %lf %lf %lf\n", i, num_node, num_node2, new_node->theta[0], new_node->theta[1], new_node->time, new_node->eta[0], new_node->eta[1]);
					num_node2++;
					new_node = new_node->succ;
				}

				fflush(subject->psub_l);
				num_node = 0;
				new_node = subject->response->succ;
				while (new_node != NULL){
					num_node++;
					new_node = new_node->succ;
				}
				/*Print subject specific parameters to c1sk.out files*/
				fprintf(subject->csub_f, "%d %d %lf %lf %lf %lf\n", i, num_node, subject->theta_f[0], subject->theta_f[1], subject->basehalf_f[0], subject->basehalf_f[1]);
				fflush(subject->csub_f);
				new_node = subject->response->succ;
				num_node2 = 0;
				/*Print pulse specific parameters to p1sk.out files*/
				while (new_node != NULL) {
					fprintf(subject->psub_f, "%d %d %d %lf %lf %lf %lf %lf\n", i, num_node, num_node2, new_node->theta[0], new_node->theta[1], new_node->time, new_node->eta[0], new_node->eta[1]);
					num_node2++;
					new_node = new_node->succ;
				}

				fflush(subject->psub_f);

				
				subject = subject->succ;


			}
			for (k = 0; k < 2; k++)
			for (j = 0; j < 2; j++)
				sigma_p[k][j] = priors->fe_precision[k][j];

			if (!cholesky_decomp(sigma_p, 2)){
				printf("not PSD matrix A\n");
				exit(0);
			}
			priors->re_var = cholesky_invert(2, sigma_p);
			/*Print common parameters to c1.out file*/
			fprintf(common_l, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf   %lf %lf\n", i, parms_l->nu,parms_l->rho,priors->fe_mean_l[0], priors->fe_mean_l[1], priors->re_var[0][0], priors->re_var[0][1], priors->fe_precision_wl, parms_l->re_precision[0], parms_l->re_precision[1],priors->meanbh_l[0], priors->meanbh_l[1], priors->varbh_l[0], priors->varbh_l[1], parms_l->sigma);
			fflush(common_l);

			fprintf(common_f, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", i, priors->fe_mean_f[0], priors->fe_mean_f[1], priors->re_var[1][1], priors->re_var[0][1],priors->fe_precision_wl, parms_f->re_precision[0], parms_f->re_precision[1], priors->meanbh_f[0], priors->meanbh_f[1], priors->varbh_f[0], priors->varbh_f[1], parms_f->sigma);
			fflush(common_f);
			fflush(stdout);
		}

		/*Print acceptance rates*/
		if (!(i%NNN)) {
			printf("\n\n");
			printf("iter = %d \n", i);
			printf("lh: mu_b %.2lf mu_h %.4lf mu_a %.2lf mu_w %.2lf  v %.4le\n",
				priors->meanbh_l[0], priors->meanbh_l[1], priors->fe_mean_l[0], priors->fe_mean_l[1], parms_l->sigma);
			printf("fsh: mu_b %.2lf mu_h %.4lf mu_a %.2lf mu_w %.2lf  v %.4le\n",
				priors->meanbh_f[0], priors->meanbh_f[1], priors->fe_mean_f[0], priors->fe_mean_f[1], parms_f->sigma);
			subject = sublist->succ;
			printf("lh: pct s_ma = %.2lf  pct A_ki = %.2lf , eta_m=%.2lf, eta_w=%.2lf\n",
				(double)afepmv_l / (double)nfepmv_l, (double)arem_l / (double)nrem_l, (double)ae_ml / (double)ne_ml, (double)ae_wl / (double)ne_wl);
			printf("pct s_mw = %.2lf pct s_w = %.2lf pct s2p_ki = %.2lf \n",
				(double)afepwv_l / (double)nfepwv_l, (double)afewv_l / (double)nfewv_l, (double)arew_l / (double)nrew_l);
			printf("pct s_b = %.2lf pct s_h = %.2lf pct B-HL = %.2lf pct time = %.2lf \n",



				(double)afebv_l / (double)nfebv_l, (double)afehv_l / (double)nfehv_l, (double)adelta_l / (double)ndelta_l, (double)atime_l / (double)ntime_l);
			printf("fsh: pct s_ma = %.2lf  pct A_ki = %.2lf , eta_m=%.2lf, eta_w=%.2lf\n",
				(double)afepmv_f / (double)nfepmv_f, (double)arem_f / (double)nrem_f, (double)ae_mf / (double)ne_mf, (double)ae_wf / (double)ne_wf);
			printf("pct s_mw = %.2lf pct s_w = %.2lf pct s2p_ki = %.2lf \n",
				(double)afepwv_f / (double)nfepwv_f, (double)afewv_f / (double)nfewv_f, (double)arew_f / (double)nrew_f);
			printf("pct s_b = %.2lf pct s_h = %.2lf pct B-HL = %.2lf pct time = %.2lf \n \n \n ",
				(double)afebv_f / (double)nfebv_f, (double)afehv_f / (double)nfehv_f, (double)adelta_f / (double)ndelta_f, (double)atime_f / (double)ntime_f);
			
			printf("%d    \n", sublist->succ->numnode_l);
			printf("%d    \n", sublist->succ->numnode_f);

			tdriver = sublist->succ->driver->succ;
			while (tdriver != NULL){
			printf("%lf  %lf  %lf   \n", tdriver->time,tdriver->theta[0],tdriver->theta[1]);
			tdriver = tdriver->succ;
		}

	/*	tresponse = tsublist->succ->response->succ;
			while (tresponse != NULL){
				printf("%lf  %lf  %lf   \n", tresponse->time, tresponse->theta[0], tresponse->theta[1]);		
				tresponse = tresponse->succ;
			}*/


			fflush(stdout);
		}

		/*Every 500th iteration until 25000, update proposal variances*/

		if (!(i % 500) && i<25000 && i >0) {

			adjust2_acceptance((double)adelta_l / (double)ndelta_l, pmd_var_l, propvar[17]);
			adjust2_acceptance((double)adelta_f / (double)ndelta_f, pmd_var_f, propvar[20]);
			adjust_acceptance((double)atime_l / (double)ntime_l, &vtime_l);
			adjust_acceptance((double)arem_l / (double)nrem_l, &vrem_l);
			adjust_acceptance((double)arew_l / (double)nrew_l, &vrew_l);
			adjust_acceptance((double)atime_f / (double)ntime_f, &vtime_f);
			adjust_acceptance((double)arem_f / (double)nrem_f, &vrem_f);
			adjust_acceptance((double)arew_f / (double)nrew_f, &vrew_f);
			adjust_acceptance((double)afepwv_l / (double)nfepwv_l, &vfepwv_l);
			adjust_acceptance((double)afepmv_l / (double)nfepmv_l, &vfepmv_l);
			adjust_acceptance((double)afepwv_f / (double)nfepwv_f, &vfepwv_f);
			adjust_acceptance((double)afepmv_f / (double)nfepmv_f, &vfepmv_f);
			adjust_acceptance((double)afepwv_l / (double)nfepwv_l, &vfepwv_l);
			adjust_acceptance((double)afepwv_f / (double)nfepwv_f, &vfepwv_f);
			adjust_acceptance((double)ae_ml / (double)ne_ml, &veta_l[0]);



			adjust2_acceptance((double)afemv / (double)nfemv, pmean_var, propvar[8]);
			adjust_acceptance((double)afewv_l / (double)nfewv_l, &vfewv_l);
			adjust_acceptance((double)afebv_l / (double)nfebv_l, &vfebv_l);
			adjust_acceptance((double)afehv_l / (double)nfehv_l, &vfehv_l);
			adjust_acceptance((double)afewv_f / (double)nfewv_f, &vfewv_f);
			adjust_acceptance((double)afebv_f / (double)nfebv_f, &vfebv_f);
			adjust_acceptance((double)afehv_f / (double)nfehv_f, &vfehv_f);


			/*After adjusting proposal variances, reset the counters for acceptance rates*/
			adelta_l = ndelta_l = 0;
			adelta_f = ndelta_f = 0;

			atime_l = ntime_l = 0;
			arem_l = nrem_l = 0;
			arew_l = nrew_l = 0;
			atime_f = ntime_f = 0;
			arem_f = nrem_f = 0;
			arew_f = nrew_f = 0;
			afepmv_l = nfepmv_l = 0;
			afepwv_l = nfepwv_l = 0;
			afepwv_f = nfepwv_f = 0;
			afepmv_f = nfepmv_f = 0;

			afemv = nfemv = 0;
			afewv_l = nfewv_l = 0;
			afebv_l = nfebv_l = 0;
			afehv_l = nfehv_l = 0;
			afewv_f = nfewv_f = 0;
			afebv_f = nfebv_f = 0;
			afehv_f = nfehv_f = 0;
			ae_wf = ne_wf = 0;
			ae_wl = ne_wl = 0;
			ae_mf = ne_mf = 0;
			ae_ml = ne_ml = 0;

			/***cholesky decompose the proposal var-covar matrix***/
			for (k = 0; k < 2; k++)
			for (l = 0; l < 2; l++)
				pmd_vch_l[k][l] = pmd_var_l[k][l];

			if (!cholesky_decomp(pmd_vch_l, 2)){
				printf("pmd not PSD matrix\n");
				exit(0);
			}
			for (k = 0; k < 2; k++)
			for (l = 0; l < 2; l++)
				pmd_vch_f[k][l] = pmd_var_f[k][l];

			if (!cholesky_decomp(pmd_vch_f, 2)){
				printf("pmd not PSD matrix\n");
				exit(0);
			}

		} /*End of adjusting acceptance if-then statement*/

	} /*End of loop through iterations*/

	for (i = 0; i < 2; i++)
		free(sigma_p[i]);
	free(sigma_p);
	for (i = 0; i < 2; i++)
		free(pmd_var_f[i]);
	free(pmd_var_f);
	for (i = 0; i < 2; i++)
		free(pmd_var_l[i]);
	free(pmd_var_l);
	for (i = 0; i < 2; i++)
		free(pmd_vch_f[i]);
	free(pmd_vch_f);
	for (i = 0; i < 2; i++)
		free(pmd_vch_l[i]);
	free(pmd_vch_l);

	fclose(common_l);



	fclose(common_f);









} /*End of MCMC */
/*draw the log of cluster size */
void mh_logalphamean(Subject_type *subject, Common_parms *parms_cox, Priors *priors_cox, unsigned long *seed, double palphamean_sd) {
	int nparents, npulse;
	double  NewMeanAlpha, NewMeanLogAlpha, lprop_ratio, lprior_ratio, lratio, alpha, integral, pulselocs;
	Node_type *tmpparent, *fnode;
	Node_type *tmplist;
	Subject_type  * subject2;

	double rnorm(double, double, unsigned long *);
	double kiss(unsigned long *);
	double phi(double, double, double);

	NewMeanLogAlpha = rnorm(parms_cox->lrho, palphamean_sd, seed);
	NewMeanAlpha = exp(NewMeanLogAlpha);

	/*   printf("newmeanalpha %lf palphamean_sd %lf\n",NewMeanAlpha,palphamean_sd);*/
	nrho++;
	lratio = 0;
	lprior_ratio = 0;
	lprop_ratio = 0;

	subject2 = subject->succ;  /****pull information from all the subject  ***/;
	while (subject2 != NULL)

	{
		nparents = subject2->numnode_l;
		npulse = subject2->numnode_f;
		tmpparent = subject2->driver->succ;
		integral = 0;
		while (tmpparent != NULL){
			pulselocs = tmpparent->time;
			integral = phi(fitend, pulselocs, sqrt(parms_cox->nu)) - phi(fitstart, pulselocs, sqrt(parms_cox->nu));
			integral = (-NewMeanAlpha + exp(parms_cox->lrho))*integral;
			lratio += integral;
			tmpparent = tmpparent->succ;
		}
		lratio += npulse*(NewMeanLogAlpha - parms_cox->lrho);

		subject2 = subject2->succ;
	}


	/*    printf("lratio pt1 %lf\n",lratio);*/

	/*log ratio of the prior distributions */
	lprior_ratio = 1 / (2 * priors_cox->rho_prior[1]) * ((parms_cox->lrho - priors_cox->rho_prior[0]) * (parms_cox->lrho - priors_cox->rho_prior[0]) - (NewMeanLogAlpha - priors_cox->rho_prior[0]) * (NewMeanLogAlpha - priors_cox->rho_prior[0]));

	lprop_ratio = lratio + lprior_ratio;

	/*  printf("lprop_ratio %lf lratio %lf lprior_ra  tio %lf\n",lprop_ratio,lratio,lprior_ratio);*/

	alpha = (0 < lprop_ratio) ? 0 : lprop_ratio;

	if (log(kiss(seed)) < alpha) {
		/*  printf("accept\n"); */
		arho++;
		parms_cox->lrho = NewMeanLogAlpha;
		parms_cox->rho = exp(parms_cox->lrho);
		/*update the conditional density function for each response pulse*/
		subject2 = subject->succ;
		while (subject2 != NULL)
		{
			fnode = subject2->response->succ;
			while (fnode != NULL)
			{
				fnode->lambda = skernel_1(fnode, subject2->driver, parms_cox);

				fnode = fnode->succ;
			}
			subject2 = subject2->succ;
		}
	}
}
void mh_logsigmean(Subject_type *subject, Common_parms *parms_cox, Priors *priors_cox, unsigned long *seed, double psigmean_sd) {
	/*declare variables */
	int i, k, npulse2, npulses, nparents;
	double tmp;
	Node_type *tmpparent, *parent, *node, *fnode;
	Node_type *ploclist;
	Subject_type  * subject2;

	double NewMeanLogOmega, NewMeanOmega, lprop_ratio, *sumlogk, sumk, temp, diff_alpha, alpha, integral, tempsum, lprior_ratio, logexpterm, kold, knew, *denomsum, lratio, *pulselocs;

	/*declare functions */
	double rnorm(double, double, unsigned long *);
	double kiss(unsigned long *);
	double phi(double, double, double);
	double *skernel(Node_type *, Node_type *, Common_parms *, int);
	double skernel_1(Node_type *, Node_type *, Common_parms *);
	

			/*draw the new log(c_omega)   ***********/
	NewMeanLogOmega = rnorm(parms_cox->lnu, psigmean_sd, seed);

	NewMeanOmega = exp(NewMeanLogOmega);

    lratio = 0.;
	integral = 0;
	nnu++;
	knew = 0;
	kold = 0; 
	tmp = parms_cox->nu;
	parms_cox->nu = NewMeanOmega;
	subject2 = subject->succ;  /****pull information from all the subject  ***/;

    while (subject2 != NULL){

		nparents = subject2->numnode_l;
		npulses = subject2->numnode_f;
		parent = subject2->driver;
		node = subject2->response;
		tmpparent = subject2->driver->succ;

		
		while (tmpparent != NULL){
			integral += exp(parms_cox->lrho)*(phi(fitend, tmpparent->time, sqrt(exp(NewMeanLogOmega))) - phi(fitstart, tmpparent->time, sqrt(exp(NewMeanLogOmega))) - phi(fitend, tmpparent->time, sqrt(exp(parms_cox->lnu))) + phi(fitstart, tmpparent->time, sqrt(exp(parms_cox->lnu))));
			tmpparent = tmpparent->succ;
		}
		fnode = node->succ;
		/*calculate the sum of old log lambda */
		while (fnode != NULL){
			kold += log(fnode->lambda);

			fnode = fnode->succ;
		}


		
		fnode = node->succ;
		while (fnode != NULL)
		{
			fnode->lambda = skernel_1(fnode, parent, parms_cox);
			knew += log(fnode->lambda);

			fnode = fnode->succ;
		}
		/*calculate the sum of old log lambda */







		subject2 = subject2->succ;
	}
	/*printf("endloop\n");*/
	/*calculate the acceptance ratio */
	lratio = knew - kold - integral;
	lprior_ratio = 1 / (2 * priors_cox->nu_prior[1]) * ((parms_cox->lnu - priors_cox->nu_prior[0]) * (parms_cox->lnu - priors_cox->nu_prior[0]) - (NewMeanLogOmega - priors_cox->nu_prior[0]) * (NewMeanLogOmega - priors_cox->nu_prior[0]));
	lprop_ratio = lratio + lprior_ratio;
	/*lprop_ratio = knew-kold;*/

	/*   printf("lprop_ratio %lf\n",lprop_ratio);*/
	/*make the decision on whether to keep the new value */
	alpha = (0 < lprop_ratio) ? 0 : lprop_ratio;

	if (log(kiss(seed)) < alpha) {
		/*         printf("accept\n");*/
		anu++;
		parms_cox->lnu = NewMeanLogOmega;
	}
	else{
		parms_cox->nu = tmp;
		subject2 = subject->succ;
		while (subject2 != NULL)
		{
			fnode = subject2->response->succ;
					parent = subject2->driver;

			while (fnode != NULL)
			{
				fnode->lambda = skernel_1(fnode, parent, parms_cox);

				fnode = fnode->succ;
			}
			subject2 = subject2->succ;
		}
	}


}


/*********************************************************************/
/*START OF draw_fe_prior_a_mean SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_fe_priors: this runs the Gibbs sampler draw for overall mean mass and width through all subjects;
	ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
	Priors *priors; the current values of the prior parameters;
	Common_parms *parms; the current values of the common parameters;
	unsigned long *seed; seed values needed for the randon number generator;
	Hyper_priors *hyper; parameters of the prior distributions;
	RETURNS: None; all updates are made internally
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 j: generic counter
 fesum: sum of subject specific mean amplitudes or widths;
 part of evaluation of gmean and gvar
 gmean: mean of the distribution used in Gibbs sampler
 gvar: variance of distribution used in Gibbs sampler
 *subject: current list of subjects and their characteristics

 SUBROUTINES USED
 rnorm: found in randgen.h; draws from the normal distribution
 kiss: found in randgen.h; draws from U(0,1) distribution
**************************************************************************/

void draw_fe_prior_a_mean(Subject_type *sublist, Priors *priors, Common_parms *parms_f, unsigned long *seed, Hyper_priors *hyper)
{
	/* declare variables */

	Subject_type *subject;

	int i, j, nsubj,flag;
	double sum_lmass[2], **fcond_var, **fcond_var_inv, temp_cond_mean[2], *cond_mean, temp[2];
	int rmvnorm(double *, double **, int, double *, unsigned long *, int);


	/*declare functions */
	double rnorm(double, double, unsigned long *);
	double kiss(unsigned long *);
	int cholesky_decomp(double **, int);
	double **cholesky_invert(int, double **);
	/* set up space for dynamically allocated variables */
	fcond_var = (double **)calloc(2, sizeof(double *));
	fcond_var_inv = (double **)calloc(2, sizeof(double *));
	cond_mean = (double *)calloc(2, sizeof(double));
	for (i = 0; i < 2; i++) {
		fcond_var[i] = (double *)calloc(2, sizeof(double));
		fcond_var_inv[i] = (double *)calloc(2, sizeof(double));
	}
	nsubj = parms_f->numsub;
	/* sum the current log masses */
	sum_lmass[0] = 0;
	sum_lmass[1] = 0;


	/*priors->fe_precision[0][0] = 5.263;
	priors->fe_precision[1][1] = 5.263;
	priors->fe_precision[0][1] = priors->fe_precision[1][0] = -4.737;*/



	subject = sublist->succ;

	/*This while loop goes through the pulses to get the sum involved*/
	while (subject != NULL){
		sum_lmass[0] += subject->theta_l[0]; /*sum of lh mean*/
		sum_lmass[1] += subject->theta_f[0];  /*sum of fsh mean*/
		subject = subject->succ;
	} /*end of loop through pulses*/
	for (i = 0; i < 2; i++)
	for (j = 0; j < 2; j++)
		fcond_var_inv[i][j] = (double)nsubj * priors->fe_precision[i][j] + hyper->prec[i][j]; /*the inverse*/

	if (!cholesky_decomp(fcond_var_inv, 2)) {
		printf("fcond_var_inv is not PSD matrix\n");
		exit(0);
	}
	/* invert to get the variance matrix.  This is needed for the mean calculation */
	/* and for the draw of the MVN distribution                                                 */
	fcond_var = cholesky_invert(2, fcond_var_inv);
	/* calculate the means of the full cond. Bivariate normal distn */
	for (i = 0; i < 2; i++) {
		temp_cond_mean[i] = 0;

		temp_cond_mean[i] += hyper->prec[i][0] * hyper->hmean_l[0] + hyper->prec[i][1] * hyper->hmean_f[0];
		temp_cond_mean[i] += priors->fe_precision[i][0] * sum_lmass[0] + priors->fe_precision[i][1] * sum_lmass[1];

	}

	for (i = 0; i < 2; i++) {
		cond_mean[i] = 0;
		for (j = 0; j < 2; j++)
			cond_mean[i] += fcond_var[i][j] * temp_cond_mean[j];
	}

	/*draw the new means*/
	flag =0;
	while (flag == 0){
		rmvnorm(temp, fcond_var, 2, cond_mean, seed, 1);
		if (temp[0]>0 && temp[1]>0)
			flag = 1;
	}
	priors->fe_mean_l[0] = temp[0];
	priors->fe_mean_f[0] = temp[1];
	for (i = 0; i < 2; i++) {
		free(fcond_var[i]);
		free(fcond_var_inv[i]);
	}
	free(fcond_var);
	free(fcond_var_inv);
	free(cond_mean);
}

/*draw_fe_prior_w_mean*/
void draw_fe_prior_w_mean(Subject_type *sublist, Priors *priors, Common_parms *parms_f, unsigned long *seed, Hyper_priors *hyper)
{
	/* declare variables */
	int j;
	double fesum[2], gmean[2], gvar[2];
	Subject_type *subject;

	/*declare functions */
	double rnorm(double, double, unsigned long *);
	double kiss(unsigned long *);



	fesum[0] = fesum[1] = 0.0;
	subject = sublist->succ;

	/*This while loop goes through the pulses to get the sum involved*/
	while (subject != NULL){
		fesum[0] += subject->theta_l[1];
		fesum[1] += subject->theta_f[1];

		subject = subject->succ;
	} /*end of loop through pulses*/

	/*Gibbs sampler mean*/

	gmean[0] = (hyper->hvar_l*fesum[0]) + (priors->fe_precision_wl*priors->fe_precision_wl*hyper->hmean_l[1]);
	gmean[0] /= (parms_f->numsub*hyper->hvar_l) + (priors->fe_precision_wl*priors->fe_precision_wl);


	gmean[1] = (hyper->hvar_f*fesum[1]) + (priors->fe_precision_wf*priors->fe_precision_wf*hyper->hmean_f[1]);
	gmean[1] /= (parms_f->numsub*hyper->hvar_f) + (priors->fe_precision_wf*priors->fe_precision_wf);

	/*Gibbs sampler variance*/
	gvar[0] = hyper->hvar_l*priors->fe_precision_wl*priors->fe_precision_wl;
	gvar[0] /= (parms_f->numsub*hyper->hvar_l) + (priors->fe_precision_wl*priors->fe_precision_wl);
	gvar[1] = hyper->hvar_f*priors->fe_precision_wf*priors->fe_precision_wf;
	gvar[1] /= (parms_f->numsub*hyper->hvar_f) + (priors->fe_precision_wf*priors->fe_precision_wf);


	/*This draws the new value for mua/muw using a Gibbs sampler*/
	priors->fe_mean_l[1] = rnorm(gmean[0], sqrt(gvar[0]), seed);
	priors->fe_mean_f[1] = rnorm(gmean[1], sqrt(gvar[1]), seed);



}

void draw_fe_prior_a_var(Subject_type *sublist, Priors *priors, Common_parms *parms_f, unsigned long *seed, Hyper_priors *hyper) {
	/*declare varibles */
	int i, j, nsubj;
	double diff[2], **sigma_a_w, **sigma_a_w_inv;
	Subject_type *subject;
	/*declare functions */
	int cholesky_decomp(double **, int);
	double **cholesky_invert(int, double **);
	void rwishart(double **, double **, int, int, unsigned long *,int);

	/* allocate memory*/
	sigma_a_w_inv = (double **)calloc(2, sizeof(double *));
	sigma_a_w = (double **)calloc(2, sizeof(double *));

	for (j = 0; j < 2; j++) {
		sigma_a_w_inv[j] = (double *)calloc(2, sizeof(double));
		sigma_a_w[j] = (double *)calloc(2, sizeof(double));
	}

	/* create the input matrix for the wishart */
	subject = sublist->succ;
	nsubj = parms_f->numsub;
	/*  priors->fe_mean_l[0] = 20.0;
	  priors->fe_mean_f[0] = 25.0;
	  subject->theta_l[0] = 21.0;
	  subject->theta_f[1] = 25.184;*/

	while (subject != NULL) {
		diff[0] = subject->theta_l[0] - priors->fe_mean_l[0];
		diff[1] = subject->theta_f[0] - priors->fe_mean_f[0];
		for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
			sigma_a_w_inv[i][j] += diff[i] * diff[j];
		subject = subject->succ;

	}

	for (i = 0; i < 2; i++)
	for (j = 0; j < 2; j++)
		sigma_a_w_inv[i][j] += hyper->sig_a_inv[i][j];

	if (!cholesky_decomp(sigma_a_w_inv, 2)) {
		printf("S_inv matrix for full cond of sigma_a_inv is not PSD \n");
		exit(0);
	}

	sigma_a_w = cholesky_invert(2, sigma_a_w_inv);
	if (!cholesky_decomp(sigma_a_w, 2)) {
		printf("S matrix for full cond of sigma_a_inv is not PSK \n");
		exit(0);
	}
	/* acutually should not be decomposed*/
	rwishart(priors->fe_precision, sigma_a_w, 2, 4 + nsubj, seed,1);
	
	for (i = 0; i < 2; i++) {
		free(sigma_a_w[i]);
		free(sigma_a_w_inv[i]);
	}
	free(sigma_a_w);
	free(sigma_a_w_inv);
}




/*********************************************************************/
/*START OF draw_fe_priors2 SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_fe_priors2: this runs the M-H draw for standard deviations of subject
				   level means of mass and width;
				   ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
				   Priors *priors; the current values of the prior parameters;
				   double v1; the proposal variance for overall st-dev of mean mass;
				   double v2; the proposal variance for overall st-dev of mean width;
				   unsigned long *seed; seed values needed for the randon number generator;
				   Hyper_priors *hyper; parameters of the prior distributions;
				   RETURNS: None; all updates are made internally; this does update the global
				   variables regarding acceptance rate: afepmv, afepwv, nfepmv, and nfepwv
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 j,t: generic counters
 *tmp1: current acceptance counters are saved to this vector so we can increment
 inside a loop
 *new_var: proposed values of st. dev. of mass and width
 prop_new: part of evaluation of prop_ratio portion of log(rho)
 prop_old: same as prop_new but with current (not proposed) value
 psum: part of evaluation of prop_ratio portion of log(rho)
 pcomp: part of evaluation of prop_ratio portion of log(rho)
 prior_old: part of evaluation of prior_ratio portion of log(rho)
 prior_new: same as prior_old, but with proposed value
 prior_ratio: comparison between prior distributions under current and proposed values
 prop_ratio: comparison between random effects distributions under current and
 proposed values
 temp: value used in evaluation of log(rho)
 alpha: minimum of 0 and log(rho); used to determine acceptance in M-H
 *subject: current list of subjects and their qualities

 SUBROUTINES USED
 kiss: found in randgen.h; draws from U(0,1) distribution
 rnorm: found in randgen.h; draws from the normal distribution
**************************************************************************/

void draw_fe_priors_w_var(Subject_type *sublist, Priors *priors, double v1, double v2, unsigned long *seed, Hyper_priors *hyper)
{
	int j, t, *tmp1;
	double *new_var, prop_new[2], prop_old[2], psum[2], pcomp[2], prior_old, prior_new;
	double prior_ratio[2], prop_ratio[2], temp, alpha;
	Subject_type *subject;
	double kiss(unsigned long *);
	double rnorm(double, double, unsigned long *);

	/*Allocate Memory*/
	new_var = (double *)calloc(2, sizeof(double));
	tmp1 = (int *)calloc(2, sizeof(int));

	/*Add 1 to the counters for acceptance rates of sigma_ma and sigma_mw*/
	nfewv_l++;
	nfewv_f++;

	/*Assign current acceptance counts to temporary vector*/
	tmp1[0] = afewv_l;
	tmp1[1] = afewv_f;

	/*Draw proposed values for sigma_ma and sigma_mw*/
	new_var[0] = rnorm(priors->fe_precision_wl, v1, seed);
	new_var[1] = rnorm(priors->fe_precision_wf, v2, seed);

	/*Accept or Reject sigma_ma, then accept or reject for sigma_mw*/

	/*t and j are both used to test what happens if we reverse the
	order of sma and smw*/

	/*Calculate ratio of Cauchy priors*/
	/*      prior_new = new_var[j]/hyper->prec[j];
		  prior_old = priors->fe_precision[j]/hyper->prec[j];
		  prior_new *= prior_new;
		  prior_old *= prior_old;
		  prior_new += 1.0;
		  prior_old += 1.0;
		  prior_ratio = prior_old/prior_new;
		  */
	/*Compute the sum included in the "likelihood"*/
	/*Also compute the old value divided by the new value raised to num_node*/
	/*This sum is the same assuming the current value of sigma_j in both*/
	/*the numerator and denominator of rho*/
	psum[0] = psum[1] = 0;
	subject = sublist->succ;
	prior_ratio[0] = prior_ratio[1] = 1;
	while (subject != NULL){
		prior_ratio[0] *= priors->fe_precision_wl / new_var[0];
		pcomp[0] = subject->theta_l[1] - priors->fe_mean_l[1];
		pcomp[0] *= pcomp[0];
		psum[0] += pcomp[0];

		prior_ratio[1] *= priors->fe_precision_wf / new_var[1];
		pcomp[1] = subject->theta_f[1] - priors->fe_mean_f[1];
		pcomp[1] *= pcomp[1];
		psum[1] += pcomp[1];

		subject = subject->succ;
	}

	/*Log prior ratio value*/
	prior_ratio[0] = log(prior_ratio[0]);
	prior_ratio[1] = log(prior_ratio[1]);


	/*Complete the "likelihood" portion of rho*/
	prop_old[1] = 1.0 / priors->fe_precision_wf;
	prop_old[1] /= priors->fe_precision_wf;
	prop_old[0] = 1.0 / priors->fe_precision_wl;
	prop_old[0] /= priors->fe_precision_wl;


	prop_new[0] = 1.0 / new_var[0];
	prop_new[0] /= new_var[0];

	prop_new[1] = 1.0 / new_var[1];
	prop_new[1] /= new_var[1];
	prop_ratio[0] = psum[0] * 0.5*(prop_old[0] - prop_new[0]);
	prop_ratio[1] = psum[1] * 0.5*(prop_old[1] - prop_new[1]);


	/*We only can accept the proposed value if it is positive*/
	if (new_var[0] > 0 && new_var[0] < hyper->precision_wl){

		/*Compute log rho, and set alpha equal to min(log rho,0)*/
		alpha = (0 < (temp = (prior_ratio[0] + prop_ratio[0]))) ? 0 : temp;

		/*If log(U) < log rho, accept the proposed value*/
		/*Increase acceptance count by 1*/
		if (log(kiss(seed))<alpha){
			tmp1[0]++;
			priors->fe_precision_wl = new_var[0];
		}
	}
	if (new_var[1]>0 && new_var[1] < hyper->precision_wf){

		/*Compute log rho, and set alpha equal to min(log rho,0)*/
		alpha = (0 < (temp = (prior_ratio[1] + prop_ratio[1]))) ? 0 : temp;

		/*If log(U) < log rho, accept the proposed value*/
		/*Increase acceptance count by 1*/
		if (log(kiss(seed)) < alpha){
			tmp1[1]++;
			priors->fe_precision_wf = new_var[1];
		}
	}

	/*Set acceptance count equal to temp vector components*/
	afewv_l = tmp1[0];
	afewv_f = tmp1[1];

	free(new_var);
	free(tmp1);

}

/*********************************************************************/
/*START OF draw_bh_mean SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_bh_mean: this runs the Gibbs sampler draw for mean baseline and half life;
	ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
	Priors *priors; the current values of the prior parameters;
	Common_parms *parms; the current values of the common parameters;
	unsigned long *seed; seed values needed for the randon number generator;
	Hyper_priors *hyper; parameters of the prior distributions;
	RETURNS: None; all updates are made internally
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 j: generic counter
 fesum: sum of subject specific mean baselines and half-lifes;
 part of evaluation of gmean and gvar
 gmean: mean of the distribution used in Gibbs sampler
 gvar: variance of distribution used in Gibbs sampler
 *subject: current list of subjects and their characteristics

 SUBROUTINES USED
 rnorm: found in randgen.h; draws from the normal distribution
 kiss: found in randgen.h; draws from U(0,1) distribution
**************************************************************************/

void draw_bh_mean_l(Subject_type *sublist, Priors *priors, Common_parms *parms, unsigned long *seed, Hyper_priors *hyper)
{
	/* declare variables */
	int j;
	double fesum, gmean, gvar;
	Subject_type *subject;

	/*declare functions */
	double rnorm(double, double, unsigned long *);
	double kiss(unsigned long *);

	for (j = 0; j < 2; j++){

		fesum = 0;
		subject = sublist->succ;

		/*This while loop goes through the pulses to get the sum involved*/
		while (subject != NULL){
			fesum += subject->basehalf_l[j];
			subject = subject->succ;
		} /*end of loop through pulses*/

		gmean = (hyper->meanvarbh_l[j] * fesum) + (priors->varbh_l[j] * priors->varbh_l[j] * hyper->meanmeanbh_l[j]);
		gmean /= (parms->numsub*hyper->meanvarbh_l[j]) + (priors->varbh_l[j] * priors->varbh_l[j]);

		gvar = hyper->meanvarbh_l[j] * priors->varbh_l[j] * priors->varbh_l[j];
		gvar /= (parms->numsub*hyper->meanvarbh_l[j]) + (priors->varbh_l[j] * priors->varbh_l[j]);

		priors->meanbh_l[j] = rnorm(gmean, sqrt(gvar), seed);
	}

}
void draw_bh_mean_f(Subject_type *sublist, Priors *priors, Common_parms *parms, unsigned long *seed, Hyper_priors *hyper)
{
	/* declare variables */
	int j;
	double fesum, gmean, gvar;
	Subject_type *subject;

	/*declare functions */
	double rnorm(double, double, unsigned long *);
	double kiss(unsigned long *);

	for (j = 0; j < 2; j++){

		fesum = 0;
		subject = sublist->succ;

		/*This while loop goes through the pulses to get the sum involved*/
		while (subject != NULL){
			fesum += subject->basehalf_f[j];
			subject = subject->succ;
		} /*end of loop through pulses*/

		gmean = (hyper->meanvarbh_f[j] * fesum) + (priors->varbh_f[j] * priors->varbh_f[j] * hyper->meanmeanbh_f[j]);
		gmean /= (parms->numsub*hyper->meanvarbh_f[j]) + (priors->varbh_f[j] * priors->varbh_f[j]);

		gvar = hyper->meanvarbh_f[j] * priors->varbh_f[j] * priors->varbh_f[j];
		gvar /= (parms->numsub*hyper->meanvarbh_f[j]) + (priors->varbh_f[j] * priors->varbh_f[j]);

		priors->meanbh_f[j] = rnorm(gmean, sqrt(gvar), seed);
	}

}
/*********************************************************************/
/*START OF draw_bh_var SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_bh_var: this runs the M-H draw for standard deviations of subject
				   level baselines and half-lifes;
				   ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
				   Priors *priors; the current values of the prior parameters;
				   double v1; the proposal variance for overall st-dev of baseline;
				   double v2; the proposal variance for overall st-dev of half-life;
				   unsigned long *seed; seed values needed for the randon number generator;
				   Hyper_priors *hyper; parameters of the prior distributions;
				   RETURNS: None; all updates are made internally; this does update the global
				   variables regarding acceptance rate: afebv, afehv, nfebv, and nfehv
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 j: generic counter
 *tmp1: current acceptance counters are saved to this vector so we can increment
 inside a loop
 *new_var: proposed values of st. dev. of baseline and half-life
 prop_new: part of evaluation of prop_ratio portion of log(rho)
 prop_old: same as prop_new but with current (not proposed) value
 psum: part of evaluation of prop_ratio portion of log(rho)
 pcomp: part of evaluation of prop_ratio portion of log(rho)
 prior_old: part of evaluation of prior_ratio portion of log(rho)
 prior_new: same as prior_old, but with proposed value
 prior_ratio: comparison between prior distributions under current and proposed values
 prop_ratio: comparison between random effects distributions under current and
 proposed values
 temp: value used in evaluation of log(rho)
 alpha: minimum of 0 and log(rho); used to determine acceptance in M-H
 *subject: current list of subjects and their qualities

 SUBROUTINES USED
 kiss: found in randgen.h; draws from U(0,1) distribution
 rnorm: found in randgen.h; draws from the normal distribution
**************************************************************************/

void draw_bh_var_l(Subject_type *sublist, Priors *priors, double v1, double v2, unsigned long *seed, Hyper_priors *hyper)
{
	int j, *tmp1;
	double *new_var, prop_new, prop_old, psum, pcomp, prior_old, prior_new;
	double prior_ratio, prop_ratio, temp, alpha;
	Subject_type *subject;
	double kiss(unsigned long *);
	double rnorm(double, double, unsigned long *);

	/*Allocate Memory*/
	new_var = (double *)calloc(2, sizeof(double));
	tmp1 = (int *)calloc(2, sizeof(int));

	/*Add 1 to the counters for acceptance rates of sigma_b and sigma_h*/
	nfebv_l++;
	nfehv_l++;

	/*Assign current acceptance counts to temporary vector*/
	tmp1[0] = afebv_l;
	tmp1[1] = afehv_l;

	/*Draw proposed values for sigma_b and sigma_h*/
	new_var[0] = rnorm(priors->varbh_l[0], v1, seed);
	new_var[1] = rnorm(priors->varbh_l[1], v2, seed);

	/*Accept or Reject sigma_b, then accept or reject for sigma_h*/
	for (j = 0; j < 2; j++){
		/*Calculate ratio of Cauchy priors*/
		/*prior_new = new_var[j]/hyper->varbh[j];
		prior_old = priors->varbh[j]/hyper->varbh[j];
		prior_new *= prior_new;
		prior_old *= prior_old;
		prior_new += 1.0;
		prior_old += 1.0;
		prior_ratio = prior_old/prior_new;
		*/
		/*Compute the sum included in the "likelihood"*/
		/*Also compute the old value divided by the new value raised to num_node*/
		/*This sum is the same assuming the current value of sigma_j in both*/
		/*the numerator and denominator of rho*/
		psum = 0;
		subject = sublist->succ;
		prior_ratio = 1;
		while (subject != NULL){
			prior_ratio *= priors->varbh_l[j] / new_var[j];
			pcomp = subject->basehalf_l[j] - priors->meanbh_l[j];
			pcomp *= pcomp;
			psum += pcomp;
			subject = subject->succ;
		}

		/*Log prior ratio value*/
		prior_ratio = log(prior_ratio);

		/*Complete the "likelihood" portion of rho*/
		prop_old = 1.0 / priors->varbh_l[j];
		prop_old /= priors->varbh_l[j];
		prop_new = 1.0 / new_var[j];
		prop_new /= new_var[j];
		prop_ratio = psum*0.5*(prop_old - prop_new);

		/*We only can accept the proposed value if it is positive*/
		/*With a uniform(0,2) prior, skip this step if we have a sd less than 2*/
		if (new_var[j] > 0 && new_var[j] < hyper->varbh_l[j]){

			/*Compute log rho, and set alpha equal to min(log rho,0)*/
			/*With a uniform prior, get rid of prior_ratio*/
			alpha = (0 < (temp = (prior_ratio + prop_ratio))) ? 0 : temp;

			/*If log(U) < log rho, accept the proposed value*/
			/*Increase acceptance count by 1*/
			if (log(kiss(seed)) < alpha){
				tmp1[j]++;
				priors->varbh_l[j] = new_var[j];
			}
		}
	}

	/*Set acceptance count equal to temp vector components*/
	afebv_l = tmp1[0];
	afehv_l = tmp1[1];

	free(new_var);
	free(tmp1);

}

void draw_bh_var_f(Subject_type *sublist, Priors *priors, double v1, double v2, unsigned long *seed, Hyper_priors *hyper)
{
	int j, *tmp1;
	double *new_var, prop_new, prop_old, psum, pcomp, prior_old, prior_new;
	double prior_ratio, prop_ratio, temp, alpha;
	Subject_type *subject;
	double kiss(unsigned long *);
	double rnorm(double, double, unsigned long *);

	/*Allocate Memory*/
	new_var = (double *)calloc(2, sizeof(double));
	tmp1 = (int *)calloc(2, sizeof(int));

	/*Add 1 to the counters for acceptance rates of sigma_b and sigma_h*/
	nfebv_f++;
	nfehv_f++;

	/*Assign current acceptance counts to temporary vector*/
	tmp1[0] = afebv_f;
	tmp1[1] = afehv_f;

	/*Draw proposed values for sigma_b and sigma_h*/
	new_var[0] = rnorm(priors->varbh_f[0], v1, seed);
	new_var[1] = rnorm(priors->varbh_f[1], v2, seed);

	/*Accept or Reject sigma_b, then accept or reject for sigma_h*/
	for (j = 0; j < 2; j++){
		/*Calculate ratio of Cauchy priors*/
		/*prior_new = new_var[j]/hyper->varbh[j];
		prior_old = priors->varbh[j]/hyper->varbh[j];
		prior_new *= prior_new;
		prior_old *= prior_old;
		prior_new += 1.0;
		prior_old += 1.0;
		prior_ratio = prior_old/prior_new;
		*/
		/*Compute the sum included in the "likelihood"*/
		/*Also compute the old value divided by the new value raised to num_node*/
		/*This sum is the same assuming the current value of sigma_j in both*/
		/*the numerator and denominator of rho*/
		psum = 0;
		subject = sublist->succ;
		prior_ratio = 1;
		while (subject != NULL){
			prior_ratio *= priors->varbh_f[j] / new_var[j];
			pcomp = subject->basehalf_f[j] - priors->meanbh_f[j];
			pcomp *= pcomp;
			psum += pcomp;
			subject = subject->succ;
		}

		/*Log prior ratio value*/
		prior_ratio = log(prior_ratio);

		/*Complete the "likelihood" portion of rho*/
		prop_old = 1.0 / priors->varbh_f[j];
		prop_old /= priors->varbh_f[j];
		prop_new = 1.0 / new_var[j];
		prop_new /= new_var[j];
		prop_ratio = psum*0.5*(prop_old - prop_new);

		/*We only can accept the proposed value if it is positive*/
		/*With a uniform(0,2) prior, skip this step if we have a sd less than 2*/
		if (new_var[j] > 0 && new_var[j] < hyper->varbh_f[j]){

			/*Compute log rho, and set alpha equal to min(log rho,0)*/
			/*With a uniform prior, get rid of prior_ratio*/
			alpha = (0 < (temp = (prior_ratio + prop_ratio))) ? 0 : temp;

			/*If log(U) < log rho, accept the proposed value*/
			/*Increase acceptance count by 1*/
			if (log(kiss(seed)) < alpha){
				tmp1[j]++;
				priors->varbh_f[j] = new_var[j];
			}
		}
	}

	/*Set acceptance count equal to temp vector components*/
	afebv_f = tmp1[0];
	afehv_f = tmp1[1];

	free(new_var);
	free(tmp1);

}



void draw_fixed_mass(Subject_type *sublist, Priors *priors, Common_parms *parms_l, Common_parms *parms_f, unsigned long *seed, double **pmean_var)
{
	int i, j, k, numnode_l, numnode_f, tmp1[2];
	double diff_mean_old[2], diff_mean_new[2], log_prior_ratio, temp_old_mean[2], theta[2], temp, alpha, *new_mean, *tmpold, *tmpnew;
	double A_inv[2][2], **re_var_inv, **re_var, stdnew[2], stdold[2], log_RE_ratio, new_sum[2], old_sum[2], **pvartmp;
	double newint, oldint, junk, iden[2][2], stdxold, stdyold, stdxnew, stdynew, corrxy, log_prop_ratio;


	Subject_type *subject;
	Node_type *subnode_f;
	Node_type *subnode_l;
	double MDBNOR(double, double, double);
	int rmvnorm(double *, double **, int, double *, unsigned long *, int);
	double rgamma(double, double, unsigned long *);
	double kiss(unsigned long *);
	void print_list(Node_type *);
	int cholesky_decomp(double **, int);
	double **cholesky_invert(int, double **);
	double Phi(double y);

	double kiss(unsigned long *);
	double rnorm(double, double, unsigned long *);
	tmpold = (double *)calloc(2, sizeof(double));
	tmpnew = (double *)calloc(2, sizeof(double));
	new_mean = (double *)calloc(2, sizeof(double));
	re_var_inv = (double **)calloc(2, sizeof(double *));
	re_var = (double **)calloc(2, sizeof(double *));
	pvartmp = (double **)calloc(2, sizeof(double *));
	/*for (i = 0; i < 2; i++) {
		re_var_inv[i] = (double *)calloc(2, sizeof(double));
		re_var[i] = (double *)calloc(2, sizeof(double));
		pvartmp[i] = (double *)calloc(2, sizeof(double));
	}*/
	/*    printf("iter %d \n",iter);*/




	/*re_var = cholesky_invert(2,priors->fe_precision);*/

	/*Go to start of subject list*/
	subject = sublist->succ;

	/*Go through each existing subject*/
	while (subject != NULL) {


		nfepm++;

		tmp1[0] = afepm;
		/*draw the new pulse pair*/
		theta[0] = subject->theta_l[0];
		theta[1] = subject->theta_f[0];
		old_sum[0] = old_sum[1] = new_sum[0] = new_sum[1] = 0.0;

		/*Go to start of pulse list for this subject*/
		rmvnorm(new_mean, pmean_var, 2, theta, seed, 1);
		if (new_mean[0] > 0 && new_mean[1] > 0) {

			/*calculate the prior ratio*/
			log_prior_ratio = 0;
			diff_mean_old[0] = subject->theta_l[0] - priors->fe_mean_l[0];
			diff_mean_old[1] = subject->theta_f[0] - priors->fe_mean_f[0];
			diff_mean_new[0] = new_mean[0] - priors->fe_mean_l[0];
			diff_mean_new[1] = new_mean[1] - priors->fe_mean_f[0];

			for (i = 0; i < 2; i++)
			for (j = 0; j < 2; j++)
				log_prior_ratio += diff_mean_old[i] * priors->fe_precision[i][j] * diff_mean_old[j] - diff_mean_new[i] * priors->fe_precision[i][j] * diff_mean_new[j];

			log_prior_ratio *= 0.5;

			/*calculate the likelihood ratio*/
			/*first the driver's pulse*/
			subnode_l = subject->driver->succ;
			numnode_l = 0;
			newint = 0.;
			oldint = 0.;

			/*Go through each existing pulse and add up log amplitude/width
			 Also count pulses*/
			while (subnode_l != NULL){

				numnode_l++;
				old_sum[0] += (subnode_l->theta[0] - subject->theta_l[0])*(subnode_l->theta[0] - subject->theta_l[0])*subnode_l->eta[0];
				new_sum[0] += (subnode_l->theta[0] - new_mean[0])*(subnode_l->theta[0] - new_mean[0])*subnode_l->eta[0];
				stdxnew = (new_mean[0]) / (parms_l->re_precision[0] / sqrt(subnode_l->eta[0]));
				stdxold = (subject->theta_l[0]) / (parms_l->re_precision[0] / sqrt(subnode_l->eta[0]));
				newint += log(Phi(stdxnew));
				oldint += log(Phi(stdxold));

				subnode_l = subnode_l->succ;

			}
			log_RE_ratio = 0.5/(parms_l->re_precision[0]* parms_l->re_precision[0])*(old_sum[0] - new_sum[0]) + oldint - newint;
			/*start of the fsh loop*/
			subnode_f = subject->response->succ;
			numnode_f = 0;
			newint = 0.;
			oldint = 0.;
			/*Go through each existing pulse and add up log amplitude/width
			 Also count pulses*/
			while (subnode_f != NULL){

				numnode_f++;
				old_sum[1] += (subnode_f->theta[0] - subject->theta_f[0])*(subnode_f->theta[0] - subject->theta_f[0])*subnode_f->eta[0];
				new_sum[1] += (subnode_f->theta[0] - new_mean[1])*(subnode_f->theta[0] - new_mean[1])*subnode_f->eta[0];
				stdxnew = (new_mean[1]) / (parms_f->re_precision[0] / sqrt(subnode_f->eta[0]));
				stdxold = (subject->theta_f[0]) / (parms_f->re_precision[0] /sqrt(subnode_f->eta[0]));
				newint += log(Phi(stdxnew));
				oldint += log(Phi(stdxold));
				subnode_f = subnode_f->succ;


				log_RE_ratio += 0.5/ (parms_f->re_precision[0]* parms_f->re_precision[0])*(old_sum[1] - new_sum[1]) + oldint - newint;

			}

			alpha = (0 < (temp = (log_prior_ratio + log_RE_ratio))) ? 0 : temp;
			/*If log(U) < log rho, accept the proposed value*/
			/*Increase acceptance count by 1*/
			if (log(kiss(seed)) < alpha){
				tmp1[0]++;
				subject->theta_l[0] = new_mean[0];
				subject->theta_f[0] = new_mean[1];

			}
		}
		afepm = tmp1[0];
		/*Go through each existing pulse and add up log amplitude/width
		 Also count pulses*/


		/*Advance to next subject*/
		subject = subject->succ;

	} /*end of loop through subjects*/
	free(re_var_inv);
	free(re_var);
	free(new_mean);
	free(tmpnew);
	free(tmpold);
}


/*********************************************************************/
/*START OF draw_fe_precision SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_fe_precision: this runs the M-H draw for standard deviations of pulse
					 masses and widths;
					 ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
					 Priors *priors; the current values of the prior parameters;
					 Common_parms *parms; the current values of the common parameters;
					 double v1; the proposal variance for overall st-dev of pulse mass;
					 double v2; the proposal variance for overall st-dev of pulse width;
					 unsigned long *seed; seed values needed for the randon number generator;
					 RETURNS: None; all updates are made internally; this does update the global
					 variables regarding acceptance rate: afemv, afewv, nfemv, and nfewv
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 j: generic counter
 *tmp1: current acceptance counters are saved to this vector so we can increment
 inside a loop
 *new_var: proposed values of st. dev. of pulse mass and width
 prop_new: part of evaluation of prop_ratio portion of log(rho)
 prop_old: same as prop_new but with current (not proposed) value
 psum: part of evaluation of prop_ratio portion of log(rho)
 pcomp: part of evaluation of prop_ratio portion of log(rho)
 prior_old: part of evaluation of prior_ratio portion of log(rho)
 prior_new: same as prior_old, but with proposed value
 prior_ratio: comparison between prior distributions under current and proposed values
 prop_ratio: comparison between random effects distributions under current and
 proposed values
 temp: value used in evaluation of log(rho)
 alpha: minimum of 0 and log(rho); used to determine acceptance in M-H
 *subnode: for a given subject, list of pulses and their characteristics
 *subject: current list of subjects and their qualities

 SUBROUTINES USED
 kiss: found in randgen.h; draws from U(0,1) distribution
 rnorm: found in randgen.h; draws from the normal distribution
**************************************************************************/

void draw_fixed_width_l(Subject_type *sublist, Priors *priors, Common_parms *parms, double v, unsigned long *seed)
{
	int j, tmp1;
	double new_mean, prop_new, prop_old, psum_old, psum_new, pcomp, old_prior, new_prior;
	double prior_ratio, prop_ratio, temp, alpha;
	double stdxold, stdxnew, newint, oldint, prior_old, prior_new;
	Node_type *subnode;
	Subject_type *subject;
	double kiss(unsigned long *);
	double Phi(double);
	double rnorm(double, double, unsigned long *);

	/*Allocate Memory*/

	/*Add 1 to the counters for acceptance rates of sigma_a and sigma_w*/
	nfepw_l++;


	/*Assign current acceptance counts to temporary vector*/
	tmp1 = afepw_l;




	/* do the simulation for each subject */

	subject = sublist->succ;

	while (subject != NULL){
		new_mean = rnorm(subject->theta_l[1], v, seed);
		if (new_mean > 0){
			old_prior = (subject->theta_l[1] - priors->fe_mean_l[1])*(subject->theta_l[1] - priors->fe_mean_l[1]);
			old_prior /= priors->fe_precision_wl;
			old_prior /= priors->fe_precision_wl;
			new_prior = (new_mean - priors->fe_mean_l[1])*(new_mean - priors->fe_mean_l[1]);
			new_prior /= priors->fe_precision_wl;
			new_prior /= priors->fe_precision_wl;

			prior_ratio = 0.5*(old_prior - new_prior);
			psum_old = 0.0;
			psum_new = 0.0;
			newint = 0.0;
			oldint = 0.0;
			subnode = subject->driver->succ;
			while (subnode != NULL){
				psum_old += (subnode->theta[1] - subject->theta_l[1])*(subnode->theta[1] - subject->theta_l[1])*subnode->eta[1];
				psum_new += (subnode->theta[1] - new_mean)*(subnode->theta[1] - new_mean)*subnode->eta[1];
				stdxnew = (new_mean) / (parms->re_precision[1] / sqrt(subnode->eta[1]));
				stdxold = (subject->theta_l[1]) / (parms->re_precision[1] / sqrt(subnode->eta[1]));
				newint += log(Phi(stdxnew));
				oldint += log(Phi(stdxold));

				subnode = subnode->succ;
			}

			prop_ratio = (psum_old - psum_new)*0.5 / (parms->re_precision[1] * parms->re_precision[1]) - newint + oldint;

			alpha = (0 < (temp = (prior_ratio + prop_ratio))) ? 0 : temp;
			/*If log(U) < log rho, accept the proposed value*/
			/*Increase acceptance count by 1*/
			if (log(kiss(seed)) < alpha){
				tmp1++;
				subject->theta_l[1] = new_mean;



			}
		}

		subject = subject->succ;
	}





	/*end of loop through s_a and s_w */

	/*Set acceptance count equal to temp vector components*/

	afewv_l = tmp1;



}
void draw_fixed_width_f(Subject_type *sublist, Priors *prior, Common_parms *parms, double v, unsigned long *seed)
{
	int j, tmp1;
	double  new_mean, prop_new, prop_old, psum_old, psum_new, pcomp, old_prior, new_prior;
	double prior_ratio, prop_ratio, temp, alpha;
	double stdxold, stdxnew, newint, oldint, prior_old, prior_new;
	Node_type *subnode;
	Subject_type *subject;
	double kiss(unsigned long *);
	double Phi(double);
	double rnorm(double, double, unsigned long *);

	/*Allocate Memory*/

	/*Add 1 to the counters for acceptance rates of sigma_a and sigma_w*/
	nfepw_f++;


	/*Assign current acceptance counts to temporary vector*/
	tmp1 = afepw_f;




	/* do the simulation for each subject */

	subject = sublist->succ;

	while (subject != NULL){
		new_mean = rnorm(subject->theta_f[1], v, seed);
		if (new_mean > 0){
			old_prior = (subject->theta_f[1] - prior->fe_mean_f[1])*(subject->theta_f[1] - prior->fe_mean_f[1]);
			old_prior /= prior->fe_precision_wf;
			old_prior /= prior->fe_precision_wf;

			new_prior = (new_mean - prior->fe_mean_f[1])*(new_mean - prior->fe_mean_f[1]);
			new_prior /= prior->fe_precision_wf;
			new_prior /= prior->fe_precision_wf;

			prior_ratio = 0.5*(old_prior - new_prior);
			newint = 0;
			oldint = 0;
			psum_old = 0.0;
			psum_new = 0.0;
			subnode = subject->response->succ;
			while (subnode != NULL){
				psum_old += (subnode->theta[1] - subject->theta_f[1])*(subnode->theta[1] - subject->theta_f[1])*subnode->eta[1];
				psum_new += (subnode->theta[1] - new_mean)*(subnode->theta[1] - new_mean)*subnode->eta[1];
				stdxnew = (new_mean) / (parms->re_precision[1] / sqrt(subnode->eta[1]));
				stdxold = (subject->theta_f[1]) / (parms->re_precision[1] / sqrt(subnode->eta[1]));
				newint += log(Phi(stdxnew));
				oldint += log(Phi(stdxold));


				subnode = subnode->succ;
			}

			prop_ratio = (psum_old - psum_new)*0.5 / (parms->re_precision[1] * parms->re_precision[1]) - newint + oldint;

			alpha = (0 < (temp = (prior_ratio + prop_ratio))) ? 0 : temp;
			/*If log(U) < log rho, accept the proposed value*/
			/*Increase acceptance count by 1*/
			if (log(kiss(seed)) < alpha){
				tmp1++;
				subject->theta_f[1] = new_mean;



			}
		}

		subject = subject->succ;
	}





	/*end of loop through s_a and s_w */

	/*Set acceptance count equal to temp vector components*/

	afewv_f = tmp1;



}
void draw_fe_precision_l(Subject_type *sublist, Priors *priors, Common_parms *parms, double v1, double v2, unsigned long *seed)
{
	int *tmp1, j;
	double *new_var, prop_new, prop_old, psum, pcomp, prior_old, prior_new;
	double prior_ratio, prop_ratio, temp, alpha;
	double stdold[2], stdnew[2], newint[2], oldint[2];
	Node_type *subnode;
	Subject_type *subject;
	double kiss(unsigned long *);
	double Phi(double);
	double rnorm(double, double, unsigned long *);

	/*Allocate Memory*/
	new_var = (double *)calloc(2, sizeof(double));
	tmp1 = (int *)calloc(2, sizeof(int));

	/*Add 1 to the counters for acceptance rates of sigma_a and sigma_w*/
	nfepmv_l++;
	nfepwv_l++;

	/*Assign current acceptance counts to temporary vector*/
	tmp1[0] = afepmv_l;
	tmp1[1] = afepwv_l;

	/*Draw proposed values for sigma_a and sigma_w*/
	new_var[0] = rnorm(parms->re_precision[0], v1, seed);
	new_var[1] = rnorm(parms->re_precision[1], v2, seed);

	/*Accept or Reject sigma_a, then accept or reject for sigma_w*/
	for (j = 0; j < 2; j++){

		/*Calculate ratio of Cauchy priors*/
		/*      prior_new = new_var[j]/priors->re_var[j];
			  prior_old = parms->re_precision[j]/priors->re_var[j];
			  prior_new *= prior_new;
			  prior_old *= prior_old;
			  prior_new += 1.0;
			  prior_old += 1.0;
			  prior_ratio = prior_old/prior_new;
			  */
		/*Compute the sum included in the "likelihood"*/
		/*Also compute the old value divided by the new value raised to num_node*/
		/*This sum is the same assuming the current value of sigma_j in both*/
		/*the numerator and denominator of rho*/
		psum = 0;
		prior_ratio = 1;
		newint[0] = 0;
		oldint[0] = 0;
		newint[1] = 0;
		oldint[1] = 0;

		subject = sublist->succ;
		while (subject != NULL){
			subnode = subject->driver->succ;
			while (subnode != NULL){
				prior_ratio *= parms->re_precision[j] / new_var[j];
				pcomp = subnode->theta[j] - subject->theta_l[j];
				pcomp *= pcomp;
				psum += pcomp*subnode->eta[j];
				stdold[j] = (subnode->theta[j]) / (parms->re_precision[j] / sqrt(subnode->eta[j]));
				stdnew[j] = (subnode->theta[j]) / (new_var[j] / sqrt(subnode->eta[j]));

				newint[j] += log(Phi(stdnew[j]));
				oldint[j] += log(Phi(stdold[j]));


				subnode = subnode->succ;
			}
			subject = subject->succ;
		}

		/*Log prior ratio value*/
		prior_ratio = log(prior_ratio) - newint[j] + oldint[j];

		/*Complete the "expo " portion of likelihood ratio*/
		prop_old = 1.0 / parms->re_precision[j];
		prop_old /= parms->re_precision[j];
		prop_new = 1.0 / new_var[j];
		prop_new /= new_var[j];
		prop_ratio = psum*0.5*(prop_old - prop_new);

		/*We only can accept the proposed value if it is positive*/
		if (new_var[j] > 0.1 && new_var[j] < priors->re_var_l[j]){

			/*Compute log rho, and set alpha equal to min(log rho,0)*/
			alpha = (0 < (temp = (prior_ratio + prop_ratio))) ? 0 : temp;
			/*If log(U) < log rho, accept the proposed value*/
			/*Increase acceptance count by 1*/
			if (log(kiss(seed)) < alpha){
				tmp1[j]++;
				parms->re_precision[j] = new_var[j];
			}
		}
	} /*end of loop through s_a and s_w */

	/*Set acceptance count equal to temp vector components*/
	afepmv_l = tmp1[0];
	afepwv_l = tmp1[1];

	free(new_var);
	free(tmp1);

}


void draw_fe_precision_f(Subject_type *sublist, Priors *priors, Common_parms *parms, double v1, double v2, unsigned long *seed)
{
	int *tmp1, j;
	double *new_var, prop_new, prop_old, psum, pcomp, prior_old, prior_new;
	double prior_ratio, prop_ratio, temp, alpha;
	double stdold[2], stdnew[2], newint[2], oldint[2];
	Node_type *subnode;
	Subject_type *subject;
	double kiss(unsigned long *);
	double Phi(double);
	double rnorm(double, double, unsigned long *);

	/*Allocate Memory*/
	new_var = (double *)calloc(2, sizeof(double));
	tmp1 = (int *)calloc(2, sizeof(int));

	/*Add 1 to the counters for acceptance rates of sigma_a and sigma_w*/
	nfepmv_f++;
	nfepwv_f++;

	/*Assign current acceptance counts to temporary vector*/
	tmp1[0] = afepmv_f;
	tmp1[1] = afepwv_f;

	/*Draw proposed values for sigma_a and sigma_w*/
	new_var[0] = rnorm(parms->re_precision[0], v1, seed);
	new_var[1] = rnorm(parms->re_precision[1], v2, seed);


	/*Accept or Reject sigma_a, then accept or reject for sigma_w*/
	for (j = 0; j < 2; j++){

		/*Calculate ratio of Cauchy priors*/
		/*      prior_new = new_var[j]/priors->re_var[j];
			  prior_old = parms->re_precision[j]/priors->re_var[j];
			  prior_new *= prior_new;
			  prior_old *= prior_old;
			  prior_new += 1.0;
			  prior_old += 1.0;
			  prior_ratio = prior_old/prior_new;
			  */
		/*Compute the sum included in the "likelihood"*/
		/*Also compute the old value divided by the new value raised to num_node*/
		/*This sum is the same assuming the current value of sigma_j in both*/
		/*the numerator and denominator of rho*/
		psum = 0;
		prior_ratio = 1;
		newint[0] = 0;
		newint[1] = 0;
		oldint[0] = 0;
		oldint[1] = 0;

		subject = sublist->succ;
		while (subject != NULL){
			subnode = subject->response->succ;
			while (subnode != NULL){
				prior_ratio *= parms->re_precision[j] / new_var[j];
				pcomp = subnode->theta[j] - subject->theta_f[j];
				pcomp *= pcomp;
				psum += pcomp*subnode->eta[j];
				stdold[j] = (subnode->theta[j]) / (parms->re_precision[j] / sqrt(subnode->eta[j]));
				stdnew[j] = (subnode->theta[j]) / (new_var[j] / sqrt(subnode->eta[j]));

				newint[j] += log(Phi(stdnew[j]));
				oldint[j] += log(Phi(stdold[j]));


				subnode = subnode->succ;
			}
			subject = subject->succ;
		}

		/*Log prior ratio value*/
		prior_ratio = log(prior_ratio) - newint[j] + oldint[j];

		/*Complete the "expo " portion of likelihood ratio*/
		prop_old = 1.0 / parms->re_precision[j];
		prop_old /= parms->re_precision[j];
		prop_new = 1.0 / new_var[j];
		prop_new /= new_var[j];
		prop_ratio = psum*0.5*(prop_old - prop_new);

		/*We only can accept the proposed value if it is positive*/
		if (new_var[j] > 0.1 && new_var[j] < priors->re_var_f[j]){

			/*Compute log rho, and set alpha equal to min(log rho,0)*/
			alpha = (0 < (temp = (prior_ratio + prop_ratio))) ? 0 : temp;
			/*If log(U) < log rho, accept the proposed value*/
			/*Increase acceptance count by 1*/
			if (log(kiss(seed)) < alpha){
				tmp1[j]++;
				parms->re_precision[j] = new_var[j];
			}
		}
	} /*end of loop through s_a and s_w */

	/*Set acceptance count equal to temp vector components*/
	afepmv_f = tmp1[0];
	afepwv_f = tmp1[1];

	free(new_var);
	free(tmp1);

}
/*********************************************************************/
/*START OF draw_bh SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_bh: this runs the M-H draw for baseline and halflife (drawn together);
	ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
	Common_parms *parms; the current values of the common parameters;
	Priors *priors; the current values of the prior parameters;
	double **ts; this is the matrix of observed data (a column of
	times and S columns of log(concentration);
	int N; the number of observations in each column of **ts;
	unsigned long *seed; seed values needed for the randon number generator;
	double **var; the proposal variance-covariance matrix for
	baseline and halflife
	RETURNS: None; all updates are made internally; this does update the global
	variables regarding acceptance rate: adelta and ndelta
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i,j,k: generic counters
 num_node: number of pulses
 *subnode: for a given subject, current list of pulses and their qualities
 *new_node: counter for number of pulses
 *subject: current list of subjects and their characteristics
 alpha: minimum of 0 and log(rho); used to determine acceptance in M-H
 plikelihood: likelihood under proposed value
 *pmd: vector of proposed baseline and halflife
 like_ratio: difference between current and proposed log-likelihoods
 priorb_old: part of evaluation of prior_ratio portion of log(rho)
 priorb_new: same as priorb_old, but with proposed baseline
 priorh_old: part of evaluation of prior_ratio portion of log(rho)
 priorh_new: same as priorb_old, but with proposed halflife
 prior_ratio: comparison between prior distributions under current and proposed values
 currentmd[2]: current baseline and halflife saved to this variable name so we can
 evaluate likelihood under proposed baseline and halflife
 **currentmc: current mean_contrib for each pulse is saved to this matrix
 so we can evaluate likelihood under proposed pulses' mean_contrib
 temp: value used in evaluation of log(rho)
 current_decay: current decay rate is saved to this variable name so we can
 evaluate likelihood under proposed decay rate
 current_like: current likelihood for a given subject; saved so we can compare

 SUBROUTINES USED
 rmvnorm: found in randgen.h; draws from the multivariate normal distribution
 kiss: found in randgen.h; draws from U(0,1) distribh_f
 ution
 mean_contribution: found in birthdeath.c; evaluates mean_contrib for a
 pulse with specific parameters
 likelihood2: found in birthdeath.c; calculates and returns likelihood for a
 given subject under current parameters*/
/**************************************************************************/

void draw_bh_l(Subject_type *sublist, Common_parms *parms, Priors *priors, double **ts,
	int N, unsigned long *seed, double **var)
{
	int i, j, k, num_node;
	Node_type *subnode, *new_node;
	Subject_type *subject;
	double alpha, plikelihood, *pmd, like_ratio, priorb_old, priorb_new, priorh_old;
	double priorh_new, prior_ratio, currentmd[2], **currentmc, temp, current_decay;
	double current_like;
	int rmvnorm(double *, double **, int, double *, unsigned long *, int);
	double kiss(unsigned long *);
	void mean_contribution(Node_type *, double **, Common_parms *, int, double);
	double likelihood2(Node_type *, double **, Common_parms *, int, Node_type *, double);

	/*Go to start of list of subjects*/
	subject = sublist->succ;

	/*this tells likelihood2 which subject's likelihood to calculate*/
	parms->subindex = 0;

	/*Loop thru all subjects*/
	while (subject != NULL){

		/*Calculate the current likelihood for this subject*/
		current_like = likelihood2(subject->driver, ts, parms, N, subject->driver, subject->basehalf_l[0]);

		/*Count number of pulses for this subject*/
		num_node = 0;
		new_node = subject->driver->succ;
		while (new_node != NULL) {
			new_node = new_node->succ;
			num_node++;
		}

		/*Allocate memory*/
		pmd = (double *)calloc(2, sizeof(double));
		currentmc = (double **)calloc(num_node, sizeof(double *));
		for (i = 0; i<num_node; i++)
			currentmc[i] = (double *)calloc(N, sizeof(double));

		/*Increase denominator of acceptance rate for b and hl*/
		ndelta_l++;

		/*Draw proposal values for b and hl*/
		rmvnorm(pmd, var, 2, subject->basehalf_l, seed, 1);

		/*Only proceed if we draw reasonable values*/
		if (pmd[0] > 0 && pmd[1] > 3) {

			/*Compute ratio of prior densities*/
			priorb_old = subject->basehalf_l[0] - priors->meanbh_l[0];
			priorb_old *= 0.5*priorb_old / (priors->varbh_l[0] * priors->varbh_l[0]);
			priorb_new = pmd[0] - priors->meanbh_l[0];
			priorb_new *= 0.5*priorb_new / (priors->varbh_l[0] * priors->varbh_l[0]);
			priorh_old = subject->basehalf_l[1] - priors->meanbh_l[1];
			priorh_old *= 0.5*priorh_old / (priors->varbh_l[1] * priors->varbh_l[1]);
			priorh_new = pmd[1] - priors->meanbh_l[1];
			priorh_new *= 0.5*priorh_new / (priors->varbh_l[1] * priors->varbh_l[1]);

			prior_ratio = priorb_old + priorh_old - priorb_new - priorh_new;

			/*Save current values of b and hl and set new current values equal to
			proposal values */
			for (k = 0; k < 2; k++) {
				currentmd[k] = subject->basehalf_l[k];
				subject->basehalf_l[k] = pmd[k];
			}

			/*Save current decay rate; calculate new current decay rate based on
			  proposal value of hl */
			current_decay = subject->decay_l;
			subject->decay_l = log(2) / subject->basehalf_l[1];

			/*Save current mean_contrib for each pulse; calculate new current
			  mean_contrib for each pulse based on proposed values of b and hl*/
			i = 0;
			subnode = subject->driver->succ;
			while (subnode != NULL) {
				for (j = 0; j < N; j++)
					currentmc[i][j] = subnode->mean_contrib[j];
				mean_contribution(subnode, ts, parms, N, subject->basehalf_l[1]);
				i++;
				subnode = subnode->succ;
			}

			/*Calculate proposed likelihood and then calculate likelihood ratio*/
			plikelihood = likelihood2(subject->driver, ts, parms, N, subject->driver, subject->basehalf_l[0]);
			like_ratio = plikelihood - current_like;

			/*Calculate log rho; set alpha equal to min(0,log rho) */
			alpha = (0 < (temp = (prior_ratio + like_ratio))) ? 0 : temp;

			/*If log U < log rho, increase acceptance rate by 1  */
			if (log(kiss(seed)) < alpha) {
				adelta_l++;
			}

			/*Otherwise, we need to revert values back to their current state */
			else {

				/*Set b and hl back equal to current values*/
				subject->basehalf_l[0] = currentmd[0];
				subject->basehalf_l[1] = currentmd[1];

				/*Set mean_contrib back equal to current values for each pulse*/
				i = 0;
				subnode = subject->driver->succ;
				while (subnode != NULL) {
					for (j = 0; j < N; j++) {
						subnode->mean_contrib[j] = currentmc[i][j];
					}
					i++;
					subnode = subnode->succ;
				}

				/*Set decay rate back equal to current value*/
				subject->decay_l = current_decay;

			} /*End of if else statement*/

		}

		/*Go to next subject*/
		subject = subject->succ;
		parms->subindex++;

		/*Free memory */
		for (i = 0; i < num_node; i++)
			free(currentmc[i]);
		free(currentmc);
		free(pmd);

	} /*End of loop through subjects*/

}

void draw_bh_f(Subject_type *sublist, Common_parms *parms, Priors *priors, double **ts,
	int N, unsigned long *seed, double **var)
{
	int i, j, k, num_node;
	Node_type *subnode, *new_node;
	Subject_type *subject;
	double alpha, plikelihood, *pmd, like_ratio, priorb_old, priorb_new, priorh_old;
	double priorh_new, prior_ratio, currentmd[2], **currentmc, temp, current_decay;
	double current_like;
	int rmvnorm(double *, double **, int, double *, unsigned long *, int);
	double kiss(unsigned long *);
	void mean_contribution(Node_type *, double **, Common_parms *, int, double);
	double likelihood2(Node_type *, double **, Common_parms *, int, Node_type *, double);

	/*Go to start of list of subjects*/
	subject = sublist->succ;

	/*this tells likelihood2 which subject's likelihood to calculate*/
	parms->subindex = 0;

	/*Loop thru all subjects*/
	while (subject != NULL){

		/*Calculate the current likelihood for this subject*/
		current_like = likelihood2(subject->response, ts, parms, N, subject->response, subject->basehalf_f[0]);

		/*Count number of pulses for this subject*/
		num_node = 0;
		new_node = subject->response->succ;
		while (new_node != NULL) {
			new_node = new_node->succ;
			num_node++;
		}

		/*Allocate memory*/
		pmd = (double *)calloc(2, sizeof(double));
		currentmc = (double **)calloc(num_node, sizeof(double *));
		for (i = 0; i<num_node; i++)
			currentmc[i] = (double *)calloc(N, sizeof(double));

		/*Increase denominator of acceptance rate for b and hl*/
		ndelta_f++;

		/*Draw proposal values for b and hl*/
		rmvnorm(pmd, var, 2, subject->basehalf_f, seed, 1);

		/*Only proceed if we draw reasonable values*/
		if (pmd[0] > 0 && pmd[1] > 3) {

			/*Compute ratio of prior densities*/
			priorb_old = subject->basehalf_f[0] - priors->meanbh_f[0];
			priorb_old *= 0.5*priorb_old / (priors->varbh_f[0] * priors->varbh_f[0]);
			priorb_new = pmd[0] - priors->meanbh_f[0];
			priorb_new *= 0.5*priorb_new / (priors->varbh_f[0] * priors->varbh_f[0]);
			priorh_old = subject->basehalf_f[1] - priors->meanbh_f[1];
			priorh_old *= 0.5*priorh_old / (priors->varbh_f[1] * priors->varbh_f[1]);
			priorh_new = pmd[1] - priors->meanbh_f[1];
			priorh_new *= 0.5*priorh_new / (priors->varbh_f[1] * priors->varbh_f[1]);

			prior_ratio = priorb_old + priorh_old - priorb_new - priorh_new;

			/*Save current values of b and hl and set new current values equal to
			proposal values */
			for (k = 0; k < 2; k++) {
				currentmd[k] = subject->basehalf_f[k];
				subject->basehalf_f[k] = pmd[k];
			}

			/*Save current decay rate; calculate new current decay rate based on
			  proposal value of hl */
			current_decay = subject->decay_f;
			subject->decay_f = log(2) / subject->basehalf_f[1];

			/*Save current mean_contrib for each pulse; calculate new current
			  mean_contrib for each pulse based on proposed values of b and hl*/
			i = 0;
			subnode = subject->response->succ;
			while (subnode != NULL) {
				for (j = 0; j < N; j++)
					currentmc[i][j] = subnode->mean_contrib[j];
				mean_contribution(subnode, ts, parms, N, subject->basehalf_f[1]);
				i++;
				subnode = subnode->succ;
			}

			/*Calculate proposed likelihood and then calculate likelihood ratio*/
			plikelihood = likelihood2(subject->response, ts, parms, N, subject->response, subject->basehalf_f[0]);
			like_ratio = plikelihood - current_like;

			/*Calculate log rho; set alpha equal to min(0,log rho) */
			alpha = (0 < (temp = (prior_ratio + like_ratio))) ? 0 : temp;

			/*If log U < log rho, increase acceptance rate by 1  */
			if (log(kiss(seed)) < alpha) {
				adelta_f++;
			}

			/*Otherwise, we need to revert values back to their current state */
			else {

				/*Set b and hl back equal to current values*/
				subject->basehalf_f[0] = currentmd[0];
				subject->basehalf_f[1] = currentmd[1];

				/*Set mean_contrib back equal to current values for each pulse*/
				i = 0;
				subnode = subject->response->succ;
				while (subnode != NULL) {
					for (j = 0; j < N; j++) {
						subnode->mean_contrib[j] = currentmc[i][j];
					}
					i++;
					subnode = subnode->succ;
				}

				/*Set decay rate back equal to current value*/
				subject->decay_f = current_decay;

			} /*End of if else statement*/

		}

		/*Go to next subject*/
		subject = subject->succ;
		parms->subindex++;

		/*Free memory */
		for (i = 0; i < num_node; i++)
			free(currentmc[i]);
		free(currentmc);
		free(pmd);

	} /*End of loop through subjects*/

}

/*********************************************************************/
/*START OF draw_times SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_times: this runs the M-H draws for each individual pulse location;
	ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
	Common_parms *parms; the current values of the common parameters;
	double **ts; this is the matrix of observed data (a column of
	times and S columns of log(concentration);
	int N; the number of observations in **ts;
	unsigned long *seed; seed values needed for the randon number generator;
	double v; the proposal variance for pulse location
	RETURNS: None; all updates are made internally; this does update the global
	variables regarding acceptance rate: atime and ntime
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i,k: generic counters
 *subject: current list of subjects and their qualities
 *subnode: for a given subject, current list of pulses and thier qualities
 alpha: minimum of 0 and log(rho); used to determine acceptance in M-H
 plikelihood: likelihood under proposed value
 ptime: proposed pulse location
 like_ratio: difference between current and proposed log-likelihoods
 prior_ratio: comparison between prior distributions under current and proposed values
 tmp: value used in evaluation of log(rho)
 current_time: current pulse location saved to this variable name so we can
 evaluate likelihood under proposed pulse location
 *tmp_mean_contrib: current mean_contrib for a pulse is saved to this vector
 so we can evaluate likelihood under proposed pulse's mean_contrib
 current_like: current likelihood for a given subject; used for comparison
 time_diff1: part of evaluation of prior_ratio portion of log(rho)
 time_diff1_new: same as time_diff1, but with proposed pulse location
 time_diff2: part of evaluation of prior_ratio portion of log(rho)
 time_diff2_new: same as time_diff1, but with proposed pulse location

 SUBROUTINES USED
 rnorm: found in randgen.h; draws from the normal distribution
 kiss: found in randgen.h; draws from U(0,1) distribution
 mean_contribution: found in birthdeath.c; evaluates mean_contrib for a
 pulse with specific parameters
 likelihood2: found in birthdeath.c; calculates and returns likelihood for a
 given subject under current parameters*/
/**************************************************************************/
void draw_times_old(Subject_type *sublist, Common_parms *parms, double **ts,
	int N, unsigned long *seed, double v) {
	int i, k;
	Subject_type *subject;
	Node_type *subnode;
	double alpha, plikelihood, ptime, like_ratio, prior_ratio, tmp, current_time, *tmp_mean_contrib;
	double current_like, time_diff1, time_diff1_new, time_diff2, time_diff2_new;
	double rnorm(double, double, unsigned long *);
	double kiss(unsigned long *);
	void mean_contribution(Node_type *, double **, Common_parms *, int, double);
	double likelihood2(Node_type *, double **, Common_parms *, int, Node_type *, double);

	/*Allocate Memory*/
	tmp_mean_contrib = (double *)calloc(N, sizeof(double));

	/*Go to start of subject list*/
	subject = sublist->succ;

	/*Index used for likelihood2*/
	parms->subindex = 0;

	/*Loop thru subjects*/
	while (subject != NULL){
		k = 0;
		/*Go to first pulse for this subject*/
		subnode = subject->driver->succ;

		/*Loop thru pulses for this subject*/
		while (subnode != NULL) {

			/*Increase denominator of acceptance rate for time*/
			ntime_l++;

			/*Find current likelihood for this subject*/
			current_like = likelihood2(subject->driver, ts, parms, N, subject->driver, subject->basehalf_l[0]);

			/*Compute proposal time*/
			ptime = rnorm(subnode->time, v, seed);

			/*Only proceed if our proposed time is reasonable */
			if (ptime <= fitend && ptime > fitstart) {

				/*If we're at the first pulse, use these as the first numerator/
				denominator of the prior ratio*/
				if (subnode->pred != NULL) {
					time_diff1_new = ptime - subnode->pred->time;
					time_diff1 = subnode->time - subnode->pred->time;
				}

				/*Otherwise, use these as the first numerator/denominator of
				the prior ratio*/
				else {
					time_diff1_new = ptime - fitstart;
					time_diff1 = subnode->time - fitstart;
				}

				/*If we're at the last pulse, use these as the second numerator/
				denominator of the prior ratio*/
				if (subnode->succ != NULL) {
					time_diff2_new = subnode->succ->time - ptime;
					time_diff2 = subnode->succ->time - subnode->time;
				}

				/*Otherwise, use these as the second numerator/denominator*/
				else {
					time_diff2_new = fitend - ptime;
					time_diff2 = fitend - subnode->time;
				}

				/*Combine it all for the prior ratio*/
				prior_ratio = (mmm - 1)*log((time_diff1_new*time_diff2_new) / (time_diff1 * time_diff2));

				/*Save current time*/
				current_time = subnode->time;

				/*Set pulse's time to the proposal value*/
				subnode->time = ptime;

				/*Save mean_contrib of this pulse*/
				for (i = 0; i<N; i++)
					tmp_mean_contrib[i] = subnode->mean_contrib[i];

				/*Recalculate mean_contrib of this pulse assuming proposed value*/
				mean_contribution(subnode, ts, parms, N, subject->basehalf_l[1]);

				/*Calculate the likelihood for this subject under the proposed value*/
				plikelihood = likelihood2(subject->driver, ts, parms, N, subject->driver, subject->basehalf_l[0]);

				/*Calculate the likelihood ratio*/
				like_ratio = plikelihood - current_like;

				/*Calculate log rho; set alpha equal to min(0,log rho)*/
				alpha = (0 < (tmp = (prior_ratio + like_ratio))) ? 0 : tmp;

				/*If log U < log rho, we accept proposed value. Increase acceptance
				count by one*/
				if (log(kiss(seed)) < alpha) {
					atime_l++;
				}

				/*Otherwise, we reject, and we have to revert back to current values*/
				else {

					/*Set pulse's time back to current time*/
					subnode->time = current_time;

					/*Set mean_contrib of this pulse back to current value*/
					for (i = 0; i<N; i++)
						subnode->mean_contrib[i] = tmp_mean_contrib[i];
				}
			} /*End of if time is feasible statement*/

			/*Advance one pulse*/
			subnode = subnode->succ;
			k++;
		}/*End of loop through pulses*/

		/*Go to next subject*/
		subject = subject->succ;
		parms->subindex++;

	}/*End of loop through subjects*/

	/*Free Memory*/
	free(tmp_mean_contrib);

}
void draw_times_old_f(Subject_type *sublist, Common_parms *parms, double **ts,
	int N, unsigned long *seed, double v) {
	int i, k;
	Subject_type *subject;
	Node_type *subnode;
	double alpha, plikelihood, ptime, like_ratio, prior_ratio, tmp, current_time, *tmp_mean_contrib;
	double current_like, time_diff1, time_diff1_new, time_diff2, time_diff2_new;
	double rnorm(double, double, unsigned long *);
	double kiss(unsigned long *);
	void mean_contribution(Node_type *, double **, Common_parms *, int, double);
	double likelihood2(Node_type *, double **, Common_parms *, int, Node_type *, double);

	/*Allocate Memory*/
	tmp_mean_contrib = (double *)calloc(N, sizeof(double));

	/*Go to start of subject list*/
	subject = sublist->succ;

	/*Index used for likelihood2*/
	parms->subindex = 0;

	/*Loop thru subjects*/
	while (subject != NULL){
		k = 0;
		/*Go to first pulse for this subject*/
		subnode = subject->response->succ;

		/*Loop thru pulses for this subject*/
		while (subnode != NULL) {

			/*Increase denominator of acceptance rate for time*/
			ntime_f++;

			/*Find current likelihood for this subject*/
			current_like = likelihood2(subject->response, ts, parms, N, subject->response, subject->basehalf_f[0]);

			/*Compute proposal time*/
			ptime = rnorm(subnode->time, v, seed);

			/*Only proceed if our proposed time is reasonable */
			if (ptime <= fitend && ptime > fitstart) {

				/*If we're at the first pulse, use these as the first numerator/
				denominator of the prior ratio*/
				if (subnode->pred != NULL) {
					time_diff1_new = ptime - subnode->pred->time;
					time_diff1 = subnode->time - subnode->pred->time;
				}

				/*Otherwise, use these as the first numerator/denominator of
				the prior ratio*/
				else {
					time_diff1_new = ptime - fitstart;
					time_diff1 = subnode->time - fitstart;
				}

				/*If we're at the last pulse, use these as the second numerator/
				denominator of the prior ratio*/
				if (subnode->succ != NULL) {
					time_diff2_new = subnode->succ->time - ptime;
					time_diff2 = subnode->succ->time - subnode->time;
				}

				/*Otherwise, use these as the second numerator/denominator*/
				else {
					time_diff2_new = fitend - ptime;
					time_diff2 = fitend - subnode->time;
				}

				/*Combine it all for the prior ratio*/
				prior_ratio = (mmm - 1)*log((time_diff1_new*time_diff2_new) / (time_diff1 * time_diff2));

				/*Save current time*/
				current_time = subnode->time;

				/*Set pulse's time to the proposal value*/
				subnode->time = ptime;

				/*Save mean_contrib of this pulse*/
				for (i = 0; i<N; i++)
					tmp_mean_contrib[i] = subnode->mean_contrib[i];

				/*Recalculate mean_contrib of this pulse assuming proposed value*/
				mean_contribution(subnode, ts, parms, N, subject->basehalf_f[1]);

				/*Calculate the likelihood for this subject under the proposed value*/
				plikelihood = likelihood2(subject->response, ts, parms, N, subject->response, subject->basehalf_f[0]);

				/*Calculate the likelihood ratio*/
				like_ratio = plikelihood - current_like;

				/*Calculate log rho; set alpha equal to min(0,log rho)*/
				alpha = (0 < (tmp = (prior_ratio + like_ratio))) ? 0 : tmp;

				/*If log U < log rho, we accept proposed value. Increase acceptance
				count by one*/
				if (log(kiss(seed)) < alpha) {
					atime_f++;
				}

				/*Otherwise, we reject, and we have to revert back to current values*/
				else {

					/*Set pulse's time back to current time*/
					subnode->time = current_time;

					/*Set mean_contrib of this pulse back to current value*/
					for (i = 0; i<N; i++)
						subnode->mean_contrib[i] = tmp_mean_contrib[i];
				}
			} /*End of if time is feasible statement*/

			/*Advance one pulse*/
			subnode = subnode->succ;
			k++;
		}/*End of loop through pulses*/

		/*Go to next subject*/
		subject = subject->succ;
		parms->subindex++;

	}/*End of loop through subjects*/

	/*Free Memory*/
	free(tmp_mean_contrib);

}

void draw_times_f(Subject_type *sublist, Common_parms *parms, Common_parms *parms_cox, double **ts,
	int N, unsigned long *seed, double v) {
	int i, ninter1,flag,k;
	Subject_type *subject;
	Node_type *subnode, *dnode,*node;
	double alpha, plikelihood, ptime, like_ratio, prior_ratio, tmp, current_time, *tmp_mean_contrib, condintnew, condintold, prop_ratio;
	double current_like, time_diff1, time_diff1_new, time_diff2, time_diff2_new;
	double rnorm(double, double, unsigned long *);
	double kiss(unsigned long *);
	void mean_contribution(Node_type *, double **, Common_parms *, int, double);
	double likelihood2(Node_type *, double **, Common_parms *, int, Node_type *, double);
	double phi(double, double, double);
	double skernel_1(Node_type *, Node_type *, Common_parms *);

	/*Allocate Memory*/
	tmp_mean_contrib = (double *)calloc(N, sizeof(double));

	/*Go to start of subject list*/
	subject = sublist->succ;

	/*Index used for likelihood2*/
	parms->subindex = 0;

	/*Loop thru subjects*/
	while (subject != NULL){
		k = 0;
		/*Go to first pulse for this subject*/
		subnode = subject->response->succ;

		/*Loop thru pulses for this subject*/
		while (subnode != NULL) {

			/*Increase denominator of acceptance rate for time*/
			ntime_f++;

			/*Find current likelihood for this subject*/
			current_like = likelihood2(subject->response, ts, parms, N, subject->response, subject->basehalf_f[0]);

			/*Compute proposal time*/
			ninter1 = 0;
			flag = 1;
			/*simulate a pulse position uniformly */
			/*if the new pulse location is too close to the old one, try it again*/
			/*like an hard core process*/
			while (flag == 1)
			{

				ptime = rnorm(subnode->time, v, seed);
				node = subject->response->succ;
				ninter1 = 0;
				while (node != NULL)
				{

					if (fabs(node->time - ptime) <= 0) {
						ninter1++;
					}
					node = node->succ;
				}
				if (ninter1 == 0) {
					flag = 0;


				}

			}

			/*Only proceed if our proposed time is reasonable */
			if (ptime <= fitend && ptime > fitstart) {

				prior_ratio = 0;

				condintnew = 0;
				condintold = 0;

				dnode = subject->driver->succ;
				while (dnode != NULL) {
					condintnew += exp(parms_cox->lrho) / sqrt(2 * 3.14159*parms_cox->nu) * exp(-1 / (2 * parms_cox->nu) * (ptime - dnode->time) * (ptime - dnode->time));
					condintold += exp(parms_cox->lrho) / sqrt(2 * 3.14159*parms_cox->nu) * exp(-1 / (2 * parms_cox->nu) * (subnode->time - dnode->time) * (subnode->time - dnode->time));
					dnode = dnode->succ;
				}

				if (condintnew == 0)
				{
					condintnew = 1e-300;
				}
				if (condintold == 0)
				{
					condintold = 1e-300;
				}

				/*Combine it all for the prior ratio*/
				condintnew = log(condintnew);
				condintold = log(condintold);

				prior_ratio = condintnew - condintold;
				prop_ratio = 0;
				/*prop_ratio = log(1 - phi(0, ptime, sqrt(v))) - log(1 - phi(0, subnode->time, sqrt(v)));*/

				/*Save current time*/
				current_time = subnode->time;

				/*Set pulse's time to the proposal value*/
				subnode->time = ptime;

				/*Save mean_contrib of this pulse*/
				for (i = 0; i < N; i++)
					tmp_mean_contrib[i] = subnode->mean_contrib[i];

				/*Recalculate mean_contrib of this pulse assuming proposed value*/
				mean_contribution(subnode, ts, parms, N, subject->basehalf_f[1]);

				/*Calculate the likelihood for this subject under the proposed value*/
				plikelihood = likelihood2(subject->response, ts, parms, N, subject->response, subject->basehalf_f[0]);

				/*Calculate the likelihood ratio*/
				like_ratio = plikelihood - current_like;

				/*Calculate log rho; set alpha equal to min(0,log rho)*/
				alpha = (0 < (tmp = (prior_ratio + like_ratio))) ? 0 : tmp;

				/*If log U < log rho, we accept proposed value. Increase acceptance
				count by one*/
				if (log(kiss(seed)) < alpha) {
					atime_f++;
					subnode->lambda = skernel_1(subnode, subject->driver, parms_cox);

				}

				/*Otherwise, we reject, and we have to revert back to current values*/
				else {

					/*Set pulse's time back to current time*/
					subnode->time = current_time;

					/*Set mean_contrib of this pulse back to current value*/
					for (i = 0; i < N; i++)
						subnode->mean_contrib[i] = tmp_mean_contrib[i];
				}
			} /*End of if time is feasible statement*/

			/*Advance one pulse*/
			subnode = subnode->succ;
			k++;
		}/*End of loop through pulses*/

		/*Go to next subject*/
		subject = subject->succ;
		parms->subindex++;

	}/*End of loop through subjects*/

	/*Free Memory*/
	free(tmp_mean_contrib);

}
void draw_times_l(Subject_type *sublist, Common_parms *parms_cox, double **ts,
	int N, unsigned long *seed, double v) {
	int i, j, k, l, ninterdiff1, ninterdiff2, ninterdiff3, num_parents, npulse;
	Subject_type *subject;
	Node_type *subnode, *node, *node2, *fnode;
	double alpha, plikelihood, new_loc, like_ratio, prior_ratio, tmp, current_time, *tmp_mean_contrib, lratio, lprop_ratio, tempdenomsum;
	double current_like, time_diff1, time_diff1_new, time_diff2, time_diff2_new, lkparent, kparent, lknew, knew;
	double rnorm(double, double, unsigned long *);
	double kiss(unsigned long *);
	void mean_contribution(Node_type *, double **, Common_parms *, int, double);
	double likelihood2(Node_type *, double **, Common_parms *, int, Node_type *, double);

	/*Allocate Memory*/

	tmp_mean_contrib = (double *)calloc(N, sizeof(double));

	/*Go to start of subject list*/
	subject = sublist->succ;
	npulse = subject->numnode_l;

	/*Index used for likelihood2*/
	parms_cox->subindex = 0;

	/*Loop thru subjects*/
	while (subject != NULL){
		k = 0;
		/*Go to first pulse for this subject*/
		subnode = subject->driver->succ;

		/*Loop thru pulses for this subject*/
		while (subnode != NULL) {

			/*Increase denominator of acceptance rate for time*/
			ntime_l++;

			/*Find current likelihood for this subject*/
			current_like = likelihood2(subject->driver, ts, parms_cox, N, subject->driver, subject->basehalf_l[0]);

			/*Compute proposal time*/
			new_loc = rnorm(subnode->time, sqrt(v), seed);

			/*Only proceed if our proposed time is reasonable */
			if (new_loc <= fitend && new_loc > fitstart) {

				ninterdiff1 = 0;
				ninterdiff2 = 0;
				ninterdiff3 = 0;
				node = subject->driver->succ;

				while (node != NULL && ninterdiff1 == 0){
					/*        printf("node->time %lf new_loc %lf\n",node->time,new_loc);*/
					if (subnode->time != node->time && fabs(node->time - new_loc) <= parms_cox->Rinter1) {
						ninterdiff1++;
					}
					else {
						if (subnode->time != node->time && fabs(node->time - new_loc) <= parms_cox->Rinter2) {
							ninterdiff2++;
						}
						if (subnode->time != node->time && fabs(subnode->time - node->time) <= parms_cox->Rinter2) {
							ninterdiff3++;
						}
					}
					node = node->succ;
				}
				/* printf("new_loc %lf node2->time %lf\n",new_loc,node2->time);
				 printf("ninterdiff1 %d\n",ninterdiff1);
				 printf("ninter2 %d\n",ninterdiff2);*/
				if (ninterdiff1 == 0) {
					/*conditional intensity for strauss/hardcore priors_cox*/
					lprop_ratio = (double)(ninterdiff2 - ninterdiff3) * log(parms_cox->gamma);

					/*component 1 of the conditional intensity for the parent given the data*/
					//lprop_ratio += exp(parms_cox->lrho) * (phi(fitend, subnode->time, sqrt(parms_cox->nu)) - phi(fitstart, subnode->time, sqrt(parms_cox->nu)) - phi(fitend, new_loc, sqrt(parms_cox->nu)) + phi(fitstart, new_loc, sqrt(parms_cox->nu)));

					/*ratio of proposal densities*/
					/*  lprop_ratio += log(phi(fitend, new_loc, pparents_sd) - phi(fitstart, new_loc, pparents_sd)) - log(phi(fitend, node2->time, pparents_sd) - phi(fitstart, node2->time, pparents_sd));*/

					/*product in the posterior cond. intensity */
					lratio = 0;
					k = 0;
					fnode = subject->response->succ;
					while (fnode != NULL){
						lkparent = -1 / (2 * parms_cox->nu) * ((fnode->time - subnode->time) * (fnode->time - subnode->time));
						kparent = exp(parms_cox->lrho) / sqrt(2 * 3.14159 * parms_cox->nu) * exp(lkparent);
						lknew = -1 / (2 * parms_cox->nu) * ((fnode->time - new_loc) * (fnode->time - new_loc));
						knew = exp(parms_cox->lrho) / sqrt(2 * 3.14159 * parms_cox->nu) * exp(lknew);
						tempdenomsum = fnode->lambda - kparent + knew;
						if (tempdenomsum == 0)
							tempdenomsum = 1e-300;
						lratio += log(tempdenomsum) - log(fnode->lambda);
						fnode = fnode->succ;
						k++;

					}            /*Save current time*/
					lprop_ratio += lratio;

					current_time = subnode->time;

					/*Set pulse's time to the proposal value*/
					subnode->time = new_loc;

					/*Save mean_contrib of this pulse*/
					for (i = 0; i < N; i++)
						tmp_mean_contrib[i] = subnode->mean_contrib[i];

					/*Recalculate mean_contrib of this pulse assuming proposed value*/
					mean_contribution(subnode, ts, parms_cox, N, subject->basehalf_l[1]);

					/*Calculate the likelihood for this subject under the proposed value*/
				plikelihood = likelihood2(subject->driver, ts, parms_cox, N, subject->driver, subject->basehalf_l[0]);
				

					like_ratio = plikelihood - current_like;

					/*Calculate log rho; set alpha equal to min(0,log rho)*/
					alpha = (0 < (tmp = (like_ratio + lprop_ratio))) ? 0 : tmp;

					/*If log U < log rho, we accept proposed value. Increase acceptance
					count by one*/
					if (log(kiss(seed)) < alpha) {
						atime_l++;
						fnode = subject->response->succ;
						while (fnode != NULL){

							fnode->lambda = skernel_1(fnode, subject->driver, parms_cox);
							fnode = fnode->succ;
						}
					}

					/*Otherwise, we reject, and we have to revert back to current values*/
					else {

						/*Set pulse's time back to current time*/
						subnode->time = current_time;

						/*Set mean_contrib of this pulse back to current value*/
						for (i = 0; i < N; i++)
							subnode->mean_contrib[i] = tmp_mean_contrib[i];
					}
				}
			}		/*End of if time is feasible statement*/

			/*Advance one pulse*/
			subnode = subnode->succ;
			k++;
		}/*End of loop through pulses*/

		/*Go to next subject*/
		subject = subject->succ;
		parms_cox->subindex++;

	}/*End of loop through subjects*/

	/*Free Memory*/
	free(tmp_mean_contrib);

}

/*********************************************************************/
/*START OF draw_random_effects SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_random_effects: this runs the M-H draw for individual pulse masses and
					  widths;
					  ARGUMENTS: double **ts; this is the matrix of observed data (a column of
					  times and a column of log(concentration);
					  Subject_type *list; this is the current list of subjects that exist;
					  Common_parms *parms; the current values of the common parameters;
					  int N; the number of observations in **ts;
					  double v1; the proposal variance for individual pulse masses;
					  double v2; the proposal variance for individual pulse widths;
					  unsigned long *seed; seed values needed for the randon number generator;
					  RETURNS: None; all updates are made internally; this does update the global
					  variables regarding acceptance rate: arem, arew, nrem, nrew
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i,j,k: generic counters
 temp: value used in evaluation of log(rho)
 prior_old: part of evaluation of prior_ratio portion of log(rho)
 prior_new: same as prior_old, but with proposed value
 old_val: current value of the random effect is saved to this variable so we
 can evaluate likelihood under the proposed value
 *tmp1: current acceptance counters are saved to this vector so we can increment
 inside a loop
 *pRE: proposed values of individual mass and width
 prior_ratio: comparison between prior distributions under current and proposed values
 like_ratio: difference between current and proposed log-likelihoods
 alpha: minimum of 0 and log(rho); used to determine acceptance in M-H
 plikelihood: likelihood under proposed value
 *old_contrib: current mean_contrib for a pulse is saved to this vector
 so we can evaluate likelihood under proposed pulse's mean_contrib
 current_like: current likelihood for a given subject; used for comparison
 *subnode: for a given subject, current list of pulses and their qualities
 *subject: current list of subjects and their qualities

 SUBROUTINES USED
 mean_contribution: found in birthdeath.c; evaluates mean_contrib for a
 pulse with specific parameters
 kiss: found in randgen.h; draws from U(0,1) distribution
 rnorm: found in randgen.h; draws from the normal distribution
 likelihood2: found in birthdeath.c; calculates and returns likelihood for a
 given subject under current parameters
**************************************************************************/

void draw_random_effects_l(double **ts, Subject_type *sublist, Common_parms *parms, int N, double v1, double v2, unsigned long *seed)
{
	int *tmp1, i, j, k;
	double temp, prior_old, prior_new, old_val, *pRE, prior_ratio, like_ratio;
	double alpha, plikelihood, *old_contrib, current_like, stdold[2], stdnew[2], newint[2], oldint[2];
	Node_type *subnode;
	Subject_type *subject;
	void mean_contribution(Node_type *, double **, Common_parms *, int, double);
	double kiss(unsigned long *);
	double rnorm(double, double, unsigned long *);
	double likelihood2(Node_type *, double **, Common_parms *, int, Node_type *, double);
	double Phi(double);

	/*Set acceptance counts equal to temporary vector*/
	tmp1 = (int *)calloc(2, sizeof(int));
	tmp1[0] = arem_l;
	tmp1[1] = arew_l;

	/*Go to first subject*/
	subject = sublist->succ;

	/*Index of subjects used for likelihood2*/
	parms->subindex = 0;

	/*Loop thru subjects*/
	while (subject != NULL){

		/*Go to start of pulses for this subject*/
		subnode = subject->driver->succ;
		k = 0;

		/*Go through each existing pulse for this subject*/
		while (subnode != NULL) {

			/*Allocate Memory*/
			old_contrib = (double *)calloc(N, sizeof(double));
			pRE = (double *)calloc(2, sizeof(double));

			/*Increase the denominators of the acceptance rates*/
			nrem_l++;
			nrew_l++;

			/*Draw proposed values of current pulse's mass and width*/
			pRE[0] = rnorm(subnode->theta[0], v1, seed);
			pRE[1] = rnorm(subnode->theta[1], v2, seed);
			/*only proceed when we draw positive value*/
			if (pRE[0] > 0.0 && pRE[1] > 0.01 && pRE[1] < 100){

				/*Determine if we accept or reject proposed pulse mass then determine*/
				/*if we accept or reject proposed pulse width*/
				for (j = 0; j < 2; j++){
					current_like = likelihood2(subject->driver, ts, parms, N, subject->driver, subject->basehalf_l[0]);

					/*Compute the log of the ratio of the priors*/
					prior_old = subnode->theta[j] - subject->theta_l[j];
					prior_old *= 0.5*prior_old;
					prior_new = pRE[j] - subject->theta_l[j];
					prior_new *= 0.5*prior_new;
					prior_ratio = (prior_old - prior_new)*subnode->eta[j];
					prior_ratio /= parms->re_precision[j];
					prior_ratio /= parms->re_precision[j];
					stdold[j] = (subnode->theta[j]) / (parms->re_precision[j] / sqrt(subnode->eta[j]));
					stdnew[j] = pRE[j] / (parms->re_precision[j] / sqrt(subnode->eta[j]));

					newint[j] = log(Phi(stdnew[j]));
					oldint[j] = log(Phi(stdold[j]));
					/*prior_ratio = prior_ratio-newint[j]+oldint[j];*/
					/*Save the current value of mass/width*/
					old_val = subnode->theta[j];

					/*Set the pulse's mass/width equal to the proposed value*/
					subnode->theta[j] = pRE[j];

					/*Save the mean_contrib for that pulse*/
					for (i = 0; i < N; i++)
						old_contrib[i] = subnode->mean_contrib[i];

					/*Recalculate that pulse's mean_contrib assuming proposed mass/width */
					mean_contribution(subnode, ts, parms, N, subject->basehalf_l[1]);

					/*Calculate likelihood assuming proposed mass/width */
					plikelihood = likelihood2(subject->driver, ts, parms, N, subject->driver, subject->basehalf_l[0]);

					like_ratio = plikelihood - current_like;

					/*Compute the log of the ratio between the two likelihoods*/

					/*Calculate log rho; set alpha equal to min(0,log rho) */
					alpha = (0 < (temp = (prior_ratio + like_ratio))) ? 0 : temp;

					/*If log U < log rho, accept the proposed value, increase acceptance*/
					/*counter */
					if (log(kiss(seed)) < alpha) {
						tmp1[j]++;
					}

					/*Otherwise, reject the proposed value, set pulse's mass/width back to */
					/*saved current value, and set pulse's mean_contrib back to saved value */
					else {
						subnode->theta[j] = old_val;
						for (i = 0; i < N; i++)
							subnode->mean_contrib[i] = old_contrib[i];
					}
				}
			}		 /*end of loop through mass & width */

			/*Advance to next pulse*/
			subnode = subnode->succ;
			k++;
			/*free memory*/
			free(pRE);
			free(old_contrib);

		} /*end of loop through pulses*/

		/*advance to next subject*/
		subject = subject->succ;
		parms->subindex++;

	} /*end of loop through subjects*/
	/*Set counter equal to temporary vector values*/
	arem_l = tmp1[0];
	arew_l = tmp1[1];
	free(tmp1);

}
void draw_random_effects_f(double **ts, Subject_type *sublist, Common_parms *parms, int N, double v1, double v2, unsigned long *seed)
{
	int *tmp1, i, j, k;
	double temp, prior_old, prior_new, old_val, *pRE, prior_ratio, like_ratio;
	double alpha, plikelihood, *old_contrib, current_like, stdold[2], stdnew[2], newint[2], oldint[2];
	Node_type *subnode;
	Subject_type *subject;
	void mean_contribution(Node_type *, double **, Common_parms *, int, double);
	double kiss(unsigned long *);
	double rnorm(double, double, unsigned long *);
	double likelihood2(Node_type *, double **, Common_parms *, int, Node_type *, double);
	double Phi(double);

	/*Set acceptance counts equal to temporary vector*/
	tmp1 = (int *)calloc(2, sizeof(int));
	tmp1[0] = arem_f;
	tmp1[1] = arew_f;

	/*Go to first subject*/
	subject = sublist->succ;

	/*Index of subjects used for likelihood2*/
	parms->subindex = 0;

	/*Loop thru subjects*/
	while (subject != NULL){

		/*Go to start of pulses for this subject*/
		subnode = subject->response->succ;
		k = 0;

		/*Go through each existing pulse for this subject*/
		while (subnode != NULL) {

			/*Allocate Memory*/
			old_contrib = (double *)calloc(N, sizeof(double));
			pRE = (double *)calloc(2, sizeof(double));

			/*Increase the denominators of the acceptance rates*/
			nrem_f++;
			nrew_f++;

			/*Draw proposed values of current pulse's mass and width*/
			pRE[0] = rnorm(subnode->theta[0], v1, seed);
			pRE[1] = rnorm(subnode->theta[1], v2, seed);
			/*only proceed when we draw positive value*/
			if (pRE[0] > 0.0 && pRE[1] > 0.01 && pRE[1] < 100){

				/*Determine if we accept or reject proposed pulse mass then determine*/
				/*if we accept or reject proposed pulse width*/
				for (j = 0; j < 2; j++){
					current_like = likelihood2(subject->response, ts, parms, N, subject->response, subject->basehalf_f[0]);

					/*Compute the log of the ratio of the priors*/
					prior_old = subnode->theta[j] - subject->theta_f[j];
					prior_old *= 0.5*prior_old;
					prior_new = pRE[j] - subject->theta_f[j];
					prior_new *= 0.5*prior_new;
					prior_ratio = (prior_old - prior_new)*subnode->eta[j];
					prior_ratio /= parms->re_precision[j];
					prior_ratio /= parms->re_precision[j];
					stdold[j] = (subnode->theta[j]) / (parms->re_precision[j] / sqrt(subnode->eta[j]));
					stdnew[j] = (pRE[j]) / (parms->re_precision[j] / sqrt(subnode->eta[j]));

					newint[j] = log(Phi(stdnew[j]));
					oldint[j] = log(Phi(stdold[j]));
					/*prior_ratio = prior_ratio-newint[j]+oldint[j];*/
					/*Save the current value of mass/width*/
					old_val = subnode->theta[j];

					/*Set the pulse's mass/width equal to the proposed value*/
					subnode->theta[j] = pRE[j];

					/*Save the mean_contrib for that pulse*/
					for (i = 0; i < N; i++)
						old_contrib[i] = subnode->mean_contrib[i];

					/*Recalculate that pulse's mean_contrib assuming proposed mass/width */
					mean_contribution(subnode, ts, parms, N, subject->basehalf_f[1]);

					/*Calculate likelihood assuming proposed mass/width */
					plikelihood = likelihood2(subject->response, ts, parms, N, subject->response, subject->basehalf_f[0]);

					like_ratio = plikelihood - current_like;

					/*Compute the log of the ratio between the two likelihoods*/

					/*Calculate log rho; set alpha equal to min(0,log rho) */
					alpha = (0 < (temp = (prior_ratio + like_ratio))) ? 0 : temp;

					/*If log U < log rho, accept the proposed value, increase acceptance*/
					/*counter */
					if (log(kiss(seed)) < alpha) {
						tmp1[j]++;
					}

					/*Otherwise, reject the proposed value, set pulse's mass/width back to */
					/*saved current value, and set pulse's mean_contrib back to saved value */
					else {
						subnode->theta[j] = old_val;
						for (i = 0; i < N; i++)
							subnode->mean_contrib[i] = old_contrib[i];
					}
				}
			}		 /*end of loop through mass & width */

			/*Advance to next pulse*/
			subnode = subnode->succ;
			k++;
			/*free memory*/
			free(pRE);
			free(old_contrib);

		} /*end of loop through pulses*/

		/*advance to next subject*/
		subject = subject->succ;
		parms->subindex++;

	} /*end of loop through subjects*/
	/*Set counter equal to temporary vector values*/
	arem_f = tmp1[0];
	arew_f = tmp1[1];
	free(tmp1);

}
/*********************************************************************/
/*START OF error_squared SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*error_squared: This subroutine calculates the sum of squared error. It uses
				the list of pulses and the common parameters to evaluate the
				likelihood and calculate the sum of the squared error between
				observed concentration and expected under the current parameters
				ARGUMENTS: double **ts; this is the matrix of observed data (a column of
				times and S columns of log(concentration)
				Subject_type *sublist; this is the current list of subjects that exist
				Common_parms *parms; the current values of the common parameters
				int N; the number of observations in **ts;
				RETURNS: double ssq; this value is needed for one of the parameters
				of the posterior distribution of model error variance
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i,j: generic counters
 *mean: mean concentration under current parameters and pulses
 ssq: sum of squared differences between observed log(concentration) and expected
 values
 *subject: current list of subjects and their qualities

 SUBROUTINES USED
 mean_concentration: found in birthdeath.c; sums across mean_contrib vectors for
 each pulse; returns this vector
 **************************************************************************/

double error_squared_l(double **ts, Subject_type *sublist, Common_parms *parms, int N) {
	int i, j;
	double *mean, ssq;
	Subject_type *subject;
	double *mean_concentration(Node_type *, Common_parms *, int, Node_type *, double);

	subject = sublist->succ;
	ssq = 0;
	for (j = 0; j < parms->numsub; j++){

		mean = mean_concentration(subject->driver, parms, N, subject->driver, subject->basehalf_l[0]);
		for (i = 0; i < N; i++)
			ssq += (ts[i][j + 1] - mean[i])*(ts[i][j + 1] - mean[i]);

		subject = subject->succ;
		free(mean);
	}

	return ssq;
}

double error_squared_f(double **ts, Subject_type *sublist, Common_parms *parms, int N) {
	int i, j;
	double *mean, ssq;
	Subject_type *subject;
	double *mean_concentration(Node_type *, Common_parms *, int, Node_type *, double);

	subject = sublist->succ;
	ssq = 0;
	for (j = 0; j < parms->numsub; j++){

		mean = mean_concentration(subject->response, parms, N, subject->response, subject->basehalf_f[0]);
		for (i = 0; i < N; i++)
			ssq += (ts[i][j + 1] - mean[i])*(ts[i][j + 1] - mean[i]);

		subject = subject->succ;
		free(mean);
	}

	return ssq;
}
void draw_eta_l(Subject_type *sublist, Common_parms *parms, unsigned long *seed, double *v)
{
	int *tmp, i, j;
	double x, peta, prior_ratio, re_ratio, old_gamma, new_gamma, stdold, stdnew, re_old, re_new, alpha, temp;
	Subject_type *subject;
	Node_type *subnode;
	double rnorm(double, double, unsigned long *);
	double kiss(unsigned long *);
	double Phi(double);


	tmp = (int *)calloc(2, sizeof(int));
	tmp[0] = ae_ml;
	tmp[1] = ae_wl;

	/*first pulse mass and then pulse width*/

	subject = sublist->succ;
	while (subject != NULL){

		subnode = subject->driver->succ;
		while (subnode != NULL) {

			/*draw the new eta*/

			for (j = 0; j<2; j++){
				peta = rnorm(subnode->eta[j], v[j], seed);
				ne_ml++;
				ne_wl++;
				if (peta > 0) {
					old_gamma = gammaPdf(subnode->eta[j], 2., 2.);
					new_gamma = gammaPdf(peta, 2., 2.);
					prior_ratio = log(new_gamma) - log(old_gamma);
					stdold = (subnode->theta[j]) / (parms->re_precision[j] / sqrt(subnode->eta[j]));
					stdnew = (subnode->theta[j]) / (parms->re_precision[j] / sqrt(peta));
					re_old = subnode->theta[j] - subject->theta_l[j];
					re_old *= 0.5*re_old*subnode->eta[j];
					re_new = subnode->theta[j] - subject->theta_l[j];
					re_new *= 0.5*re_new*peta;
					re_ratio = re_old - re_new;
					re_ratio /= parms->re_precision[j];
					re_ratio /= parms->re_precision[j];
					re_ratio += log(Phi(stdold)) - log(Phi(stdnew)) - 0.5*log(subnode->eta[j]) + 0.5*log(peta); /*the 1/2pi term in normal distirbution*/
					alpha = (0 < (temp = (prior_ratio + re_ratio))) ? 0 : temp;

					/*If log U < log rho, accept the proposed value, increase acceptance*/
					/*counter */
					if (log(kiss(seed)) < alpha) {
						tmp[j]++;
						subnode->eta[j] = peta;
					}


				}
				/*subnode->eta[j] = 1;*/

			}

			subnode = subnode->succ;

		}
		subject = subject->succ;

	}

	free(tmp);
	ae_ml = tmp[0];
	ae_wl = tmp[1];
}

void draw_eta_f(Subject_type *sublist, Common_parms *parms, unsigned long *seed, double *v)
{
	int *tmp, i, j;
	double x, peta, prior_ratio, re_ratio, old_gamma, new_gamma, stdold, stdnew, re_old, re_new, alpha, temp;
	Subject_type *subject;
	Node_type *subnode;
	double rnorm(double, double, unsigned long *);
	double kiss(unsigned long *);
	double Phi(double);


	tmp = (int *)calloc(2, sizeof(int));
	tmp[0] = ae_mf;
	tmp[1] = ae_wf;

	/*first pulse mass and then pulse width*/

	subject = sublist->succ;
	while (subject != NULL){

		subnode = subject->response->succ;
		while (subnode != NULL) {

			/*draw the new eta*/

			for (j = 0; j<2; j++){
				peta = rnorm(subnode->eta[j], v[j], seed);
				ne_mf++;
				ne_wf++;
				if (peta > 0) {
					old_gamma = gammaPdf(subnode->eta[j], 2., 2.);
					new_gamma = gammaPdf(peta, 2., 2.);
					prior_ratio = log(new_gamma) - log(old_gamma);
					stdold = (subnode->theta[j]) / (parms->re_precision[j] / sqrt(subnode->eta[j]));
					stdnew = (subnode->theta[j]) / (parms->re_precision[j] / sqrt(peta));
					re_old = subnode->theta[j] - subject->theta_f[j];
					re_old *= 0.5*re_old*subnode->eta[j];
					re_new = subnode->theta[j] - subject->theta_f[j];
					re_new *= 0.5*re_new*peta;
					re_ratio = re_old - re_new;
					re_ratio /= parms->re_precision[j];
					re_ratio /= parms->re_precision[j];
					re_ratio += log(Phi(stdold)) - log(Phi(stdnew)) - 0.5*log(subnode->eta[j]) + 0.5*log(peta);
					alpha = (0 < (temp = (prior_ratio + re_ratio))) ? 0 : temp;

					/*If log U < log rho, accept the proposed value, increase acceptance*/
					/*counter */
					if (log(kiss(seed)) < alpha) {
						tmp[j]++;
						subnode->eta[j] = peta;
					}

					/*Otherwise, reject the proposed value, set pulse's mass/width back to */
					/*saved current value, and set pulse's mean_contrib back to saved value */

				}
				//subnode->eta[j] = 1;

			}

			subnode = subnode->succ;

		}
		subject = subject->succ;

	}
    ae_mf = tmp[0];
	ae_wf = tmp[1];
	free(tmp);
	
}
double Phi(double y)
{
	/* Returns standard normal distribution function evaluated at y */
	/* Has absolute error of order 10^(-7) */
	/* Uses approximation to erfc from Numerical Recipes in C */
	double t, z, ans, x, inv_sqrt2 = 0.7071067811865475;
	x = y*inv_sqrt2;
	z = fabs(x);
	t = 1.0 / (1.0 + .5*z);
	ans = t*exp(-z*z - 1.26551223 + t*(1.00002368 + t*(0.37409196 + t*(0.09678418 +
		t*(-0.18628806 + t*(0.27886807 + t*(-1.13520398 + t*(1.48851587 +
		t*(-0.82215223 + t*0.17087277)))))))));
	if (!(x < 0))
		return 1.0 - 0.5*ans;
	else
		return 0.5*ans;
}
/*********************************************************************/
/*START OF adjust_acceptance SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*adjust_acceptance: this adjusts the proposal variances based on the inputted
					acceptance rate and proposal variance. If the acceptance rate is
					too high or too low, proposal variance will be adjusted.
					ARGUMENTS: double x; the inputted acceptance rate; usually inputted as the
					acceptance counter divided by the attempt counter
					double *X; the current proposal variance
					RETURNS: None; update to the proposal variance is made internally
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 y: new proposal variance based on inputs

 SUBROUTINES USED
 None
 **************************************************************************/

void adjust_acceptance(double x, double *X)
{
	double y;

	y = 1. + 1000.*(x - .35)*(x - .35)*(x - .35);
	if (y < .9)
		y = .9;
	if (y > 1.1)
		y = 1.1;

	*X *= y;
}

/*********************************************************************/
/*START OF adjust2_acceptance SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*adjust2_acceptance: this adjusts the proposal variance-covriance matrix for
					 baseline and halflife based on the inputted acceptance rate
					 and proposal matrix. If the acceptance rate is too high or
					 too low, proposal variance will be adjusted.
					 ARGUMENTS: double x; the inputted acceptance rate; usually inputted as the
					 acceptance counter divided by the attempt counter
					 double **X; the current proposal variance-covariance matrix
					 double corr; the correlation between the two proposal variances
					 RETURNS: None; update to the proposal variance is made internally
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 y: new diagonal elements of proposal variance-covariance matrix based on inputs

 SUBROUTINES USED
 None
**************************************************************************/

void adjust2_acceptance(double x, double **X, double corr)
{
	double y;

	y = 1. + 1000.*(x - .25)*(x - .25)*(x - .25);
	if (y < .90) {
		y = .90;
		X[0][0] *= y;
		X[1][1] *= y;
		X[0][1] = X[1][0] = corr * sqrt(X[0][0] * X[1][1]);
	}
	if (y > 1.1) {
		y = 1.1;
		X[0][0] *= y;
		X[1][1] *= y;
		X[0][1] = X[1][0] = corr * sqrt(X[0][0] * X[1][1]);
	}
}

// the gamma distribution PDF
double gamm(double x)
{
	double ret = (1.000000000190015 +
		76.18009172947146 / (x + 1) +
		-86.50532032941677 / (x + 2) +
		24.01409824083091 / (x + 3) +
		-1.231739572450155 / (x + 4) +
		1.208650973866179e-3 / (x + 5) +
		-5.395239384953e-6 / (x + 6));

	return ret * sqrt(2 * M_PI) / x * pow(x + 5.5, x + .5) * exp(-x - 5.5);
}




// the gamma distribution PDF
double gammaPdf(double x, double a, double b)
{
	if (x <= 0 || a <= 0 || b <= 0)
		return 0.0;
	else
		return exp(-x * b) * pow(x, a - 1.0) * pow(b, a) / gamm(a);
}




double phi(double y, double mu, double s)
{
	/* Returns normal distribution function evaluated at y */
	/* Returns standard normal distribution function evaluated at y */
	/* Has absolute error of order 10^(-7) */
	/* Uses approximation to erfc from Numerical Recipes in C */
	double t, z, ans, x;

	x = (y - mu) / (1.4142135623730951*s);
	z = fabs(x);
	t = 1.0 / (1.0 + .5*z);
	ans = t*exp(-z*z - 1.26551223 + t*(1.00002368 + t*(0.37409196 + t*(0.09678418 +
		t*(-0.18628806 + t*(0.27886807 + t*(-1.13520398 + t*(1.48851587 +
		t*(-0.82215223 + t*0.17087277)))))))));
	if (x >= 0){
		return 1.0 - 0.5*ans;
	}
	else{
		return 0.5*ans;
	}
}