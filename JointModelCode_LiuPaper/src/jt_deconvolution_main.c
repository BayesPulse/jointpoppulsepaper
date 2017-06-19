/*******************************************************************/
/********************MAIN BDMCMC PROGRAM***************************/
/*******************************************************************/

#include "jt_deconvolution_main.h"

/*******************************************************************
*******************GLOBAL VARIABLE DEFINITIONS**********************

 * M: Sets the environment precision and used for all random number generation
    * Needed specfically for the KISS subroutine
 * mmm: Order statistic used for distribution of pulse locations.
    * This is inputted by the user and is typically 3
 * fitstart: The first time in hours that a pulse may occur
 * fitend: The last time in hours that a pulse may occur

*********************************************************************/

double M;
int mmm;
double fitstart;
double fitend;

#if _WIN32 || _WIN64
#if _WIN64
#define ENVIRONMENT 64
#else
#define ENVIRONMENT 32
#endif
#endif

// Check GCC
#if __GNUC__
#if __x86_64__ || __ppc64__
#define ENVIRONMENT 64
#else
#define ENVIRONMENT 32
#endif
#endif

/********************************************************************/
/*SUBROUTINES THAT EXIST IN THIS PROGRAM

 start_position: This subroutine initializes all the parameters in the study.
    Initial baseline and halflife are user-inputted; all others are set in
    the subroutine. This subroutine also creates our list of nodes and inserts
    one new node with parameters equal to 0.
    ARGUMENTS: Common_parms *parms; The values in this data structure are
               modified in this routine;
               double pmean; starting value for baseline inputted by user in
               main; must be greater than 0
               double pdelta; starting value for halflife inputted by user in
               main; must be greater than 0
    RETURNS: Node_type *list; returns the initialized node list
 **********************************************************************/

/**********************************************************************
 MAIN PROGRAM STARTS HERE

 VARIABLE DEFINITIONS
    i: generic counter
    *N: number of observations in inputted file; determined by read_data_file
        subroutine
    iter: number of iterations to run the MCMC subroutine called in main;
          this variable is inputted by the user
    a: generic counter
    *seed: 3 element vector with seed values for the random number generator;
           read in from an external file
    **ts: matrix containing column of time points and column of hormone
          concentrations
    **pmd_var: Variance-Covariance matrix of the proposal distribution of
               baseline and halflife; matrix is defined in main
    pmean: Starting value for baseline; inputted by the user
    pdelta: Starting value for halflife; inputted by the user
    *fseed: The pointer to the file that contains the random numbers for the
            random number generator
    *parms: Contains the parameters that are common throughout the model
    *priors: Contains all parameters of the prior distributions used
    *hyper: Contains all the hyperpriors (though this is not used)
    *list: Contains the list of nodes and their characteristics

 SUBROUTINES USED
    **read_data_file: Found in format_data.c; scans inputted data file and
                      returns the matrix of data (time and concentration)
    rnorm: Found in randgen.h; draws from the normal distribution
    rgamma: Found in randgen.h; draws from the gamma distribution
    destroy_list: Found in hash.h; frees all memory used by the nodes
    mcmc: Found in linklistv2.c; This runs birth-death mcmc routine
    *start_position: Found in this file; discussed above 
 ************************************************************************/

int main(int argc,char *argv[])
{
  int i,*N,iter,a,subj;
  char datafile_l[90];
  char datafile_f[100];

  char commonf[100];
    char commonl[100];

  char parml[100];
    char parmf[100];

  char tmp[5];
  unsigned long *seed;
  double **ts_l,**ts_f,propvar[33] ,**temp2;
 double mprior1, mprior2, mprior3, mprior4;
  double psprior1, psprior2, psprior3, psprior4,psprior5, psprior6, psprior7, psprior8,psprior9;
  double bprior1, bprior2, bprior3, bprior4, bprior5, bprior6;
  double hprior1, hprior2, hprior3,hprior4, hprior5, hprior6;
  double sprior1, sprior2,sprior3,sprior4;
  /*int nprior1, nprior2;*/
  double gprior,betaprior;
  double rhoprior1,rhoprior2;
  double nuprior1,nuprior2;
  
  double mstart1, mstart2, mstart3, mstart4;
  double pstart1, pstart2, pstart3,pstart4, pstart5, pstart6,pstart7, pstart8, pstart9,psprior10, psprior11;
  double bstart1, bstart2,bstart3, bstart4;
  double hstart1, hstart2,hstart3, hstart4;
  double sstart1,sstart2;
  double rhostart,nustart;

  FILE *fseed, *finput;
  Common_parms *parms_l,*parms_f;
  Priors *priors;
  Hyper_priors *hyper;
  Subject_type *subject, *sublist;

  double **read_data_file(char *,int *,int);
  double rnorm(double,double,unsigned long *);
  double rgamma(double,double,unsigned long *);
  void destroy_sublist(Subject_type *);
  Subject_type *initialize_subject(void);
  void insert_subject(Subject_type *,Subject_type *);
  void mcmc(Subject_type *, Common_parms *, Common_parms *, double **, double **,long, int, Priors *, unsigned long *,
	  char *,  char *, double[], Hyper_priors *);
  char *itoa(int , char *, int );
  int cholesky_decomp(double **,int);
  double **cholesky_invert(int, double **);



  Subject_type *initialize_subject(void);

  if (argc != 2 ) exit(0);

/************************************/
/* M is used in the KISS random number generator */
    M = exp(-ENVIRONMENT*log(2.0));
/*************************************************************/
    
  for (a=0;a<1;a++) {
  
    seed = (unsigned long *)calloc(3,sizeof(unsigned long));
    fseed = fopen("seed.dat","r");
    fscanf(fseed,"%lu %lu %lu\n",&seed[0],&seed[1],&seed[2]);
    fclose(fseed);
    
    finput = fopen(argv[1],"r");

    fscanf(finput,"%s %s \n", datafile_l, datafile_f);  /*First line: data for lh and fsh*/
    fscanf(finput,"%s %s %s %s \n", commonl, commonf, parml,parmf); /*second line: output file name for result*/
    fscanf(finput,"%d %d \n", &subj, &iter);  /*third line: number of subject, number of iteration*/

    /* read in the hormonal time series */
    N = (int *)calloc(1,sizeof(int));
    ts_l = read_data_file(datafile_l,N,subj);
	    ts_f = read_data_file(datafile_f,N,subj);

    mmm = 3;
    
    /* perform the simulation */
    
    parms_l = (Common_parms *)calloc(1,sizeof(Common_parms));
	    parms_f = (Common_parms *)calloc(1,sizeof(Common_parms));
    parms_l->re_precision = (double *)calloc(2,sizeof(double ));
        parms_f->re_precision = (double *)calloc(2,sizeof(double ));



    fitend = ts_l[*N-1][0]+ ts_l[1][0] * 4;  /*search 2 units farther in time*/
    fitstart = -ts_l[1][0] * 4;  /*search 4 units in the past*/

  


  
    priors = (Priors *)calloc(1,sizeof(Priors));
    hyper = (Hyper_priors *)calloc(1,sizeof(Hyper_priors));
    
    priors->re_var_f = (double *)calloc(2,sizeof(double));
    priors->re_var_l = (double *)calloc(2,sizeof(double));
	
	priors->fe_mean_f = (double *)calloc(2,sizeof(double));
    priors->fe_mean_l = (double *)calloc(2,sizeof(double));

    priors->fe_precision = (double **)calloc(2,sizeof(double *));
    for (i=0;i<2;i++)
    priors->fe_precision[i] = (double *)calloc(2,sizeof(double));

	 priors->re_var = (double **)calloc(2,sizeof(double *));
    for (i=0;i<2;i++)
    priors->re_var[i] = (double *)calloc(2,sizeof(double));
    
	hyper->prec = (double **)calloc(2,sizeof(double *));
    for (i=0;i<2;i++)
    hyper->prec[i] = (double *)calloc(2,sizeof(double));
	hyper->sig_a_inv = (double **)calloc(2, sizeof(double *));
	for (i = 0; i<2; i++)
		hyper->sig_a_inv[i] = (double *)calloc(2, sizeof(double));

		temp2 = (double **)calloc(2,sizeof(double *));
    for (i=0;i<2;i++)
    temp2[i] = (double *)calloc(2,sizeof(double));
	hyper->sig_a_inv[0][0] = 0.01;
	hyper->sig_a_inv[1][1] = 0.01;
	hyper->sig_a_inv[1][0] = hyper->sig_a_inv[1][0]= 0.0;

    fscanf(finput,"%lf %lf %lf %lf\n", &mprior1, &mprior2, &mprior3, &mprior4); /*fourth line: lh pulse mass mean, fsh pulse mass mean, lh pulse width mean, fsh pulse width mean*/
    hyper->hmean_l[0] = mprior1;
    hyper->hmean_f[0] = mprior2;
    hyper->hmean_l[1] = mprior3;
    hyper->hmean_f[1] = mprior4;

	fscanf(finput, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &psprior1, &psprior2, &psprior3, &psprior4, &psprior5, &psprior6, &psprior7, &psprior8, &psprior9,&psprior10, &psprior11); /*5th line: first 3 is var-cor for pulse mass variance, 4-5 is for pulse width */
    
	temp2[0][0] = psprior1;
    temp2[1][1]= psprior2;
	temp2[1][0]= temp2[0][1]= sqrt(psprior2)*sqrt(psprior1)*psprior3;
		  if (!cholesky_decomp(temp2,2)) {
     printf("prior variance matrix not PSD matrix\n");
     exit(0);
   }
      hyper->prec = cholesky_invert(2, temp2);

    hyper->precision_wl=psprior4; /*prior for population variance of pulse width */
    hyper->precision_wf=psprior5;

	hyper->hvar_l=psprior6;    /*variance term for population mean of puse width*/
	hyper->hvar_f=psprior7;
	priors->re_var_l[0]=psprior8;  /*prioor for subject variance */
	priors->re_var_l[1]=psprior9;
	priors->re_var_f[0]=psprior10;
	priors->re_var_f[1]=psprior11;


	fscanf(finput, "%lf %lf %lf %lf %lf %lf\n", &bprior1, &bprior2, &bprior3, &bprior4, &bprior5, &bprior6); /*6th baseline*/
    hyper->meanmeanbh_l[0] = bprior1;
    hyper->meanvarbh_l[0] = bprior2;
    hyper->varbh_l[0] = bprior3;
	hyper->meanmeanbh_f[0] = bprior4;
    hyper->meanvarbh_f[0] = bprior5;
    hyper->varbh_f[0] = bprior6;

    fscanf(finput,"%lf %lf %lf %lf %lf %lf\n", &hprior1, &hprior2, &hprior3,&hprior4, &hprior5, &hprior6); /*7th halflife*/
    hyper->meanmeanbh_l[1] = hprior1;
    hyper->meanvarbh_l[1] = hprior2;
    hyper->varbh_l[1] = hprior3;
	 hyper->meanmeanbh_f[1] = hprior4;
    hyper->meanvarbh_f[1] = hprior5;
    hyper->varbh_f[1] = hprior6;

    fscanf(finput,"%lf %lf %lf %lf\n", &sprior1, &sprior2, &sprior3, &sprior4); /*8th error*/
    priors->alpha_l = sprior1; 
    priors->beta_l= sprior2;
	 priors->alpha_f = sprior3;
    priors->beta_f = sprior4;
	
	/*fscanf(finput,"%d %d\n", &nprior1, &nprior2);*/ /*9th number of pulse*/
/*	parms_l->nprior = nprior1;
	parms_f->nprior = nprior2;*/
	priors->rho_prior[0] = -0.105;/*cluster size*/
	priors->rho_prior[1] = 1;
	priors->nu_prior[0] = 1.5;
	priors->nu_prior[1] = 10;
	parms_l->Rinter1 = 36;
	parms_l->Rinter2=96;



 /***********saved for future use **************/
/*********************************************/ 
  /*input prior for rho and nu here*/
/*************************************************/
    fscanf(finput,"%lf %lf %lf %lf\n", &mstart1, &mstart2, &mstart3,&mstart4); /*10th line, starting value for mean pulse features*/
    priors->fe_mean_l[0] = mstart1;
	priors->fe_mean_f[0]=mstart2;
	priors->fe_mean_l[1] = mstart3;
	priors->fe_mean_f[1]=mstart4;
  
    
    fscanf(finput,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &pstart1, &pstart2, &pstart3,&pstart4,&pstart5,&pstart6,&pstart7,&pstart8,&pstart9 );/*11th line starting value for variance pulse features*/
    /*first three is the variance-covariance for population variance and cor of pulse mass*/
	/*next two (4-5) is the population variance for pulse width */
	/*next four (6-9) is the subject-level variance for pulse mass and width*/
	priors->re_var[0][0] = pstart1;
	priors->re_var[1][1] = pstart2;
	priors->re_var[0][1]=priors->re_var[1][0] = sqrt(pstart2)*sqrt(pstart1)*pstart3;
	  if (!cholesky_decomp(priors->re_var,2)) {
     printf("proposal variance matrix not PSD matrix\n");
     exit(0);
   }
   priors->fe_precision = cholesky_invert(2,priors->re_var);
   priors->fe_precision_wl = pstart4;
   priors->fe_precision_wf = pstart5;  
   parms_l->re_precision[0]= pstart6;
   parms_l->re_precision[1]= pstart7;
   parms_f->re_precision[0]= pstart8;
   parms_f->re_precision[1]= pstart9;

   
   

   
    fscanf(finput,"%lf %lf %lf %lf\n", &bstart1, &bstart2,&bstart3, &bstart4); /*12th line baseline*/
    priors->meanbh_l[0] = bstart1;
    priors->varbh_l[0] = bstart2;
	priors->meanbh_f[0] = bstart3;
    priors->varbh_f[0] = bstart4;
    
    fscanf(finput,"%lf %lf %lf %lf\n", &hstart1, &hstart2,&hstart3, &hstart4); /*13th line half life */
    priors->meanbh_l[1] = hstart1;
    priors->varbh_l[1] = hstart2;
	priors->meanbh_f[1] = hstart3;
    priors->varbh_f[1] = hstart4;

    fscanf(finput,"%lf %lf\n", &sstart1, &sstart2);  /*14th line error*/
    parms_l->sigma = sstart1;
    parms_l->lsigma = log(parms_l->sigma);
	parms_f->sigma = sstart2;
    parms_f->lsigma = log(parms_f->sigma);
	parms_l -> rho= 1;
    parms_l ->lrho = log(parms_l->rho);
	
	parms_l ->nu= 4;
	parms_l -> lnu=log(parms_l ->nu);
	parms_l->gamma = 0.1;
	parms_l->beta =4;
	


    parms_l->numsub =parms_f->numsub = subj;

    fscanf(finput,"%lf %lf \n", &propvar[0], &propvar[1]);  /*15th lineï¼?prop for population width variance */
    fscanf(finput,"%lf %lf %lf %lf\n", &propvar[2],&propvar[3], &propvar[4], &propvar[5]); /*16 th line: prop for four variance term of population md*/
    fscanf(finput,"%lf %lf %lf\n", &propvar[6], &propvar[7],&propvar[8]); /*17 th line: prop for subject mass, third is the cor*/
    fscanf(finput,"%lf %lf\n", &propvar[9], &propvar[10]); /*18th line: uwi*/
    fscanf(finput,"%lf %lf %lf %lf \n", &propvar[11],&propvar[12],&propvar[13],&propvar[14]); /*19th line v*/
	fscanf(finput,"%lf %lf %lf %lf %lf %lf\n", &propvar[15],&propvar[16],&propvar[17],&propvar[18],&propvar[19],&propvar[20]); /*20th line bi, hi,3rd and 6th is the cor*/
    fscanf(finput,"%lf %lf %lf %lf\n", &propvar[21],&propvar[22],&propvar[23],&propvar[24]); /*21th line v, random effect, lh first and fsh*/
    fscanf(finput,"%lf %lf \n", &propvar[25],&propvar[26]); /*22th line v, time, lh first and fsh*/
	fscanf(finput, "%lf %lf %lf %lf  \n", &propvar[27], &propvar[28], &propvar[29], &propvar[30]); /*23th line v, eta, lh first and fsh*/
	fscanf(finput, "%lf %lf \n", &propvar[31], &propvar[32]); /*24th line rho and niu*/

    fclose(finput);

     /*Create the list of subjects*/
      sublist = initialize_subject();
          
      i=0;
      for(i=0;i<subj;i++){
          subject=initialize_subject();
          subject->common_l = (char *)calloc(30,sizeof(char *));
          subject->pulse_l = (char *)calloc(30,sizeof(char *));
          subject->common_f = (char *)calloc(30,sizeof(char *));
          subject->pulse_f = (char *)calloc(30,sizeof(char *));

          sprintf(tmp,"%d",subj-i);

          strcpy(subject->common_l,commonl);
          strcat(subject->common_l,"s");
          strcat(subject->common_l,tmp);
          strcat(subject->common_l,".out");
		  
		  strcpy(subject->common_f,commonf);
          strcat(subject->common_f,"s");
          strcat(subject->common_f,tmp);
          strcat(subject->common_f,".out");
          
          strcpy(subject->pulse_l,parml);
          strcat(subject->pulse_l,"s");
          strcat(subject->pulse_l,tmp);
          strcat(subject->pulse_l,".out");

		  strcpy(subject->pulse_f,parmf);
          strcat(subject->pulse_f,"s");
          strcat(subject->pulse_f,tmp);
          strcat(subject->pulse_f,".out");
		  
          subject->basehalf_l[0]=bstart1;
		  subject->basehalf_f[0]=bstart3;
          subject->basehalf_l[1]=hstart1;
		  subject->basehalf_l[1]=hstart3;


          subject->theta_l[0]=mstart1;
		  subject->theta_f[0]=mstart2;
          subject->theta_l[1]=mstart3;
		  subject->theta_f[1]=mstart4;

          insert_subject(subject,sublist);
      }
  
      fflush(stdout);

      mcmc(sublist,parms_l,parms_f,ts_l,ts_f,iter,*N,priors,seed,commonl,commonf,propvar,hyper);
  
    destroy_sublist(sublist);
    /**************************/
    
    /* save the current random number as the seed for the next simulation */
    fseed = fopen("seed.dat","w");
    fprintf(fseed,"%lu %lu %lu\n",seed[0],seed[1],seed[2]);
    fclose(fseed);
    /**********************************************************************/

    /* deallocate resources */

    free(seed);
    for (i=0;i<*N;i++)
      free(ts_f[i]);
    free(ts_f);
	 for (i=0;i<*N;i++)
      free(ts_l[i]);
    free(ts_l);
	
    free(N);
    free(priors->re_var_f);
    free(priors->re_var_l);
	free(priors->fe_mean_l);
    free(priors->fe_mean_f);
	free(parms_f->re_precision);
    free(parms_l->re_precision);
	
	 for (i=0;i<2;i++)
      free(hyper->prec[i]);
    free(hyper->prec);

	 for (i=0;i<2;i++)
      free(priors->fe_precision[i]);
    free(priors->fe_precision);
	
	for (i=0;i<2;i++)
    free(priors->re_var[i]);
 free(priors->re_var);
	

   

    free(hyper);
    free(priors);
    free(parms_l);
	free(parms_f);
	free(temp2);
/************************/
  }
  return 0;
}


