/*******************************************************************/
/*********************MAIN BDMCMC PROGRAM***************************/
/*******************************************************************/

#include "deconvolution_main.h"
#define _USE_MATH_DEFINES

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
    rnorm: Found in randgen.c; draws from the normal distribution
    rgamma: Found in randgen.c; draws from the gamma distribution
    destroy_list: Found in hash.c; frees all memory used by the nodes
    mcmc: Found in linklistv2.c; This runs birth-death mcmc routine
    *start_position: Found in this file; discussed above 
 ************************************************************************/

int main(int argc,char *argv[])
{
  int i,*N,iter,a,subj;
  char datafile[100];
  char common1[100];
  char parm1[100];
  char tmp[5];
  unsigned long *seed;
  double **ts,propvar[15];
  double mprior1, mprior2, mprior3, mprior4;
  double wprior1, wprior2, wprior3, wprior4;
  double bprior1, bprior2, bprior3;
  double hprior1, hprior2, hprior3;
  double sprior1, sprior2;
  double nprior;
  double mstart1, mstart2, mstart3;
  double wstart1, wstart2, wstart3;
  double bstart1, bstart2;
  double hstart1, hstart2;
  double sstart1;

  FILE *fseed, *finput;
  Common_parms *parms;
  Priors *priors;
  Hyper_priors *hyper;
  Subject_type *subject, *sublist;

  double **read_data_file(char *,int *,int);
  double rnorm(double,double,unsigned long *);
  double rgamma(double,double,unsigned long *);
  void destroy_sublist(Subject_type *);
  Subject_type *initialize_subject(void);
  void insert_subject(Subject_type *,Subject_type *);
  void mcmc(Subject_type *,Common_parms *,double **,long,int,Priors *,unsigned long *,
          char *,double [], Hyper_priors *);
  char *itoa(int , char *, int );

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

    fscanf(finput,"%s \n", datafile);
    fscanf(finput,"%s %s \n", common1, parm1);
    fscanf(finput,"%d %d \n", &subj, &iter);

    /* read in the hormonal time series */
    N = (int *)calloc(1,sizeof(int));
    ts = read_data_file(datafile,N,subj);
    mmm = 3;

     
    /* perform the simulation */
    
    parms = (Common_parms *)calloc(1,sizeof(Common_parms));

    fitend = ts[*N-1][0]+ ts[0][0] * 2;  /*search 2 units farther in time*/
    fitstart = -ts[0][0] * 4;  /*search 4 units in the past*/
  
    priors = (Priors *)calloc(1,sizeof(Priors));
    hyper = (Hyper_priors *)calloc(1,sizeof(Hyper_priors));
    
    priors->re_var = (double *)calloc(2,sizeof(double));
    
    fscanf(finput,"%lf %lf %lf %lf\n", &mprior1, &mprior2, &mprior3, &mprior4);
    hyper->hmean[0] = mprior1;
    hyper->hvar[0] = mprior2;
    priors->re_var[0] = mprior3;
    hyper->prec[0] = mprior4;

    fscanf(finput,"%lf %lf %lf %lf\n", &wprior1, &wprior2, &wprior3, &wprior4);
    hyper->hmean[1] = wprior1;
    hyper->hvar[1] = wprior2;
    priors->re_var[1] = wprior3;
    hyper->prec[1] = wprior4;

    fscanf(finput,"%lf %lf %lf\n", &bprior1, &bprior2, &bprior3);
    hyper->meanmeanbh[0] = bprior1;
    hyper->meanvarbh[0] = bprior2;
    hyper->varbh[0] = bprior3;

    fscanf(finput,"%lf %lf %lf\n", &hprior1, &hprior2, &hprior3);
    hyper->meanmeanbh[1] = hprior1;
    hyper->meanvarbh[1] = hprior2;
    hyper->varbh[1] = hprior3;

    fscanf(finput,"%lf %lf\n", &sprior1, &sprior2);
    priors->alpha = sprior1;
    priors->beta = sprior2;

    fscanf(finput,"%lf\n", &nprior);
    parms->nprior = nprior;

    priors->fe_mean = (double *)calloc(2,sizeof(double));
    parms->re_precision = (double *)calloc(2,sizeof(double));
    priors->fe_precision = (double *)calloc(2,sizeof(double));

    fscanf(finput,"%lf %lf %lf\n", &mstart1, &mstart2, &mstart3);
    priors->fe_mean[0] = mstart1;
    parms->re_precision[0] = mstart2;
    priors->fe_precision[0] = mstart3;
    
    fscanf(finput,"%lf %lf %lf\n", &wstart1, &wstart2, &wstart3);
    priors->fe_mean[1] = wstart1;
    parms->re_precision[1] = wstart2;
    priors->fe_precision[1] = wstart3;
    
    fscanf(finput,"%lf %lf\n", &bstart1, &bstart2);
    priors->meanbh[0] = bstart1;
    priors->varbh[0] = bstart2;
    
    fscanf(finput,"%lf %lf\n", &hstart1, &hstart2);
    priors->meanbh[1] = hstart1;
    priors->varbh[1] = hstart2;

    fscanf(finput,"%lf\n", &sstart1);
    parms->sigma = sstart1;
    parms->lsigma = log(parms->sigma);
	

    parms->numsub = subj;

    fscanf(finput,"%lf %lf %lf\n", &propvar[0], &propvar[1], &propvar[2]);
    fscanf(finput,"%lf %lf %lf\n", &propvar[3], &propvar[4], &propvar[5]);
    fscanf(finput,"%lf %lf\n", &propvar[6], &propvar[7]);
    fscanf(finput,"%lf %lf\n", &propvar[8], &propvar[9]);
    fscanf(finput,"%lf\n", &propvar[10]);
	fscanf(finput, "%lf %lf %lf %lf\n", &propvar[11], &propvar[12], &propvar[13], &propvar[14]);

    fclose(finput);

     /*Create the list of subjects*/
      sublist = initialize_subject();
          
      i=0;
      for(i=0;i<subj;i++){
          subject=initialize_subject();
          subject->common = (char *)calloc(30,sizeof(char *));
          subject->pulse = (char *)calloc(30,sizeof(char *));
		
          sprintf(tmp,"%d",subj-i);

          strcpy(subject->common,common1);
          strcat(subject->common,"s");
          strcat(subject->common,tmp);
          strcat(subject->common,".out");
          
          strcpy(subject->pulse,parm1);
          strcat(subject->pulse,"s");
          strcat(subject->pulse,tmp);
          strcat(subject->pulse,".out");

          subject->basehalf[0]=bstart1;
          subject->basehalf[1]=hstart1;

          subject->theta[0]=mstart1;
          subject->theta[1]=wstart1;

          insert_subject(subject,sublist);
      }
  
      fflush(stdout);

      mcmc(sublist,parms,ts,iter,*N,priors,seed,common1,propvar,hyper);
  
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
      free(ts[i]);
    free(ts);
    free(N);
    free(priors->fe_mean);

    free(priors->fe_precision);

    free(parms->re_precision);
    free(priors->re_var);

    free(hyper);
    free(priors);
    free(parms);
/************************/
  }
  return 0;
}


