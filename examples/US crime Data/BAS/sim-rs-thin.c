/*gcc -g -O2 -o sim-rs-thin.exe sim-rs-thin.c -llapack -lblas -lm */

/* sim-rs-thin.c */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include <time.h>
#include <limits.h>
#include "arrays.hpp" /*dynamic memory allocation routines */
#define pi M_PI	/* math.h includes important consts, renamed here */
#define clocks_per_second <integer constant expression>0>



main()
{
  FILE *fpt;         // define a pointer to a file for reading*/
  FILE *fmodelmat;   // define a pointer to a file for writing the array modelmat
  int i,j,jseed,k,l,nrow = 47, p=15, ncolumn,  dimgam, prop, itno = p*3276, g=47, *indices, ok, *varin, *varout ;               
  double Xfull[nrow][p+1], XTfull[p+1][nrow], XTXfull[p+1][p+1], y[nrow], XTyfull[p+1];
  double *XTXg, *XTyg;
  double meany,TSS,SSRgamma,Rsqgamma, logMHratio,logmarggammaold,logmarggammanew,dif;
  int modelmat[p], tempvec[p]; //counter[(int)(pow(2,p))], totalunique, modelunique[(int)(pow(2,p))], tempint;
  time_t start, end;
  struct timeval start_time, end_time;
  double total_usecs;
 
  char filename[50];
 
  
 
for (jseed=1; jseed<101; jseed++)
   {
 srand(jseed); 
 sprintf (filename,"sim-rs-thin.%d.dat",jseed);
 //printf ("Filename is %s \n",filename);
 fmodelmat =   fopen(filename,"w");
 
  /* First, call gettimeofday() to get start time */
   gettimeofday(&start_time, (struct timeval*)0);


fpt = fopen("simcen-x.txt","r");  // open file for reading only
    for (i=0; i<nrow; i++)
       {
	  for (j=0; j<p; j++)
	    {	      
               fscanf(fpt,"%lf",& Xfull[i][j+1]);
            }
         }
    fclose(fpt);                             // close the data file  

   fpt = fopen("simcen-y.txt","r");  // open file for reading only
    for (i=0; i<nrow; i++)
       {	      
          fscanf(fpt,"%lf",& y[i]);
       }
    fclose(fpt);                             // close the data file  

    for (k=0; k<nrow; k++)   {Xfull[k][0] = 1; }
  
   
     
      for (j=0; j<(p+1); j++)
        {
	  for(k=0; k<nrow; k++) {XTfull[j][k] = Xfull[k][j];}	     // 1. find XTfull
        }

     for (j=0; j<(p+1); j++)
        {
	for (k=0; k<(p+1); k++)
	  {
	     XTXfull[j][k] = 0.0;
	     for(l=0; l<nrow; l++)
	      {
                XTXfull[j][k] += XTfull[j][l]*Xfull[l][k];
		                                              // 2. find XTXfull
	      }
	  }
      }

       for (j=0; j<(p+1); j++)
      {
          XTyfull[j] = 0.0;
	  for (k=0; k<nrow; k++)
	   {
	     XTyfull[j] += XTfull[j][k]*y[k];            // 3. find XTyfull
	   }
      }
          double sumy = 0.0;
          for (j=0; j<nrow; j++) {sumy += y[j];}
          meany = sumy/(double)(nrow);
          double sumysq = 0.0;
          for (j=0; j<nrow; j++) {sumysq += y[j]*y[j];}
          TSS = sumysq - (double)(nrow)*(meany*meany); 
	  // printf("Total Sum of squares is %lf \n",TSS);
      

	  

 for (j=0; j<p; j++)
  {
    modelmat[j] = 0 ;// initializing modelmat
  }

  for (j=0; j<p; j++)
  {
    fprintf(fmodelmat,"%d \t",modelmat[j]) ;
  }
 
 logmarggammaold = 0.0; 
 fprintf(fmodelmat,"%lf \t",logmarggammaold);
    fprintf(fmodelmat,"\n");
   
 for (i=1; i<itno; i++)       // starting iterations for MH Algorithm
   {
      dimgam = 0;
      for (j=0; j<p; j++) {dimgam += modelmat[j];} // dimension of current model 
    
      for (j=0; j<p; j++)    // tempvec initially is the model from the previous iteration
       {
	   tempvec[j] = modelmat[j];
       }

 if (dimgam==0) {prop = 0;} else if (dimgam==p)  
        {prop = 0;} else
        {prop = (rand())%2;} // prop: picks randomly from the i)add/delete or ii)swap proposal
 // printf("prop is %d \n",prop);
       
      if (prop==0)      
        {
            int  index = (rand())%p; // index :corresponds to the indicator gamma_j chosen by the proposal q in MH
            tempvec[index] = 1 - tempvec[index]; // now tempvec is changed by one bit as in MC^3 for the proposed move
        } else if (prop==1)
      
        {
          varin = ivector(dimgam) ;//vector of included variables in current model (don't consider intercept)
     
      j = 0; k = 0;
      while (j<p && k<dimgam)
        {
             if (modelmat[j]==1) {varin[k] = j; k++;}
	       j++ ;
        }

    varout = ivector(p-dimgam) ;//vector of excluded variables in current model (don't consider intercept)
     
      j = 0; k = 0;
      while (j<p && k<(p-dimgam))
        {
             if (modelmat[j]==0) {varout[k] = j; k++;}
	       j++ ;
        }

          int swapin = (rand())%dimgam;    // swapin :corresponds to position of randomly chosen included variable
          int swapout = (rand())%(p-dimgam);// swapout :corresponds to position of randomly chosen excluded variable
          tempvec[varin[swapin]] = 0;
          tempvec[varout[swapout]] =1; 
	}

     
      ncolumn = 0;
      for (j=0; j<p; j++) {ncolumn += tempvec[j];}  
      // printf("Dimesion of model is %d \n",ncolumn);
      
     indices = ivector(ncolumn) ;//indices  indicates the position of the nonzero gamma_j's, it always has 0 as the first argument
     
      j = 0; k = 0;
      while (j<p && k<ncolumn)
        {
             if (tempvec[j]==1) {indices[k] = j+1; k++;}
	       j++ ;
        }
       // need to compute R^2_gamma and eventually the marg lik under g-prior

      // First trying to compute betahat_gamma using lapack
       // 5. Compute Rsqgamma as follows:
      //                           a. Transpose (Xgamma^T*y ) to get y^T* Xgamma
      //                           b. Multiply  y^T* Xgamma by betagammahat
      // 6. Calculate marggamma (marginal likelihood under model gamma)
       

       XTXg = vector(ncolumn*ncolumn); XTyg = vector(ncolumn);  // allocate memory for XTXg,XTyg
       for (j=0; j<(ncolumn); j++)
        {
	  XTyg[j] = XTyfull[indices[j]];  
	  //printf("XTyg[%d] is %lf \n",j,XTyg[j]);
	 for (k=0; k<(ncolumn); k++)
	  {
	     XTXg[j*ncolumn+k] = XTXfull[indices[k]][indices[j]];
	  }
      }
      	  
       int c2 = 1; 
       int *pivot;
       
     
       if (ncolumn == 0)  
             { 
	       logmarggammanew = 0.0;
             } else
	 {
       pivot = ivector(ncolumn);  //allocate memory
       //printf("%d \n",ncolumn);
       dgesv_(&ncolumn, &c2, XTXg, &ncolumn, pivot, XTyg, &ncolumn, &ok); // replaces XTyg by the soln i.e.(XTXg)^(-1)*(XTyg)
       // 5b. Multiply yTX by betagammahat to get SSRgamma(SS due to Regression)
	 
          SSRgamma = 0.0; for (j=0; j<(ncolumn); j++){SSRgamma += XTyg[j]*XTyfull[indices[j]];}
	   Rsqgamma = SSRgamma/ TSS;
	   
	  logmarggammanew = .5*((double)(log(1 + g)) * (double) (nrow - ncolumn - 1)
				- log(1.0 + (double)(g)*(1.0 - Rsqgamma)) * (double)(nrow-1));
              freevector(pivot); 
            }  
          // MH step
       if (dimgam==0 | dimgam==p) {logMHratio = (double)log(0.5) + logmarggammanew - logmarggammaold;} 
       else if (ncolumn==0 | ncolumn==p)  {logMHratio = (double)log(2.0) + logmarggammanew - logmarggammaold;} 
       else {logMHratio = logmarggammanew - logmarggammaold;} 
	 
	 
          double randnum = (double)(rand())/RAND_MAX;
            if (log(randnum) <= logMHratio)
             {
                for (j=0; j<p; j++)
                 {
		   modelmat[j] = tempvec[j];  // accepting the move with prob MHratio
                 }
               logmarggammaold = logmarggammanew;
             }

	    if(i%p == 0) // storing every pth
	      {
            for (j=0; j<p; j++)
             {
                fprintf(fmodelmat,"%d \t",modelmat[j]);
             }
            fprintf(fmodelmat,"%lf \t",logmarggammaold);
           
            fprintf(fmodelmat,"\n") ;
	    if (i%(p*10000)==0){printf("%d \n",i/p);}
	       }
	    
          
	    //printf("%d \n",bintoint(p,modelmat));
            
             freevector(varin);  
             freevector(varout);  
             freevector(indices);                        //   free memory
             freevector(XTXg);                           
             freevector(XTyg);                            
	  ////////////////////////////////////////////////////////////////////////////
  
   }    // end of the MCMC i-iterations loop for the MH Algorithm
  
 

   /* Now call gettimeofday() to get end time */
   gettimeofday(&end_time, (struct timeval*)0);  /* after time */
   
  
  /* Print the execution time */
   total_usecs = (end_time.tv_sec-start_time.tv_sec)  + (end_time.tv_usec-start_time.tv_usec)/1000000.0;

   //printf("Total time was %lf Sec.\n", total_usecs);
fclose(fmodelmat); // finished writing the file
printf("Finished replicate %d \n",jseed) ;   
} // end jseed
 
}



  
