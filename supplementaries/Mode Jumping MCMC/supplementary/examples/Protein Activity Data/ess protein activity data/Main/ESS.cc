/* This file is part of ESS++.
 *      Copyright (c) Marc Chadeau (m.chadeau@imperial.ac.uk)
 *                    Leonardo Bottolo (l.bottolo@imperial.ac.uk)
 *                    David Hastie (d.hastie@imperial.ac.uk)
 *      2010
 *
 * ESS++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESS++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ESS++.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <iostream>
#include <ctime>
#include "xml_file_write.h"
#include "xml_file_read.h"
#include <vector>
#include <algorithm>
#include <memory>
#include <fstream>
#include <../Routines/struc.h>
#include <../Routines/dyn_name.h>
#include <../Routines/matrix_handling.h>
#include <../Routines/rand.h>
#include <../Routines/moves.h>
#include <../Routines/regression.h>
#include <../Routines/cond_post.h>
#include <../Routines/post_processing.h>
#include <../Routines/xml_file_read.h>
#include <../Classes/String_Matrices.h>
#include <../Classes/Int_Matrices.h>
#include <../Classes/Int_Matrices_var_dim.h>
#include <../Classes/Double_Matrices.h>
#include <../Classes/Double_Matrices_cont.h>
#include <../Classes/Prior_param.h>
#include <../Classes/Temperatures.h>
#include <../Classes/Move_monitor.h>
#include <../Classes/g_AdMH.h>
#include <../Classes/DR.h>
#include <../Classes/CM.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#define DEBUG 0
#define Missing_data -1
using namespace std;

//******************************************************************
//*main
//******************************************************************

int main(int argc, char *  argv[])
{

  std::clock_t startTime = std::clock();
  std::clock_t tmpTime,endTime;
  double setupTime,mainLoopTime,postProcessTime;

  char filename_in_mat_X[1000];
  char filename_in_mat_Y[1000];
  char filename_par[1000];

  char path_name_out[1000];
  int na=0;
  long MY_SEED=-1;
  unsigned int n_sweeps=0;
  unsigned int burn_in=0;
  unsigned int n_top_models_from_read=0;
  // If g_opt=1 we are using g-priors otherwise indep priors
  // In this version g_opt set to 1 here and no option to change
  int g_opt=1;
  int g_sample=1;
  double g_init=1.0; 

  bool HistoryFlag=false;
  bool Time_monitorFlag=false;
  bool X_Flag=false;
  bool Y_Flag=false;
  bool Par_Flag=false;
  bool Out_full_Flag=false;
  bool Log_Flag=false;
  bool iso_T_Flag=false;
  // This is not changed from false in this version
  bool cudaFlag=false;
  na++;
  while(na < argc)
    {
      if ( 0 == strcmp(argv[na],"-X") )
	{
	  X_Flag=true;
	  strcpy(filename_in_mat_X,argv[++na]);
	  if (na+1==argc) break;
	  na++;
	}
      else if ( 0 == strcmp(argv[na],"-Y") )
	{
	  Y_Flag=true;
	  strcpy(filename_in_mat_Y,argv[++na]);
	  if (na+1==argc) break;
	  na++;
	}
      else if ( 0 == strcmp(argv[na],"-par") )
	{
	  Par_Flag=true;
	  strcpy(filename_par,argv[++na]);
	  if (na+1==argc) break;
	  na++;
	}
      else if ( 0 == strcmp(argv[na],"-nsweep") )
	{
	  n_sweeps=atoi(argv[++na]);
	  if (na+1==argc) break;
	  na++;
	}
      else if ( 0 == strcmp(argv[na],"-burn_in") )
	{
	  burn_in=atoi(argv[++na]);
	  if (na+1==argc) break;
	  na++;
	}
      else if ( 0 == strcmp(argv[na],"-seed") ){
	MY_SEED=(long)((atoi(argv[++na])));
	if (na+1==argc) break;
	na++;
      }
      else if ( 0 == strcmp(argv[na],"-out") )
	{
	  strcpy(path_name_out,argv[++na]);
	  if (na+1==argc) break;
	  na++;
	}
      else if ( 0 == strcmp(argv[na],"-out_full") )
	{
	  Out_full_Flag=true;
	  strcpy(path_name_out,argv[++na]);
	  if (na+1==argc) break;
	  na++;
	}
      else if ( 0 == strcmp(argv[na],"-history") )
	{
	  HistoryFlag=true;
	  if (na+1==argc) break;
	  na++;
	}
      else if ( 0 == strcmp(argv[na],"-time") )
	{
	  Time_monitorFlag=true;
	  if (na+1==argc) break;
	  na++;
	}
      else if ( 0 == strcmp(argv[na],"-top") )
	{
	  n_top_models_from_read=atoi(argv[++na]);
	  if (na+1==argc) break;
	  na++;
	}   
      else if ( 0 == strcmp(argv[na],"-log") )
	{
	  Log_Flag=true;
	  na++;
	}
      else if ( 0 == strcmp(argv[na],"-isoT") )
	{
	  iso_T_Flag=true;
	  if (na+1==argc) break;
	  na++;
	}
      else if ( 0 == strcmp(argv[na],"-g_set") )
	{
	  g_sample=0;
	  g_init=(double)(atof(argv[++na]));
	  if (na+1==argc) break;
	  na++;
	}
      else
        {
          cout << "Unknown option: " << argv[na] << endl;
	  exit(1);
        }
  }
  
  if(!X_Flag){
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	 << "   The predictor matrix X has not been specified, RUN STOPPED" << endl
	 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    exit(1);
  }

  if(!Y_Flag){
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	 << "   The outcome matrix Y has not been specified, RUN STOPPED" << endl
	 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    exit(1);
  }

  if(!Par_Flag){
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	 << "   The parameters matrix has not been specified, RUN STOPPED" << endl
	 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    exit(1);
  }

  if(MY_SEED<0){
    MY_SEED=(long)time(0);
  }

  smyrand((long)(MY_SEED));

 //Reading X matrix.
  Double_Matrices mat_X;
  mat_X.Read_from_file(filename_in_mat_X);
  unsigned int nX=mat_X.nb_rows;// # observations
  unsigned int pX=mat_X.nb_columns;// # SNPs
 
 
  //Reading Y matrix.
  Double_Matrices mat_Y;

  mat_Y.Read_from_file(filename_in_mat_Y);
  unsigned int pY=mat_Y.nb_columns;// # outcomes

  /////Testing the number of variables in Y
  if(nX<pY){
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	 << "There are too many outcomes ("<< pY
	 << ") compared to the number of observations (" << nX
	 << "), run stopped" << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    exit(0);
  }


  //////////////////////////////////
  //  Running options
  //////////////////////////////////
  
  if(n_sweeps==0 || burn_in==0){
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	 << "The Number of sweeps and/or the burn-in has not been specified" << endl
	 << "Use -iter and/or -burn-in option(s) in the command line" << endl
	 << "Run stopped" << endl
	 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    exit(1);
  }
  if(n_sweeps <= burn_in){
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	 << "The Number of sweeps: " << n_sweeps << " is lower than " << endl
	 << "(or equal to) the burn-in: " << burn_in << " -- Run stopped" << endl
	 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    exit(1);
  }
   //Setting up modelling parameters.
 
  double E_p_gam_from_read=2;
  double Sd_p_gam_from_read=1;
  unsigned int nb_chains=3;
  unsigned int n_top_models=n_sweeps;
  double g=0;
  if(g_sample==0){
    g=g_init;
  }
  else{
    g=pow((double)(pX),2);
  }
  double P_mutation_from_read=0.5;
  double P_sel_from_read=0.5;
  double P_csvr_r_from_read=0.375;
  double P_DR_from_read=0.5;

  //Reading the parameter File
  FILE *fparameter=NULL;
  char str[256]; // string used to read parameters

  fparameter = fopen(filename_par,"r");
  MaXmlTagRead XML_PAR(fparameter); // assign class
  //Reading elements of the par file
  if(XML_PAR.ReadTag("E_P_GAM", 0, 0, str,256)){
    sscanf(str,"%lf",&E_p_gam_from_read);
  }
  if(XML_PAR.ReadTag("SD_P_GAM", 0, 0, str,256)){
    sscanf(str,"%lf",&Sd_p_gam_from_read);
  }
  if(XML_PAR.ReadTag("NB_CHAINS", 0, 0, str,256)){
    sscanf(str,"%u",&nb_chains);
  }
  if(XML_PAR.ReadTag("P_MUTATION", 0, 0, str,256)){
    sscanf(str,"%lf",&P_mutation_from_read);
  }
  if(XML_PAR.ReadTag("P_SEL", 0, 0, str,256)){
    sscanf(str,"%lf",&P_sel_from_read);
  }
  if(XML_PAR.ReadTag("P_CSRV_R", 0, 0, str,256)){
    sscanf(str,"%lf",&P_csvr_r_from_read);
  }
  if(XML_PAR.ReadTag("P_DR", 0, 0, str,256)){
    sscanf(str,"%lf",&P_DR_from_read);
  }

 if(n_top_models_from_read!=0){
    n_top_models=n_top_models_from_read;
  }

  //Setting up regression parameters.
  double n_Pvalue_enter=0.01;
  double n_Pvalue_remove=0.01;

  if(XML_PAR.ReadTag("N_P_VALUE_ENTER", 0, 0, str,256)){
    sscanf(str,"%lf",&n_Pvalue_enter);
  }
  if(XML_PAR.ReadTag("N_P_VALUE_REMOVE", 0, 0, str,256)){
    sscanf(str,"%lf",&n_Pvalue_remove);
  } 

  
  //Regression Setting up parameters.
  double Pvalue_enter = 1.0 - pow((1.0 - n_Pvalue_enter),(1.0/(double)(mat_X.nb_columns)));
  double Pvalue_remove = 1.0 - pow((1.0 - n_Pvalue_remove),(1.0/(double)(mat_X.nb_columns)));
  

  //Moves Parameters
  //g Adaptative M-H
  unsigned int g_n_batch_from_read=100;
  double g_AdMH_optimal_from_read=0.44;
  double  g_AdMH_ls_from_read=0.0;
  double g_M_min_input=0.0;
  double g_M_max_input=0.0;
  if(XML_PAR.ReadTag("G_N_BATCH", 0, 0, str,256)){
    sscanf(str,"%u",&g_n_batch_from_read);
  }
  if(XML_PAR.ReadTag("G_ADMH_OPTIMAL", 0, 0, str,256)){
    sscanf(str,"%lf",&g_AdMH_optimal_from_read);
  }
  if(XML_PAR.ReadTag("G_ADMH_LS", 0, 0, str,256)){
    sscanf(str,"%lf",&g_AdMH_ls_from_read);
  }
  if(XML_PAR.ReadTag("G_M_MIN", 0, 0, str,256)){
    sscanf(str,"%lf",&g_M_min_input);
  }
  if(XML_PAR.ReadTag("G_M_MAX", 0, 0, str,256)){
    sscanf(str,"%lf",&g_M_max_input);
  }
  //Crossover Move
  unsigned int k_max_from_read=2;
  if(XML_PAR.ReadTag("K_MAX", 0, 0, str,256)){
    sscanf(str,"%u",&k_max_from_read);
  }
  //Gibbs Move
  unsigned int Gibbs_n_batch=500;
  if(XML_PAR.ReadTag("GIBBS_N_BATCH", 0, 0, str,256)){
    sscanf(str,"%u",&Gibbs_n_batch);
  }

  //Temperature Placement
  double b_t_input=2.0;
  double a_t_den_inf_5k=2.0;
  double a_t_den_5_10k=4.0;
  double a_t_den_sup_10k=2.0;
  unsigned int temp_n_batch=50;
  double temp_optimal_input=0.5;
  vector < double > M_input;
  M_input.resize(2);
  double M_min_input=1.0;
  double M_max_input=4.0;

  

  if(XML_PAR.ReadTag("B_T", 0, 0, str,256)){
    sscanf(str,"%lf",&b_t_input);
  }
  if(XML_PAR.ReadTag("A_T_DEN_INF_5K", 0, 0, str,256)){
    sscanf(str,"%lf",&a_t_den_inf_5k);
  }
  if(XML_PAR.ReadTag("A_T_DEN_5_10K", 0, 0, str,256)){
    sscanf(str,"%lf",&a_t_den_5_10k);
  }
  if(XML_PAR.ReadTag("A_T_DEN_SUP_10K", 0, 0, str,256)){
    sscanf(str,"%lf",&a_t_den_sup_10k);
  }
  if(XML_PAR.ReadTag("TEMP_N_BATCH", 0, 0, str,256)){
    sscanf(str,"%u",&temp_n_batch);
  }
  if(XML_PAR.ReadTag("TEMP_OPTIMAL", 0, 0, str,256)){
    sscanf(str,"%lf",&temp_optimal_input);
  }
  if(XML_PAR.ReadTag("M_MIN", 0, 0, str,256)){
    sscanf(str,"%lf",&M_min_input);
  }
  if(XML_PAR.ReadTag("M_MAX", 0, 0, str,256)){
    sscanf(str,"%lf",&M_max_input);
  }
 
  M_input[0]=M_min_input;
  M_input[1]=M_max_input;

  fclose(fparameter);


  cout << "**********************************************************" << endl
       << "***************** Setup options **********************" << endl
       << "Random seed " << MY_SEED << endl
       << "nb_chains " << nb_chains << endl
       << "n_sweeps " << n_sweeps << endl
       << "n_top_models " << n_top_models << endl
       << "burn_in " << burn_in << endl
       << "E_p_gam_input " << E_p_gam_from_read << endl
       << "Sd_p_gam_input " << Sd_p_gam_from_read << endl
       << "g_sample " << g_sample << endl;
  if(g_sample==1){
    cout << "g " << g << endl;
  }
  cout << "**********************************************************" << endl
       << "**********************************************************" << endl << endl;


  cout << "**********************************************************" << endl
       << "*************** Regression parameters ********************" << endl
       << "n_Pvalue_enter " << n_Pvalue_enter << endl
       << "n_Pvalue_enter " << n_Pvalue_enter << endl
       << "Pvalue_enter stepwise " << Pvalue_enter << endl
       << "Pvalue_remove stepwise " << Pvalue_remove << endl
       << "**********************************************************" << endl
       << "**********************************************************" << endl << endl;
  

 
 
  //Testing PATH-file names for output
  
  string Extension_out=".txt";
  string Name_number="sweeps";
  string Main_Output_name="output_models_history";
  string OutputName=Get_stddzed_name(path_name_out,
				     n_sweeps,
				     Name_number,
				     Main_Output_name,
				     Extension_out);
  ofstream f_out;

  string OutputName_n_vars_in=Get_stddzed_name(path_name_out,
					       n_sweeps,
					       Name_number,
					       "output_model_size_history",
					       Extension_out);
  ofstream f_out_n_vars_in;

  string OutputName_log_cond_post_per_chain=Get_stddzed_name(path_name_out,
							     n_sweeps,
							     Name_number,
							     "output_log_cond_post_prob_history",
							     Extension_out);
  ofstream f_out_log_cond_post_per_chain;

  if(HistoryFlag){

    f_out.open(OutputName.c_str(),ios::out);
    if(f_out.fail()){
      cout << "Invalid Path and/or permission rights for " << OutputName << " -- run stopped." << endl;
      exit(1);
    }
    else{
      f_out << "Sweep\tModel_size\tlog_marg\tlog_cond_post\tModel"<<endl;
    }
    

    f_out_n_vars_in.open(OutputName_n_vars_in.c_str(),ios::out);
    if(f_out_n_vars_in.fail()){
      cout << "Invalid Path and/or permission rights for " << OutputName_n_vars_in << " -- run stopped." << endl;
      exit(1);
    }
    else{
      f_out_n_vars_in << "Sweep\t";
      for(unsigned int tmp_chain=0;tmp_chain<nb_chains;tmp_chain++){
	f_out_n_vars_in << "Chain_"<< tmp_chain+1 << "\t";
      }
      f_out_n_vars_in << endl;
    }
    

    f_out_log_cond_post_per_chain.open(OutputName_log_cond_post_per_chain.c_str(),ios::out);
    if(f_out_log_cond_post_per_chain.fail()){
      cout << "Invalid Path and/or permission rights for " << OutputName_log_cond_post_per_chain << " -- run stopped." << endl;
      exit(1);
    }
    else{
      f_out_log_cond_post_per_chain << "Sweep"<< "\t";
      for(unsigned int tmp_chain=0;tmp_chain<nb_chains;tmp_chain++){
	f_out_log_cond_post_per_chain << "Chain_"<< tmp_chain+1 << "\t";
      }
      f_out_log_cond_post_per_chain << endl;
    }
  }


  string OutputName_time; 
  ofstream f_out_time;
  if(Time_monitorFlag){
    OutputName_time=Get_stddzed_name(path_name_out,
					    n_sweeps,
					    Name_number,
					    "output_time_monitor",
					    Extension_out);
    f_out_time.open(OutputName_time.c_str(),ios::out);
    f_out_time << "Sweep\tTime\tTime_per_eval_model"<<endl;
  }


  string OutputName_best_models=Get_stddzed_name(path_name_out,
						 n_sweeps,
						 Name_number,
						 "output_best_visited_models",
						 Extension_out);
  ofstream f_out_best_models;
  f_out_best_models.open(OutputName_best_models.c_str(),ios::out);
  if(f_out_best_models.fail()){
    cout << "Invalid Path and/or permission rights for " << OutputName_best_models << " -- run stopped." << endl;
    exit(1);
  }
  else{
    if(Out_full_Flag){
      f_out_best_models << "Rank\t#Visits\tSweep_1st_visit\t#models_eval_before_1st_visit\tModel_size\tlog_Post_Prob\tModel_Post_Prob\tJeffreys_scale\tModel"<<endl;
    }
    else{
      f_out_best_models << "Rank\t#Visits\tModel_size\tlog_Post_Prob\tModel_Post_Prob\tJeffreys_scale\tModel"<<endl;
    }
  }
  string OutputName_marg_gam=Get_stddzed_name(path_name_out,
						 n_sweeps,
						 Name_number,
						 "output_marg_prob_incl",
						 Extension_out);
  ofstream f_out_marg_gam;
  f_out_marg_gam.open(OutputName_marg_gam.c_str(),ios::out);
  if(f_out_marg_gam.fail()){
    cout << "Invalid Path and/or permission rights for " << OutputName_marg_gam << " -- run stopped." << endl;
    exit(1);
  }
  else{
    f_out_marg_gam << "Predictor\tMarg_Prob_Incl"<< endl;
  }


  bool scaling_flag=true;


  
  
  
  Double_Matrices_cont mat_X_copy;
  mat_X_copy.Copy_from_double_matrix(mat_X);
  
  if(scaling_flag){
    standardize_matrix_cont(mat_X_copy);
  }
  
  gsl_matrix *mat_X_work=Double_matrices_cont_2_gsl_matrix(mat_X_copy);

  vector < vector <unsigned int> > Gam_step_regr;
  Gam_step_regr.resize(pY);
  Double_Matrices Gam_step_regr_pvals;
  Gam_step_regr_pvals.Alloc_double_matrix(pY,
					  pX);
  Double_Matrices Gam_step_regr_SE;
  Gam_step_regr_SE.Alloc_double_matrix(pY,
				       pX);

  Double_Matrices Gam_step_regr_beta;
  Gam_step_regr_beta.Alloc_double_matrix(pY,
					 pX);

  double tolerance= 6.6835e-14;
  gsl_vector *vect_RMSE = gsl_vector_calloc(pY);
  Double_Matrices_cont mat_Y_copy;
  mat_Y_copy.Copy_from_double_matrix(mat_Y);
  
  gsl_matrix *mat_Y_work=Double_matrices_cont_2_gsl_matrix(mat_Y_copy);
  gsl_vector *current_outcome =gsl_vector_calloc(nX);
  gsl_vector *vect_residuals=gsl_vector_calloc(nX);
  gsl_vector *vect_p_value=gsl_vector_calloc(mat_X_work->size2);
  gsl_vector *vect_beta_full=gsl_vector_calloc(mat_X_work->size2);
  gsl_vector *vect_SE_full=gsl_vector_calloc(mat_X_work->size2);

  cout << "***************************************************" << endl
       << "*************  Stepwise regression   *************" << endl
       << "***************************************************" << endl << endl;

  for(unsigned int outcome=0;outcome<pY;outcome++){ 
    vector < unsigned int > list_columns_X_gam;
    vector < unsigned int > list_columns_X_gam_bar;
    vector < unsigned int > is_var_in;
    
    is_var_in.resize(pX);
    gsl_matrix_get_col(current_outcome,
		       mat_Y_work,
		       outcome);
    int stop=0;
    int loop=0;
    int n_loop_max=100;

    while(stop==0){
      get_list_var_in_and_out(list_columns_X_gam,
			      list_columns_X_gam_bar,
			      is_var_in);
      
   
      gsl_matrix * mat_X_gam=get_X_reduced_and_constant(list_columns_X_gam,
							mat_X_work);
      
      gsl_matrix * mat_X_gam_bar=get_X_reduced(list_columns_X_gam_bar,
					       mat_X_work);
      //Calculating p-values for all variables
      get_full_pvalues_and_beta(vect_p_value,
				vect_beta_full,
				vect_SE_full,
				mat_X_gam,
				mat_X_gam_bar,
				current_outcome,
				vect_residuals,
				vect_RMSE,
				tolerance,
				list_columns_X_gam,
				list_columns_X_gam_bar,
				outcome); 
      //Updating the list of var in
      stop=update_is_var_in(is_var_in,
			    list_columns_X_gam,
			    list_columns_X_gam_bar,
			    vect_p_value,
			    Pvalue_enter,
			    Pvalue_remove,
			    loop,
			    n_loop_max);
      gsl_matrix_free(mat_X_gam);
      gsl_matrix_free(mat_X_gam_bar);
 
      loop++;
      list_columns_X_gam.clear();
      list_columns_X_gam_bar.clear();
      
      if(stop==1){
	get_list_var_in_and_out(list_columns_X_gam,
				list_columns_X_gam_bar,
				is_var_in);
      }
      
    }//end of while
    //Filling the output file


 
    store_model_per_outcome(Gam_step_regr,
			    list_columns_X_gam,
			    vect_p_value,
			    vect_beta_full,
			    vect_SE_full,
			    Gam_step_regr_pvals,
			    Gam_step_regr_SE,
			    Gam_step_regr_beta,
			    outcome);
    
    
    list_columns_X_gam.clear();
    is_var_in.clear();
  }

  cout << "Result From Step-Wise Regression" << endl;
  display_matrix_var_dim(Gam_step_regr);

  cout << endl;


  cout << "***************************************************" << endl
       << "**********  End of Stepwise regression  **********" << endl
       << "***************************************************" << endl << endl;



  gsl_vector_free(vect_residuals);
  gsl_vector_free(current_outcome);
  gsl_vector_free(vect_p_value);
  gsl_vector_free(vect_beta_full);
  gsl_vector_free(vect_SE_full);

  Gam_step_regr_pvals.Free_double_matrix();
  Gam_step_regr_SE.Free_double_matrix();
  Gam_step_regr_beta.Free_double_matrix();

  gsl_matrix_free(mat_X_work);
  gsl_matrix_free(mat_Y_work);
  mat_X_copy.Free_double_matrix_cont();
 
 cout << "**********************************************************" << endl
       << "****************** MOVES parameters **********************" << endl
       << "g-adaptative M-H" << endl
       << "\t-g_n_batch: " << g_n_batch_from_read << endl
       << "\t-g_AdMH_optimal: " << g_AdMH_optimal_from_read << endl
       << "\t-g_AdMH_ls: " << g_AdMH_ls_from_read << endl
       << "Crossover Move" << endl
       << "\tk_max: " <<k_max_from_read << endl
       << "Gibbs Move" << endl
       << "\tGibbs_n_batch: " <<Gibbs_n_batch << endl
       << "**********************************************************" << endl
       << "**********************************************************" << endl << endl;
 

//   cout << "**********************************************************" << endl
//        << "****************** TEMP parameters **********************" << endl
//        << "b_t " << b_t_input << endl
//        << "a_t_den_inf_5k " << a_t_den_inf_5k << endl
//        << "a_t_den_5_10k " << a_t_den_5_10k << endl
//        << "a_t_den_sup_10k " << a_t_den_sup_10k << endl
//        << "temp_n_batch " << temp_n_batch << endl
//        << "temp_optimal " << temp_optimal_input << endl
//        << " M= [" << M_input[0] << " - " << M_input[1] << "]" << endl
//        << "**********************************************************" << endl
//        << "**********************************************************" << endl << endl;


  //////////////////////////////////
  //  Setting up prior parameters
  //////////////////////////////////
  
  Prior_param PR;
  PR.set_PR_param(g_opt,
		  E_p_gam_from_read,
		  Sd_p_gam_from_read,
		  pX,
		  pY,
		  nX,
		  vect_RMSE,
		  P_mutation_from_read,
		  P_sel_from_read,
		  P_csvr_r_from_read,
		  P_DR_from_read);
  
  PR.display_prior_param();
  
  //////////////////////////////////
  //  Setting up g_AdMH parameters
  //////////////////////////////////
  
  g_AdMH *My_g_AdMH=new g_AdMH;
  (*My_g_AdMH).set_g_AdMH(g_sample,
			  g_n_batch_from_read,
			  g_AdMH_optimal_from_read,
			  g_AdMH_ls_from_read,
			  pX,
			  burn_in,
			  g_M_min_input,
			  g_M_max_input);
  
  (*My_g_AdMH).display_g_AdMH();
  
  
  //////////////////////////////////
  //  Setting up DR parameters
  //////////////////////////////////

  DR *My_DR=new DR();
  (*My_DR).set_DR(nb_chains);
  (*My_DR).display_DR();
  
  
  //////////////////////////////////
  //  Setting up CM parameters
  //////////////////////////////////
  vector <unsigned int > list_CM_moves_enabled_from_read;
  list_CM_moves_enabled_from_read.push_back(1);
  list_CM_moves_enabled_from_read.push_back(2);
  
  CM *My_CM=new CM;
  
  (*My_CM).set_CM(k_max_from_read,
		  list_CM_moves_enabled_from_read);
  
  (*My_CM).display_CM();
  
  list_CM_moves_enabled_from_read.clear();

  //////////////////////////////////
  //  Defining Move Monitor Object
  //////////////////////////////////

 
  Move_monitor *My_Move_monitor=new Move_monitor(nb_chains,
						 (*My_CM).n_possible_CM_moves);
  //////////////////////////////////
  //   Initializing the chains
  //////////////////////////////////
 

  vector < vector <unsigned int> > vect_gam;
  vect_gam.resize(nb_chains);
  for(unsigned int chain=0;chain<nb_chains;chain++){
    vect_gam[chain].resize(pX);
  }
  
  get_vect_gam_init(vect_gam,
		    Gam_step_regr,
		    iso_T_Flag);
  
  
  Double_Matrices_cont mat_X_copy2;
  mat_X_copy2.Copy_from_double_matrix(mat_X);
  mat_X.Free_double_matrix();
  mat_Y.Free_double_matrix();
  
  //display_matrix_var_dim(vect_gam);
  gsl_matrix *mat_X_work2=Double_matrices_cont_2_gsl_matrix(mat_X_copy2);
  gsl_matrix *mat_Y_work2=Double_matrices_cont_2_gsl_matrix(mat_Y_copy);
  
  
  
  //Centering X and Y (if g_opt!=1, stddize X as well)
  center_matrix_gsl(mat_X_work2);
  center_matrix_gsl(mat_Y_work2);
  
  Gam_step_regr.clear();


  //////////////////////////////////
  // Intializing chain temparature
  //////////////////////////////////

  Temperatures *t_tun=new Temperatures();
  (*t_tun).set_Temp_param(nb_chains,
			  pX,
			  b_t_input,
			  a_t_den_inf_5k,
			  a_t_den_5_10k,
			  a_t_den_sup_10k,
			  temp_n_batch,
			  M_input,
			  burn_in,
			  temp_optimal_input,
			  iso_T_Flag);
  M_input.clear();
  (*t_tun).display_Temp_param();

  ///////////////////////////////////
  //   Declaring output matrices
  ///////////////////////////////////

  Double_Matrices mat_log_marg;
  mat_log_marg.Alloc_double_matrix(nb_chains,
				   n_sweeps+1);
  
  Double_Matrices mat_log_cond_post;
  mat_log_cond_post.Alloc_double_matrix(nb_chains,
					n_sweeps+1);


  gsl_permutation *MyPerm=gsl_permutation_calloc(pX);
  My_Permut_unsigned_int(MyPerm);
  
  unsigned int sweep=0;

  //Defining the chain ID:
  //chain_idx[i], says on which line chain number i is located
  vector < unsigned int > chain_idx;
  intialize_chain_idx(chain_idx,
		      nb_chains);
  //description_exchange_move: each move is presented in columns
  // -line 1: first chain
  // -line 2: second chain
  // -line 3: Pbty of the move

  gsl_matrix *description_exchange_move=description_exch_moves(nb_chains);

  endTime = std::clock();
  tmpTime = startTime;
  setupTime = (endTime-tmpTime)/(double)(CLOCKS_PER_SEC);


  cout << "***************************************" << endl
       << "             FIRST SWEEP" << endl
       <<  "***************************************" << endl << endl;
  for(unsigned int chain=0;chain<nb_chains;chain++){
  
    //Initial Calculation of the logMarg and log_cond_post;
    vector < unsigned int > list_columns_X_gam;
    vector < unsigned int > list_columns_X_gam_bar;
    
    //Step 1 Getting X_gamma  
    get_list_var_in_and_out(list_columns_X_gam,
			    list_columns_X_gam_bar,
			    vect_gam[chain]);
    
    cout << "List_columns_X" << endl;
    for(unsigned int col=0;col<list_columns_X_gam.size();col++){
      cout << list_columns_X_gam[col] << " ";
    }
    cout << endl;
    unsigned int n_vars_in=list_columns_X_gam.size();
    cout << "n_vars_in " << n_vars_in << endl;


    if(n_vars_in>0){
      gsl_matrix * mat_X_gam=get_X_reduced(list_columns_X_gam,
					   mat_X_work2);
      
      //Step 2: Calculate log-marginal
      
      full_log_marg_computation(mat_log_marg,
				mat_log_cond_post,
				mat_X_gam,
				mat_Y_work2,
				PR,
				g_opt,
				g_sample,
				g,
				pX,
				chain,
				sweep,
				cudaFlag);
      gsl_matrix_free(mat_X_gam);  

    }
    else{
      full_log_marg_computation_no_var_in(mat_log_marg,
					  mat_log_cond_post,
					  mat_Y_work2,
					  PR,
					  g_opt,
					  g_sample,
					  g,
					  pX,
					  chain,
					  sweep);
      
    }
    
    
    cout << endl << "**************Results, chain " << chain << " -- sweep " << sweep << "***************" << endl;
    cout << "\tlog_marg " << mat_log_marg.matrix[chain][sweep] << endl
	 << "\tlog_cond_post " << mat_log_cond_post.matrix[chain][sweep] << endl;
    cout << "********************************************************" << endl;
    
    
    
    list_columns_X_gam.clear();
    list_columns_X_gam_bar.clear();
    
  }//end of for chain.
  
  cout << "***************************************" << endl
       << "          END OF FIRST SWEEP" << endl
       <<  "***************************************" << endl << endl;
 
  print_main_results_per_sweep(f_out,
			      vect_gam,
			      chain_idx,
			      mat_log_marg,
			      mat_log_cond_post,
			      0);

  double cum_g=0.0;
  unsigned int count_g=0.0;


  ////////////////////////////////////
  ////////////////////////////////////
  //   START OF THE ITERATIVE
  //   ALGORITHM
  ////////////////////////////////////
  ////////////////////////////////////

  vector <unsigned int> n_Models_visited;
  vector < vector <unsigned int> > List_Models;
  vector < vector <unsigned int> > Unique_List_Models;
  
  n_Models_visited.assign(n_sweeps+1,0);
    
  for(sweep=1;sweep<n_sweeps+1;sweep++){
    clock_t Time_start,Time_end;
    Time_start=clock();
 
    n_Models_visited[sweep]=n_Models_visited[sweep-1];
    if(Log_Flag){
      cout << "***************************************" << endl
	   << "             SWEEP #" << sweep << "/" << n_sweeps << endl
	   << "***************************************" << endl << endl; 
    }
    ////////////////////////////////////////////////////
    /////////////////// Local Moves  ///////////////////
    ////////////////////////////////////////////////////
    for(unsigned int chain=0;chain<nb_chains;chain++){
      mat_log_marg.matrix[chain][sweep]=mat_log_marg.matrix[chain][sweep-1];
      mat_log_cond_post.matrix[chain][sweep]=mat_log_cond_post.matrix[chain][sweep-1];
    }

    //Gibbs Moves
    if(sweep%Gibbs_n_batch==0){
      if(Log_Flag){
	cout << "Gibbs" << endl;
      }
      Gibbs_move(mat_log_marg,
		 mat_log_cond_post,
		 mat_X_work2,
		 mat_Y_work2,
		 t_tun,
		 MyPerm,
		 vect_gam,
		 sweep,
		 g_opt,
		 g_sample,
		 g,
		 PR,
		 My_Move_monitor,
		 chain_idx,
		 n_Models_visited,
		 cudaFlag);


    }


    

    ////Fast Scan Metropolis Hastings (FSMH)
    
    double local_move_rand=myrand();
    if(nb_chains==1){
      local_move_rand=0.0;
    }
    if(local_move_rand<PR.Prob_mut){
      if(Log_Flag){
	cout << "FSMH" << endl;
      }
      FSMH_move(mat_log_marg,
		mat_log_cond_post,
		mat_X_work2,
		mat_Y_work2,
		t_tun,
		MyPerm,
		vect_gam,
		sweep,
		g_opt,
		g_sample,
		g,
		PR,
		My_Move_monitor,
		chain_idx,
		n_Models_visited,
		cudaFlag);
 
    }
    else{
      ///////////////////////////////////////////////////////
      /////////////////////// CM Moves  /////////////////////
      ///////////////////////////////////////////////////////

      //Defining the number of CO move to simulate
      if(Log_Flag){
	cout << "CM" << endl;
      }
      Crossover_move(mat_log_marg,
		     mat_log_cond_post,
		     mat_X_work2,
		     mat_Y_work2,
		     t_tun,
		     vect_gam,
		     g,
		     g_opt,
		     g_sample,
		     PR,
		     My_Move_monitor,
		     chain_idx,
		     sweep,
		     My_CM,
		     n_Models_visited,
		     cudaFlag);
    }

    ///////////////////////////////////////////////////////
    ///////////////////// Global Moves  ///////////////////
    ///////////////////////////////////////////////////////
    if(nb_chains>1){
      double global_move_rand=myrand();
      if(sweep<=burn_in){
	global_move_rand=0.0;//during burn-in, DR is the only global move
      }
      
      if(global_move_rand<PR.Prob_DR){
	
	if(Log_Flag){
	  cout << "DR" << endl;
	}
	
	DR_move(My_DR,
		chain_idx,
		mat_log_cond_post,
		t_tun,
		sweep,
		My_Move_monitor);
	
      }
      else{ 
	if(Log_Flag){
	  cout << "AE" << endl;
	}
	All_exchange_move(description_exchange_move,
			  chain_idx,
			  mat_log_cond_post,
			  t_tun,
			  sweep,
			  My_Move_monitor);
      }

    }
    
    /////////////////////////////////////////////////////
    ///////////////////// Sampling g  ///////////////////
    /////////////////////////////////////////////////////

    if(g_sample==1){
  
      sample_g(mat_log_marg,
	       mat_log_cond_post,
	       My_g_AdMH,
	       g,
	       t_tun,
	       vect_gam,
	       mat_X_work2,
	       mat_Y_work2,
	       g_opt,
	       g_sample,
	       PR,
	       sweep,
	       My_Move_monitor,
	       chain_idx,
	       n_Models_visited,
	       cudaFlag);
      if(Log_Flag){
	cout << "g_sample= " << g << endl;
      }

    }
    
    
    ///////////////////////////////////////////////////
    ///////////////// Temp placement //////////////////
    ///////////////////////////////////////////////////
    if(nb_chains>1){
      if(sweep<=burn_in){
	if(((*My_DR).nb_calls_adj==(*t_tun).nbatch) || ((*My_DR).nb_calls==5*(*t_tun).nbatch)){
	  unsigned int n_vars_in_last_chain=sum_line_std_mat(vect_gam,
							     chain_idx[nb_chains-1]);
	  

          //cout << "t_placement" << endl;
	  
	  temp_placement(t_tun,
			 My_DR,
			 My_Move_monitor,
			 sweep,
			 n_vars_in_last_chain,
			 nX,
			 iso_T_Flag);
	  //cout << "end -- t_placement" << endl;
	  
	}
      }
    }
 
    if(Log_Flag){
      
      display_summary_result_per_sweep(vect_gam,
				      chain_idx,
				      mat_log_marg,
				      mat_log_cond_post,
				      sweep,
				      t_tun);  
    }
    
    
    
    print_and_save_main_results_per_sweep(f_out,
					 f_out_n_vars_in,
					 f_out_log_cond_post_per_chain,
					 vect_gam,
					 List_Models,
					 chain_idx,
					 mat_log_marg,
					 mat_log_cond_post,
					 sweep,
					 burn_in,
					 HistoryFlag);
    
    //     cout << "Sweep\tTime\t\t#models\t#cum models\tTime/model" << endl;
    
    Time_end=clock();
    double time_taken=((double)Time_end-Time_start)/CLOCKS_PER_SEC;
    if(Log_Flag){
      cout << "Sweep\tTime\t#models\tTime/model" << endl;
      cout << sweep << "\t"
          << time_taken << "\t"
          << n_Models_visited[sweep]-n_Models_visited[sweep-1] << "\t"
          << time_taken/(double)(n_Models_visited[sweep]-n_Models_visited[sweep-1])
          <<  endl;

      cout << "***************************************" << endl
          << "           END OF SWEEP #" << sweep+1 << "/" << n_sweeps << endl
          <<  "***************************************" << endl << endl;
    }
    
    if(Time_monitorFlag){
      f_out_time << sweep << "\t" << time_taken << "\t" << time_taken/(double)(n_Models_visited[sweep]-n_Models_visited[sweep-1]) << endl;
    }
    
    if(sweep>burn_in){
      cum_g+=g;
      count_g++;
    }
 


  }//end of for sweep


  if(HistoryFlag){
    f_out.close();
    f_out_n_vars_in.close();
    f_out_log_cond_post_per_chain.close();
  }
  if(Time_monitorFlag){
    f_out_time.close();
  }
 
    cout << "***************************************" << endl
	 << "***************************************" << endl
	 << "           END OF SWEEPS" << endl
	 << "       # models evaluated: " << n_Models_visited[n_sweeps-1] << endl
	 << "***************************************" << endl
	 <<  "***************************************" << endl << endl;

  tmpTime = endTime;
  endTime = std::clock();
  mainLoopTime = (endTime-tmpTime)/(double)(CLOCKS_PER_SEC);


  //////////////////
  //Post-Processing
  //////////////////

  //Step1: Calculating E(g);
  double mean_g=cum_g/(double)(count_g);
  //Step2: Get the Unique list of models and integrate log posterior
  
  
  vector< vector<double> >vect_marg_log_post;
  int pos_null_model;
  pos_null_model=getUniqueList(Unique_List_Models,
                                           List_Models,
                                           n_Models_visited,
                                           burn_in,
                                           pX);

  getLogPost(vect_marg_log_post,
            Unique_List_Models,
            mat_X_work2,
            mat_Y_work2,
            PR,
            g_opt,
            g_sample,
            mean_g,
            cudaFlag);


  //Step4: Get the posterior of the model
  unsigned int n_unique = Unique_List_Models.size();
  gsl_permutation *idx_post_gam_sort=gsl_permutation_calloc(vect_marg_log_post.size());
  gsl_vector *vect_post_gam = gsl_vector_calloc(vect_marg_log_post.size());

  unsigned int n_retained_models=min(n_unique,n_top_models);


  getAndSortPostGam(vect_post_gam,
                   idx_post_gam_sort,
                   vect_marg_log_post);

  combineAndPrintBestModel(f_out_best_models,
				   idx_post_gam_sort,
				   vect_post_gam,
				   vect_marg_log_post,
				   Unique_List_Models,
				   n_retained_models,
				   pos_null_model,
				   Out_full_Flag);
  

  getAndPrintMargGam(f_out_marg_gam,
			 Unique_List_Models,
			 vect_post_gam,
			 pX);
  
  tmpTime = endTime;
  endTime = std::clock();
  postProcessTime = (endTime-tmpTime)/(double)(CLOCKS_PER_SEC);
  cout << "Setup Time: " << setupTime  << endl;
  cout << "Main Loop Time: " << mainLoopTime  << endl;
  cout << "Post Processing Time: " << postProcessTime  << endl;

  f_out_marg_gam.close();
  
  
  f_out_best_models.close();
  
  
  gsl_permutation_free(idx_post_gam_sort);
  gsl_vector_free(vect_post_gam);
    
  for(unsigned int line=0;line<List_Models.size();line++){
    List_Models[line].clear();
  }
  List_Models.clear();
  
  for(unsigned int line=0;line<Unique_List_Models.size();line++){
    Unique_List_Models[line].clear();
    vect_marg_log_post[line].clear();
  }
  
  Unique_List_Models.clear();
    
  if(HistoryFlag){  
    
    string OutputName_FSMH=Get_stddzed_name(path_name_out,
					    n_sweeps,
					    Name_number,
					    "output_fast_scan_history",
					    Extension_out);
    
    ofstream f_out_FSMH;
    f_out_FSMH.open(OutputName_FSMH.c_str(),ios::out);
 

    string OutputName_CM=Get_stddzed_name(path_name_out,
					    n_sweeps,
					    Name_number,
					    "output_cross_over_history",
					    Extension_out);
    
    ofstream f_out_CM;
    f_out_CM.open(OutputName_CM.c_str(),ios::out);
    
    

    string OutputName_AE=Get_stddzed_name(path_name_out,
					  n_sweeps,
					  Name_number,
					  "output_all_exchange_history",
					  Extension_out);
    
    ofstream f_out_AE;
    f_out_AE.open(OutputName_AE.c_str(),ios::out);
    
    
    string OutputName_DR=Get_stddzed_name(path_name_out,
					  n_sweeps,
					  Name_number,
					  "output_delayed_rejection_history",
					  Extension_out);
    ofstream f_out_DR;
    f_out_DR.open(OutputName_DR.c_str(),ios::out);
    string OutputName_g=Get_stddzed_name(path_name_out,
					 n_sweeps,
					 Name_number,
					 "output_g_history",
					 Extension_out);
    ofstream f_out_g;
    if(g_sample!=0){
      f_out_g.open(OutputName_g.c_str(),ios::out);
    }
    string OutputName_g_adapt=Get_stddzed_name(path_name_out,
					       n_sweeps,
					       Name_number,
					       "output_g_adaptation_history",
					       Extension_out);
    ofstream f_out_g_adapt;
    if(g_sample!=0){
      f_out_g_adapt.open(OutputName_g_adapt.c_str(),ios::out);
    }
    
    
    string OutputName_Gibbs=Get_stddzed_name(path_name_out,
					     n_sweeps,
					     Name_number,
					     "output_gibbs_history",
					     Extension_out);
    ofstream f_out_Gibbs;
    f_out_Gibbs.open(OutputName_Gibbs.c_str(),ios::out);
    
    string OutputName_t_tun=Get_stddzed_name(path_name_out,
					     n_sweeps,
					     Name_number,
					     "output_temperature_history",
					     Extension_out);  
    ofstream f_out_t_tun;
    if(iso_T_Flag==false){   
      f_out_t_tun.open(OutputName_t_tun.c_str(),ios::out);
    }
    
    (*My_Move_monitor).print_move_monitor_full(f_out_FSMH,
					       f_out_CM,
					       f_out_AE,
					       f_out_DR,
					       f_out_g,
					       f_out_g_adapt,
					       f_out_Gibbs,
					       f_out_t_tun,
					       g_sample,
					       iso_T_Flag);
    
    
    f_out_FSMH.close();
    f_out_CM.close();
    f_out_AE.close();
    f_out_DR.close();
    if(g_sample!=0){
      f_out_g.close();
      f_out_g_adapt.close();
    }
    f_out_Gibbs.close();
    if(iso_T_Flag==false){   
      f_out_t_tun.close();
    }
  }
  
  n_Models_visited.clear();

  
  for(unsigned int chain=0;chain<nb_chains;chain++){  
    vect_gam[chain].clear();  
  }
  vect_gam.clear();  
  mat_Y_copy.Free_double_matrix_cont();
  gsl_matrix_free(mat_X_work2);
  gsl_matrix_free(mat_Y_work2);
  mat_X_copy2.Free_double_matrix_cont();
  
  gsl_vector_free(vect_RMSE);
  mat_log_marg.Free_double_matrix();
  mat_log_cond_post.Free_double_matrix();
  gsl_permutation_free(MyPerm);
  chain_idx.clear();
  gsl_matrix_free(description_exchange_move);
  delete (My_g_AdMH);
  delete (My_Move_monitor);
  delete (My_CM);
  delete (My_DR);
  delete (t_tun);

  return 0; 
}

