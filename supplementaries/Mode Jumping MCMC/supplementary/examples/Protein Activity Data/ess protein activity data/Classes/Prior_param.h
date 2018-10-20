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

#ifndef PRIOR_PARAM_H
#define PRIOR_PARAM_H 1

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>
#include <math.h>
#include "../Routines/matrix_handling.h"
#include <gsl/gsl_statistics_double.h>

using namespace std;

class Prior_param
{
public:
  Prior_param();
  ~Prior_param(){};
  double E_p_gam;
  double Sd_p_gam;
  double a;
  double b;
  double a_pi;
  double b_pi;
  vector<double> w;
  double delta;
  double k;
  double alpha;
  double beta;
  double mu;
  double Prob_mut;
  double Prob_sel;
  double Prob_crsv_r;
  double Prob_DR;

 

  void set_PR_param(int g_opt,
		    double E_p_gam_from_read,
		    double Sd_p_gam_from_read,
		    unsigned int pX,
		    unsigned int pY,
		    unsigned int nX,
		    gsl_vector *vect_RMSE,
		    double P_mutation_from_read,
		    double Prob_sel_from_read,
		    double P_crsv_r_from_read,
		    double P_DR_from_read);

  void display_prior_param();

 private:
  
  void get_prior_parameters_g_prior(double E_p_gam_from_read,
				    double Sd_p_gam_from_read,
				    unsigned int pX,
				    unsigned int pY,
				    unsigned int nX,
				    gsl_vector *vect_RMSE,
				    double P_mutation_from_read,
				    double Prob_sel_from_read,
				    double P_crsv_r_from_read,
				    double P_DR_from_read);
  
  void get_prior_parameters_ind_prior(double E_p_gam_from_read,
				      double Sd_p_gam_from_read,
				      unsigned int pX,
				      unsigned int pY,
				      gsl_vector *vect_RMSE,
				      double P_mutation_from_read,
				      double Prob_sel_from_read,
				      double P_crsv_r_from_read,
				      double P_DR_from_read);
};

#endif /* !defined PRIOR_PARAM_H */
