////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//    Copyright 2023 Amol Upadhye                                             //
//                                                                            //
//    This file is part of Cosmic-Enu.                                        //
//                                                                            //
//    Cosmic-Enu is free software: you can redistribute it and/or modify      //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation, either version 3 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    Cosmic-Enu is distributed in the hope that it will be useful,           //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with Cosmic-Enu.  If not, see <http://www.gnu.org/licenses/>.     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "AU_ncint.h"
#include "AU_fftgrid.h"
#include "AU_D2_D2ehNorm_deciles.h"
#include "AU_fftm_emu_data.h"

int main(int argn, char *args[]){

  //emulator inputs: assume everything after first 9 is a redshift
  //model_string, om, ob, s8, h0, ns, w0, wa, on, z0, z1, ...
  char *model = args[1];
  double om = atof(args[2]);
  double ob = atof(args[3]);
  double s8 = atof(args[4]);
  double h0 = atof(args[5]);
  double ns = atof(args[6]);
  double w0 = atof(args[7]);
  double wa = atof(args[8]);
  double on = atof(args[9]);

  //initialize Eisenstein-Hu power spectrum
  init_D2eh(om,ob,on,h0,ns,s8);

  int Nzi = argn-10;
  double *zi = malloc(Nzi*sizeof(double));
  for(int i=0; i<Nzi; i++){
    zi[i] = atof(args[i+10]);
    if(zi[i]<ENU_z_i[26] || zi[i]>ENU_z_i[0]){
      printf("ERROR: redshift %g out of bounds.\n", zi[i]);
      fflush(stdout);
      abort();
    }
  } 

  //remap input cosmic parameters to theta
  const double xmin[] = {0.12,  0.0215, 0.7, 0.55, 0.85, -1.3, 0.3,  0.0};
  const double xmax[] = {0.155, 0.0235, 0.9, 0.85, 1.05, -0.7, 1.29, 0.01};
  double x[] = {om, ob, s8, h0, ns, w0, pow(-w0-wa,0.25), on};
  double *t = malloc(ENU_Npar*sizeof(double));
  for(int ell=0; ell<ENU_Npar; ell++)
    t[ell] = (x[ell]-xmin[ell]) / (xmax[ell]-xmin[ell]);

  //get weights
  double *w = malloc(10*ENU_Npc*sizeof(double));
#pragma omp parallel for schedule(dynamic) collapse(2)
  for(int L=0; L<10; L++){
    for(int j=0; j<ENU_Npc; j++){
      w[L*ENU_Npc+j] = 0;
      for(int m=0; m<ENU_Nsim; m++){
        double lnRm = 0;
        for(int ell=0; ell<ENU_Npar; ell++)
	  lnRm += -ENU_beta_jl[L][j*ENU_Npar+ell]
	    * sq(t[ell] - ENU_thetaStar_ml[m*ENU_Npar+ell]);
        w[L*ENU_Npc+j] += exp(lnRm)*ENU_Rw_jm[L][j*ENU_Nsim+m] / ENU_lam[L][j];
      }
    }
  }

  //emulate D^2 / D^2_{EH,L}, but reverse order of k and z for interpolation
  double *ykz = malloc(10*ENU_Nkz*sizeof(double));
#pragma omp parallel for schedule(dynamic) collapse(2)
  for(int L=0; L<10; L++){
    for(int i=0; i<ENU_Nkz; i++){
      int iz=i/78, ik=i%78;
      double lny = ENU_mean[L][i];
      for(int j=0; j<ENU_Npc; j++) 
        lny += ENU_stddev[L] * ENU_PhiP_ji[L][j*ENU_Nkz+i] * w[L*ENU_Npc+j];
      ykz[L*ENU_Nkz + ik*27 + iz] = exp(lny);
    }
  }

  //interpolate to desired redshift
  double *yout = malloc(10*78*Nzi*sizeof(double));
  for(int iy=0; iy<10*78*Nzi; iy++) yout[iy] = 0;
  
  double a_list[27];
  for(int iz=0; iz<27; iz++) a_list[iz] = 1.0 / (1.0 + ENU_z_i[iz]);

  //precompute Eisenstein-Hu power spectrum
  //Note: This computes all 128 k values, so offset by 25.
  double *D2 = malloc(10*Nzi*NK*sizeof(double));
#pragma omp parallel for schedule(dynamic) collapse(2)
  for(int L=0; L<10; L++){
    for(int iz=0; iz<Nzi; iz++){
      D2eh(L, zi[iz], D2 + L*Nzi*NK + iz*NK);
    }
  }

  //emulate power in deciles
#pragma omp parallel for schedule(dynamic) collapse(2)
  for(int L=0; L<10; L++){
    for(int ik=0; ik<78; ik++){
      gsl_interp_accel *ac = gsl_interp_accel_alloc();
      gsl_spline *sp = gsl_spline_alloc(gsl_interp_cspline, 27);
      gsl_spline_init(sp, a_list, ykz + L*ENU_Nkz + ik*27, 27);
      for(int iz=0; iz<Nzi; iz++)
        yout[L*78*Nzi + iz*78 + ik] = 
          D2[L*Nzi*NK + iz*NK + ik+25] *gsl_spline_eval(sp,1.0/(1.0+zi[iz]),ac);
      gsl_spline_free(sp);
      gsl_interp_accel_free(ac);
    }
  }

  //output y (deciles and total) vs. k, one file per z
  for(int iz=0; iz<Nzi; iz++){
    char filename[200];
    sprintf(filename, "enu_%s_z%g.dat", model, zi[iz]);
    FILE *fp = fopen(filename,"w");

    for(int ik=0; ik<78; ik++){
      double dnutot=0; //total nu density perturbation
      fprintf(fp, "%g", ENU_k_i[ik]);
      for(int L=0; L<10; L++){
        fprintf(fp, " %g", yout[L*78*Nzi + iz*78 + ik]);
        dnutot += 0.1*sqrt(yout[L*78*Nzi + iz*78 + ik]);
      }
      fprintf(fp, " %g\n", dnutot*dnutot);
    }
    fclose(fp);
  }

  free(zi);
  free(t);
  free(w);
  free(ykz);
  free(yout);
  free(D2);
  return 0;
}
