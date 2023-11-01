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

/*******************************************************************************
composite Newton-Cotes integration
Use closed 4th-degree method in steps until fewer than 4 points remain,
then switch to a lower-degree method.
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double ncint_boole_coeffs[5] =
  {0.311111111111111, 1.42222222222222, 0.533333333333333, 1.42222222222222,
   0.311111111111111};

double ncint_boole_step(double dx, const double *f){
  double res = 0;
  for(int i=0; i<5; i++) res += ncint_boole_coeffs[i] * f[i];
  return res * dx;
}

const double ncint_simpson38_coeffs[4] = {0.375, 1.125, 1.125, 0.375};

double ncint_simpson38_step(double dx, const double *f){
  double res = 0;
  for(int i=0; i<4; i++) res += ncint_simpson38_coeffs[i] * f[i];
  return res * dx;
}

const double ncint_simpson_coeffs[3] =
  {0.333333333333333, 1.33333333333333, 0.333333333333333};

double ncint_simpson_step(double dx, const double *f){
  double res = 0;
  for(int i=0; i<3; i++) res += ncint_simpson_coeffs[i] * f[i];
  return res * dx;
}

const double ncint_trapezoid_coeffs[2] = {0.5, 0.5};

double ncint_trapezoid_step(double dx, const double *f){
  double res = 0;
  for(int i=0; i<2; i++) res += ncint_trapezoid_coeffs[i] * f[i];
  return res * dx;
}

//composite fixed-step-size Newton-Cotes integrator
double ncint_cf(int N, double dx, const double *f){
  if(N<2) return 0;
  int Nint = N-1; //number of intervals
  double result = 0;
  
  while(Nint >= 6){
    Nint -= 4;
    result += ncint_boole_step(dx, &f[Nint]);
  }

  switch(Nint){
  case 5:
    result += ncint_simpson38_step(dx,&f[2]) + ncint_simpson_step(dx,f);
    break;
  case 4:
    result += ncint_boole_step(dx,f);
    break;
  case 3:
    result += ncint_simpson38_step(dx,f);
    break;
  case 2:
    result += ncint_simpson_step(dx,f);
    break;
  case 1: //should only get here if N=2, otherwise don't use trapezoid rule
    result += ncint_trapezoid_step(dx,f);
    break;
    //case 0: //shouldn't get here
    //break;
  default:
    printf("ERROR in ncint: Nint=%i invalid. Quitting.\n",Nint);
    fflush(stdout);
    abort();
    break;
  }

  return result;   
}

