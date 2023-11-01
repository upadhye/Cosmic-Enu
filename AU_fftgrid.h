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

////////////////////////////////////////////////////////////////////////////////
//grid parameters

#define NK (128)
#define KMIN (1e-4)
#define KMAX (1e1)

////////////////////////////////////////////////////////////////////////////////
//k grid and windowing: NK is number of points in "real" k grid,                 
//I extend this by a factor of four for fast-pt.  The extended grid: 
//   -3NK/2  <= i < -17NK/16 : P=0
//  -17NK/16 <= i < -NK      : P tapered smoothly from 0 to P[0]
//     -NK   <= i < 0        : Power spectrum extrapolated to left
//       0   <= i < NK       : Power spectrum region of interest
//      NK   <= i < 2NK      : Power spectrum extrapolated to right
//     2NK   <= i < 33NK/16  : P tapered smoothly from P[NK-1] to 0
//  33NK/16  <= i < 5NK/2    : P=0

#define NKP (4*NK)

#define NSHIFT ((NKP-NK)/2)
#define DLNK (log(KMAX/KMIN)/(NK-1.0+1e-300))
#define LNKMIN (log(KMIN))

//split up interval between zero-padding and tapering:
//these are measured in units of NK/16, and must add to (NKP/NK-1)*16
#define S_PADL (7)
#define S_TAPL (1)
#define S_EXTL (8*(NKP/NK-2))
#define S_EXTR (S_EXTL)
#define S_TAPR (1)
#define S_PADR (7)

#define LNK_PAD_MIN   (LNKMIN        - DLNK*NSHIFT)
#define LNK_PAD_WINLO (LNK_PAD_MIN   + DLNK*NK*S_PADL/16)
#define LNK_PAD_WINLI (LNK_PAD_WINLO + DLNK*NK*S_TAPL/16)
#define LNK_PAD_WINRI (LNK_PAD_WINLI + DLNK*(NK*(16+S_EXTL+S_EXTR)/16 -1))
#define LNK_PAD_WINRO (LNK_PAD_WINRI + DLNK*NK*S_TAPR/16)
#define LNK_PAD_MAX   (LNK_PAD_WINRO + DLNK*NK*S_PADR/16)

////////////////////////////////////////////////////////////////////////////////
//FAST-PT constants and functions

//useful functions: square, cube, quadrature add
double sq(double x){ return x*x; }
double cu(double x){ return x*x*x; }
inline double qu(double x){ return sq(sq(x)); }
inline double qadd(double x, double y){ return sqrt(sq(x)+sq(y)); }

//tilt applied to power spectrum before fft
#define FASTPT_TILT_NU (-2)
const int nu_int = FASTPT_TILT_NU;
const double nu = FASTPT_TILT_NU;

//power spectrum and Fourier coefficient window functions
inline double W_edge(double x){
  if(x<0) return 0;
  if(x>1) return 1;
  return x - sin(2.0*M_PI*x)/(2.0*M_PI);
}

double WP(double lnk){
  if(lnk <= LNK_PAD_WINLO) return 0;
  else if(lnk < LNK_PAD_WINLI) 
    return W_edge((lnk-LNK_PAD_WINLO)/(LNK_PAD_WINLI-LNK_PAD_WINLO));
  else if(lnk < LNK_PAD_WINRI) return 1;
  else if(lnk < LNK_PAD_WINRO) 
    return W_edge((LNK_PAD_WINRO-lnk)/(LNK_PAD_WINRO-LNK_PAD_WINRI));
  return 0;
}

double WC(int n){ //coeffs in GSL halfcomplex notation

  const int nl = NKP/4, nc = NKP/2, nr = NKP-nl, Dn = NKP/4;

  if(n<=nl || n>=nr) return 1;
  else if(n<nc) return W_edge((double)(Dn+nl-n)/Dn);
  else if(n<nr) return W_edge((double)(n-nr+Dn)/Dn);
  return 1;
}

