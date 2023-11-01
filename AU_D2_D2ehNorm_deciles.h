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

//un-normalized Eisenstein-Hu power spectrum; Ldec=0-9 for nu, 10 for m
int D2ehUz0(const double *p, int Ldec, double z, double *D2u, double nrm){

  //only need to call with param vector once
  static int init = 0;
  static double om, ob, on, h0, ns;
  if(!init){ om=p[0]; ob=p[1]; on=p[2]; h0=p[3]; ns=p[4]; init=1; }
  
  const double Nn=3.044, TcmbK=2.726, T27=TcmbK/2.7;
  const double tau_Ldec_eV[] = {0.0001501274, 0.00024775, 0.0003170484,
				0.0003806534, 0.0004446064, 0.0005131268,
				0.0005912368, 0.0006876718, 0.0008242774,
				0.0011337452};
  const double fb = ob / om;
  const double fn = on / om;
  const double fcb = 1.0 - fn, fc = fcb - fb;

  const double a1 = pow(46.9*om,0.670)
    * (1.0 + pow(32.1*om,-0.532));
  const double a2 = pow(12.0*om,0.424)
    * (1.0 + pow(45.0*om,-0.582));
  const double alc = pow(a1,-fb/fcb) * pow(a2,-pow(fb/fcb,3));
  const double b1 = 0.944 / (1.0 + pow(458.0*om,-0.708));
  const double b2 = pow(0.395*om,-0.0266);
  const double bec = 1.0 / (1.0 + b1 * ( pow(fc/fcb,b2) - 1.0 ) );
  const double c1 = 0.313 * pow(om,-0.419)
    * (1.0 + 0.607*pow(om,0.674));
  const double c2 = 0.238 * pow(om,0.223);
  
  const double zd = 1291.0 * pow(om,0.251)
    * (1.0 + c1*pow(ob,c2)) / (1.0 + 0.659*pow(om,0.828));
  const double zeq = 2.50e4 * om * pow(T27,-4);
  const double keq = 0.0746 * om / h0 / pow(T27,2);
  const double Rd = 31.5e3 * ob * pow(T27,-4) / zd;
  const double Req = 31.5e3 * ob * pow(T27,-4) / zeq;
  const double s = (2.0 / (3.0*keq)) * sqrt(6.0/Req)
    * log( (sqrt(1.0+Rd) + sqrt(Rd + Req)) / (1.0 + sqrt(Req)) );
  
  const double pc = 0.25*(5.0-sqrt(1.0+24.0*fc));
  const double pcb = 0.25*(5.0-sqrt(1.0+24.0*fcb));

  const double m_nu_eV = 93.25 * on / 3.0;
  const double T_nu_eV = 0.000168;
  const double zeta_3 = 1.20205690315959;
  const double H0h = 100.0 / 299792.458;
  const double k_fs_2_0 = 1.5 * om/(h0*h0)
    * H0h*H0h * pow(m_nu_eV/tau_Ldec_eV[Ldec],2);

  const double y = (1.0+zeq)/(1.0+zd);
  const double G = y * (-6.0*sqrt(1.0+y)
			+ (2.0+3.0*y)*log((sqrt(1.0+y)+1.0)
					  /(sqrt(1.0+y)-1.0)) );
  const double alb = 2.07 * keq * s * pow(1.0+Rd,-0.75) * G;
  const double beb = 0.5 + fb
    + (3.0 - 2.0*fb)*sqrt( pow(17.2*om,2) + 1.0 );
  const double benode = 8.41 * pow(om,0.435);
  const double kSilk = 1.6 / h0 * pow(ob,0.52)
    * pow(om,0.73) * (1.0 + pow(10.4*om,-0.95));
  
  double g2h2	= om*pow(1.0+z,3)
    + (h0*h0-om);
  double Om_z	= om*pow(1.0+z,3) / g2h2;
  double OL_z	= (h0*h0-om) / g2h2;
  double D1 =	(1.0+zeq)/(1.0+z) * 2.5*Om_z
    / (pow(Om_z,(4.0/7.0)) - OL_z + (1.0+0.5*Om_z)/(1.0+OL_z/70.0));
    
  for(int i=0; i<NK; i++){
    double lnk = LNKMIN + DLNK*i, k = exp(lnk);
    double q = k*h0 * T27*T27 / om;
    double C1 = 14.0/1.0 + 386.0/(1 + 69.9*pow(q,1.08));
    double Ca = 14.0/alc + 386.0/(1 + 69.9*pow(q,1.08));
    double T0t_11 = log(exp(1.0) + 1.8*1.0*q)
      / ( log(exp(1.0) + 1.8*1.0*q) + C1*q*q );
    double T0t_1b = log(exp(1.0) + 1.8*bec*q)
      / ( log(exp(1.0) + 1.8*bec*q) + C1*q*q );
    double T0t_ab = log(exp(1.0) + 1.8*bec*q)
	/ ( log(exp(1.0) + 1.8*bec*q) + Ca*q*q );
    double fTc = 1.0 / (1.0 + pow(k*s/5.4,4));
    double Tc = fTc*T0t_1b + (1.0-fTc)*T0t_ab;
    double stil = s / pow(1.0 + pow(benode/k/s,3),1.0/3.0);
    double j0 = sin(k*stil) / (k*stil);
    double Tb = (T0t_11 / (1.0 + pow(k*s/5.2,2))
		 + alb * exp(-pow(k/kSilk,1.4)) / (1.0 + pow(beb/k/s,3))) * j0;
    double qn = 3.92 * q * sqrt(Nn/fn);
    double Bnu = 1.0 + 1.24*pow(fn,0.64)*pow(Nn,0.3+0.6*fn)
      / (pow(qn,-1.6) + pow(qn,0.8));
    double yfs_n = 17.2*fn * (1.0+0.488*pow(fn,(-7.0/6.0))) * pow(Nn*q/fn,2);
    double Dcb_D1 = pow(1.0 + pow(D1/(1.0+yfs_n),0.7),(pcb/0.7)) / pow(D1,pcb);
    double Dcbn_D1 = pow( pow(fcb,(0.7/pcb))
			  + pow(D1/(1.0+yfs_n),0.7), (pcb/0.7) ) / pow(D1,pcb);
    double Tcb = Bnu*Dcb_D1 * ( (fc/fcb) * Tc + (fb/fcb) * Tb );
    double Tm = Bnu*Dcbn_D1 * ( (fc/fcb) * Tc + (fb/fcb) * Tb );

    double Pdec_Pm = 1.0; //super-easy fit
    if(Ldec<10) Pdec_Pm = pow(1.0 + k*sqrt((1.0+z)/k_fs_2_0), -4);
    
    D2u[i] = nrm * Pdec_Pm * D1*D1 * Tm*Tm * pow(k,ns+3.0) / (2.0*M_PI*M_PI);
  }

  return 0;
}

double norm(const double *p){

  //only need to call with param vector once
  static double nrm=-1;
  if(nrm<0){
    double integrand[NK];
    double om=p[0], ob=p[1], on=p[2], h0=p[3], ns=p[4], s8=p[5];
    D2ehUz0(p,10,0,integrand,1);
    
    for(int i=0; i<NK; i++){
      double lnk = LNKMIN + DLNK*i, k = exp(lnk), kR=8.0*k, kR2=kR*kR;
      double W = ( kR<1e-2 ? 1.0 - 0.1*kR2
		   : 3.0 * ( sin(kR)/(kR2*kR) - cos(kR)/kR2 ) );
      integrand[i] *= W*W;
    }
    nrm = s8*s8 / ncint_cf(NK,DLNK,integrand); 
  }

  return nrm;
}

//initialize D2eh norm
int init_D2eh(double om,double ob,double on,double h0,double ns,double s8){
  double p[6] = {om, ob, on, h0, ns, s8};
  return norm(p)<0;
}
		    
//normalized Eisenstein-Hu power spectrum for nu or total matter
//NOTE: Make sure to initialize norm first!
int D2eh(int Ldec, double z, double *D2){
  double dum, nrm=norm(&dum);
  D2ehUz0(&dum,Ldec,z,D2,nrm);
  return 0;
}
  
