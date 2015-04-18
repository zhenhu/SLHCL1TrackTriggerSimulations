/* --- An Analytical Track Fitter for rapid least-chi2 fitting ---

Algorithm described in CMS Detector Note DN-2015/004
cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=DN2015_004_v2.pdf

Souvik Das, University of Florida
15 October 2014
--- ---  --- --- --- --- --- --- --- --- --- --- --- --- --- --- */

#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TrackFitterAlgoATF.h"
#include <iostream>

int stubUncertaintyStripWidth(double r, double &d_phi, double &d_z)
{
  
  double d_w_PS=0.01/sqrt(12);
  double d_w_SS=0.009/sqrt(12);
  double d_z_PS=0.15/sqrt(12);
  double d_z_SS=5/sqrt(12);
  
  if (fabs(r-22)<5)
  {
    d_phi=d_w_PS/r;
    d_z=d_z_PS;
  }
  else if (fabs(r-34)<5)
  {
    d_phi=d_w_PS/r;
    d_z=d_z_PS;
  }
  else if (fabs(r-50)<5)
  {
    d_phi=d_w_PS/r;
    d_z=d_z_PS;
  }
  else if (fabs(r-67)<5)
  {
    d_phi=d_w_SS/r;
    d_z=d_z_SS;
  }
  else if (fabs(r-88)<5)
  {
    d_phi=d_w_SS/r;
    d_z=d_z_SS;
  }
  else if (fabs(r-107)<5)
  {
    d_phi=d_w_SS/r;
    d_z=d_z_SS;
  }
  else
  {
    d_phi=-1;
    d_z=-1;
    return 0;
  }
  return 1;
}

int stubUncertaintyClusterWidth(double r, double &d_phi, double &d_z)
{
  
  double d_z_PS=0.15/sqrt(12);
  double d_z_SS=5/sqrt(12);
  
  if (fabs(r-22)<5)
  {
    d_phi=0.0086/r;
    d_z=d_z_PS;
  }
  else if (fabs(r-34)<5)
  {
    d_phi=0.0082/r;
    d_z=d_z_PS;
  }
  else if (fabs(r-50)<5)
  {
    d_phi=0.0115/r;
    d_z=d_z_PS;
  }
  else if (fabs(r-67)<5)
  {
    d_phi=0.0176/r;
    d_z=d_z_SS;
  }
  else if (fabs(r-88)<5)
  {
    d_phi=0.0228/r;
    d_z=d_z_SS;
  }
  else if (fabs(r-107)<5)
  {
    d_phi=0.0280/r;
    d_z=d_z_SS;
  }
  else
  {
    d_phi=-1;
    d_z=-1;
    return 0;
  }
  return 1;
}

int stubUncertaintyFountain(double r, double &d_phi, double &d_z)
{
  
  double d_z_PS=0.15/sqrt(12);
  double d_z_SS=5/sqrt(12);
  
  if (fabs(r-22)<5)
  {
    d_phi=0.00381;
    d_z=d_z_PS;
  }
  else if (fabs(r-34)<5)
  {
    d_phi=0.00439;
    d_z=d_z_PS;
  }
  else if (fabs(r-50)<5)
  {
    d_phi=0.00459;
    d_z=d_z_PS;
  }
  else if (fabs(r-67)<5)
  {
    d_phi=0.00485;
    d_z=d_z_SS;
  }
  else if (fabs(r-88)<5)
  {
    d_phi=0.00523;
    d_z=d_z_SS;
  }
  else if (fabs(r-107)<5)
  {
    d_phi=0.00575;
    d_z=d_z_SS;
  }
  else
  {
    d_phi=-1;
    d_z=-1;
    return 0;
  }
  
  d_phi/=2.0;
  
  return 1;
}

int stubUncertaintyResidual(double r, double &d_phi, double &d_z)
{
  
  double d_z_PS=0.15/sqrt(12);
  double d_z_SS=5/sqrt(12);
  
  if (fabs(r-22)<5)
  {
    d_phi=0.0002;
    d_z=d_z_PS;
  }
  else if (fabs(r-34)<5)
  {
    d_phi=0.00009;
    d_z=d_z_PS;
  }
  else if (fabs(r-50)<5)
  {
    d_phi=0.0001;
    d_z=d_z_PS;
  }
  else if (fabs(r-67)<5)
  {
    d_phi=0.0003;
    d_z=d_z_SS;
  }
  else if (fabs(r-88)<5)
  {
    d_phi=0.0005;
    d_z=d_z_SS;
  }
  else if (fabs(r-107)<5)
  {
    d_phi=0.0007;
    d_z=d_z_SS;
  }
  else
  {
    d_phi=-1;
    d_z=-1;
    return 0;
  }
  
  return 1;
}

int TrackFitterAlgoATF::fit(const TTRoadComb& acomb, TTTrack2& track)
{
  double A=0, B=0, C=0, D=0, E=0, F=0, G=0, H=0, I=0, J=0, K=0, L=0, M=0, P=0, Q=0;
  for (unsigned int i=0; i<acomb.stubs_phi.size(); ++i)
  {
    if (!acomb.stubs_bool.at(i))
        continue;

    double r=acomb.stubs_r.at(i);
    double phi=acomb.stubs_phi.at(i);
    double z=acomb.stubs_z.at(i);
    
    double d_phi=-1, d_z=-1;
    if (!stubUncertaintyClusterWidth(r, d_phi, d_z)) {std::cout<<"ERROR: Unrecognized layer at r = "<<r<<std::endl; continue;}

    double d_phi_sq=d_phi*d_phi;
    double d_z_sq=d_z*d_z;

    A+=phi/d_phi_sq;
    B+=1./d_phi_sq;
    C+=r/d_phi_sq;
    D+=1./(r*d_phi_sq);
    E+=r*phi/d_phi_sq;
    F+=r*r/d_phi_sq;
    G+=z*r/d_z_sq;
    H+=r/d_z_sq;
    I+=r*r/d_z_sq;
    J+=phi/(r*d_phi_sq);
    K+=1./(r*r*d_phi_sq);
    L+=z/d_z_sq;
    M+=1./d_z_sq;
    
    P+= phi*phi/d_phi_sq;
    Q+= z*z/d_z_sq;
  }
  
  if (fiveParameters_)
  {  
    double denom1 = B*B*B - 2*B*C*D + D*D*F + C*C*K - B*F*K;
    double denom2 = I*M - H*H;
  
    double phi0      = (A*B*B - B*D*E - B*C*J + D*F*J + C*E*K - A*F*K)/denom1;
    double rinv      = 2.*(A*B*D - D*D*E - B*B*J + C*D*J - A*C*K + B*E*K)/denom1;
    double d0        = (A*B*C - B*B*E + C*D*E - A*D*F - C*C*J + B*F*J)/denom1;
    double z0        = (I*L - G*H)/denom2;
    double cottheta0 = ((G*M - H*L)/denom2)/(1 - d0*rinv/2.);
  
    double chi2 = P + pow(phi0, 2)*B + F*pow(rinv/2., 2) + pow(d0, 2)*K
                - 2.*A*phi0 + E*rinv + 2.*d0*J
                - C*phi0*rinv - 2.*phi0*d0*D
                + d0*rinv*B
                + Q + M*pow(z0, 2) + I*pow(cottheta0*(1.-(d0*rinv/2.)), 2) - 2.*L*z0 - 2.*G*cottheta0*(1.-(d0*rinv/2.))
                + 2.*H*z0*cottheta0*(1.-(d0*rinv/2.));
                 
    int    ndof = (acomb.nstubs-5);

    double chi2_phi = P + pow(phi0, 2)*B + F*pow(rinv/2., 2) + pow(d0, 2)*K
                    - 2.*A*phi0 + E*rinv + 2.*d0*J
                    - C*phi0*rinv - 2.*phi0*d0*D
                    + d0*rinv*B;

    double chi2_z   = Q + M*pow(z0, 2) + I*pow(cottheta0*(1.-(d0*rinv/2.)), 2) - 2.*L*z0 - 2.*G*cottheta0*(1.-(d0*rinv/2.))
                    + 2.*H*z0*cottheta0*(1.-(d0*rinv/2.));

    track.setTrackParams(rinv, phi0, cottheta0, z0, d0, chi2, ndof, chi2_phi, chi2_z);
  }
  else
  {
    double denom1 = B*F - C*C;
    double denom2 = I*M - H*H;
    
    double phi0      = (A*F - C*E)/denom1;    
    double rinv      = 2.*(A*C - B*E)/denom1; 
    double z0        = (I*L - G*H)/denom2;    
    double cottheta0 = (G*M-H*L)/denom2;
    
    double chi2 = P + pow(phi0, 2)*B + F*pow(rinv/2., 2)
                - 2.*A*phi0 + E*rinv 
                - C*phi0*rinv
                + Q + M*pow(z0, 2) + I*pow(cottheta0, 2) - 2.*L*z0 - 2.*G*cottheta0 
                + 2.*H*z0*cottheta0;
                 
    int    ndof = (acomb.nstubs-4);

    double chi2_phi = P + E*rinv/2. - phi0*A;

    double chi2_z = Q - L*z0 - G*cottheta0;
    
    track.setTrackParams(rinv, phi0, cottheta0, z0, 0, chi2, ndof, chi2_phi, chi2_z);
  }

  return 0;
}
