function G = BRE_KERNELS(rho,zeta,delta,a,mu,kap)
%BRE_KERNELS Calculates the kernels for a ring dislocation br.
%   G = BRE_KERNELS(rho,zeta,delta,a,mu,kap) returns the 
%   kernels G for a dislocation burgers vector eBr as a 
%   function of the normalized variables zeta and delta. The dislocation
%   has a radius 'a', is put at (rho,zeta) and evaluated at a depth
%   'delta'.
%
%   The modulus of rigidity is 'mu' and 'kappa' the Kolosov's constant.
%
%   The dislocation path cut is cylindrical, to the EXTERIOR of the
%   material. ALPHA is -pi/2.
%
%   University of Oxford 
%   Department of Engineering Science
%   Jhonatan Da Ponte Lopes, PhD
%   June, 2017; Last revision: 2017-06-13


%-------------------------------------------------------------------
%                         INTIALIZATION
%-------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arguments and initial variables

alpha=-pi/2;

%-------------------------------------------------------------------
%                   PREVIOUS SOLUTIONS KERNELS
%-------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inside cut Kernel

[iGrz,iGrr,iGzz]=BRI_KERNELS(rho,zeta,delta,alpha,a,mu,kap);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outside cut Kernel

[oGrz,oGrr,oGzz]=BRO_KERNELS(rho,zeta,delta,alpha,a,mu,kap);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kelvin solution Kernel

[kGrz,kGrr,kGzz]=KEL_KERNELS(rho,zeta,delta,alpha,a,kap);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boussinesq correction Kernel

[bGrz,bGrr,bGzz]=BOU_KERNELS(rho,zeta,a,kap);


%-------------------------------------------------------------------
%                         KERNELS
%-------------------------------------------------------------------

Grr=(kap-1)./(kap+1).*iGrr+2/(kap+1).*oGrr-...
    2.*mu.*(3-kap)./(a.*(kap+1)).*(kGrr-bGrr);

Grz=(kap-1)./(kap+1).*iGrz+2./(kap+1).*oGrz-...
    2.*mu.*(3-kap)./(a.*(kap+1)).*(kGrz-bGrz);

Gzz=(kap-1)./(kap+1).*iGzz+2./(kap+1).*oGzz-...
    2.*mu.*(3-kap)./(a.*(kap+1)).*(kGzz-bGzz);

G.rr = Grr;
G.rz = Grz;
G.zz = Gzz;

end

