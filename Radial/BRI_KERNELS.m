function G = BRI_KERNELS(rho,zeta,delta,alpha,a,mu,kap)
%BRI_KERNELS Calculates the kernels for a ring dislocation br.
%   G = BRI_KERNELS(rho,zeta,delta,alpha,a,mu,kap) returns 
%   the kernels G for a dislocation burgers vector iBr as a 
%   function of the normalized variables zeta and delta. The dislocation
%   has a radius 'a', is put at (rho,zeta) and evaluated at a depth
%   'delta'.
%
%   The modulus of rigidity is 'mu' and 'kappa' the Kolosov's constant.
%
%   The dislocation path cut is INSIDE the dislocation ring.
%
%   University of Oxford 
%   Department of Engineering Science
%   Jhonatan Da Ponte Lopes, PhD 
%   May, 2017; Last revision: 2017-05-09


%-------------------------------------------------------------------
%                         INTIALIZATION
%-------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arguments and initial variables

r=rho;

z1=zeta-delta;
z2=-(zeta+delta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elliptic integrals

par=a;
k2m=(4.*par.*r)./((par+r).^2+z1.^2);
[Km,Em]=ellipke(k2m);

k2p=(4.*par.*r)./((par+r).^2+z2.^2);
[Kp,Ep]=ellipke(k2p);

%-------------------------------------------------------------------
%                         LH INTEGRALS
%-------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J TYPE (ZETA - DELTA)

J201=LH_INTEGRALS(a,r,z1,alpha,Km,Em,'201');
J202=LH_INTEGRALS(a,r,z1,0,Km,Em,'202');
J210byr=LH_INTEGRALS(a,r,z1,alpha,Km,Em,'210byr');
J211=LH_INTEGRALS(a,r,z1,0,Km,Em,'211');
J211byr=LH_INTEGRALS(a,r,z1,0,Km,Em,'211byr');
J212=LH_INTEGRALS(a,r,z1,0,Km,Em,'212');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I TYPE  -(ZETA + DELTA)

I201=LH_INTEGRALS(a,r,z2,abs(alpha),Kp,Ep,'201');
I202=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'202');
I203=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'203');
I210byr=LH_INTEGRALS(a,r,z2,abs(alpha),Kp,Ep,'210byr');
I211byr=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'211byr');
I211=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'211');
I212=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'212');
I212byr=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'212byr');
I213=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'213');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernels

cst=2.*mu.*a./(kap+1);

Grr=cst.*(2.*J201+2.*I201-z1.*J202+(zeta+3.*delta).*I202+...
    2.*zeta.*delta.*I203-(kap+1)./2.*J210byr-(kap+1)./2.*I210byr+...
    z1.*J211byr-(zeta+kap.*delta).*I211byr-2.*zeta.*delta.*I212byr);

Grz=cst.*(-J211+I211+z1.*J212+(zeta+delta).*I212+2.*zeta.*delta.*I213); %OK

Gzz=cst.*(z1.*J202-z1.*I202-2.*zeta.*delta.*I203);

G.rr = Grr;
G.rz = Grz;
G.zz = Gzz;

end