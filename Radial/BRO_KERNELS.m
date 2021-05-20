function G = BRO_KERNELS(rho,zeta,delta,alpha,a,mu,kap)
%BRO_KERNELS Calculates the kernels for a ring dislocation br.
%   G = BRO_KERNELS(rho,zeta,delta,alpha,a,mu,kap) returns 
%   the kernels G for a dislocation burgers vector iBr as a 
%   function of the normalized variables zeta and delta. The dislocation
%   has a radius 'a', is put at (rho,zeta) and evaluated at a depth
%   'delta'.
%
%   The modulus of rigidity is 'mu' and 'kappa' the Kolosov's constant.
%
%   The dislocation path cut is OUTSIDE the dislocation ring.
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
[Km,Em] = ellipke(k2m);

k2p=(4.*par.*r)./((par+r).^2+z2.^2);
[Kp,Ep] = ellipke(k2p);

%-------------------------------------------------------------------
%                         LH INTEGRALS
%-------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J TYPE (ZETA - DELTA)

J001=LH_INTEGRALS(a,r,z1,0,Km,Em,'001');
J002=LH_INTEGRALS(a,r,z1,0,Km,Em,'002');
J010byr=LH_INTEGRALS(a,r,z1,alpha,Km,Em,'010byr');
J011=LH_INTEGRALS(a,r,z1,0,Km,Em,'011');
J011byr=LH_INTEGRALS(a,r,z1,0,Km,Em,'011byr');
J012=LH_INTEGRALS(a,r,z1,0,Km,Em,'012');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I TYPE  -(ZETA + DELTA)

I001=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'001');
I002=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'002');
I003=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'003');
I010byr=LH_INTEGRALS(a,r,z2,abs(alpha),Kp,Ep,'010byr');
I011=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'011');
I011byr=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'011byr');
I012=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'012');
I012byr=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'012byr');
I013=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'013');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernels

cst=2.*mu.*a./(kap+1);

Grr=cst.*(-2.*J001-2.*I001+z1.*J002-(zeta+3.*delta).*I002-...
    2.*zeta.*delta.*I003+(kap+1)./2.*J010byr+(kap+1)./2.*I010byr-...
    z1.*J011byr+(zeta+kap.*delta).*I011byr+2.*zeta.*delta.*I012byr);

Grz=cst.*(J011-I011-z1.*J012+z2.*I012-2.*zeta.*delta.*I013);

Gzz=cst.*(-z1.*J002+z1.*I002+2.*zeta.*delta.*I003);

G.rr = Grr;
G.rz = Grz;
G.zz = Gzz;

end