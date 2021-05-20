function G = KEL_KERNELS(rho,zeta,delta,alpha,a,kap)
%KEL_KERNELS Calculates the kernels for a ring dislocation br.
%   G = KEL_KERNELS(rho,zeta,delta,alpha,a,kap) returns 
%   the kernels G for the KELVIN correction term
%   necessary for cylindrical cut paths.  The dislocation has a radius 'a', 
%   is put at (rho,zeta) and evaluated at a depth 'delta'.
%
%   The modulus of rigidity is 'mu' and 'kappa' the Kolosov's constant.
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

J100=LH_INTEGRALS(a,r,z1,alpha,Km,Em,'100');
J101=LH_INTEGRALS(a,r,z1,0,Km,Em,'101');
J110=LH_INTEGRALS(a,r,z1,0,Km,Em,'110');
J110byr=LH_INTEGRALS(a,r,z1,0,Km,Em,'110byr');
J111=LH_INTEGRALS(a,r,z1,0,Km,Em,'111');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I TYPE  -(ZETA + DELTA)

I100=LH_INTEGRALS(a,r,z2,abs(alpha),Kp,Ep,'100');
I101=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'101');
I102=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'102');
I11m1byr=LH_INTEGRALS(a,r,z2,abs(alpha),Kp,Ep,'11m1byr');
I110=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'110');
I110byr=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'110byr');
I111=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'111');
I111byr=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'111byr');
I112=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'112');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernels

cst=a./(kap+1);

Grr=cst.*(-(3-kap)./2.*J100+(3.*kap-1)./2.*I100+...
    z1.*J101+(kap.*zeta-3.*delta).*I101-2.*zeta.*delta.*I102-...
    z1.*J110byr-z1.*kap.*I110byr+2.*zeta.*delta.*I111byr+...
    (1-kap^2)./2.*I11m1byr);

Grz=cst.*(-(kap-1)./2.*J110+(kap-1)./2.*I110-...
    z1.*J111+(kap.*zeta-delta).*I111-2.*zeta.*delta.*I112);

Gzz=cst.*(-(kap+1)./2.*J100+(kap+1)./2.*I100-...
    z1.*J101-(delta+zeta.*kap).*I101+2.*zeta.*delta.*I102);

G.rr = Grr;
G.rz = Grz;
G.zz = Gzz;

end