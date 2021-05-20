function G = BZ_KERNELS(rho,zeta,delta,a,mu,kap)
%BZ_KERNELS Calculates the kernels for a ring dislocation bz.
%   G = BZ_KERNELS(rho,zeta,delta,a,mu,kap) returns the kernel structure 
%   G for a dislocation burgers vector Bz as a function of the
%   normalized variables zeta and delta. The dislocation has a radius 'a',
%   is put at (rho,zeta) and evaluated at a depth 'delta'.
%
%   The modulus of rigidity is 'mu' and 'kappa' the Kolosov's constant.
%
%   University of Oxford 
%   Department of Engineering Science
%   Jhonatan Da Ponte Lopes, PhD
%   May, 2017; Last revision: 2017-05-09


%-------------------------------------------------------------------
%                         KERNELS
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

J101=LH_INTEGRALS(a,r,z1,0,Km,Em,'101');
J102=LH_INTEGRALS(a,r,z1,0,Km,Em,'102');
J110byr=LH_INTEGRALS(a,r,z1,0,Km,Em,'110byr');
J111byr=LH_INTEGRALS(a,r,z1,0,Km,Em,'111byr');
J112=LH_INTEGRALS(a,r,z1,0,Km,Em,'112');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I TYPE  -(ZETA + DELTA)

I101=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'101');
I102=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'102');
I103=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'103');
I110byr=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'110byr');
I111byr=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'111byr');
I112=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'112');
I113=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'113');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernels

cst=2.*mu.*a./(kap+1);

G.rr=cst.*(I101-J101+(zeta-3.*delta).*I102+z1.*J102-2.*zeta.*delta.*I103-...
    (kap-1)./(2).*I110byr+(kap-1)./(2).*J110byr-...
    (zeta-kap.*delta).*I111byr-z1.*J111byr+2.*zeta.*delta./rho.*I112); % OK!

G.rz=cst.*(-z1.*J112+z1.*I112-2.*zeta.*delta.*I113); % OK!

G.zz=cst.*(I101-J101+(-delta-zeta).*I102-z1.*J102+2.*delta.*zeta.*I103);


end