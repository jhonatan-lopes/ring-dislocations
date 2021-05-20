function G = BP_KERNELS_RTZ(r,t,z,d,a,mu,kap)
%BP_KERNELS_RTZ Calculates the kernels for a circular PRISMATIC GLIDE loop.
%   G = BP_KERNELS_RTZ(rho,zeta,delta,a,mu,kap) returns the kernels 
%   G for a dislocation burgers vector Bz as a function of the
%   normalized variables zeta and delta. The dislocation has a radius 'a',
%   is put at (r,t,z) and evaluated at a depth 'd'.
%
%   CYLINDRICAL COORDINATES R, THETA, Z
%
%   The modulus of rigidity is 'mu' and 'kappa' the Kolosov's constant.
%
%   University of Oxford 
%   Department of Engineering Science
%   Jhonatan Da Ponte Lopes, PhD 
%   November, 2018; Last revision: 2018-11-28


%-------------------------------------------------------------------
%                         KERNELS
%-------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arguments and initial variables

% Angle cosines and sines
c=cos(t);
c2=cos(t).^2;
s=sin(t);
s2=sin(t).^2;

% Sign ( + for z-d>0, - for z-d<0);
signa=sign(z-d);

z1=abs(z-d);
z2=-abs(z+d);

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
% J TYPE |ZETA - DELTA|

J101=LH_INTEGRALS(a,r,z1,0,Km,Em,'101');
J102=LH_INTEGRALS(a,r,z1,0,Km,Em,'102');

J111=LH_INTEGRALS(a,r,z1,0,Km,Em,'111');
J112=LH_INTEGRALS(a,r,z1,0,Km,Em,'112');

J110byr=LH_INTEGRALS(a,r,z1,0,Km,Em,'110byr');
J111byr=LH_INTEGRALS(a,r,z1,0,Km,Em,'111byr');

J121=LH_INTEGRALS(a,r,z1,0,Km,Em,'121');
J122=LH_INTEGRALS(a,r,z1,0,Km,Em,'122');

J121byr=LH_INTEGRALS(a,r,z1,0,Km,Em,'121byr');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I TYPE  -|ZETA + DELTA|

I101=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'101');
I102=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'102');
I103=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'103');

I111=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'111');
I112=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'112');
I113=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'113');

I110byr=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'110byr');
I111byr=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'111byr');

I121=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'121');
I122=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'122');
I123=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'123');

I121byr=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'121byr');
I122byr=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'122byr');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernels

cst=mu.*a./(kap+1);

% Normal components

Grr=cst.*c.*(4.*I111-4.*J111.*signa+...
    2.*I112.*(-3.*d-z)+2.*J121byr.*(d-z)+...
    4.*d.*I113.*z-4.*d.*I122byr.*z+2.*J112.*(-d+z)+...
    2.*I121byr.*(z+d.*kap));

Gtt=cst.*c.*(4.*d.*I122byr.*z+J121byr.*(-2.*d+2.*z)+...
    I111.*(3-kap)+2.*d.*I112.*(-3+kap)+...
    J111.*signa.*(-3+kap)+I121byr.*(-2.*z-2.*d.*kap));

Gzz=cst.*c.*(-2.*(z-d).*J112+2.*(z-d).*I112-...
    4.*z.*d.*I113);

% Shear components

Grt=cst.*s./2.*(J121byr.*(4.*d-4.*z)-8.*d.*I122byr.*z+...
    I111.*(-1-kap)+J111.*signa.*(1+kap)+I121byr.*(4.*z+4.*d.*kap));

Grz=cst.*c./2.*(8.*d.*I103.*z.*c2+4.*J102.*abs(d-z).*c2+...
    I102.*(-4.*d.*c2-4.*z.*c2)+...
    I101.*(1+kap+(3-kap).*c2)+...
    J101.*(-1-kap+(-3+kap).*c2)-...
    8.*d.*I112./r.*z.*cos(2.*t)+...
    I110byr.*(-3+kap).*cos(2.*t)-...
    4.*J111byr.*abs(d-z).*cos(2.*t)+...
    I111byr.*(4.*d.*cos(2.*t)+4.*z.*cos(2.*t))+...
    J110byr.*(3.*cos(2.*t)-kap.*cos(2.*t))-...
    8.*d.*I123.*z.*s2+I121.*(-3+kap).*s2-...
    4.*J122.*abs(d-z).*s2+...
    I122.*(4.*d.*s2+4.*z.*s2)+...
    J121.*(3.*s2-kap.*s2));

Gzt=cst.*1/2.*s.*(-(I101-J101).*(1+kap)+(4.*d.*(I102+I122)+...
    4.*(I102+I122-2.*d.*(I103+I123)).*z+...
    I101.*(-3+kap)+(I121-J101-J121).*(-3+kap)-...
    4.*(J102+J122).*abs(d-z)).*c2+(-4.*d.*I111byr-...
    3.*J110byr-4.*I111byr.*z+8.*d.*1./r.*I112.*z-...
    I110byr.*(-3+kap)+J110byr.*kap+...
    4.*J111byr.*abs(d-z)).*cos(2.*t));

% Kernel strucutre

G.rr=Grr;
G.tt=Gtt;
G.zz=Gzz;
G.rz=Grz;
G.rt=Grt;
G.zt=Gzt;

end