function G = BP_KERNELS_XYZ(x,y,z,d,a,mu,kap)
%BP_KERNELS Calculates the kernels for a circular PRISMATIC GLIDE loop.
%   G = BP_KERNELS(rho,zeta,delta,a,mu,kap) returns the kernels 
%   G for a dislocation burgers vector Bz as a function of the
%   normalized variables zeta and delta. The dislocation has a radius 'a',
%   is put at (x,y,z) and evaluated at a depth 'd'.
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

% Angle
t=atan(y./x);
c=cos(t);
c2=cos(t).^2;
s=sin(t);
s2=sin(t).^2;

% Sign ( + for z-d>0, - for z-d<0);
signa=sign(z-d);

r=sqrt(x.^2+y.^2);

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

Gxx=cst.*c.*(-signa.*4.*J111+4.*I111+...
    2.*(z-d).*c2.*J112-2.*((3-kap).*d+(z+kap.*d).*c2).*I112+...
    4.*z.*d.*c2.*I113+...
    2.*(z-d).*(3-4.*c2).*J121byr-2.*(z+kap.*d).*(3-4.*c2).*I121byr+...
    (3-4.*c2).*4.*z.*d.*I122byr);

Gyy=cst.*c.*(-signa.*(3-kap).*J111+(3-kap).*I111+...
    2.*(z-d).*s2.*J112-2.*((3-kap).*d+(z+kap.*d).*s2).*I112+...
    4.*z.*d.*s2.*I113-...
    2.*(z-d).*(3-4.*c2).*J121byr+2.*(z+kap.*d).*(3-4.*c2).*I121byr-...
    (3-4.*c2).*4.*z.*d.*I122byr);

Gzz=cst.*c.*(-2.*(z-d).*J112+2.*(z-d).*I112-...
    4.*z.*d.*I113);

% Shear components

Gxy=cst.*s.*(-signa.*1/2.*(1+kap).*J111+1/2.*(1+kap).*I111+...
    2.*(z-d).*c2.*J112-2.*(z+kap.*d).*c2.*I112+...
    4.*z.*d.*c2.*I113-...
    2.*(z-d).*(3-4.*s2).*J121byr+2.*(z+kap.*d).*(3-4.*s2).*I121byr-...
    (3-4.*s2).*4.*z.*d.*I122byr);


Gxz=cst.*(-1/2.*(kap+1+(3-kap).*c2).*J101+1/2.*(kap+1+(3-kap).*c2).*I101+...
    2.*c2.*abs(z-d).*J102-2.*c2.*(z+d).*I102+...
    4.*c2.*z.*d.*I103+...
    (3-kap)/2.*cos(2.*t).*J110byr-(3-kap)/2.*cos(2.*t).*I110byr-...
    2.*cos(2.*t).*abs(z-d).*J111byr+2.*cos(2.*t).*(z+d).*I111byr-...
    4.*z.*d.*cos(2.*t).*1./r.*I112);

Gyz=cst/2.*sin(2.*t).*((3-kap)./2.*J121-(3-kap)./2.*I121-...
    2.*abs(z-d).*J122+2.*(z+d).*I122-...
    4.*z.*d.*I123);

% Kernel strucutre

G.xx=Gxx;
G.yy=Gyy;
G.zz=Gzz;
G.xy=Gxy;
G.xz=Gxz;
G.yz=Gyz;



end