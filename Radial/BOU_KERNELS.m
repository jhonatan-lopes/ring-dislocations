function G = BOU_KERNELS(rho,zeta,a,kap)
%BOU_KERNELS Calculates the kernels for a ring dislocation br.
%   G = BOU_KERNELS(rho,zeta,a,kap) returns the kernel structure
%   G for the BOUSSINESQ correction term necessary for 
%   cylindrical cut paths.  The dislocation has a radius 'a', 
%   is put at (rho,zeta) and evaluated at a depth 'delta'=0.
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
[rl,rc]=size(r);

alpha=-pi/2;

delta=0;
z1=(zeta-delta).*ones(rl,rc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elliptic integrals

par=a;
k2m=(4.*par.*r)./((par+r).^2+z1.^2);
[Km,Em] = ellipke(k2m);


%-------------------------------------------------------------------
%                         LH INTEGRALS
%-------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J TYPE (ZETA - DELTA)

J100=LH_INTEGRALS(a,r,z1,alpha,Km,Em,'100');
J101=LH_INTEGRALS(a,r,z1,0,Km,Em,'101');
J11m1byr=LH_INTEGRALS(a,r,z1,alpha,Km,Em,'11m1byr');
J110byr=LH_INTEGRALS(a,r,z1,0,Km,Em,'110byr');
J111=LH_INTEGRALS(a,r,z1,0,Km,Em,'111');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernels

cst=a;

Grr=cst.*(-J100+zeta.*J101-zeta.*J110byr+(kap-1)./2.*J11m1byr);

Grz=cst.*(-zeta.*J111);

Gzz=cst.*(-J100-z1.*J101);

G.rr = Grr;
G.rz = Grz;
G.zz = Gzz;


end