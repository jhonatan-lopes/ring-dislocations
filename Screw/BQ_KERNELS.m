function G = BQ_KERNELS(rho,zeta,delta,a,mu,kap)
%BQ_KERNELS Calculates the kernels for a screw ring dislocation bq.
%   G = BQ_KERNELS(rho,zeta,delta,a,mu,kap) returns the kernels 
%   Grq,Gzq for a dislocation burgers vector Bq as a function of the
%   normalized variables zeta and delta. The dislocation has a radius 'a',
%   is put at (rho,zeta) and evaluated at a depth 'delta'.
%
%   The modulus of rigidity is 'mu' and 'kappa' the Kolosov's constant.
%
%   University of Oxford 
%   Department of Engineering Science
%   Jhonatan Da Ponte Lopes, PhD 
%   May, 2019; Last revision: 2019-05-13


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

J211=LH_INTEGRALS(a,r,z1,0,Km,Em,'211');
J221=LH_INTEGRALS(a,r,z1,0,Km,Em,'221');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I TYPE  -(ZETA + DELTA)

I211=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'211');
I221=LH_INTEGRALS(a,r,z2,0,Kp,Ep,'221');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernels

cst=mu.*a./(2);

G.rq=cst*(J221-I221);

G.zq=cst*(J211-I211);

end