%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_discretize_all.m     PIETOOLS 2021d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform discretization of everything  
%
% Inputs:
% PIE
% uinput
% psize
%
% Outputs:
% Dop
% coeff
% ngrid
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 9_18_2021
function [Dop, coeff, grid]=PIESIM_discretize_all(PIE, uinput, psize);

% Define local variables

nx=psize.nx;
nf=psize.nf;
nw=psize.nw;
nu=psize.nu;
n0=0;
n1=psize.n(1);
n2=0;
if(size(psize.n,2)>1)
n1=psize.n(2);
end
if(size(psize.n,2)>2)
n2=psize.n(3);
end
a=uinput.a;
b=uinput.b;
N=psize.N;


% Define "Degree of smoothness for each PDE state". Currently we keep track of it in a discretization of PIE
% for a consistent representation of the original states. This is needed to properly
% account for the number of Chebyshev coefficients that are being solved
% for for each state: N+1 coefficients for n0 states (polynomilas of degree N), N coefficients for n1
% states (polynomials of degree N-1), N-1 coefficients for n2 states
% (polynomials of degree N-2)

p(1:n0)=0;
p(n0+1:n0+n1)=1;
p(n0+n1+1:n0+n1+n2)=2;

% Define grids in the physical space as vector arrays. Spatial grid points
% are colocated with the Chebyshev nodes.

% Define primary state Chebyshev grid in computational domain
% The same grid is used for n0 fundamental states
n0grid_comp = cos(pi*(0:N)/N)';

% Define Chebyshev grid in computational domain for n1 fundamental states
n1grid_comp = cos(pi*(0:N-1)/(N-1))';

% Define Chebyshev grid in computational domain for n2 fundamental states
n2grid_comp = cos(pi*(0:N-2)/(N-2))';

% Define Chebyshev grid in physical domain
n0grid=0.5*(b-a)*n0grid_comp+0.5*(b+a);
n1grid=0.5*(b-a)*n1grid_comp+0.5*(b+a);
n2grid=0.5*(b-a)*n2grid_comp+0.5*(b+a);

% Define initial conditions on states on the physcal grid for states.
% var_f denotes the value of the solution variable of the fundamental
% states.
% Define Chebyshev coefficients of initial conditions in the same loop
% acheb_f denotes the Chebyshev coefficients of the fundamental states
% Each state vector array coefficients are arranged into a global column
% vector


% For i=1:nf Evaluate Chebyshev coefficients of the sptaital components of the forcing function
ic=uinput.ic.PDE;
if (nf>0)
force=uinput.force;
end

acheb_fn0=double.empty(N+1,0);
acheb_forcen0=double.empty(N+1,0);
for n=1:n0;
    var_fn0(:,n)=double(subs(ic(n),n0grid));
    acheb_fn0(:,n)=fcht(var_fn0(:,n));
    for i=1:nf 
    var_forcen0(:,n)=double(subs(force(n,i),n0grid));
    acheb_forcen0(:,n,i)=fcht(var_forcen0(:,n));
    end
end
acheb_fn0_glob=reshape(acheb_fn0, [], 1);
acheb_forcen0_glob=reshape(acheb_forcen0, [], nf);


acheb_fn1=double.empty(N,0);
acheb_forcen1=double.empty(N,0);
for n=1:n1;
    var_fn1(:,n)=double(subs(ic(n0+n)*0.5*(b-a),n1grid));
    acheb_fn1(:,n)=fcht(var_fn1(:,n));
    for i=1:nf 
    var_forcen1(:,n)=double(subs(force(n0+n,i),n1grid));
    acheb_forcen1(:,n,i)=fcht(var_forcen1(:,n));
    end
end

acheb_fn1_glob= reshape(acheb_fn1, [], 1);
acheb_forcen1_glob= reshape(acheb_forcen1, [], nf);

acheb_fn2=double.empty(N-1,0);
acheb_forcen2=double.empty(N-1,0);
for n=1:n2;
    var_fn2(:,n)=double(subs(ic(n0+n1+n)*(0.5*(b-a))^2,n2grid));
    acheb_fn2(:,n)=fcht(var_fn2(:,n));
    for i=1:nf 
    var_forcen2(:,n)=double(subs(force(n0+n1+n,i),n2grid));
    acheb_forcen2(:,n,i)=fcht(var_forcen2(:,n));
    end
end

acheb_fn2_glob= reshape(acheb_fn2, [], 1);
acheb_forcen2_glob= reshape(acheb_forcen2, [], nf);


% Concatenate coefficients of all states into a single column vector

acheb_f0=cat(1, acheb_fn0_glob, acheb_fn1_glob, acheb_fn2_glob);

% Add initial conditions on ODE states to the front of initial conditions
if (nx>0)
acheb_f0=cat(1,uinput.ic.ODE',acheb_f0);
end

for i=1:nf
acheb_force_PDE(:,i)=cat(1, acheb_forcen0_glob(:,i), acheb_forcen1_glob(:,i), acheb_forcen2_glob(:,i));
end
% Currently assume that ODE states do no contribute into forcing function
% defined that way
ODE.forcing(1:nx,1:nf)=0;
for i=1:nf
acheb_force(:,i)=cat(1, ODE.forcing(:,i), acheb_force_PDE(:,i));
end

coeff.acheb_f0=acheb_f0;
if (nf>0)
    coeff.acheb_force=acheb_force;
end


disp('Setting up Chebyshev matrices for the PIE system');

% Discretize 4PI operators

% Set the last entry to PIESIM_4PI2Mat_cheb to 1 if a structure is a full 4PI operator (for A and T)

% Set the last entry to PIESIM_4PI2Mat_cheb to 0 if a structure has an empty right side (for Tw, Tu, B1 and B2
% operators)

 Dop.Twcheb=PIESIM_4PI2Mat_cheb(N,nw,PIE.Tw,p,0);
 Dop.Tucheb=PIESIM_4PI2Mat_cheb(N,nu,PIE.Tu,p,0);
 Dop.B1cheb=PIESIM_4PI2Mat_cheb(N,nw,PIE.B1,p,0);
 Dop.B2cheb=PIESIM_4PI2Mat_cheb(N,nu,PIE.B2,p,0);
 Dop.Acheb=PIESIM_4PI2Mat_cheb(N,nx,PIE.A,p,1);
[Mcheb, Dop.Mcheb_nonsquare]=PIESIM_4PI2Mat_cheb(N,nx,PIE.T,p,1);


 
%  
  Dop.Mcheb_inv=inv(Mcheb);
  Dop.Atotal=Dop.Mcheb_inv*Dop.Acheb;
  
  grid.phys=n0grid;
  grid.comp=n0grid_comp;



