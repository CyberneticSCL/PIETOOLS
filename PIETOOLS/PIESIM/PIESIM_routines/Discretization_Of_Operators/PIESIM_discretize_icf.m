%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_discretize_icf.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform discretization of initial conditions and forcing functions 
%
% Inputs:
% 1) uinput - user-defined boundary inputs, forcing and initial conditions
% 2) psize - size of the problem
% 3) gridall - cell array containing grids of different degrees of differentiability
%
% Outputs:
% 1) coeff - Chebyshev coefficients for initial conditions and forcing functions
% 2) B1_nonpol - contribution to PIE B1 operator arising from non-polynomial forcing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  11_05_2021
% YP 6_16_2022 Renamed uinput.B21_nonpol to uinput.Bpw_nonpol, updated description of
% outputs
% YP - added functionality to support infinite-dimensional disturbances
% through parser - 6_1_2025
% YP 1/7/2026 - changed the notaiton of ic.PDE to ic.PIE for PIE simulation


function [coeff,B1_nonpol]=PIESIM_discretize_icf(uinput,psize,gridall);

% Define local variables

no=psize.no;
nw=psize.nw;
nu=psize.nu;
ns=sum(psize.n);
a=uinput.a;
b=uinput.b;
N=psize.N;

% Define degree of smoothness p

p = repelem(0:length(psize.n)-1, psize.n);

% Define initial conditions on states on the physical gridall for states.
% var_f denotes the value of the solution variable of the fundamental
% states.
% Define Chebyshev coefficients of initial conditions in the same loop
% acheb_f denotes the Chebyshev coefficients of the fundamental states
% Each state vector array coefficients are arranged into a global column
% vector

ic=uinput.ic.PIE;

for i=1:ns
     acheb=fcht(double(subs(ic(i),gridall{p(i)+1})));
     acheb_glob{i}=reshape(acheb, [], 1);
     clear('acheb');
end


% Concatenate coefficients of all states into a single column vector

 acheb_f0=cat(1, acheb_glob{:});

% Add initial conditions on ODE states to the front of initial conditions
if (no>0)
acheb_f0=cat(1,double(uinput.ic.ODE)',acheb_f0);
end

coeff.acheb_f0=acheb_f0;


% Discretize matrix operator for non-polynomial in space forcing

if isfield(uinput,'Bpw_nonpol') 
B1_nonpol = PIESIM_NonPoly2Mat_cheb(N, uinput.Bpw_nonpol, p, gridall);
else
B1_nonpol=[];   
end

coeff.w=1;
coeff.u=1;

% Discretize spatial contribution of u and w disturbances entered in parser
% format as spatially-varying disturbances

if isfield(uinput,'wspace')
     Nforce=psize.nw0+(N+1)*psize.nwx;
     coeff.w=zeros(Nforce,nw);
         k=0;
         index=1;
         for kk=1:psize.nw0
         k=k+1;
             coeff.w(index,k)=1;
             index=index+1;
         end
         for kk=1:psize.nwx
             k=k+1;
             coeff.w(index:index+N,k)=PIESIM_NonPoly2Mat_cheb(N, uinput.wspace(k), 0, gridall);
             index=index+N+1;
         end 
 end

 if isfield(uinput,'uspace')
     Nforce=psize.nu0+(N+1)*psize.nux;
     coeff.u=zeros(Nforce,nu);
         k=0;
         index=1;
         for kk=1:psize.nu0
         k=k+1;
             coeff.u(index,k)=1;
             index=index+1;
         end
         for kk=1:psize.nux
             k=k+1;
             coeff.u(index:index+N,k)=PIESIM_NonPoly2Mat_cheb(N, uinput.uspace(k), 0, gridall);
             index=index+N+1;
         end 
 end












  



