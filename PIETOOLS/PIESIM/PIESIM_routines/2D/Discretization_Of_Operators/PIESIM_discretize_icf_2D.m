%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_discretize_icf_2D.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform discretization of initial conditions and forcing functions in 2D
%
% Inputs:
% 1) uinput - user-defined boundary inputs, forcing and initial conditions
% 2) psize - size of the PIE problem: all variables defining the size of the PIE problem
% 3) gridall - cell array of size dmax containing physical grid for all states
% depending on their degree of differentiability; dmax corresponds to the
% maximum degree of differentiability among all the states
%  gridall.x - grids in x direction
%  gridall.y - grids in y direction
%
% Output:
% coeff - Chebyshev coefficients for initial conditions and forcing functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  04_16_2024

function coeff=PIESIM_discretize_icf_2D(uinput,psize,gridall)

syms sx sy st;

% Define local variables

no=psize.no;
nw=psize.nw;
a=uinput.dom(1,1);
b=uinput.dom(1,2);
c=uinput.dom(2,1);
d=uinput.dom(2,2);
N=psize.N;



psize_aux0=[0 psize.nx psize.ny psize.n];
nsump0=cumsum(psize_aux0);

% Define degree of smoothness p of 2D-1var states
psize_aux1=[1 psize.nx];
nsum=cumsum(psize.nx);
nsump1=cumsum(psize_aux1);
for i=1:length(psize.nx);
px(nsump1(i):nsum(i))=i-1;
end

psize_aux1=[1 psize.ny];
psize_aux0=[0 psize.ny];
nsum=cumsum(psize.ny);
nsump1=cumsum(psize_aux1);
for i=1:length(psize.ny);
py(nsump1(i):nsum(i))=i-1;
end

% Define degree of smoothness p of 2D-2var states
psize_aux1=[1 psize.n];
psize_aux0=[0 psize.n];
nsum=cumsum(psize.n);
nsump1=cumsum(psize_aux1);
for i=1:length(psize.n);
p(nsump1(i):nsum(i))=i-1;
end


% Define initial conditions on states on the physcal grid for states.
% var_f denotes the value of the solution variable of the fundamental
% states.
% Define Chebyshev coefficients of initial conditions in the same loop
% acheb_f denotes the Chebyshev coefficients of the fundamental states
% Each state vector array coefficients are arranged into a global column
% vector

ic=uinput.ic.PDE;

acheb_glob_x{1}=[];
acheb_glob_y{1}=[];
acheb_glob{1}=[];

for i=1:length(psize.nx)
 acheb=double.empty(N-i+2,0);
 for n=1:psize.nx(i);
     var(:,n)=double(subs(ic(n+nsump0(i)),gridall.x{i}));
     acheb(:,n)=fcht(var(:,n));
 end
 acheb_glob_x{i}=reshape(acheb, [], 1);
 clear('acheb');
 clear('var');
end

for i=1:length(psize.ny)
    ii=i+length(psize.nx);
 acheb=double.empty(N-i+2,0);
 for n=1:psize.ny(i);
     var(:,n)=double(subs(ic(n+nsump0(ii)),gridall.y{i}));
     acheb(:,n)=fcht(var(:,n));
 end
 acheb_glob_y{i}=reshape(acheb, [], 1);
 clear('acheb');
 clear('var');
 end


 for i=1:length(psize.n)
    ii=i+length(psize.nx)+length(psize.ny);
 acheb=double.empty((N-i+2)*(N-i+2),0);
 for n=1:psize.n(i)
     var=double(subs(subs(ic(n+nsump0(ii)),sx,gridall.x{i}),sy,gridall.y{i}'));
     aacheb=fcgltran2d(var,1);
     acheb(:,n)=aacheb(:);
 end
 acheb_glob{i}=reshape(acheb, [], 1);
 clear('acheb');
 clear('var');
 end

% Concatenate coefficients of all states into a single column vector

 acheb_f0=cat(1, acheb_glob_x{:}, acheb_glob_y{:}, acheb_glob{:});

% Add initial conditions on ODE states to the front of initial conditions
if (no>0)
acheb_f0=cat(1,uinput.ic.ODE',acheb_f0);
end

coeff.acheb_f0=acheb_f0;

% Discretize spatial contribution of u and w disturbnaces

coeff.u=1;
coeff.w=1;

 if isfield(uinput,'wspace')
     Nforce=psize.nw0+(N+1)*(psize.nwx+psize.nwy)+(N+1)^2*psize.nw2;
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
             coeff.w(index:index+N,k)=PIESIM_NonPoly2Mat_cheb(N, uinput.wspace(k), 0, gridall.x);
             index=index+N+1;
         end 
         for kk=1:psize.nwy
             k=k+1;
             coeff.w(index:index+N,k)=PIESIM_NonPoly2Mat_cheb(N, uinput.wspace(k), 0, gridall.y);
             index=index+N+1;
         end
         for kk=1:psize.nw2
             k=k+1;
             coeff.w(index:index+(N+1)^2-1,k)=PIESIM_NonPoly2Mat_cheb_2D(N, uinput.wspace(k), 0, gridall);
             index=index+(N+1)^2;
         end
 end

 if isfield(uinput,'uspace')
     Nforce=psize.nu0+(N+1)*(psize.nux+psize.nuy)+(N+1)^2^psize.nu2;
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
             coeff.u(index:index+N,k)=PIESIM_NonPoly2Mat_cheb(N, uinput.uspace(k), 0, gridall.x);
             index=index+N+1;
         end 
         for kk=1:psize.nwy
             k=k+1;
             coeff.u(index:index+N,k)=PIESIM_NonPoly2Mat_cheb(N, uinput.uspace(k), 0, gridall.y);
             index=index+N+1;
         end
         for kk=1:psize.nw2
             k=k+1;
             coeff.u(index:index+(N+1)^2-1,k)=PIESIM_NonPoly2Mat_cheb_2D(N, uinput.uspace(k), 0, gridall);
             index=index+(N+1)^2;
         end
 end
