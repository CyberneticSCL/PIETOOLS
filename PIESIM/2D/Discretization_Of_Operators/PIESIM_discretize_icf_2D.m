%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_discretize_icf_2D.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform discretization of initial conditions and forcing functions in 2D
%
% Inputs:
% 1) uinput - user-defined boundary inputs, forcing and initial conditions
% 2) psize - size of the PIE problem: all variables defining the size of the PIE problem
% 3) grid - physical and computational grid for states differentiable up to order zero (corresponding to a primary = PDE state discretization)
% 4) gridall - cell array of size dmax containing physical grid for all states
% depending on their degree of differentiability; dmax corresponds to the
% maximum degree of differentiability among all the states
%  gridall.x - grids in x direction
%  gridall.y - grids in y direction
% 5) B1_neq0 - nx x nw array of type 'logical', specifying for each row
%    i in (1:nx) and each column j in (1:nw) whether the element B1(i,j)
%    of the PI operator B1 is nonzero, i.e. whether the disturbance w
%    contributes to the state equation i.
%
% Outputs:
% 1) coeff - Chebyshev coefficients for initial conditions and forcing functions
% 2) B1_nonpol - contribution to the PIE B1 operator arising from non-polynomial in space forcing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  04_16_2024
% DJ, 07/16/2024 - add input "B1_neq0", to keep track of which state
%                   equation each disturbance contributes to.
%                  Also, add support for forcing to multiple state
%                  components in multiple variables.

function [coeff,B1_nonpol]=PIESIM_discretize_icf_2D(uinput,psize,grid,gridall,B1_neq0)

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
    ii=i+sum(psize.nx);
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
    ii=i+sum(psize.nx)+sum(psize.ny);
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


% Discretize matrix operator for non-polynomial in space forcing
% Only assume forcing on one variable at this point 
if isfield(uinput,'wspace')
    ns_list = [psize.no; sum(psize.nx); sum(psize.ny); sum(psize.n)];
    ns_list_cum = cumsum(ns_list);
    Nonpol=sx*sy*zeros(ns_list_cum(end),nw);
    % Set the non-polynomial spatial contribution of each disturbance to
    % each state variable.
    for k=1:nw
        Nonpol(B1_neq0(:,k),k)=uinput.wspace(k);
    end
    % Express contribution in terms of Chebyshev coefficients.
    % First set contribution to ODE states.
    B1_nonpol = double(Nonpol(1:ns_list_cum(1),:));
    % Then contribution to 1D PDE states in just x.
    if ns_list_cum(1)<ns_list_cum(2)
        B1_nonpol = [B1_nonpol;
                    PIESIM_NonPoly2Mat_cheb(N, nw, Nonpol(ns_list_cum(1)+1:ns_list_cum(2),:), px, gridall.x)];
    end
    % Then contribution to 1D PDE states in just y.
    if ns_list_cum(2)<ns_list_cum(3)
        B1_nonpol = [B1_nonpol; 
                    PIESIM_NonPoly2Mat_cheb(N, nw, Nonpol(ns_list_cum(2)+1:ns_list_cum(3),:), py, gridall.y)];
    end
    % And finally contribution to 2D PDE states.
    if ns_list_cum(3)<ns_list_cum(4)
        B1_nonpol = [B1_nonpol; 
                    PIESIM_NonPoly2Mat_cheb_2D(N, nw, Nonpol(ns_list_cum(3)+1:ns_list_cum(4),:), p, gridall)];
    end
else
    B1_nonpol=[];   
end