%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_PI2Mat_cheb_opint_discretize_2to0.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A, A_nonsquare - discrete matrix representations of a 2D to 0D
% operator (R02)
%
% Called by 'PIESIM_fullPI2Mat_cheb_2D.m'
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% R -  polymonial block corresponding to R20 operator
% var1: 2x1 pvar array specifying the primary spatial variables (x,y)
% p - vector with the degrees of differentiability for the 2D states
%
%
% Outputs:
% A - Chebyshev discretizaiton of R02 block that represents a block of a
% square total matrix operator for time-advancement of the
% spatially-discretized PIE solution (square ODE system matrix)
% A_nonsquare - Chebyshev
% discretizaiton of R02 block that represents a block of a
% nonsquare total matrix operator for reconstruction of the primary (PDE)
% solution (nonsquare transformation matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024
% DJ, 12/16/2024: Remove hard-coded variables. Instead, pass variables
%                   defining R as additional inputs.

function A=PIESIM_PI2Mat_cheb_opint_discretize_2to0(N, R, var1, p)

% Extract second spatial variable (y) defining the operator R               % DJ, 12/16/2024
s2 = var1(2);

ns=size(p,2);
no=size(R,1);

for i=1:no
    Astring=zeros(1,0);
    for j=1:ns
        Rloc=R(i,j);
        int=0;
    for k=1:Rloc.nterms
        Rstrip=polynomial(Rloc.coeff(k),Rloc.degmat(k,:),Rloc.varname,Rloc.matdim);
            Rs1=subs(Rstrip,s2,1);
            % Integrate in s1 direction over the interval [-1,1]
            int_s1=PIESIM_PI2Mat_cheb_opint_discretize(N, Rs1, p(j));
            index = find(strcmp(Rstrip.varname, s2.varname{1}));
            Rs2=polynomial(1,Rstrip.degmat(1,index),Rstrip.varname(index),Rstrip.matdim);
            % Integrate in s2 direction over the interval [-1,1]
            int_s2=PIESIM_PI2Mat_cheb_opint_discretize(N, Rs2, p(j));
            int=int+kron(int_s2,int_s1);
    end
            Astring=[Astring,int];
    end
    A(i,:)=Astring;
end




