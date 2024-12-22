%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_PI2Mat_cheb_opmult_discretize_2D.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A, A_nonsquare - discrete matrix representations of a polymonial
% multiplicative block (R00) of 9PI operator in 2D
%
% Called by 'PIESIM_9PI2Mat_cheb_2D.m' and 'PIESIM_9PI2Mat_cheb_opint_discretize_line.m'
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% R -  polymonial multiplicative block (R00) of 9PI operator in 2D
% corresponding to a single solution state
% rsize - the number of rows in the resulting A matrix block (number of the PIE state Chebyshev coefficients per dimension in the solution state on the left-hand side of the ODE matrix system corresponding to the discrete block in question)
% csize - the number of columns in the resulting A matrix block (number of the PIE state Chebyshev coefficients per dimension in the right-hand side of the ODE matrix system corresponding to the discrete block in question)
% var1: 2x1 pvar array specifying the primary spatial variables (x,y)
%
% Outputs:
% A - Chebyshev
% discretizaiton of a R00 block that represents a block of a
% square total matrix operator for time-advancement of the
% spatially-discretized PIE solution (square ODE system matrix)
% A_nonsquare - Chebyshev
% discretizaiton of a R00 block that represents a block of a
% nonsquare total matrix operator for reconstruction of the primary (PDE)
% solution (nonsquare transformation matrix)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024
% DJ, 07/17/24  - bugfix in case "R" is not polynomial;
% DJ, 12/16/2024: Remove hard-coded variables. Instead, pass variables
%                   defining R as additional inputs.

function [A, A_nonsquare]=PIESIM_PI2Mat_cheb_opmult_discretize_2D(N, R, rsize, csize, var1)

% Extract second spatial variable (y) defining the operator R               % DJ, 12/16/2024
s2 = var1(2);
p=N+1-csize;

A=0;
A_nonsquare=0;
for k=1:R.nterms
    Rstrip=polynomial(R.coeff(k),R.degmat(k,:),R.varname,R.matdim);
        Rs1=subs(Rstrip,s2,1);
        % Multiplicative operator in s1 direction
        [mult_s1,mult_s1_nonsquare]=PIESIM_PI2Mat_cheb_opmult_discretize(N, rsize, Rs1, p);
        index = find(strcmp(Rstrip.varname, s2.varname{1}));
        Rs2=polynomial(1,Rstrip.degmat(1,index),Rstrip.varname(index),Rstrip.matdim);
        % Multiplicative operator in s2 direction
        [mult_s2,mult_s2_nonsquare]=PIESIM_PI2Mat_cheb_opmult_discretize(N, rsize, Rs2, p);
        A=A+kron(mult_s2,mult_s1);
        A_nonsquare=A_nonsquare+kron(mult_s2_nonsquare,mult_s1_nonsquare);
end