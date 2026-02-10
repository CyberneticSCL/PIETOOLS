%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_9PI_opmult_discretize_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A, A_nonsquare - discrete matrix representations of a polymonial
% multiplicative block (R00) of 9PI operator in 2D
% Discretizes an operator component (not a full operator)
% Called by 'PIESIM_9PI2Mat_cheb_2D.m' and 'PIESIM_9PI2Mat_cheb_opint_discretize_line.m'
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% R -  polymonial multiplicative block (R00) of 9PI operator in 2D
% corresponding to a single solution state
% var1: 2x1 pvar array specifying the primary spatial variables (x,y)
%
% Outputs:
% A - Chebyshev
% discretizaiton of a R00 block that represents a block of a
% spatially-discretized total matrix operator - size (N+1)x(N+1)
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
% YP, 1/31/2026 - changed operator constuction to square matrices, truncation will be done later

function A=PIESIM_9PI_opmult_discretize_2D(N, R, var1)

% Extract second spatial variable (y) defining the operator R               % DJ, 12/16/2024
s2 = var1(2);

A=0;
for k=1:R.nterms
    Rstrip=polynomial(R.coeff(k),R.degmat(k,:),R.varname,R.matdim);
        Rs1=subs(Rstrip,s2,1);
        % Multiplicative operator in s1 direction
        mult_s1=PIESIM_3PI_opmult_discretize(N(1), Rs1);
        index = find(strcmp(Rstrip.varname, s2.varname{1}));
        Rs2=polynomial(1,Rstrip.degmat(1,index),Rstrip.varname(index),Rstrip.matdim);
        % Multiplicative operator in s2 direction
        mult_s2=PIESIM_3PI_opmult_discretize(N(2), Rs2);
        A=A+kron(mult_s2,mult_s1);
end