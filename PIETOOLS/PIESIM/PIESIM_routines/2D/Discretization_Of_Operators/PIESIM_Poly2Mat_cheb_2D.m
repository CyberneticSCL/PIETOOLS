%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_Poly2Mat_cheb_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representation of a polymonial operator in
% 2D
% Discretizes a FULL OPERATOR
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% Rop -  polynomial matrix operator 
% var1: 2x1 pvar array specifying the primary spatial variables (x,y)
% p - scalar - a "degree of smoothness" vector
%
% Outputs:
% A - block of a discrete matrix that represents a Chebyshev
% discretizaiton of the polynomial matrix operator. Represents a block of a
% total matrix operator for PIE solution
% A_2PDEstate - block of a discrete matrix that represents a Chebyshev
% discretizaiton of the polynomial matrix operator. Represents a block of a
% 2PDEstate total matrix operator for reconstruction of the primary (PDE)
% solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024
% DJ, 07/17/2024: bugfix in case "Rop" is not polynomial;
% DJ, 12/16/2024: Correct limits of for loop j=1:size(acheb);
% DJ, 12/16/2024: Remove hard-coded variables. Instead, pass variables
%                   defining R as additional inputs.
% YP, 1/25/2026: Fixed construction of a non-squre matrix for the 2D state
% reconstruction

function [A, A_2PDEstate]=PIESIM_Poly2Mat_cheb_2D(N, Rop, var1, p)

% Extract spatial variables (x,y) defining the operator R                   % DJ, 12/16/2024
s1 = var1(1);   s2 = var1(2);

% Evaluate everything on maximum grid and then compress

ns=size(p,2);
no=size(Rop,2);

if isa(Rop,'polynomial')
    deg=Rop.maxdeg;
else
    deg=1;
end

chebgrid=cos(pi*(0:deg+1)/(deg+1));

for m=1:ns
Alocal(1:prod(N+1),1:no)=0;

for k=1:no
    if (isempty(Rop))
        Rop_local=0;
    else
        Rop_local=Rop(m,k);
    end
    if isa(Rop_local,'polynomial')
    Reval=subs(subs(Rop_local,s1,chebgrid)',s2,chebgrid);
    else
    Reval=Rop_local*ones(deg+2,deg+2);
    end

    % Initiailize square array of size (N(1)+1)x(N(2)+1)
    acheb2D_all=zeros(N(1)+1,N(2)+1);
    acheb2D=fcgltran2d(double(Reval),1);
    acheb2D_all(1:deg+2, 1:deg+2) = acheb2D;

    acheb=reshape(acheb2D_all,[],1);

for j=1:length(acheb)                                                     
Alocal(j,k)=acheb(j);
end

% YP 1/25/2026

end % k loop

% Truncate time propagator matrix to the correct dimensions
A_propagator=PIESIM_Mat_Truncate(Alocal,N,p(:,m),'row');
A_cell{m,1}=A_propagator;
A_cell_2PDEstate{m,1}=Alocal;

end % m loop


A_2PDEstate=cell2mat(A_cell_2PDEstate);
A=cell2mat(A_cell);
