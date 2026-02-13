%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_NonPoly2Mat_cheb.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representation of a non-polymonial operator
% Discretizes FULL OPERATOR but does not provide a PIE-2-PDE state map 
%
% Inputs:
% 1) N   - polynomial order of Chebyshev discretization polynomial
% 2) Rop -  non-polynomial matrix operator 
% 3) p - scalar - a "degree of smoothness" vector
% 4) grid - field containing the following sub-fields:
% grid.phys - physical grid for states differentiable up to order zero (corresponding to a primary = PDE state discretization)
% grid - cell array containing grids of different degrees of differentiability
%
% Outputs:
% 1) A - matrix that represents a Chebyshev
% discretizaiton of a non-polynomial matrix operator. 
% 
% Requires: fcht
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 11_27_2021
% DJ, 12/16/2024: Remove redundant variables.

function A=PIESIM_NonPoly2Mat_cheb(N, Rop, p, grid)

ns=size(p,2);
no=size(Rop,2);

for m=ns:-1:1
    rsize=N-p(m)+1;

    Alocal(1:rsize,1:no)=0;
    for i=1:no 
        var_force=double(subs(Rop(m,i),grid{p(m)+1}));
        Alocal(:,i)=fcht(var_force);
    end % i loop

    A_cell{m,1}=Alocal;

end % m loop

A=cell2mat(A_cell);
