%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_NonPoly2Mat_cheb.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representation of a non-polymonial operator
% Discretizes FULL OPERATOR but does not provide a PIE-2-PDE state map 
%
% Inputs:
% 1) N   - polynomial order of Chebyshev discretization polynomial
% 2) Rop -  non-polynomial matrix operator size (n0+n1+n2) x no
% 3) p - scalar - a "degree of smoothness" vector
% 4) gridall - cell array of size 3 containing physical grid for n0, n1 and
% n2 states
%
% Outputs:
% 1) A - matrix of size n0(N+1) n1N n2(N-1) x no that represents a Chebyshev
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

function A=PIESIM_NonPoly2Mat_cheb(N, Rop, p, gridall)

ns=size(p,2);
no=size(Rop,2);

for m=ns:-1:1
    rsize=N-p(m)+1;

    Alocal(1:rsize,1:no)=0;
    for i=1:no 
        var_force=double(subs(Rop(m,i),gridall{p(m)+1}));
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_NonPoly2Mat_cheb.m
        A(:,i)=fcht(var_force);
=======
        Alocal(:,i)=fcht(var_force);
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_NonPoly2Mat_cheb.m
    end % i loop

    A_cell{m,1}=Alocal;

end % m loop

A=cell2mat(A_cell);
