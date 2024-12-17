%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_NonPoly2Mat_cheb.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representation of a non-polymonial operator
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
function A=PIESIM_NonPOLy2Mat_cheb(N, Rop, p, gridall)
pvar s theta
syms sym_s sym_theta
ns=size(p,2);
no=size(Rop,2);


for m=ns:-1:1
    rsize=N-p(m)+1;

    A(1:rsize,1:no)=0;
    for i=1:no 
    var_force=double(subs(Rop(m,i),gridall{p(m)+1}));
    A(:,i)=fcht(var_force);
    end % i loop

    A_cell{m,1}=A;

end % m loop

A=cell2mat(A_cell);
