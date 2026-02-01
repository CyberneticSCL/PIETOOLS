%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_NonPoly2Mat_cheb.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representation of a non-polymonial
% operator in 2D
%
% Inputs:
% 1) N   - polynomial order of Chebyshev discretization polynomial
% 2) Rop -  non-polynomial matrix operator 
% 3) p - scalar (non-polynomial) - a "degree of smoothness" vector
% 4) gridall - cell array of size dmax containing physical grid for all states
% depending on their degree of differentiability; dmax corresponds to the
% maximum degree of differentiability among all the states
%  gridall.x - grids in x direction
%  gridall.y - grids in y direction
%
%
% Outputs:
% A - block of a discrete matrix that represents a Chebyshev
% discretizaiton of the polynomial matrix operator. Represents a block of a
% square total matrix operator for PIE solution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024
function A=PIESIM_NonPoly2Mat_cheb_2D(N, Rop, p, gridall)
syms sx sy
ns=size(p,2);
no=size(Rop,2);

for m=ns:-1:1
    rsize=N-p(m)+1;

    A(1:rsize^2,1:no)=0;
    for i=1:no 
    var_force=double(subs(subs(Rop(m,i),sx,gridall.x{p(m)+1}),sy,gridall.y{p(m)+1}'));
    A(:,i)=reshape(fcgltran2d(var_force,1),[],1);
    end % i loop

    A_cell{m,1}=A;

end % m loop

A=cell2mat(A_cell);
