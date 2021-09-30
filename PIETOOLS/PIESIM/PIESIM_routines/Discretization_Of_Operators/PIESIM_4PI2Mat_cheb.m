%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_4PI2Mat_cheb.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs a discrete version of the multi-dimensional 4PI operator  
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% nx - number of ODE states if flag=1
% nx - number of disturbances or control inputs if flag=0
% Rop - 4PI operator 
% p - vector of dimension 1xns -
% a "degree of smoothness" structure for 3PI, see Peet & Peet 2021 paper
% flag = 1 if a structure is a full 4PI operator (for A and T)
% flag = 0 if a structure has an empty right side (for Tw, Tu, B1 and B2
% operators)
%
% Outputs:
% A - square discretization matrix for PIE solution
% A_nonsquare -non-square discretization matrix for transform between
% funamental and primary solution
% NOTE: nonsquare matrix is only required for Mcheb matrix
%
% Requires: chebyshevT 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 5_31_2021
function [A, A_nonsquare]=PIESIM_4PI2Mat_cheb(N, nx, Rop, p, flag)

% Pblock: size nx x nx
Pblock=double(Rop.P);

% Q1block: size nx x n0(N+1) n1N n2(N-1)

if (flag==1)||(flag==2)
    Q1block=PIESIM_PI2Mat_cheb_intop_discretize(N, nx, Rop.Q1, p);
end


% Q2block: size n0(N+1) n1N n2(N-1) x nx
% Rblock: size  n0(N+1) n1N n2(N-1) x n0(N+1) n1N n2(N-1)

if (nargout==1)
    if flag~=2
        Q2block= PIESIM_Poly2Mat_cheb(N, nx, Rop.Q2, p);
    end
if (flag==1)
    Rblock=PIESIM_3PI2Mat_cheb(N, Rop, p);
end
else
% Q2block_nonsquare: size n0n1n2(N+1)^3 x nx
% Rblock_nonsquare: size  n0n1n2(N+1)^3 x n0(N+1) n1N n2(N-1)
if flag~=2
[Q2block, Q2block_nonsquare]= PIESIM_Poly2Mat_cheb(N, nx, Rop.Q2, p);
end
if (flag==1)
[Rblock, Rblock_nonsquare]=PIESIM_3PI2Mat_cheb(N, Rop, p);
end
end

if (flag==1)
A=[Pblock Q1block; Q2block Rblock];

if (nargout>1)
    A_nonsquare=[Pblock Q1block; Q2block_nonsquare Rblock_nonsquare];
end

elseif (flag==2)
    A=[Pblock, Q1block];
else
    A=[Pblock; Q2block];
end
    




