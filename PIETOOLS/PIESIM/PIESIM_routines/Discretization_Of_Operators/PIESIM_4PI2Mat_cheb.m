%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_4PI2Mat_cheb.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs a discrete version of the multi-dimensional 4PI operator  
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% Rop - 4PI operator 
% p - vector of dimension 1xns -
% a "degree of smoothness" structure for 3PI, see Peet & Peet 2021 paper
%
% flag = 0 if a structure has only the first row and is acting on disturbances or control inputs (for D11, D12,
% D21, D22 operators)
% flag = 1 if a structure is a full operator acting on disturbances or control inputs (for Tw, Tu, B1 and B2
% operators)
% flag = 2 if a structure has only the first row and acts on the PDE + ODE states (for C1, C2 operators)
% flag = 3 if a structure is a full operator acting on the PDE+ODE states (for A and T)
%
%
% Outputs:
% A - square discretization matrix for PIE solution
% A_nonsquare -non-square discretization matrix for transform between
% funamental and primary solution
% NOTE: nonsquare matrix is only required for Mcheb matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 5_31_2021
% YP - added functionality to support infinite-dimensional disturbances
% through parser - 6_1_2025

function [A, A_nonsquare]=PIESIM_4PI2Mat_cheb(N, Rop, p, flag)

% Iniatilize the matrix blocks of the full PI operator as empty
 Pblock=[];
 Q1block=[];
 Q2block=[];
 Rblock=[];

 Q2block_nonsquare=[];
 Rblock_nonsquare=[];

% Pblock: size nscalar1 x nscalar2
Pblock=double(Rop.P);

% Q1block: size nscalar1 x nspat2


if (size(Rop.Q1)~=0)

if(flag<=1)
pcol=zeros(size(Rop.Q1,2),1);
else
pcol=p;
end

Q1block=PIESIM_PI2Mat_cheb_opint_discretize(N, Rop.Q1, pcol);
end

if (mod(flag,2)==1)

% Q2block: size nspat1 x nscalr2

% Q2block_nonsquare: nspat1_full x nscalar2


if (size(Rop.Q2)~=0)

if (nargout==1)
 Q2block= PIESIM_Poly2Mat_cheb(N, Rop.Q2, p);
else
[Q2block, Q2block_nonsquare]= PIESIM_Poly2Mat_cheb(N, Rop.Q2, p);
end

end % size(Rop.Q2)~=0


% Rblock: size  nspat1 x nspat2


if (size(Rop.R.R0,2)~=0)

if(flag<=1)
pcol=zeros(size(Rop.R.R0,2),1);
else
pcol=p;
end

if (nargout==1)
Rblock=PIESIM_3PI2Mat_cheb(N, Rop.R, p, pcol);
else
[Rblock, Rblock_nonsquare]=PIESIM_3PI2Mat_cheb(N, Rop.R, p, pcol);
end

end % size(Rop.R.R0,2)~=0

end % mod(flag,2)==1



if (mod(flag,2)==1)
A=[Pblock Q1block; Q2block Rblock];
if (nargout>1)
    A_nonsquare=[Pblock Q1block; Q2block_nonsquare Rblock_nonsquare];
end
else
    A=[Pblock, Q1block];
end
    




