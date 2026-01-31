%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_4PI2Mat_cheb.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs a discrete version of the multi-dimensional 4PI operator  
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% Rop - 4PI operator 
% p - vector of dimension 1xns -
% a "degree of smoothness" structure for 3PI, see Peet & Peet 2021 paper
%
% flag = 0 if a structure maps disturbances or control inputs to regulated or observed outputs (for D11, D12,
% D21, D22 operators)
% flag = 1 if a structure maps solution states (ODE+PDE) to regulated or observed outputs (for C1, C2 operators)
% flag = 2 if a structure maps disturbances or control inputs to solution states (ODE+PDE) (for Tw, Tu, B1 and B2
% operators)
% flag = 3 if a structure maps solution states (ODE+PDE) to solution states (ODE+PDE) (for A and T)
%
%
%
% Outputs:
% A - discretization matrix of the time-dependent PIE propagator
% A_2PDEstate -discretization matrix for transform between
% fundamental and primary solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 5_31_2021
% YP - added functionality to support infinite-dimensional disturbances and
% regulated/observed outputs - 6_2_2025

function [A, A_2PDEstate]=PIESIM_4PI2Mat_cheb(N, Rop, p, flag)

% Iniatilize the matrix blocks of the full PI operator as empty
 Pblock=[];
 Q1block=[];
 Q2block=[];
 Rblock=[];

 Q2block_2PDEstate=[];
 Rblock_2PDEstate=[];

% Pblock: size nscalar1 x nscalar2
Pblock=double(Rop.P);

% Q1block: size nscalar1 x nspat2


if (size(Rop.Q1)~=0)

if(mod(flag,2)==0)
pcol=zeros(size(Rop.Q1,2),1);
else
pcol=p;
end

Q1block=PIESIM_PI2Mat_cheb_opint_discretize(N, Rop.Q1, pcol);
end


% Q2block: size nspat1 x nscalr2

% Q2block_2PDEstate: nspat1_full x nscalar2


if (size(Rop.Q2)~=0)

    if (flag<=1)
        prow=zeros(size(Rop.Q2,1),1);
    else
        prow=p;
    end


if (nargout==1)
 Q2block= PIESIM_Poly2Mat_cheb(N, Rop.Q2, prow);
else
[Q2block, Q2block_2PDEstate]= PIESIM_Poly2Mat_cheb(N, Rop.Q2, prow);
end

end % size(Rop.Q2)~=0


% Rblock: size  nspat1 x nspat2


if (size(Rop.R.R0)~=0)

switch flag
    case 0
        prow=zeros(size(Rop.R.R0,1),1);
        pcol=zeros(size(Rop.R.R0,2),1);
    case 1
        prow=zeros(size(Rop.R.R0,1),1);
        pcol=p;
    case 2
        prow=p;
        pcol=zeros(size(Rop.R.R0,2),1);
    case 3
        prow=p;
        pcol=p;

end
      

if (nargout==1)
Rblock=PIESIM_3PI2Mat_cheb(N, Rop.R, prow, pcol);
else
[Rblock, Rblock_2PDEstate]=PIESIM_3PI2Mat_cheb(N, Rop.R, prow, pcol);
end

end % size(Rop.R.R0,2)~=0


A=[Pblock Q1block; Q2block Rblock];
if (nargout>1)
    A_2PDEstate=[Pblock Q1block; Q2block_2PDEstate Rblock_2PDEstate];
end






