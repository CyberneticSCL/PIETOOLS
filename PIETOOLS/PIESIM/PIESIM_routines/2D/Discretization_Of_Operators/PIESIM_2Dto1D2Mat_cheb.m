%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_2Dto1D2Mat_cheb.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A, A_nonsquare - discrete matrix representations of 2D to 1D
% operators (Rx2, Ry2)
%
% Called by 'PIESIM_fullPI2Mat_cheb_2D.m'
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% Rop -  a 1D to 2D component of the full 2D PI operator, corresponding to
% either Rx2 or Ry2
% p - vector with the degrees of differentiability for the 2D states
% p1 - vector with the degrees of differentiability for the 1D states
% (differentiable in x for Rx2 and in y for Ry2)
% dir - direction ('x' for Rx2 and 'y' for Ry2)
%
% Outputs:
% A - Chebyshev
% discretizaiton of a Rx2/Ry2 block that represents a block of a
% square total matrix operator for time-advancement of the
% spatially-discretized PIE solution (square ODE system matrix)
% A_nonsquare - Chebyshev
% discretizaiton of a Rx2/Ry2 block that represents a block of a
% nonsquare total matrix operator for reconstruction of the primary (PDE)
% solution (nonsquare transformation matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024
function [A, A_nonsquare]=PIESIM_2Dto1D2Mat_cheb(N, Rop, p, p1, dir)

  R0=Rop{1};
  R1=Rop{2};
  R2=Rop{3};

% Treatment of mutliplicative operator (opmult stands for multiplicative)

   if (~isempty(R0)) 
       [A_opmult, A_opmult_nonsquare]=PIESIM_PI2Mat_cheb_opmixed_discretize_2to1(N, R0, p, p1, dir, 0);
   end

% Treatment of integrative operators 

   if (~isempty(R1)) 
       [A_opint_block_1,A_opint_block_nonsquare_1]=PIESIM_PI2Mat_cheb_opmixed_discretize_2to1(N, R1, p, p1, dir, [1 0]);
   end

   if (~isempty(R2)) 
       [A_opint_block_2,A_opint_block_nonsquare_2]=PIESIM_PI2Mat_cheb_opmixed_discretize_2to1(N, R2, p, p1, dir, [0 1]);
   end

   A=A_opmult+A_opint_block_1+A_opint_block_2;
   A_nonsquare=A_opmult_nonsquare+A_opint_block_nonsquare_1+A_opint_block_nonsquare_2;