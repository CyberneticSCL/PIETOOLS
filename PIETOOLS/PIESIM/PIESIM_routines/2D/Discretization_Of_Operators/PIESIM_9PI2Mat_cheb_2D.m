%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_9PI2Mat_cheb_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A and A_2PDEstate matrices for the 9PI R22 (2D-2D) operator  
% Discretized FULL 9PI OPERATOR
%   (R u)(x,y) = R00 u(x,y) + int_{-1}^{x} R10 u(theta,y) dtheta
%     + int_{x}^{1} R20 u(theta,y) dtheta + int_{-1}^{y} R01 u(x,nu) dnu
%     + int_{y}^{1} R02 u(x,nu) dnu + int_{-1}^{x} int_{-1}^{y} R11 u dnu dtheta
%     + int_{x}^{1} int_{-1}^{y} R21 u + int_{-1}^{x} int_{y}^{1} R12 u
%     + int_{x}^{1} int_{y}^{1} R22 u,   (x,y) in [-1,1]^2.
% Block (i,j) maps the coefficients of column state j to those of row
% state i, each restricted to the coefficients resolved for its degree of
% smoothness (the rows not restricted in A_2PDEstate).
% Inputs:
% Required:
% - N:      1x2 vector [N_x, N_y], polynomial orders of the Chebyshev
%           discretization in x and y
% - Rop:    9PI operator 
%           (3x3 cell {R00 R01 R02; R10 R11 R12; R20 R21 R22} of
%           'polynomial' or 'double' matrices: first index the x part,
%           second the y part, each 1 multiplier, 2 integral from -1, 3
%           integral to 1)
% - var1: 2x1 pvar array specifying the primary spatial variables (x,y)
% - var2: 2x1 pvar array specifying the dummy variables for integration
%           (theta,nu)
% - prow:   array of dimension 2 x number of rows of the entries of Rop
%           (2 x size(Rop{},1))- a "degree of smoothness" structure for
%           variables corresponding to the rows of Rop, [p_x; p_y] per row;
% Optional:
% - pcol:   array of dimension 2 x number of columns of the entries of Rop
%           - a "degree of smoothness" structure for variables 
%           corresponding to the columns of Rop (2 x size(Rop{},2)).
%           If not supplied, it is assumed that pcol=prow;
%
% Outputs:
%
% A - block of a square discretization matrix that transforms spatio-temporal PIE into
% a square system of temporal ODEs 
% A_2PDEstate - block of a 2PDEstate discretization matrix for transform between
% PDE and PIE solution states
% NOTE: 2PDEstate matrix is only required for the discrete forms of the
% operators T, Tu, Tw (and Tmap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_9PI2Mat_cheb_2D
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or Y. Peet at ypeet@asu.edu
% Initial coding YP  - 4_16_2024
% DJ, 12/16/2024: Remove redundant variables, and add var1 and var2 as the
%                   actual variables in which Rop is defined, to pass to
%                   e.g. PIESIM_9PI2Mat_cheb_opint_discretize
% YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: output A_nonsquare
%   renamed A_2PDEstate, input p renamed prow; N is [N_x, N_y], and prow,
%   pcol are 2 x (number of states) arrays [p_x; p_y]; each block is built
%   on the full grid of prod(N+1) coefficients, then its columns are
%   truncated by pcol (A_2PDEstate) and its rows by prow (A); a term with
%   an empty kernel adds zero.
function [A, A_2PDEstate]=PIESIM_9PI2Mat_cheb_2D(N,Rop,var1,var2,prow,pcol)

% pcol = prow if not given
if nargin<=5
    pcol=prow;
end

% Symbolic x, y for the variable integration limits
syms sx sy

% Numbers of row and column states
ns_row=size(prow,2);
ns_col=size(pcol,2);

  % Kernels of the nine terms
  R00=Rop{1,1};
  R01=Rop{1,2};
  R02=Rop{1,3};
  R10=Rop{2,1};
  R11=Rop{2,2};
  R12=Rop{2,3};
  R20=Rop{3,1};
  R21=Rop{3,2};
  R22=Rop{3,3};

%         % Go over the rows of the 9PI operator row by row and column by column

%        rsize is the number of rows in each block 
%        csize is the number of columns in each block

for i=1:ns_row
     for j=1:ns_col

         % YP, MMP, 11/21/2026: BEGIN each term zero unless its kernel is given
         [A_opmult,A_opint_block_x1,A_opint_block_x2,A_opint_block_y1, ...
          A_opint_block_y2,A_opint_block_11,A_opint_block_12, ...
          A_opint_block_21,A_opint_block_22]=deal(zeros(prod(N+1)));
         % YP, MMP, 11/21/2026: END each term zero unless its kernel is given
% Treatment of mutliplicative operator (opmult stands for multiplicative)
        if (~isempty(R00))
         A_opmult=PIESIM_9PI_opmult_discretize_2D(N,R00(i,j),var1);
         end
 % Treatment of integrative operator (opint stands for integrative)
         % Line integrals in x: from -1 to x (R10), from x to 1 (R20)
         if (~isempty(R10))
         A_opint_block_x1=PIESIM_9PI_opint_line_discretize_2D(N, R10(i,j),var1,var2,'x',-1);
         end 
         if (~isempty(R20))
         A_opint_block_x2=PIESIM_9PI_opint_line_discretize_2D(N, R20(i,j), var1,var2,'x',1);
         end 
         % Line integrals in y: from -1 to y (R01), from y to 1 (R02)
         if (~isempty(R01))
         A_opint_block_y1=PIESIM_9PI_opint_line_discretize_2D(N, R01(i,j),var1,var2,'y',-1); 
         end
         if (~isempty(R02))
         A_opint_block_y2=PIESIM_9PI_opint_line_discretize_2D(N, R02(i,j),var1,var2,'y',1); 
         end
         % Area integrals over [-1,x]x[-1,y] (R11), [x,1]x[-1,y] (R21),
         % [-1,x]x[y,1] (R12), [x,1]x[y,1] (R22)
         if (~isempty(R11))
         A_opint_block_11=PIESIM_9PI_opint_area_discretize_2D(N, R11(i,j),var1,var2,[-1 sx;-1 sy]); 
         end
          if (~isempty(R21))
          A_opint_block_21=PIESIM_9PI_opint_area_discretize_2D(N, R21(i,j),var1,var2,[sx 1;-1 sy]); 
          end
          if (~isempty(R12))
          A_opint_block_12=PIESIM_9PI_opint_area_discretize_2D(N, R12(i,j),var1,var2,[-1 sx;sy 1]); 
          end
          if (~isempty(R22))
          A_opint_block_22=PIESIM_9PI_opint_area_discretize_2D(N, R22(i,j),var1,var2,[sx 1;sy 1]); 
          end

         % Sum of the nine terms, prod(N+1) x prod(N+1)
         A_sum=A_opmult+A_opint_block_x1+A_opint_block_x2+ ...
         A_opint_block_y1+A_opint_block_y2+ ...
         A_opint_block_11+A_opint_block_12+...
         A_opint_block_21+A_opint_block_22;

% Truncate the total operator to the correct dimensions

% First, PIE-to-PDE state matrix
A_2PDEstate_local=PIESIM_Mat_Truncate(A_sum,N,pcol(:,j),'col');
% Next, PIE time propagator matrix
A_local=PIESIM_Mat_Truncate(A_2PDEstate_local,N,prow(:,i),'row');

 % Create a cell array for the full operator
         A_cell{i,j}=A_local;
         A_cell_2PDEstate{i,j}=A_2PDEstate_local;
     end
end

% Concatenating from cellular to a global matrix structure 


 A=cell2mat(A_cell);
 A_2PDEstate=cell2mat(A_cell_2PDEstate);


