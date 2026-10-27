%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_2Dto1D2Mat_cheb_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A, A_2PDEstate - discrete matrix representations of 2D to 1D
% operators (Rx2, Ry2)
%
% Called by 'PIESIM_fullPI2Mat_cheb_2D.m'
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% Rop -  a 1D to 2D component of the full 2D PI operator, corresponding to
% either Rx2 or Ry2
% var1: 2x1 pvar array specifying the primary spatial variables (x,y)
% p - vector with the degrees of differentiability for the 2D states
% p1 - vector with the degrees of differentiability for the 1D states
% (differentiable in x for Rx2 and in y for Ry2)
% dir - direction ('x' for Rx2 and 'y' for Ry2)
%
% Outputs:
% A - Chebyshev
% discretizaiton of a Rx2/Ry2 block that represents a block of a
% total matrix operator for time-advancement of the
% spatially-discretized PIE solution 
% A_2PDEstate - Chebyshev
% discretizaiton of a Rx2/Ry2 block that represents a block of a
% total matrix operator for reconstruction of the primary (PDE)
% solution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024
% DJ, 12/16/2024: Add var1 as the actual variables in which Rop is defined, 
%                   to pass to the discretizer.
% YP, 2/2/2026: made global assembly here as opposed to component codes

function [A, A_2PDEstate]=PIESIM_2Dto1D2Mat_cheb_2D(N, Rop, var1, prow, pcol, dir)

ns_row=size(prow,2);
ns_col=size(pcol,2);

  R0=Rop{1};
  R1=Rop{2};
  R2=Rop{3};

if dir=='x'
    Nrow=N(1);
else
    Nrow=N(2);
end


%         % Go over the rows of the PI operator row by row and column by column

  for i=1:ns_row
     for j=1:ns_col       

% Treatment of mutliplicative operator (opmult stands for multiplicative)

   if (~isempty(R0)) 
       A_opmult=PIESIM_3PI_op_line_int_discretize_2D(N, R0(i,j), var1,dir, 0);
   end

% Treatment of integrative operators 

   if (~isempty(R1)) 
       A_opint_block_1=PIESIM_3PI_op_line_int_discretize_2D(N, R1(i,j), var1, dir, [1 0]);
   end

   if (~isempty(R2)) 
       A_opint_block_2=PIESIM_3PI_op_line_int_discretize_2D(N, R2(i,j), var1, dir, [0 1]);
   end

   A_sum=A_opmult+A_opint_block_1+A_opint_block_2;

% Truncate the total operator to the correct dimensions

% First, PIE-to-PDE state matrix
A_2PDEstate_local=PIESIM_Mat_Truncate(A_sum,N,pcol(:,j),'col');
% Next, PIE time propagator matrix
A_local=A_2PDEstate_local(1:Nrow-prow(i)+1,:);

 % Create a cell array for the full operator
         A_cell{i,j}=A_local;
         A_cell_2PDEstate{i,j}=A_2PDEstate_local;
     end
end

% Concatenating from cellular to a global matrix structure 

 A=cell2mat(A_cell);
 A_2PDEstate=cell2mat(A_cell_2PDEstate);