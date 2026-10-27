%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_3PI2Mat_cheb.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A and A_2PDEstate matrices for the entire multi-dimensional 3PI operator  
% Discretizes FULL OPERATOR
% Inputs:
% varargin - variable number of arguments (between 3 and 4)
% Required:
% 1) varargin(1): N   - polynomial order of Chebyshev discretization polynomial
% 2) varargin(2): Rop - 3PI operator 
% 3) varargin(3): p - vector of dimension 1 x number of rows of the entries of Rop (1 x size(Rop{},1))-
% a "degree of smoothness" structure for vriables corresponding to the
% rows of Rop
% Optional:
% 4) varargin(4): pcol - vector of dimension 1 x number of columns of the entries of Rop  -
% a "degree of smoothness" structure for variables corresponding to the
% columns of Rop (1 x size(Rop{},2)). If not supplied, it is assumed that pcol=p;
%
% Outputs:
% A - discretization matrix of the time-dependent PIE propagator
% A_2PDEstate -discretization matrix for transform between
% funsamental and primary solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 6_28_2022
function [A, A_2PDEstate]=PIESIM_3PI2Mat_cheb(varargin)

N=varargin{1};
Rop=varargin{2};
p=varargin{3};
if nargin==4
    pcol=varargin{4};
else
    pcol=p;
end

ns_row=length(p);
ns_col=length(pcol);


  R0=Rop.R0;
  R1=Rop.R1;
  R2=Rop.R2;


for i=1:ns_row
     for j=1:ns_col
%         % Go over the rows of the matrices Tau and A column by column

%        rsize is the number of rows in each block 
%        csize is the number of columns in each block
%         
         rsize=N-p(i)+1;
         csize=N-pcol(j)+1;
%        
% Treatment of mutliplicative operator (opmult stands for multiplicative)
         A_opmult_square=PIESIM_3PI_opmult_discretize(N, R0(i,j));
         A_opmult_2PDEstate=A_opmult_square(:,1:csize);
         A_opmult=A_opmult_2PDEstate(1:rsize,1:csize);
 % Treatment of integrative operator (opint stands for integrative)
         A_opint_square=PIESIM_3PI_opint_discretize(N,R1(i,j),R2(i,j)); 
         A_opint_2PDEstate=A_opint_square(:,1:csize);
         A_opint=A_opint_2PDEstate(1:rsize,1:csize);

 % Summing multiplicative and integrative together
         A_cell{i,j}=double(A_opmult+A_opint);
         A_cell_2PDEstate{i,j}=double(A_opmult_2PDEstate+A_opint_2PDEstate);
%         
     end
end

% Concatenating from cellular to a global matrix structure 

%  A_2PDEstate is a 2PDEstate matrix for the transfer between fundamental and primary
%  state coefficients

 A_2PDEstate=cell2mat(A_cell_2PDEstate);
% 
% %  A is a square matrix for the PIE (fundamental state
% %  coefficients)
% 
 A=cell2mat(A_cell);

