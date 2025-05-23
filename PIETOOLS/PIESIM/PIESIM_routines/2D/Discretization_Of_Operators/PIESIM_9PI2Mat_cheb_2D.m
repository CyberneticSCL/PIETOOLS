%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_9PI2Mat_cheb_2D.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A and A_nonsquare matrices for the 9PI R22 (2D-2D) operator  

% Inputs:
% Required:
% - N:      polynomial order of Chebyshev discretization polynomial
% - Rop:    9PI operator 
% - var1: 2x1 pvar array specifying the primary spatial variables (x,y)
% - var2: 2x1 pvar array specifying the dummy variables for integration
%           (theta,nu)
% - p:      vector of dimension 1 x number of rows of the entries of Rop
%           (1 x size(Rop{},1))- a "degree of smoothness" structure for 
%           variables corresponding to the rows of Rop;
% Optional:
% - pcol:   vector of dimension 1 x number of columns of the entries of Rop
%           - a "degree of smoothness" structure for variables 
%           corresponding to the columns of Rop (1 x size(Rop{},2)). 
%           If not supplied, it is assumed that pcol=p;
%
% Outputs:
%
% A - block of a square discretization matrix that transforms spatio-temporal PIE into
% a square system of temporal ODEs 
% A_nonsquare - block of a nonsquare discretization matrix for transform between
% PDE and PIE solution states
% NOTE: nonsquare matrix is only required for Mcheb matrix (discrete form
% of T operator)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024
% DJ, 12/16/2024: Remove redundant variables, and add var1 and var2 as the
%                   actual variables in which Rop is defined, to pass to
%                   e.g. PIESIM_9PI2Mat_cheb_opint_discretize
function [A, A_nonsquare]=PIESIM_9PI2Mat_cheb_2D(N,Rop,var1,var2,p,pcol)

if nargin<=5
    pcol=p;
end

syms sx sy

ns_row=size(p,2);
ns_col=size(pcol,2);

  R00=Rop{1,1};
  R01=Rop{1,2};
  R02=Rop{1,3};
  R10=Rop{2,1};
  R11=Rop{2,2};
  R12=Rop{2,3};
  R20=Rop{3,1};
  R21=Rop{3,2};
  R22=Rop{3,3};

for i=1:ns_row
     for j=1:ns_col
%         % Go over the rows of the matrices T and A column by column

%        rsize is the number of rows in each block 
%        csize is the number of columns in each block
%         
         rsize=N-p(i)+1;
         csize=N-pcol(j)+1;

% Treatment of mutliplicative operator (opmult stands for multiplicative)
        if (~isempty(R00))
        [A_opmult, A_opmult_nonsquare]=PIESIM_PI2Mat_cheb_opmult_discretize_2D(N,R00(i,j),rsize, csize,var1);
         end
 % Treatment of integrative operator (opint stands for integrative)
         if (~isempty(R10))
         [A_opint_block_x1,A_opint_block_nonsquare_x1]=PIESIM_9PI2Mat_cheb_opint_discretize_line(N, R10(i,j), rsize,csize,var1,var2,'x',-1);
         end 
         if (~isempty(R20))
         [A_opint_block_x2,A_opint_block_nonsquare_x2]=PIESIM_9PI2Mat_cheb_opint_discretize_line(N, R20(i,j), rsize,csize,var1,var2,'x',1);
         end 
         if (~isempty(R01))
         [A_opint_block_y1,A_opint_block_nonsquare_y1]=PIESIM_9PI2Mat_cheb_opint_discretize_line(N, R01(i,j),rsize,csize,var1,var2,'y',-1); 
         end
         if (~isempty(R02))
         [A_opint_block_y2,A_opint_block_nonsquare_y2]=PIESIM_9PI2Mat_cheb_opint_discretize_line(N, R02(i,j),rsize,csize,var1,var2,'y',1); 
         end
         if (~isempty(R11))
         [A_opint_block_11,A_opint_block_nonsquare_11]=PIESIM_9PI2Mat_cheb_opint_discretize_area(N, R11(i,j),rsize,csize,var1,var2,[-1 sx;-1 sy]); 
         end
          if (~isempty(R21))
          [A_opint_block_21,A_opint_block_nonsquare_21]=PIESIM_9PI2Mat_cheb_opint_discretize_area(N, R21(i,j),rsize,csize,var1,var2,[sx 1;-1 sy]); 
          end
          if (~isempty(R12))
          [A_opint_block_12,A_opint_block_nonsquare_12]=PIESIM_9PI2Mat_cheb_opint_discretize_area(N, R12(i,j),rsize,csize,var1,var2,[-1 sx;sy 1]); 
          end
          if (~isempty(R22))
          [A_opint_block_22,A_opint_block_nonsquare_22]=PIESIM_9PI2Mat_cheb_opint_discretize_area(N, R22(i,j),rsize,csize,var1,var2,[sx 1;sy 1]); 
          end

         A_opint=A_opint_block_x1+A_opint_block_x2+ ...
         A_opint_block_y1+A_opint_block_y2+ ...
         A_opint_block_11+A_opint_block_12+...
         A_opint_block_21+A_opint_block_22;
         A_opint_nonsquare=A_opint_block_nonsquare_x1+A_opint_block_nonsquare_x2+...
         A_opint_block_nonsquare_y1+A_opint_block_nonsquare_y2+...
         A_opint_block_nonsquare_11+A_opint_block_nonsquare_12+...
         A_opint_block_nonsquare_21+A_opint_block_nonsquare_22;
 % Summing multiplicative and integrative together
         A_cell{i,j}=double(A_opmult+A_opint);
         A_cell_nonsquare{i,j}=double(A_opmult_nonsquare+A_opint_nonsquare);
%         
     end
end

% Concatenating from cellular to a global matrix structure 

%  A_nonsquare is a nonsquare matrix for the transfer between fundamental and primary
%  state coefficients

 A_nonsquare=cell2mat(A_cell_nonsquare);
% 
% %  A is a square matrix for the PIE (fundamental state
% %  coefficients)
% 
 A=cell2mat(A_cell);

