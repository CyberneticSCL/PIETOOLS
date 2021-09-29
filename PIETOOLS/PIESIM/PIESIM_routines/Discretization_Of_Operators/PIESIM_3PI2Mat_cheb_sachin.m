%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_3PI2Mat_cheb.m     PIETOOLS 2021a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A and A_nonsquare matrices for the entire multi-dimensional 3PI operator  
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% Rop - 3PI operator 
% p - vector of dimension 1xns -
% a "degree of smoothness" structure for 3PI, see Peet & Peet 2021 paper
%
% Outputs:
% A - square discretization matrix for PIE solution
% A_nonsquare -non-square discretization matrix for transform between
% funamental and primary solution
%
% Requires: chebyshevT 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 5_29_2021
function [A, A_nonsquare]=PIESIM_3PI2Mat_cheb_new(N, Rop, p)
syms sym_s sym_theta

ns=size(p,2);
[nso,nsi]=size(Rop); % need to add size method

  R0=Rop.R.R0;
  R1=multipoly2sym(Rop.R.R1);
  R2=multipoly2sym(Rop.R.R2);

for i=1:nso
     for j=1:nsi
%         % Go over the rows of the matrices Tau and A column by column

%        rsize is the number of rows in each block 
%        csize is the number of columns in each block
%         
         rsize=N-p(i)+1;
         csize=N-p(j)+1;
%        
% Treatment of mutliplicative operator (opmult stands for multiplicative)
         A_opmult_cell{i,j}=R0(i,j)*eye(rsize,csize);
         A_opmult_cell_nonsquare{i,j}=R0(i,j)*eye(N+1,csize);
 % Treatment of integrative operator (opint stands for integrative)
%          A_opint_block_nonsquare(:,:,i,j)=PIESIM_3PI2Mat_cheb_intop_discretize(N,R1(i,j),R2(i,j),p(j)); 
%          A_opint_block(:,:,i,j)=A_opint_block_nonsquare(1:rsize,1:csize,i,j);
%          A_opint_cell{i,j}=A_opint_block(:,:,i,j);
%          A_opint_cell_nonsquare{i,j}=A_opint_block_nonsquare(:,:,i,j);
%             
         A_opint_block_nonsquare(1:N+1,1:csize,i,j)=PIESIM_3PI2Mat_cheb_intop_discretize(N,R1(i,j),R2(i,j),p(j)); 
         A_opint_block(1:rsize,1:csize,i,j)=A_opint_block_nonsquare(1:rsize,1:csize,i,j);
         A_opint_cell{i,j}=A_opint_block(1:rsize,1:csize,i,j);
         A_opint_cell_nonsquare{i,j}=A_opint_block_nonsquare(1:N+1,1:csize,i,j);
 % Summing multiplicative and integrative together
         A_cell{i,j}=double(A_opmult_cell{i,j}+A_opint_cell{i,j});
         A_cell_nonsquare{i,j}=double(A_opmult_cell_nonsquare{i,j}+A_opint_cell_nonsquare{i,j});
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

