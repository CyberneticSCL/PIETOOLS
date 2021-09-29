%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_PI2Mat_cheb_intop_discretize.m     PIETOOLS 2021a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representation for the integrative Q1 operator  
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% nx - number of ODE states
% Rop - Q1 operator of dimension nx x ns
% p - vector of dimension 1xns -
% a "degree of smoothness" structure for PDE, see Peet & Peet 2021 paper
%
% Outputs:
% A - discretization matrix of the Q1 operator: 
% dimension nx x n0(N+1) n1N n2(N-1)
%
% Requires: multipoly2sym, chebyshevT 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 5_31_2021
function A=PIESIM_PI2Mat_cheb_intop_discretize(N, nx, Rop, p);
syms sym_s sym_theta

ns=size(p,2);

  Q1=multipoly2sym(Rop);
  
  for j=1:ns
        Norder=N-p(j);
       
  A=zeros(nx,Norder);
  
        for k=0:Norder
            for i=1:nx
        ak=int(Q1(i,j)*chebyshevT(k,sym_s),sym_s,[-1,1]);
       
        A(i,k+1)=ak;
            end
        end
        
        A_cell{j}=double(A);
  end
        
        A=cell2mat(A_cell);
  



