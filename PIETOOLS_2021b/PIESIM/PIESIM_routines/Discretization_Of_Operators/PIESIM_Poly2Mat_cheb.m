%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_Poly2Mat_cheb.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representation of a polymonial operator
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% nx - number of ODE states
% Rop -  polynomial matrix operator size (n0+n1+n2) x nx
% p - scalar - a "degree of smoothness" vector
%
% Outputs:
% A - matrix of size n0(N+1) n1N n2(N-1) x nx that represents a Chebyshev
% discretizaiton of the polynomial matrix operator. Represents a block of a
% square total matrix operator for PIE solution
% A_nonsquare - matrix of size n0n1n2(N+1)^3 x nx that represents a Chebyshev
% discretizaiton of the polynomial matrix operator. Represents a block of a
% non-square total matrix operator for reconstruction of the primary
% solution
%
%
% Requires: multipoly2sym 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 5_31_2021
function [A, A_nonsquare]=PIESIM_POLy2Mat_cheb(N, nx, Rop, p)
pvar s theta
syms sym_s sym_theta
ns=size(p,2);
for m=1:ns
        
A(1:N+1,1:nx,m)=0;

Rsym=multipoly2sym(Rop);


for k=1:nx
    if (isempty(Rsym))
        Rop_local=0;
    else
        Rop_local=Rsym(m,k);
    end
        
    if (isa(Rop_local,'double')==1)
        b=Rop_local;
    else
    b=coeffs(Rop_local,'all');
    end

Tpoly = [b];
M = size(b,2);
if (M<2) 
    Tpoly = [0 b];
    if isempty(b) 
        Tpoly =[0 0];
    end
    M=2;
end
C = zeros(M,M);
C(1,1) = 1;
C(2,2) = 1;
for n = 3:M
  C(n,:) = [0,2*C(n-1,1:M-1)] - C(n-2,:);
end
C = fliplr( C );
TpolyLC = Tpoly/C;


for j=1:size(TpolyLC,2);
A(j,k,m)=TpolyLC(1,j);
end

end % k loop

A_block_nonsquare(1:N+1,1:nx,m)=double(A(1:N+1,1:nx,m));

rsize=N-p(m)+1;

A_block(1:rsize,1:nx,m)=A_block_nonsquare(1:rsize,1:nx,m);

A_cell{m,1}=A_block(1:rsize,1:nx,m);
A_cell_nonsquare{m,1}=A_block_nonsquare(1:N+1,1:nx,m);


end % m loop


A_nonsquare=cell2mat(A_cell_nonsquare);
A=cell2mat(A_cell);