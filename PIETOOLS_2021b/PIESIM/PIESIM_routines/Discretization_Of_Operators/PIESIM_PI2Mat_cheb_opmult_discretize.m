%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_PI2Mat_cheb_opmult_discretize.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A, A_nonsquare - discrete matrix representations of a polymonial
% multiplicative operator 
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% rsize - a number of rows for A matrix (component of a square
% assembled operator, defined in PESIM_3PI2Mat_cheb)
% Rop -  polynomial multplicative operator R0 of size 1x1
% p - scalar - a "degree of smoothness" for the function on which
% R0 operator acts
%
% Outputs:
% A - matrix of size N+1 x rsize that represents a Chebyshev
% discretizaiton of the polynomial multiplicative matrix operator. Represents a block of a
% square total matrix operator for PIE solution
% A_nonsquare - matrix of size N+1 x N-p+1 that represents a Chebyshev
% discretizaiton of the polynomial multiplicative matrix operator. Represents a block of a
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
% Initial coding YP  - 12_21_2021
function [A, A_nonsquare]=PIESIM_PI2Mat_cheb_opmult_discretize(N, rsize, Rop, p)

A_nonsquare(1:N+1,1:N-p+1)=0;

Rsym=multipoly2sym(Rop);

    if (isempty(Rsym))
        Rop_local=0;
    else
        Rop_local=Rsym;
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

for i=1:N-p+1
    id=i-1;
    for j=1:length(TpolyLC)
        jd=j-1;
        if (jd+id<=N)
        A_nonsquare(jd+id+1,i)=A_nonsquare(jd+id+1,i)+0.5*TpolyLC(j);
        end
        if (abs(jd-id)<=N)
        A_nonsquare(abs(jd-id)+1,i)=A_nonsquare(abs(jd-id)+1,i)+0.5*TpolyLC(j);
        end
    end
end
A=A_nonsquare(1:rsize,:);



