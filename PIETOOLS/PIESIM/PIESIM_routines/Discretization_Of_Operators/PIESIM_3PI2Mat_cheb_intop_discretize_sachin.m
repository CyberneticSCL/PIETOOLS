%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_3PI2Mat_cheb_intop_discretize.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - a dicrete matrix representation of an integrative portion
% of the 3PI operator
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% R1 - integrative operator R1 of size 1x1
% R2 - integrative operator R2 of size 1x1
% p - scalar - a "degree of smoothness" for the function on which
% R1, R2 operators act
%
% Output:
% A - a non-square matrix of size (N+1)x(N-p+1) that represents a Chebyshev
% discretizaiton of the integrative portion of the PI operator
%
% Requires: chebyshevT 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 5_29_2021

function A=PIESIM_3PI2Mat_cheb_intop_discretize_new(N, R1, R2, p)
syms sym_s sym_theta

Norder=N-p;

A(1:N+1,1:Norder+1)=0;

for k=0:Norder
% Coefficient of ak
ak=int(R1*chebyshevT(k,sym_theta),sym_theta,[-1,sym_s])+int(R2*chebyshevT(k,sym_theta),sym_theta,[sym_s,1]);

b=coeffs(ak,'all');

Tpoly = [b];
N = size(b,2);
if (N<2) 
    Tpoly = [0 b];
    if isempty(b) 
        Tpoly =[0 0];
    end
    N=2;
end
C = zeros(N,N);
C(1,1) = 1;
C(2,2) = 1;
for n = 3:N
  C(n,:) = [0,2*C(n-1,1:N-1)] - C(n-2,:);
end
C = fliplr( C );
TpolyLC = Tpoly/C;


for j=1:size(TpolyLC,2)
A(j,k+1)=TpolyLC(1,j);
end

end

A=double(A);