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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 06_28_2022
function [A, A_nonsquare]=PIESIM_PI2Mat_cheb_opmult_discretize(N, rsize, Rop, p)

pvar s;

A_nonsquare(1:N+1,1:N-p+1)=0;

if isa(Rop,'polynomial')
    deg=Rop.maxdeg;
else
    deg=1;
end

chebgrid=cos(pi*(0:deg+1)/(deg+1));
if isa(Rop,'polynomial')
    Reval=subs(Rop,s,chebgrid);
else
    Reval=Rop*ones(1,deg+2);
end

acheb=fcht(double(Reval));

for i=1:N-p+1
    id=i-1;
    for j=1:length(acheb)
        jd=j-1;
        if (jd+id<=N)
        A_nonsquare(jd+id+1,i)=A_nonsquare(jd+id+1,i)+0.5*acheb(j);
        end
        if (abs(jd-id)<=N)
        A_nonsquare(abs(jd-id)+1,i)=A_nonsquare(abs(jd-id)+1,i)+0.5*acheb(j);
        end
    end
end
A=A_nonsquare(1:rsize,:);



