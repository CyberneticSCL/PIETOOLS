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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 6_30_2022
function [A, A_nonsquare]=PIESIM_POLy2Mat_cheb(N, nx, Rop, p)
pvar s 
ns=size(p,2);

if isa(Rop,'polynomial')
    deg=Rop.maxdeg;
else
    deg=1;
end

chebgrid=cos(pi*(0:deg+1)/(deg+1));

for m=1:ns
        
A(1:N+1,1:nx,m)=0;


for k=1:nx
    if (isempty(Rop))
        Rop_local=0;
    else
        Rop_local=Rop(m,k);
    end
    if isa(Rop_local,'polynomial')
    Reval=subs(Rop_local,s,chebgrid);
    else
    Reval=Rop_local*ones(1,deg+2);
    end

    acheb=fcht(double(Reval));

for j=1:size(acheb);
A(j,k,m)=acheb(j);
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
