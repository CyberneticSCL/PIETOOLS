%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_Poly2Mat_cheb.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representation of a polymonial operator in
% 2D
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% Rop -  polynomial matrix operator 
% var1: 2x1 pvar array specifying the primary spatial variables (x,y)
% p - scalar - a "degree of smoothness" vector
%
% Outputs:
% A - block of a discrete matrix that represents a Chebyshev
% discretizaiton of the polynomial matrix operator. Represents a block of a
% square total matrix operator for PIE solution
% A_nonsquare - block of a discrete matrix that represents a Chebyshev
% discretizaiton of the polynomial matrix operator. Represents a block of a
% nonsquare total matrix operator for reconstruction of the primary (PDE)
% solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024
% DJ, 07/17/2024: bugfix in case "Rop" is not polynomial;
% DJ, 12/16/2024: Correct limits of for loop j=1:size(acheb);
% DJ, 12/16/2024: Remove hard-coded variables. Instead, pass variables
%                   defining R as additional inputs.

function [A, A_nonsquare]=PIESIM_Poly2Mat_cheb_2D(N, Rop, var1, p)

% Extract spatial variables (x,y) defining the operator R                   % DJ, 12/16/2024
s1 = var1(1);   s2 = var1(2);

ns=size(p,2);
no=size(Rop,2);

for m=1:ns

rsize=N-p(m)+1;

chebgrid=cos(pi*(0:rsize-1)/(rsize-1));
        
A_block(1:rsize^2,1:no,m)=0;
A_block_nonsquare(1:(N+1)^2,1:no,m)=0;


for k=1:size(Rop,2)
    if (isempty(Rop))
        Rop_local=0;
    else
        Rop_local=Rop(m,k);
    end
    if isa(Rop_local,'polynomial')
    Reval=subs(subs(Rop_local,s1,chebgrid)',s2,chebgrid);
    else
    Reval=Rop_local*ones(rsize,rsize);
    end

    acheb=reshape(fcgltran2d(double(Reval),1),[],1);

for j=1:length(acheb)                                                       % DJ, 12/16/2024
A_block(j,k,m)=acheb(j);
A_block_nonsquare(j,k,m)=acheb(j);
end

end % k loop

A_cell{m,1}=A_block(1:rsize^2,1:no,m);
A_cell_nonsquare{m,1}=A_block_nonsquare(1:(N+1)^2,1:no,m);


end % m loop


A_nonsquare=cell2mat(A_cell_nonsquare);
A=cell2mat(A_cell);
