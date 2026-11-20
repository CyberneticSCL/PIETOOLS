%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_Poly2Mat_cheb.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representation of a polymonial operator
% Discretizes FULL OPERATOR
% A maps v in R^no to the Chebyshev coefficients of R(s) v, s in [-1,1]:
% column k holds, for each row m of R, the coefficients a_j of
% R(m,k)(s) = sum_{j<=N} a_j T_j(s), the first N-p(m)+1 of them in A.
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% Rop -  polynomial matrix operator size (n0+n1+n2) x no
%       ('polynomial' in one variable, or 'double'; a polynomial in more
%       than one variable throws an error)
% p - 1 x (n0+n1+n2) vector - a "degree of smoothness" vector, p(m) for
%     row m of Rop
%
% Outputs:
% A - matrix of size n0(N+1) n1N n2(N-1) x no that represents a Chebyshev
% discretizaiton of the polynomial matrix operator. Represents a block of
% the total time-propagator matrix operator for the PIE solution
% A_2PDEstate - matrix of size (n0+n1+n2)(N+1) x no that represents a Chebyshev
% discretizaiton of the polynomial matrix operator. Represents a block of
% the total matrix operator for reconstruction of the primary solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_Poly2Mat_cheb
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or Y. Peet at ypeet@asu.edu
% Initial coding YP  - 6_30_2022
% YP -  added support for an arbitary variable name - 04_16_2024
% DJ, 12/16/2024: Make sure variable to substitute matches that of Rop;
% MMP, 11/19/2026: Release comment pass (comment-only).

function [A, A_2PDEstate]=PIESIM_Poly2Mat_cheb(N, Rop, p)

% ns row states, no columns of R
ns=length(p);
no=size(Rop,2);

% deg+2 Chebyshev-Gauss-Lobatto nodes, exact for the maximal degree deg of
% R
if isa(Rop,'polynomial')
    deg=Rop.maxdeg;
else
    deg=1;
end

chebgrid=cos(pi*(0:deg+1)/(deg+1));

% Block m (row state m): column k holds the coefficients of R(m,k)
for m=1:ns
        
Alocal(1:N+1,1:no)=0;

for k=1:no
    % Values of R(m,k) at the nodes
    if (isempty(Rop))
        Rop_local=0;
    else
        Rop_local=Rop(m,k);
    end
    if isa(Rop_local,'polynomial') && ~isdouble(Rop_local)
        if numel(Rop_local.varname)>1                                       % DJ, 12/16/2024
            error("The input polynomial depends on multiple variables; discretization not supported...")
        else
            var = polynomial(Rop_local.varname);
            Reval=subs(Rop_local,var,chebgrid);
        end        
    else
        Reval=Rop_local*ones(1,deg+2);
    end

    % Chebyshev coefficients, as column k
    acheb=fcht(double(Reval));

for j=1:length(acheb)
    Alocal(j,k)=acheb(j);
end

end % k loop

% A keeps the N-p(m)+1 coefficients resolved for row state m, A_2PDEstate
% all N+1
rsize=N-p(m)+1;
A_cell_2PDEstate{m,1}=Alocal(1:N+1,1:no);
A_cell{m,1}=Alocal(1:rsize,1:no);

end % m loop


A_2PDEstate=cell2mat(A_cell_2PDEstate);
A=cell2mat(A_cell);
