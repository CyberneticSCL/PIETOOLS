function [A_out,B_out] = common_basis(A_in,B_in)
% [A_OUT,B_OUT] = COMMON_BASIS(A_IN,B_IN) takes two distributed polynomials
% and expresses them in terms of the same basis of distributed monomials.
%
% INPUTS
% - A_in, B_in:     m x n 'polyopvar' objects;
%
% OUTPUTS
% - A_out, B_out:   m x n 'polyopvar' objects representing the same
%                   distributed polynomials as A_in, B_in, respectively,
%                   but expressed in terms of the same basis of
%                   monomials, so that
%                       A_out.varname = B_out.varname;
%                       A_out.Z = B_out.Z;
%                       A_out.vars = B_out.vars;
%                       A_out.depmat = B_out.depmat;
%
% NOTES
% - If a polynomial has several terms with the same monomial (equal rows
%   of Z), the operators acting on them are added, so that A_out, B_out
%   have unique monomials.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - common_basis
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
% DJ, 10/15/2026: Initial coding;

% Initialize the polynomials
if ~A_in.chkval
    A_in = initialize(A_in);
end
if ~B_in.chkval
    B_in = initialize(B_in);
end

% First, make sure the polynomials are expressed in terms of the same state
% variables and independent variables
[A,B] = common_vars(A_in,B_in);
Z_A = A.Z;      Z_B = B.Z;
C_A = A.C;      C_B = B.C;

% Combine the monomial bases into a unique basis
nZ1 = size(Z_A,1);
[nrA,ncA] = size(Z_A);
[nrB,ncB] = size(Z_B);
Z_A = [Z_A,zeros(nrA,max(ncB-ncA,0))];
Z_B = [Z_B,zeros(nrB,max(ncA-ncB,0))];
Z_full = [Z_A;Z_B];
% Monomials with the same Z but different numbers of factors are different monomials: e.g.
% x o c, with c a multiplier from R, Z=[1,0], and x, Z=[1], which are the same row once padded
% to the same width, see `polyopvar`
[~,nops_A,empty_A] = num_factors(A);        [~,nops_B,empty_B] = num_factors(B);
nops_full = factor_key(Z_full,[nops_A;nops_B],[empty_A;empty_B]);   % empty monomials merge with those of the same Z
[~,old2new_Cdcs,new2old_Cdcs] = unique([Z_full,nops_full],'rows');
Z_new = Z_full(old2new_Cdcs,:);
nZ = size(Z_new,1);

% Update the coefficient operators to match the new basis
dims_in_full = [C_A.dims.in; C_B.dims.in];
deps_in_full = [C_A.depmat.in; C_B.depmat.in];
new2old_Cdcs1 = new2old_Cdcs(1:nZ1);    new2old_Cdcs2 = new2old_Cdcs(nZ1+1:end);
C_A.ops = place_ops(C_A.ops,new2old_Cdcs1,nZ);
C_B.ops = place_ops(C_B.ops,new2old_Cdcs2,nZ);
C_A.dims.in = dims_in_full(old2new_Cdcs);
C_B.dims.in = dims_in_full(old2new_Cdcs);
C_A.depmat.in = deps_in_full(old2new_Cdcs,:);
C_B.depmat.in = deps_in_full(old2new_Cdcs,:);
C_A.chkval = true;                      C_B.chkval = true;

% Build polynomials in terms of the shared basis
A_out = A;              B_out = B;
A_out.C = C_A;          B_out.C = C_B;
A_out.Z = Z_new;        B_out.Z = Z_new;
A_out.chkval = true;    B_out.chkval = true;

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
function ops_new = place_ops(ops,idcs,nZ)
% Place column k of the cell array OPS in column IDCS(k) of an array with NZ columns. If the monomials
% of several columns are equal (duplicates in Z), the operators acting on them are added.
ops_new = cell(size(ops,1),nZ);
ops_new(:,idcs) = ops;
[u,~,ic] = unique(idcs);
for k=u(accumarray(ic(:),1)>1)'
    cols = find(idcs==k);
    for i=1:size(ops,1)
        op = ops{i,cols(1)};
        for m=2:numel(cols)
            if isempty(op)
                op = ops{i,cols(m)};
            elseif ~isempty(ops{i,cols(m)})
                op = op + ops{i,cols(m)};
            end
        end
        ops_new{i,k} = op;
    end
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %