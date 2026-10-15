function [A_out,B_out] = common_states(A_in,B_in)
% [A_OUT,B_OUT] = COMMON_STATES(A_IN,B_IN) takes two distributed
% polynomials and expresses them in terms of the same STATE variables
% (the monomial factors x1,...,xp), without touching either object's own
% output space (C, vars.out, dom.out, depmat.out, dims.out).
%
% INPUTS
% - A_in, B_in:     m x n 'polyopvar' objects;
%
% OUTPUTS
% - A_out, B_out:   m x n 'polyopvar' objects representing the same
%                   distributed polynomials as A_in, B_in, respectively,
%                   but expressed in terms of the same state variables,
%                   so that
%                   A_out.varname = B_out.varname;
%                   A_out.vars.in = B_out.vars.in;
%                   A_out.C, .vars.out, .dom.out, .depmat.out, .dims.out
%                   are the same as for A_in, B_in;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - common_states
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

% Extract the parameters representing each polynomial's states
var_A = A_in.varname;               var_B = B_in.varname;
dims_A = A_in.dims;                 dims_B = B_in.dims;
pvars_A = A_in.vars;                pvars_B = B_in.vars;
dom_A = A_in.dom;                   dom_B = B_in.dom;
Z_A = A_in.Z;                       Z_B = B_in.Z;
depmat_A = A_in.depmat;             depmat_B = B_in.depmat;
N1 = numel(pvars_A.in);             N2 = numel(pvars_B.in);

% First, combine the list of state variables
p1 = numel(var_A);      p2 = numel(var_B);
if isequal(var_A,var_B)
    var_new = var_A(:);
    new2old_vdcs = [1:p1,1:p1]';
    old2new_vdcs = 1:p1;
else
    var_full = [var_A(:);var_B(:)];
    [var_new,old2new_vdcs,new2old_vdcs] = unique(var_full);   % var_full = var_new(new2old_vdcs);
end
p = numel(var_new);
% Keep track of the size of each variable
dims_in_full = [dims_A.in; dims_B.in];
dims_in = dims_in_full(old2new_vdcs);
if any(dims_in(new2old_vdcs)~=dims_in_full)
    error("Dimensions of state vectors in polynomials must match.")
end
dims_A.in = dims_in;                dims_B.in = dims_in;

% Combine the list of spatial variables underlying the states (input
% variables must be merged manually)
vars_full = [pvars_A.in, pvars_B.in];
[vars2,old2new_pdcs,new2old_pdcs] = unique(vars_full);
dom_full = [dom_A.in; dom_B.in];
dom2 = dom_full(old2new_pdcs,:);
if any(any(dom2(new2old_pdcs,:)~=dom_full))
    error("Identical spatial variables in the polynomials must have the same domain.")
end
pvars_A.in = vars2;                 pvars_B.in = vars2;
dom_A.in = dom2;                    dom_B.in = dom2;
N = numel(vars2);

% Build a depmat indicating for each of the new state variables on which
% spatial variable it depends.
depmat2 = false(p,N);
new2old_vdcs1 = new2old_vdcs(1:p1);
new2old_vdcs2 = new2old_vdcs(p1+(1:p2));
depmat2(new2old_vdcs1,new2old_pdcs(1:N1)) = depmat_A.in;
depmat2(new2old_vdcs2,new2old_pdcs(N1+(1:N2))) = depmat_B.in;
depmat_A.in = depmat2;              depmat_B.in = depmat2;

% Express the monomial bases in terms of the new set of state variables
[ridcsA,cidcsA,ZvalsA] = find(Z_A);
ZvalsA_new = new2old_vdcs1(ZvalsA);
Z_A = sparse(ridcsA,cidcsA,ZvalsA_new,size(Z_A,1),size(Z_A,2));
[ridcsB,cidcsB,ZvalsB] = find(Z_B);
ZvalsB_new = new2old_vdcs2(ZvalsB);
Z_B = sparse(ridcsB,cidcsB,ZvalsB_new,size(Z_B,1),size(Z_B,2));

% Build polynomials in terms of the shared state variables; output space
% (C, vars.out, dom.out, depmat.out, dims.out) is carried through as-is
A_out = polyopvar();            B_out = polyopvar();
A_out.C = A_in.C;                B_out.C = B_in.C;
A_out.Z = Z_A;                  B_out.Z = Z_B;
A_out.varname = var_new;        B_out.varname = var_new;
A_out.dims = dims_A;            B_out.dims = dims_B;
A_out.vars = pvars_A;           B_out.vars = pvars_B;
A_out.dom = dom_A;              B_out.dom = dom_B;
A_out.depmat = depmat_A;        B_out.depmat = depmat_B;

end
