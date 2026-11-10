function [A_out,B_out] = common_outputs(A_in,B_in)
% [A_OUT,B_OUT] = COMMON_OUTPUTS(A_IN,B_IN) takes two distributed
% polynomials and lifts them to the same output space, so that they can be
% added. The output variables are those of A_in and B_in together, and the
% degree of each row-block in each variable is the largest of the two. The
% output of a polynomial is constant along the variables it was not defined
% in, see 'lift_output'. A scalar-valued polynomial is also lifted to the
% dimension of a vector-valued one, repeating the scalar in each row.
%
% INPUTS
% - A_in, B_in:     'polyopvar' objects, of the same size, or one of them
%                   scalar-valued;
%
% OUTPUTS
% - A_out, B_out:   'polyopvar' objects representing the same polynomials
%                   as A_in, B_in (with a scalar repeated for each row, if
%                   the other is vector-valued), but with the same output
%                   space, so that
%                       A_out.dims.out = B_out.dims.out;
%                       A_out.vars.out = B_out.vars.out;
%                       A_out.depmat.out = B_out.depmat.out;
%
% NOTES
% - The polynomials are returned unchanged if their output spaces already
%   match, or if the dimensions of their row-blocks differ (other than for
%   a scalar and a vector), in which case adding the polynomials throws an
%   error.
% - Output variables that A_in and B_in share must have the same domain.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - common_outputs
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
% DJ, 11/03/2026: Initial coding;

A_out = A_in;       B_out = B_in;
if ~A_out.chkval
    A_out = initialize(A_out);
end
if ~B_out.chkval
    B_out = initialize(B_out);
end

% A scalar is added to each row of a vector: lift it to the dimensions of
% the row-blocks of the vector
dims_A = A_out.dims.out;        dims_B = B_out.dims.out;
if isequal(size(A_out),[1,1]) && size(B_out,1)>1 && size(B_out,2)==1
    dims_A = dims_B;
elseif isequal(size(B_out),[1,1]) && size(A_out,1)>1 && size(A_out,2)==1
    dims_B = dims_A;
end
if ~isequal(dims_A,dims_B)
    return
end
if isequal(dims_A,A_out.dims.out) && isequal(dims_B,B_out.dims.out) ...
        && isequal(A_out.vars.out,B_out.vars.out) && isequal(A_out.depmat.out,B_out.depmat.out)
    return
end

% Establish the output variables of A and B together, which must have the
% same domain in both
[vars,idcs] = unique([A_out.vars.out,B_out.vars.out]);
[~,ia,ib] = intersect(A_out.vars.out,B_out.vars.out);
if ~isempty(ia) && any(any(A_out.dom.out(ia,:)~=B_out.dom.out(ib,:)))
    error("Output variables of the polynomials to add must have the same domain.")
end
dom_full = [A_out.dom.out; B_out.dom.out];
dom = dom_full(idcs,:);

% The degree of each row-block in each variable is the largest of the two,
% where the single row of degrees of a scalar is expanded to each row-block
% of the vector
M = numel(vars);
[~,idcs_A] = ismember(A_out.vars.out,vars);
[~,idcs_B] = ismember(B_out.vars.out,vars);
deps_A = zeros(size(A_out.depmat.out,1),M);     deps_A(:,idcs_A) = A_out.depmat.out;
deps_B = zeros(size(B_out.depmat.out,1),M);     deps_B(:,idcs_B) = B_out.depmat.out;
deps = max(deps_A,deps_B);

A_out = lift_output(A_out,vars,dom,dims_A,deps);
B_out = lift_output(B_out,vars,dom,dims_B,deps);

end
