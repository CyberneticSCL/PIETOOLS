function Y = apply_poly(P,xvals)
% Y = APPLY_POLY(P,XVALS) evaluates the matrix-valued distributed
% polynomial P at the state values XVALS, returning P(x) = C*Z(x).
%
% INPUTS
% - P:      'polyopvar' object;
% - xvals:  p x 1 cell, one entry per state variable (matching
%           P.varname, in order), each a 'polynomial'/'double' column
%           vector of size P.dims.in(i). States with an empty entry, or
%           beyond numel(xvals), are not substituted, and must not appear
%           in P;
%
% OUTPUTS
% - Y:      polynomial array representing P(x) evaluated at xvals;
%
% NOTES
% - The values are substituted for the states by 'subs', which replaces
%   the states by constants (see subs_const for its restrictions, e.g. a
%   polynomial value may only depend on the spatial variables of its
%   state). This gives a degree-0 polyopvar, of which the value is read
%   off by applying the operators to the dummy input 1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - apply_poly
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

if ~isa(P,'polyopvar')
    error("First argument must be of type 'polyopvar'.")
end
if ~P.chkval
    P = initialize(P);
end

% Replace the states by their values: stack the values of the states that are given
nstates = numel(P.varname);
xvals = xvals(:);
xvals(end+1:nstates) = {[]};
xvals = xvals(1:nstates);
has_val = ~cellfun(@isempty,xvals);
if any(any(P.depmat.out>1))
    % An output that is a tensor product of its factors, a function of several copies of a
    % variable (s_1,s_2,...), cannot be represented by the constant polyopvar that subs returns
    % (see subs_const): apply the operators to the values directly.
    Z = full(P.Z);
    nfac = max(num_factors(P),1);
    Xcell = cell(1,size(Z,1));
    for k=1:size(Z,1)
        Xk = num2cell(ones(1,nfac(k)));         % multipliers from R take the dummy input 1
        Xk(find(Z(k,:))) = reshape(xvals(Z(k,Z(k,:)~=0)),1,[]);
        Xcell{k} = Xk;
    end
    Y = apply_op(P.C,Xcell);
    return
end
if any(has_val)
    vals = xvals(has_val);
    if ~any(cellfun(@(a) isa(a,'dpvar'),vals))
        vals = cellfun(@polynomial,vals,'UniformOutput',false);
    end
    P = subs(P,P.varname(has_val),vertcat(vals{:}));
end
if any(any(P.Z))
    error("Values must be specified for all states appearing in the polynomial.")
end

% Read off the constant: each factor of each monomial is a multiplier from R, taking the dummy input 1
nfac = max(num_factors(P),1);
Xcell = cell(1,size(P.Z,1));
for k=1:numel(Xcell)
    Xcell{k} = num2cell(ones(1,nfac(k)));
end
Y = apply_op(P.C,Xcell);
% apply_op renames output variables s to s_1 for operators with a Kronecker (otype(1)=2) output;
% as depmat.out<=1 here, rename them back, as subs_const does, so that the result is in the same
% variables with or without the substitution
if isa(Y,'polynomial')
    for k=1:numel(P.vars.out)
        Y.varname(strcmp(Y.varname,[P.vars.out{k},'_1'])) = P.vars.out(k);
    end
end

end
