function Q = subs(P,invar,inval,opts)
% Q = SUBS(P,INVAR,INVAL,OPTS) evaluates the output of the PI operator P (a
% 'copvar' object) at a boundary value of one of the variables of its output
% spaces,
%       (P*u)_i(invar=inval) = (Q*u)_i,
% for each block row i of the output space that depends on the variables
% invar. The block rows that do not depend on invar are left unchanged, so
% that Q is a container with the same inputs as P, and with output spaces
% L_2^{q_i}[s^i], for s^i the variables of the i-th block row of P without
% invar, if it contained invar.
% See the 'subs' method of 'sopvar' for the evaluation of each block, and
% for the operators that are supported:
%     * int_{a}^{s} K(s,theta) u(theta) dtheta
%       --> int_{a}^{b} K(b,theta) u(theta) dtheta (s=b) or 0 (s=a);
%     * int_{s}^{b} K(s,theta) u(theta) dtheta
%       --> int_{a}^{b} K(a,theta) u(theta) dtheta (s=a) or 0 (s=b);
%     * multiplier K(s) u(s): UNSUPPORTED unless K(inval)=0.
%
% INPUT
% - P:      'copvar' object;
% - invar:  char, string, cell of char, or 'polynomial' object (pvar),
%           specifying variables of the output spaces of P;
% - inval:  numeric array, with one boundary value for each variable in
%           invar;
% - opts:   (optional) tolerance for multipliers to vanish, see the 'subs'
%           method of 'sopvar';
%
% OUTPUT
% - Q:      'copvar' object with blocks Q.C{i,j}=subs(P.C{i,j},...) in the
%           block rows depending on invar, and Q.C{i,j}=P.C{i,j} in the
%           others. Zero blocks remain zero.
%
% NOTES
% - A variable must be a variable of the output space of at least one block
%   row of P. The block rows are evaluated separately: a variable in the
%   input space of a block column becomes an integrated variable of the
%   blocks in rows that depended on it, as it already is in the blocks of
%   the other rows, so the input spaces of the columns do not change.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - subs (copvar)
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
% DJ, 11/01/2026: Initial coding;

% % % Check the inputs
if ~isa(P,'copvar')
    error("Operator must be specified as object of type 'copvar'.")
end
if nargin<4
    opts = [];
end
if ischar(invar) || isstring(invar)
    invar = cellstr(invar);
elseif isa(invar,'polynomial')
    invar = invar.varname;
end
invar = invar(:);
if isa(inval,'polynomial')
    % A variable evaluated at itself (x=x) is left unchanged
    same = false(numel(invar),1);
    for i=1:numel(invar)
        same(i) = ispvar(inval(i)) && strcmp(inval(i).varname{1},invar{i});
    end
    invar = invar(~same);   inval = double(inval(~same));
end
inval = inval(:);
if numel(invar)~=numel(inval)
    error("The number of variables must match the number of values at which to evaluate them.")
end

% % % Evaluate along each variable in turn, in the block rows whose output
% % % space contains it
Q = P;
for ii=1:numel(invar)
    name = invar{ii};       val = inval(ii);
    ridcs = find(any(Q.space_out(:,strcmp(Q.vars,name)),2))';
    if isempty(ridcs)
        error("'%s' is not a variable of the output space of the operator.",name)
    end
    % Evaluate each nonzero block of those rows; the copvar constructor
    % then rebuilds the output spaces (without name in these rows)
    C_new = Q.C;
    for ridx=ridcs
        nz_blocks = ~cellfun(@isempty,C_new(ridx,:));
        C_new(ridx,nz_blocks) = cellfun(@(Cij) subs(Cij,name,val,opts),C_new(ridx,nz_blocks),'UniformOutput',false);
    end
    Q = copvar(C_new);
end

end
