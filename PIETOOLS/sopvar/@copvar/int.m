function Pint = int(Pop,vars,dom)
% PINT = INT(POP,VARS,DOM) takes the composition of the PI operator POP
% with an integral operator along the variables VARS, over the domain DOM,
%   (Pint*u)_i(s) = int_L^U (Pop*u)_i(s,v) dv,
% where [L,U] = DOM(k,:) for v = VARS{k}, for each of the output spaces i
% of Pop. The integral operator acts as the identity along all other
% variables.
%
% INPUTS
% - Pop:    M x N 'copvar' object representing a PI operator;
% - vars:   r x 1 or 1 x r 'cellstr' or 'polynomial' (pvar) array
%           specifying the names of the unique spatial variables to
%           integrate over;
% - dom:    (optional) r x 2 'double', 'polynomial' or 'cell' array
%           specifying the lower and upper limits of integration for each
%           variable. Each limit is either the variable itself, or the
%           corresponding limit of the domain of the variable, which must
%           match that in Pop.dom, e.g. [a,b], [a,s] or [s,b] for s in
%           [a,b]. Defaults to the full domain, [a,b], of each variable;
%
% OUTPUTS
% - Pint:   M x N 'copvar' object representing the composition of the
%           integral operator with the PI operator Pop;
%
% NOTES
% - The operator is declared as an M x M 'copvar' object, with block (i,i)
%   the integral operator on the output space i of Pop, declared with
%   'int2sopvar', and composed with Pop using 'mtimes'. For each variable
%   s, the block acts on space i as
%       s in space i, [a,b]:    full integral, s leaves the output space;
%       s in space i, [a,s]:    int_a^s, s remains in the output space;
%       s in space i, [s,b]:    int_s^b, s remains in the output space;
%       s in space i, [s,s]:    multiplier (identity) along s, as if s
%                               were not in vars;
%       s not in space i:       multiplication with U-L, so the limits must
%                               be numbers for every space lacking s;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - int(copvar)
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
% DJ, 10/29/2026: Initial coding;

if ~isa(Pop,'copvar')
    error("PI operator must be specified as object of type 'copvar'.")
end
if nargin<=2
    dom = [];
end

% Variables and limits; the registry holds the domain of every variable
[vars,lims,is_var] = parse_int_inputs(vars,dom,Pop.vars,Pop.dom);
if isempty(vars)
    % Integration along no variables is the identity
    Pint = Pop;
    return
end

% A variable limit adds the variable to the output, so the variable must
% be in every output space of Pop
is_all = ismember(vars,Pop.vars(all(Pop.space_out,1)));
k = find(~is_all(:) & any(is_var,2),1);
if ~isempty(k)
    error("Variable limit of integration in '%s' not allowed, as '%s' is not in every output space of the operator.",vars{k},vars{k})
end

% No variable in any output space: each space integrates a constant, so the
% integral is U-L. Scaling Pop skips mtimes of the containers
if ~any(ismember(vars,Pop.vars(any(Pop.space_out,1))))
    Pint = prod(lims(:,2)-lims(:,1))*Pop;
    return
end

% Declare the integral operator on each output space of Pop, with kernel
% coef*I, where coef is U-L for the variables not in the space. The limits
% of integration are [s,s], the identity, for the variables of the space not
% in vars, and those of vars otherwise.
lim_cell = num2cell(lims);          % numbers, or the name where the limit is the variable
names = repmat(vars(:),1,2);        lim_cell(is_var) = names(is_var);
M = size(Pop.C,1);
C = cell(M,M);          % off-diagonal blocks are structurally zero
for i=1:M
    vout = Pop.vars(Pop.space_out(i,:));    dout = Pop.dom(Pop.space_out(i,:),:);
    is_out = ismember(vars,vout);
    [~,idcs] = ismember(vars(is_out),vout);
    int_dom = [vout(:),vout(:)];            int_dom(idcs,:) = lim_cell(is_out,:);
    is_S1 = false(1,numel(vout));           is_S1(idcs(~any(is_var(is_out,:),2))) = true;
    coef = prod(lims(~is_out,2)-lims(~is_out,1));
    C{i,i} = int2sopvar(coef*speye(Pop.dim_out(i)),struct('in',{vout},'out',{vout(~is_S1)}),...
                struct('in',dout,'out',dout(~is_S1,:)),int_dom);
end
Iop = copvar(C);
Pint = Iop*Pop;

end
