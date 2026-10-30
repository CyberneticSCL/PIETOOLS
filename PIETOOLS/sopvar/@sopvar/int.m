function Pint = int(Pop,vars,dom)
% PINT = INT(POP,VARS,DOM) takes the composition of the PI operator POP
% with an integral operator along the variables VARS, over the domain DOM,
%   (Pint*u)(s) = int_L^U (Pop*u)(s,v) dv,   [L,U] = DOM(i,:), v = VARS{i}
% The integral operator acts as the identity along all other variables.
%
% INPUTS
% - Pop:    m x n 'sopvar' object representing a PI operator;
% - vars:   N x 1 or 1 x N 'cellstr' or 'polynomial' (pvar) array
%           specifying the names of the unique spatial variables to
%           integrate over;
% - dom:    (optional) N x 2 'double', 'polynomial' or 'cell' array
%           specifying the lower and upper limits of integration for each
%           variable. Each limit is either the variable itself, or the
%           corresponding limit of the domain of the variable, which must
%           match that in Pop.dom, e.g. [a,b], [a,s] or [s,b] for s in
%           [a,b]. Defaults to the full domain, [a,b] (also if s is only an
%           input variable of Pop, in which case the domain is Pop.dom.in);
%
% OUTPUTS
% - Pint:   m x n 'sopvar' object representing the composition of the
%           integral operator with the PI operator Pop;
%
% NOTES
% - The operator is declared as a 'sopvar' object with 'int2sopvar', and
%   composed with Pop using 'mtimes'. For each variable s, it acts as
%       s in Pop.vars.out, [a,b]:    full integral, s is not an output
%                                    variable of Pint;
%       s in Pop.vars.out, [a,s]:    int_a^s, s remains an output variable;
%       s in Pop.vars.out, [s,b]:    int_s^b, s remains an output variable;
%       s in Pop.vars.out, [s,s]:    multiplier (identity) along s, as if s
%                                    were not in vars;
%       s not in Pop.vars.out:       Pop*u does not depend on s, so the
%                                    integral is multiplication with U-L.
%                                    Limits must then be numbers;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - int
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

if ~isa(Pop,'sopvar')
    error("PI operator must be specified as object of type 'sopvar'.")
end
if nargin<=2
    dom = [];
end

% Variables and limits; a variable of Pop has a known domain, in or out
[vars,lims,is_var] = parse_int_inputs(vars,dom,[Pop.vars.out(:)',Pop.vars.in(:)'],[Pop.dom.out;Pop.dom.in]);
if isempty(vars)
    % Integration along no variables is the identity
    Pint = Pop;
    return
end

% A variable limit adds the variable to the output, so the variable must
% be an output variable of Pop
is_out = ismember(vars,Pop.vars.out);
k = find(~is_out(:) & any(is_var,2),1);
if ~isempty(k)
    error("Variable limit of integration in '%s' not allowed, as '%s' is not an output variable of the operator.",vars{k},vars{k})
end

% Pop*u does not depend on the other variables, so the integral is U-L. If
% there are no other variables, scaling Pop skips mtimes
coef = prod(lims(~is_out,2)-lims(~is_out,1));
if ~any(is_out)
    Pint = coef*Pop;
    return
end

% Declare the integral operator on the output space of Pop, with kernel
% coef*I. The limits of integration are [s,s], the identity, for the output
% variables not in vars, and those of vars otherwise.
vout = Pop.vars.out;                dout = Pop.dom.out;
[~,idcs] = ismember(vars(is_out),vout);
lim_cell = num2cell(lims);          % numbers, or the name where the limit is the variable
names = repmat(vars(:),1,2);        lim_cell(is_var) = names(is_var);
int_dom = [vout(:),vout(:)];        int_dom(idcs,:) = lim_cell(is_out,:);
is_S1 = false(1,numel(vout));       is_S1(idcs(~any(is_var(is_out,:),2))) = true;
Iop = int2sopvar(coef*speye(Pop.dims(1)),struct('in',{vout},'out',{vout(~is_S1)}),...
        struct('in',dout,'out',dout(~is_S1,:)),int_dom);
Pint = Iop*Pop;

end
