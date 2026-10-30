function Pop = int2sopvar(Kmat,vars,dom,int_dom)
% POP = INT2SOPVAR(KMAT,VARS,DOM,INT_DOM) takes a kernel KMAT and returns a
% 'sopvar' object POP representing the associated integral operator between
% the function spaces defined by DOM,
%   (Pop*x)(s) = int_{INT_DOM} Kmat(s,t) x(t) dt,
% with s the variables vars.out, and t the variables vars.in.
%
% INPUTS
% - Kmat:   m x n 'double', 'polynomial' or 'quadPoly' array, representing
%           the kernel K(s,t). A 'polynomial' cannot distinguish the
%           variable s from the dummy variable t of the same name, so a
%           variable in vars.out is read as s, and a variable only in
%           vars.in as t; use 'quadPoly' for a kernel depending on both.
%           For a 'quadPoly', the left variables Kmat.ns must be in
%           vars.out and the right variables Kmat.nt in vars.in, except
%           for variables with only degree 0 in the kernel;
% - vars:   struct with fields
%     in:   1 x N cellstr array specifying the names of the input variables
%           to the operator, the variables t;
%    out:   1 x M cellstr array specifying the names of the output
%           variables to the operator, the variables s;
%           Can also be specified as a 1 x N cellstr, in which case
%           vars.in = vars.out = vars, as for 'copvar';
% - dom:    struct with fields
%    in:    N x 2 'double' array specifying the intervals on which the
%           input variables are defined;
%   out:    M x 2 'double' array specifying the intervals on which the
%           output variables are defined;
%           Can also be specified as an N x 2 'double' array, in which
%           case dom.in = dom.out = dom, as for 'copvar';
% - int_dom:    (optional) N x 2 'double', 'polynomial' or 'cell' array,
%           int_dom(i,:) = [L,U] specifying the lower and upper limit of
%           integration in the variable vars.in{i}. Each limit is either
%           the limit of dom.in(i,:), or the variable vars.in{i} itself (a
%           pvar, or its name in a 'cell'), the latter only if the variable
%           is also in vars.out. Defaults to the full domain, dom.in, for
%           all of the input variables;
%
% OUTPUTS
% - Pop:    m x n 'sopvar' object representing the integral operator
%           defined by Kmat;
%
% NOTES
% For a variable s in vars.in, with domain [a,b], and limits [L,U]:
%   [a,b]:  int_a^b K(s,t) x(t) dt. If s is also in vars.out, the output
%           still depends on s, and this is the sum of the integrals
%           int_a^s and int_s^b, parameters gamma = 2 and 3;
%   [a,s]:  int_a^s K(s,t) x(t) dt, parameter gamma = 2;
%   [s,b]:  int_s^b K(s,t) x(t) dt, parameter gamma = 3;
%   [s,s]:  multiplier along s, K(s) x(s), parameter gamma = 1. The kernel
%           may not depend on the dummy variable t of s, as a multiplier
%           acts pointwise;
% where the labels 1, 2, 3 for gamma are those of the sopvar document,
% Sec. 4 (see 'gamma_of_cell'). A variable only in vars.in must have the
% limits [a,b], as the integral is then taken over all of the interval.
% Other limits are refused, as the integral between other limits is not a
% PI operator in the variable.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - int2sopvar
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

% Check that the kernel is properly specified
if ~isa(Kmat,'double') && ~isa(Kmat,'polynomial') && ~isa(Kmat,'quadPoly')
    error("Kernel must be specified as object of type 'double', 'polynomial' or 'quadPoly'.")
end

% Check that the variables are properly specified. A cellstr is the same
% list of variables for the input and the output space, as for 'copvar'
if iscellstr(vars)
    vars = struct('in',{vars},'out',{vars});
end
if ~isa(vars,'struct') || ~isfield(vars,'in') || ~isfield(vars,'out')
    error("Variables must be specified as 'cellstr', or as struct with fields 'in' and 'out'.")
end
if ~iscellstr(vars.in) || ~iscellstr(vars.out)
    error("Input and output variable names should be specified as 'cellstr' objects.")
end
vars.in = reshape(vars.in,1,[]);        N = numel(vars.in);
vars.out = reshape(vars.out,1,[]);      M = numel(vars.out);

% Check that the domains are properly specified. An array is the same domain
% for the input and the output variable, as for 'copvar'
if isnumeric(dom)
    dom = struct('in',dom,'out',dom);
end
if ~isa(dom,'struct') || ~isfield(dom,'in') || ~isfield(dom,'out')
    error("Domains must be specified as numeric array, or as struct with fields 'in' and 'out'.")
end
if ~isnumeric(dom.in) || ~all(size(dom.in)==[N,2])
    error("Number of input domains should match number of input variables.")
end
if ~isnumeric(dom.out) || ~all(size(dom.out)==[M,2])
    error("Number of output domains should match number of output variables.")
end

% Check that the domain of integration is properly specified, defaulting to
% the full domain, and that each limit is one that has a PI form. A variable
% only in vars.in has no variable limit, so it is integrated over its full
% interval
if nargin<=3
    int_dom = [];
end
[~,lims,is_var] = parse_int_inputs(vars.in,int_dom,vars.in,dom.in);
is_out = ismember(vars.in,vars.out);        % input variables that are also output variables
is_var_bad = any(is_var,2) & ~is_out(:);
if any(is_var_bad)
    idx = find(is_var_bad,1);
    error("Variable limit of integration in '%s' not allowed, as '%s' is not an output variable of the operator.",vars.in{idx},vars.in{idx})
end
is_lim_bad = any(~(is_var | lims==dom.in),2);
if any(is_lim_bad)
    idx = find(is_lim_bad,1);
    error("Limits of integration in '%s' must each be '%s' or the limit of the domain [%g,%g].",vars.in{idx},vars.in{idx},dom.in(idx,1),dom.in(idx,2))
end

% Express the kernel as quadPoly, K(s,t) = (I_m x Z(s)') C (I_n x Z(t))
if isa(Kmat,'double')
    Kmat = quadPoly(Kmat);
elseif isa(Kmat,'polynomial')
    % The kernel can only depend on variables of the operator, and
    % polynomial2quadPoly would drop any other variable without an error
    vname_Kmat = Kmat.varname(:);
    has_coef = full(any(Kmat.coefficient,2));
    is_bad = any(full(Kmat.degmat(has_coef,:)),1)' & ~ismember(vname_Kmat,[vars.out,vars.in]);
    if any(is_bad)
        error("Kernel depends on variable(s) '%s', which are not variables of the operator.",...
            strjoin(vname_Kmat(is_bad)',"', '"))
    end
    % Variables in vars.out are the left variables s, the others are t
    Kmat = quadPoly.polynomial2quadPoly(Kmat,vars.out,vars.in(~is_out));
end

% Bind the monomial bases to the variables. The monomial index of the
% coefficients follows the variables of the kernel, so these are placed
% first, the others getting degree 0; the constructor puts the variables in
% canonical order, permuting the coefficients accordingly (sopvar.m)
[vars_out,ZL] = bind_basis(Kmat.ns,Kmat.Zs,vars.out,'output');
[vars_in,ZR] = bind_basis(Kmat.nt,Kmat.Zt,vars.in,'input');
[~,idcs_out] = ismember(vars_out,vars.out);
[~,idcs_in] = ismember(vars_in,vars.in);

% A multiplier acts pointwise, K(s) x(s), so the kernel cannot depend on
% the dummy variable t of a variable with the limits [s,s]
[~,idcs_mult] = ismember(vars.in(is_var(:,1) & is_var(:,2)),vars_in);
idx = find(cellfun(@(z) any(z~=0),ZR(idcs_mult)),1);
if ~isempty(idx)
    error("Kernel depends on the dummy variable of '%s', but the limits [s,s] define a multiplier along it, which acts pointwise.",...
        vars_in{idcs_mult(idx)})
end

% Declare the parameters, with the kernel in each integral of the sum that
% the limits define. The cells run over the shared variables, in sorted
% order, with gamma = 1 (multiplier), 2 (int_a^s) or 3 (int_s^b)
vars_S3 = intersect(vars.in,vars.out);
[~,idcs_S3] = ismember(vars_S3,vars.in);
gamma_all = zeros(1,0);
for k=1:numel(vars_S3)
    idx = idcs_S3(k);
    if is_var(idx,1) && is_var(idx,2)
        gamma_k = 1;
    elseif is_var(idx,1)
        gamma_k = 3;
    elseif is_var(idx,2)
        gamma_k = 2;
    else
        gamma_k = [2;3];
    end
    % All combinations with the options of the earlier variables
    gamma_all = [repelem(gamma_all,numel(gamma_k),1),repmat(gamma_k,size(gamma_all,1),1)];
end
params = repmat({sparse(size(Kmat.C,1),size(Kmat.C,2))},[3*ones(1,numel(vars_S3)),1,1]);
params(cell_of_gamma(gamma_all)) = {Kmat.C};

% Declare the operator
Pop = sopvar(params,struct('in',{vars_in},'out',{vars_out}),ZL,ZR,...
        struct('in',dom.in(idcs_in,:),'out',dom.out(idcs_out,:)),Kmat.dim);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [vars_op,Z_op] = bind_basis(vars_K,Z_K,vars_side,side)
% [VARS_OP,Z_OP] = BIND_BASIS(VARS_K,Z_K,VARS_SIDE,SIDE) gives each
% variable of one side of the operator a monomial basis, such that the
% coefficient matrix of the kernel can be used as it is. The variables of
% the kernel come first, with the bases of the kernel, which define the
% monomial index of its coefficients. The other variables of VARS_SIDE
% follow with the basis 0, which does not change that index. A variable of
% the kernel that is not in VARS_SIDE is dropped if the kernel has degree 0
% in it, and an error otherwise.
%
% INPUTS
% - vars_K:     1 x k 'cellstr', the variables of the kernel on this side
%               (ns or nt of the 'quadPoly');
% - Z_K:        1 x k cell of column vectors, the monomial degrees of the
%               kernel in each of these variables (Zs or Zt);
% - vars_side:  1 x p 'cellstr', the variables of the operator on this
%               side;
% - side:       'char', 'output' or 'input', for the error message;
%
% OUTPUTS
% - vars_op:    1 x p 'cellstr', the variables of vars_K in vars_side,
%               followed by the other variables of vars_side;
% - Z_op:       1 x p cell of column vectors, the monomial degrees of each
%               of the variables in vars_op;
%
vars_K = reshape(vars_K,1,[]);      Z_K = reshape(Z_K,1,[]);
is_side = ismember(vars_K,vars_side);
% A variable that the operator lacks is only allowed with degree 0, in which
% case it does not affect the monomial index
is_const = cellfun(@(z) all(z(:)==0),Z_K);
if any(~is_side & ~is_const)
    i = find(~is_side & ~is_const,1);
    error("Kernel depends on variable '%s', which is not an %s variable of the operator.",vars_K{i},side)
end
vars_K = vars_K(is_side);       Z_K = Z_K(is_side);
vars_rest = vars_side(~ismember(vars_side,vars_K));
vars_op = [vars_K,vars_rest];
Z_op = [cellfun(@(z) z(:),Z_K,'UniformOutput',false),repmat({0},1,numel(vars_rest))];

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
