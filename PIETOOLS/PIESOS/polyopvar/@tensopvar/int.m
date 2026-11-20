function KCop = int(Cop,W,vars,dom)                                         % DJ, 11/19/2026
% KCOP = INT(COP,W,VARS,DOM) takes the composition of an integral
% (functional) with a tensor-PI operator COP,
%   KCop * x = int_{Omega} W*Cop*x dOmega
% where Omega is the output domain of COP, or the domain DOM of the
% variables VARS.
%
% INPUTS
% - Cop:    p x n 'tensopvar' object representing a tensor-PI operator
%           acting on a single degree-d distributed monomial;
% - W:      (optional) q x p 'double', 'polynomial' or 'dpvar' weight
%           inside the integral, defaults to the identity. If W depends on
%           the output variables of Cop, it is integrated along with Cop,
%           and may depend only on the variables
%             'si_j'  for the jth copy of vars.out{i}, j=1,...,deps.out(i),
%                     as named by `apply_op` for a tensor product of the
%                     outputs (Cop.otype(1)==2, or degree above 1);
%             'si'    also for the only copy, if the ith output variable
%                     vars.out{i} has degree 1 (deps.out(i)==1) and
%                     Cop.otype(1)==1 (a pointwise product);
%           where si is the name vars.out{i}. Other variables are rejected,
%           as is a weight depending on both 'si' and 'si_1', unless they
%           are in VARS;
% - vars:   (optional) 1 x N 'cellstr' or 'polynomial' (pvar) array, the
%           variables to integrate over: each copy of the output variables
%           of Cop, named as for W, once, and any other variables. For the
%           other variables, W is simply multiplied by the length of its 
%           domain. Defaults to the copies of the output variables;
% - dom:    (optional) N x 2 array, the domain of each variable in VARS:
%           that of Cop for its output variables (the full domain).
%           Variable limits are not supported. Defaults to the domain of
%           Cop, so it must be given if VARS has other variables;
%
% OUTPUTS
% - KCop:   q x n 'intop' object defining a functional representing the
%           composition of an integral operator and Cop;
%
% NOTES
% - The variables "vars" must include all variables on which the output of
% Cop depends, integration with respect to only a subset of the variables
% is not supported. 
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
% DJ, 10/15/2026: Initial coding;
% DJ, 11/10/2026: Allow a weight in the copies si_j of the outputs, for
%                   tensor products;
% DJ, 11/19/2026: Add optional inputs vars, dom, integrating W over the
%                   variables on which the output of Cop does not depend;

if nargin<2
    W = [];
end
% Check that the tensor-PI operator is properly declared
if ~Cop.chkval
    Cop = initialize(Cop);
end

% Variables of integration on which the output of Cop does not depend:      % DJ, 11/19/2026
% integrate W over them, so that KCop integrates over all of VARS. The
% copies of the output variables are all integrated over below, over the
% domain of Cop
if nargin>=3
    if nargin<4
        dom = [];
    end
    % Copies of the output variables, named as for W (see below), with the
    % index of the copy each name refers to and its domain
    names_cp = cell(1,0);       idcs_cp = zeros(1,0);       dom_cp = zeros(0,2);
    ncp = 0;
    for i=1:numel(Cop.vars.out)
        d_i = Cop.deps.out(i);
        names_i = get_tensorvars(Cop.vars.out{i},d_i,'full');
        idcs_i = ncp+(1:d_i);
        if d_i==1 && Cop.otype(1)==1
            names_i = [names_i,get_tensorvars(Cop.vars.out{i},1,{'bare','full'})];
            idcs_i = [idcs_i,ncp+1];
        end
        names_cp = [names_cp,names_i];
        idcs_cp = [idcs_cp,idcs_i];
        dom_cp = [dom_cp;repmat(Cop.dom.out(i,:),numel(idcs_i),1)];
        ncp = ncp+d_i;
    end
    [vars,lims,is_var] = parse_int_inputs(vars,dom,names_cp,dom_cp);
    if any(is_var(:))
        error("Variable limits of integration are currently not supported.")
    end
    % KCop is a functional, no partial integral: each copy of the output
    % must be integrated over, once ('si' and 'si_1' are the same copy),
    % over its full domain
    [is_out,loc] = ismember(vars,names_cp);
    if ~isequal(sort(idcs_cp(loc(is_out))),1:ncp)
        error("Partial integrals are not supported: variables of integration must include each copy of the output variables of the operator once.")
    elseif any(any(lims(is_out,:)~=dom_cp(loc(is_out),:)))
        error("Output variables of the operator must be integrated over their full domain.")
    end
    for k=find(~is_out)
        if isempty(W)
            W = (lims(k,2)-lims(k,1))*speye(size(Cop,1));
        elseif isa(W,'double')
            W = (lims(k,2)-lims(k,1))*W;
        else
            W = int(W,polynomial(vars(k)),lims(k,1),lims(k,2));
        end
    end
end

% Extract the output function space of Cop
m = size(Cop,1);
M = numel(Cop.vars.out);
vars_out = Cop.vars.out;
dom_out = Cop.dom.out;
degs_out = Cop.deps.out;

% Declare an integral on the output function space
intmat = cell(1,M);
for i=1:M
    intmat{i} = repmat([0,degs_out(i)+1],1,degs_out(i));
end
params = {speye(m)};
Kop = intop(params,intmat,vars_out,dom_out,true);

% For weights depending on variables: replace the identity by W, with its
% variables renamed to the dummy variables of integration
has_spatial_W = ~isempty(W) && ~isa(W,'double') && ~isempty(W.varname);
if has_spatial_W
    % Variables allowed in W, and the dummy variables of Kop they are
    % renamed to: the copies <name>_<j>, as named by apply_op for a tensor
    % product, and, for the only copy of a pointwise output, also the bare
    % name
    dum_vars = Kop.varnames.in_dum;
    names_W = cell(1,0);        names_dum = cell(1,0);
    for i=1:M
        names_i = get_tensorvars(vars_out{i},degs_out(i),'full');
        dum_i = reshape(dum_vars{i},1,[]);
        if degs_out(i)==1 && Cop.otype(1)==1
            names_i = [names_i,get_tensorvars(vars_out{i},1,{'bare','full'})];
            dum_i = [dum_i,dum_vars{i}(1)];
        end
        names_W = [names_W,names_i];
        names_dum = [names_dum,dum_i];
    end
    [is_ok,idcs] = ismember(W.varname,names_W);
    if ~all(is_ok)
        error("The weight may only depend on the variables 'si_j' (copy j of the output variable si) or, for one copy of a pointwise output, 'si'.")
    end
    % Rename all variables at once: a copy may not be renamed twice
    new_names = names_dum(idcs);
    if numel(unique(new_names))<numel(new_names)
        error("The weight may not depend on both 'si' and 'si_1', the same variable.")
    end
    W.varname = reshape(new_names,size(W.varname));
    Kop = intop({W},intmat,vars_out,dom_out,true);
end

% Take the composition of the integral with Cop
KCop = Kop*Cop;
if ~isempty(W) && ~has_spatial_W
    KCop = W*KCop;
end

end