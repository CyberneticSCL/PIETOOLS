function KCop = int(Cop,W)
% KCOP = INT(COP,W) takes the composition of an integral (functional)
% with a tensor-PI operator COP,
%   KCop * x = int_{Omega} W*Cop*x dOmega
% where Omega is the output domain of COP.
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
%           as is a weight depending on both 'si' and 'si_1';
%
% OUTPUTS
% - KCop:   q x n 'intop' object defining a functional representing the
%           composition of an integral operator and Cop;
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

if nargin<2
    W = [];
end  
% Check that the tensor-PI operator is properly declared
if ~Cop.chkval
    Cop = initialize(Cop);
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
        names_i = arrayfun(@(j) [vars_out{i},'_',num2str(j)],1:degs_out(i),'UniformOutput',false);
        dum_i = reshape(dum_vars{i},1,[]);
        if degs_out(i)==1 && Cop.otype(1)==1
            names_i = [names_i,vars_out(i)];
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