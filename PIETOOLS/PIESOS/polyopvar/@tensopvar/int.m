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
%           and all output variables must then have degree 1 (deps.out==1);
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
    if any(degs_out~=1)
        error("A weight depending on the output variables is only supported if each output variable has degree 1.")
    elseif ~all(ismember(W.varname,vars_out)) 
        error("The weight may only depend on the output variables of the operator.") 
    end
    dum_vars = Kop.varnames.in_dum;
    for i=1:M
        W.varname(strcmp(W.varname,vars_out{i})) = dum_vars{i}(1);
    end
    Kop = intop({W},intmat,vars_out,dom_out,true);
end

% Take the composition of the integral with Cop
KCop = Kop*Cop;
if ~isempty(W) && ~has_spatial_W
    KCop = W*KCop;
end

end