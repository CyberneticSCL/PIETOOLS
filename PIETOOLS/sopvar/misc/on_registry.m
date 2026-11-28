function X = on_registry(X,R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X = ON_REGISTRY(X,R) restates the container X on the variable
% registry (vars, dom) of the container R, which must contain the registry
% of X. The blocks are untouched; only the space masks are widened.
% Operators built on different registries, e.g. -gam*I_w on R -> R and
% Tw'*(P*Bw) on the variable of the state, can then be added on one
% registry.
%
% INPUTS
% - X:  'opvar' or 'dopvar' object;
% - R:  'opvar' or 'dopvar' object whose registry contains that of X;
%
% OUTPUTS
% - X:  the same operator on the registry of R, of the class of the input.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - on_registry
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
% Initial coding MMP, 10/08/2026

% Already on the registry of R: nothing to do
if isequal(X.vars,R.vars) && isequal(X.dom,R.dom),  return,     end
% Positions of the variables of X in the registry of R, which must hold
% them all
meta = metadata(X);
[tf,loc] = ismember(X.vars,R.vars);
if ~all(tf)
    error('on_registry:notSubset','X''s registry is not contained in R''s.')
end
% The space masks widened to the registry of R, and the container rebuilt
% in its class on the new metadata
so = false(size(X.C,1),numel(R.vars));  so(:,loc) = X.space_out;
si = false(size(X.C,2),numel(R.vars));  si(:,loc) = X.space_in;
meta.vars = R.vars;     meta.dom = R.dom;
meta.space_out = so;    meta.space_in = si;
if isa(X,'dopvar'),    X = dopvar(X.C,meta);
else,                   X = opvar(X.C,meta);
end
end
