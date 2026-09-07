function Pd = sopvar2sdopvar(P,Zd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pd = sopvar2sdopvar(P,Zd) takes a fixed sopvar operator and returns the
% sdopvar object representing the same operator as a decision operator that
% happens not to depend on any decision variable.
%
% INPUT
% P:    'sopvar' object;
% Zd:   (optional) cell array naming the decision variables of the result,
%       in the order they are to be indexed. Defaults to none;
%
% OUTPUT
% Pd:   'sdopvar' object representing the same operator as P, with
%           vec(C_gamma(d)) = A_gamma + B_gamma'*d,   B_gamma = 0.
%
% NOTES:
% An sopvar stores each parameter as a coefficient MATRIX, an sdopvar stores
% the same content as vec of that matrix plus a q x nC block of decision
% variable coefficients (sopvar document Sec. 8). Promotion is therefore a
% vec and a zero block; no arithmetic is involved and the kernels are
% unchanged.
%
% Pass Zd when the result is to be combined with an existing sdopvar, so the
% two share a decision variable list and 'CombineDecisionBasis' takes its
% fast path instead of a setdiff over every name.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026 PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, 09/07/2026

if ~isa(P,'sopvar')
    error("Input must be of type 'sopvar'.")
end
if nargin<2 || isempty(Zd)
    Zd = cell(0,1);
end

q = numel(Zd);
params = struct('A',{cell(size(P.params))},'B',{cell(size(P.params))});
for k = 1:numel(P.params)
    params.A{k} = sparse(P.params{k}(:));
    % Fixed content, so no dependence on the decision variables. Allocated
    % with 'sparse' so the q axis costs only a column pointer array.
    params.B{k} = sparse(q,numel(params.A{k}));
end

Pd = sdopvar(params,P.vars,Zd,P.ZL,P.ZR,P.dom,P.dims);

end
