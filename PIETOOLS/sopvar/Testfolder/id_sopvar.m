function I = id_sopvar(m,vars,dom)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I = id_sopvar(m,vars,dom) returns the identity operator on L2^m[S] as a
% 'sopvar' object, for any number N of spatial variables,
%
%   (I*x)(s) = x(s),   s in [a1,b1] x ... x [aN,bN].
%
% I = id_sopvar(A) and I = id_sopvar(m,A) take the variables and domain of the
% output space of a given 'sopvar' operator A, so that I*A = A.
%
% INPUT
% - m:    number of components, I: L2^m[S] -> L2^m[S]. Default 1;
% - vars: 1 x N cell array of variable names, e.g. {'s1','s2'}, or a
%         single name. Default {'s1'};
% - dom:  N x 2 array, dom(k,:) = [ak,bk] the domain of vars{k}.
%         Default [0,1] in every variable. Instead, the struct A.dom of a
%         'sopvar' A can be passed, in which case A.dom.out is used;
%         Instead, the struct A.vars with fields 'in' and 'out' of a 'sopvar'
%         A can be passed, in which case the output variables A.vars.out
%         are used;
% - A:    'sopvar' operator, in place of vars (and dom): the identity then
%         acts on the output space of A, vars = A.vars.out and
%         dom = A.dom.out, and m defaults to A.dims(1). For an operator
%         with the same input and output space this is also the identity
%         for A*I = A.
%
% OUTPUT
% - I:    m x m 'sopvar' object, with all variables common to input and
%         output (S1 = S2 = empty, S3 = vars).
%
% NOTES
% The identity is the multiplier operator delta(s-s') in every direction,
% which in the 3^N cell of parameters is the all-multiplier cell (index
% 1 in each direction, the first entry of P.params). That cell is
% eye(m), with the constant monomial as the only basis element in each
% variable (ZL = ZR = {0,...,0}); all other parameters are zero. This is
% the form recognised by is_id_op, and it satisfies the canonical
% multiplier form (no degree on the dummy variable).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - id_sopvar
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
% Initial coding AT - 10/15/2026
% AT, 10/15/2026: Optional sopvar input A for the variables and domain.
% AT, 10/15/2026: Accept the structs A.vars and A.dom, id_sopvar(m,A.vars,A.dom),
%                 and use their output fields.

if nargin<3,  dom = [];  end                       % so isempty(dom) means 'not given'% AT, 10/15/2026
if nargin>=1 && isa(m,'sopvar')                                             % AT, 10/15/2026
    I = id_sopvar([],m);  return       % id_sopvar(A): m defaults to A.dims(1)% AT, 10/15/2026
end                                                                         % AT, 10/15/2026
if nargin>=2 && isa(vars,'sopvar')                                          % AT, 10/15/2026
    A = vars;                                                               % AT, 10/15/2026
    if nargin<1 || isempty(m),  m = A.dims(1);  end    % components of A      % AT, 10/15/2026
    vars = A.vars.out;                                                      % AT, 10/15/2026
    if isempty(dom),  dom = A.dom.out;  end                                 % AT, 10/15/2026
end                                                                         % AT, 10/15/2026
if nargin<1 || isempty(m),     m = 1;         end
% A.vars, A.dom structs: use the output variables and domain.            % AT, 10/15/2026
if nargin>=2 && isstruct(vars),  vars = vars.out;  end                      % AT, 10/15/2026
if nargin>=3 && isstruct(dom),   dom = dom.out;    end                      % AT, 10/15/2026
if nargin<2 || isempty(vars),  vars = {'s1'}; end
if ischar(vars) || isstring(vars),  vars = cellstr(vars);  end
vars = vars(:).';
N = numel(vars);
%if nargin<3 || isempty(dom),   dom = repmat([0,1],N,1);  end              % AT, 10/15/2026 (was)
if isempty(dom),   dom = repmat([0,1],N,1);  end                            % AT, 10/15/2026
if ~isscalar(m) || m<1 || m~=round(m)
    error('m should be a positive integer.')
end
if ~isequal(size(dom),[N,2])
    error('dom should be an N x 2 array, one row per variable.')
end

v = struct();    v.in = vars;    v.out = vars;
d = struct();    d.in = dom;     d.out = dom;
Z = repmat({0},1,N);                    % constant monomial in each variable
prm = cell([repmat(3,1,N),1]);          % 3^N parameters
for i = 1:numel(prm),  prm{i} = sparse(m,m);  end
prm{1} = speye(m);                      % all-multiplier cell: delta(s-s') in every direction
I = sopvar(prm,v,Z,Z,d,[m,m]);
end
