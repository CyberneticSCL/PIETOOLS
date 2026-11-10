function [prog,Pmat,Zop] = piesos_poslpivar(prog,dims,vars,dom,depmat,degs,opts)
% [PROG,PMAT,ZOP] = PIESOS_POSLPIVAR(PROG,DIMS,VARS,DOM,DEPMAT,DEGS,OPTS)
% declares a positive semidefinite matrix variable PMAT(s), SOS in the
% primary variables s, parameterizing the positive tensor-PI operator
%   Pop = Zop'*Pmat*Zop
% on the function space
%   L2^{n1}[Omega1] x ... x L2^{nq}[Omegaq]
% defined by DIMS, VARS, DOM, and DEPMAT. (The operator Pop itself is not
% returned, since it cannot be represented by a 'tensopmat'.)
%
% INPUTS
% - prog:   'struct' specifying a PIESOS program in which to declare the
%           operator variable (see `piesos_program`);
% - dims:   q x 1 array specifying the number of coupled state variables
%           (ni) on which each of the blocks in Zop acts. Can also be
%           specified as struct with field 'in' containing this array;
% - vars:   1 x N cellstr specifying the names of the spatial variables
%           in terms of which the states are defined. Can also be
%           specified as struct with fields 'in' and 'out', which must
%           be equal;
% - dom:    N x 2 array specifying the intervals on which the spatial
%           variables are defined. Can also be specified as struct with
%           fields 'in' and 'out', which must be equal;
% - depmat: q x N integer array specifying the degree to which each of the
%           states depends on each of the spatial variables, i.e. the
%           distributed monomial degree of the state. Can also be
%           specified as struct with fields 'in' and 'out', which must
%           be equal;
% - degs:   (optional) struct with fields
%       'in':   1 x q cell, specifying the degrees of the input variables
%               in each operator block. Each element degs.in{i} can be a
%               di x N cell, where di depends on depmat(i,:), specifying
%               the degrees of the monomial basis in each of the factors
%               of each of the input variables in block i. Element
%               degs.in{i} can also be a di x N array, specifying the
%               maximal monomial degree in each of the factors of each
%               variable. If degs.in{i} is scalar, the same degree is
%               used for all variables;
%       'out':  1 x q cell, specifying the degrees of the primary
%               variables in each operator block. Each element
%               degs.out{i} can be a 1 x N cell of column vectors,
%               specifying the exponents of the monomial basis in each of
%               the primary variables, or a 1 x N array specifying the
%               maximal monomial degree in each of the variables (a
%               scalar is used for all variables).
%       The degrees in degs.in define the monomials in the dummy
%       variables, in terms of which the kernels of Zop are defined. The
%       degrees in degs.out define the monomials in the primary
%       variables, defining the basis of the SOS matrix Pmat. If degs is
%       an array or cell, the same degrees are used for both. Defaults
%       to 0;
% - opts:   (optional) struct with fields
%       'psatz':  1 x N array of integers, specifying the maximal order
%               of the psatz multiplier in each of the spatial
%               variables. Here, psatz(i)=qi implies inclusion of all
%               multipliers
%                   gj = (bi-si)^{qi-j}*(si-ai)^{j},   j in {0,...,qi}
%               This means a total number of prod(psatz+1) terms. Can
%               also be set to a q1 x ... x qN logical array, where now
%               the dimensions of psatz determine the maximal order, and
%               only the multipliers associated with true elements of
%               psatz are included. Element k along dimension i of psatz
%               corresponds to gj with j=k-1, so psatz(i)=qi is
%               equivalent to a logical array of size qi+1 along
%               dimension i. Terms are ordered by linear index in
%               psatz. A scalar psatz is used for all variables on which
%               some state depends. For gj of order qi, the degree in si
%               of the basis monomials of Pmat is reduced by ceil(qi/2),
%               so that the degree of Pmat in si does not exceed twice
%               its degree in degs.out;
%       'exclude':  1 x q cell, specifying for each of the operator
%               blocks whether any terms need to be excluded. Element
%               opts.exclude{i} can be a 3^N3(i) logical array, where
%               N3(i) is the number of 3-PI variables (appearing both in
%               the input and output) of the operator block i, with
%               element k along direction j set to true to exclude the
%               multiplier term (k=1), lower integral term (k=2), or
%               upper integral term (k=3), along spatial variable j in
%               the list of 3-PI variables (intersect(vars.in,vars.out)).
%               If exclude is a 3 x 1 array, it is applied in every
%               variable on which a state depends, and a term is excluded
%               if its type is excluded in any of these variables (e.g.
%               exclude = [0;0;1] excludes every term with an upper
%               integral). Variables a state does not depend on are
%               ignored, see `monopmat`;
%       'sep':  1 x q cell, specifying for each of the operator blocks
%               whether to enforce separability along each of the
%               variables. Element opts.sep{i} is a 1 x N3(i) logical
%               array, specifying for each of the 3-PI variables whether
%               the operator is separable along that variable, meaning
%               that the kernel in the upper and lower integrals is the
%               same;
%       'otype':  integer in {1,2} specifying whether to take the proper
%               tensor product (otype = 2) or the pointwise product
%               (otype = 1) of the 3-PI operators in constructing the
%               tensor-PI operator basis. If otype==1, then d_out =
%               depmat.out. If otype==2, then d_out = d*depmat.out,
%               where d = max(depmat.in(i,:)). The case otype==2 is not
%               supported;
%
% OUTPUTS
% - prog:   'struct' specifying the PIESOS program with the positive
%           semidefinite matrix variable added;
% - Pmat:   q x q cell, with element Pmat{i,j} an mi x mj 'dpvar' object
%           in the primary variables, where mi is the number of rows of
%           the ith block row of Zop.ops. The cell represents the SOS
%           matrix Pmat(s) parameterizing Pop = Zop'*Pmat*Zop. With
%           psatz multipliers g_k, Pmat = sum_k g_k*Pmat_k, where Pmat_k
%           is of the form Zs'*Qk*Zs, with Qk >= 0 and Zs a basis of
%           monomials in the primary variables;
% - Zop:    m x n 'tensopmat' object representing a basis of operators of
%           specified degrees on the declared function space, where n =
%           sum(dims) and where m depends on the selected monomial
%           degrees. The kernels of the operators depend only on the
%           dummy variables;
%
% NOTES
% - The monomials in the primary variables, s1 through sN, are split from
%   those in the dummy variables, s1_dum through sN_dum: the kernels of
%   Zop are defined only in terms of the dummy variables, and the
%   dependence of Pop on the primary variables is parameterized by the
%   SOS matrix Pmat.
% - For states depending on different variables, the operators Zop map
%   to a function space in all of the variables appearing in the state.
%   For example, for state x(s1) and y(s1,s2), Zop maps x to a function
%   of both s1 and s2.
% - The operator is declared by `piesos_quadopvar`, with opts.pos=true,
%   for the same function space and operators on the left and right, see
%   there for more details on the options.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - piesos_poslpivar
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
% DJ, 10/21/2026: Describe 3 x 1 'exclude' as OR over variables of a state;

% Check that the inputs are appropriate
if isstruct(dims) && isfield(dims,'in')
    dims = dims.in;     % monopmat needs only the input dimensions
end
if isa(vars,'struct')
    if isfield(vars,'out') && ~isequal(vars.in,vars.out)
        error("Input and output variables must match for self-adjoint positive semidefinite operator.")
    end
else
    vars = struct('in',{vars});
end
if isa(dom,'struct')
    if isfield(dom,'out') && ~isequal(dom.in,dom.out)
        error("Input and output domains must match for self-adjoint positive semidefinite operator.")
    end
else
    dom = struct('in',dom);
end
if isa(depmat,'struct')
    if isfield(depmat,'out') && ~isequal(depmat.in,depmat.out)
        error("Input and output function spaces must match for self-adjoint positive semidefinite operator.")
    end
else
    depmat = struct('in',depmat);
end

% Set default monomial degrees of 0 in the independent variables
if nargin<=5
    degs = 0;
end
if nargin<=6
    opts = struct();    % defaults of monopmat
end
if ~isstruct(opts)
    error("Options must be specified as a 'struct'.")
end

% Use piesos_quadopvar with 'pos' option to declare positive semidefinite
% operator
opts.pos = true;
[prog,Pmat,ZopL,ZopR] = piesos_quadopvar(prog,dims,vars.in,dom.in,depmat.in,degs,opts);

% Zop'*Pmat*Zop is positive only if the operators on the left and right are
% the same, which opts.L, opts.R, degs.L, degs.R may violate
if ~isequal(ZopL,ZopR)
    error("A positive semidefinite operator requires the same operators on the left and right: opts.L and opts.R, and degs.L and degs.R, must be equal.")
end
Zop = ZopL;

end