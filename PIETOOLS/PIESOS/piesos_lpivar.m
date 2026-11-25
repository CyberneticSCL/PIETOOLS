function [prog,Pop,Pmat,Zop] = piesos_lpivar(prog,dims,vars,dom,depmat,degs,opts)
% [PROG,POP,PMAT,ZOP] = PIESOS_LPIVAR(PROG,DIMS,VARS,DOM,DEPMAT,DEGS,OPTS)
% declares an indefinite operator variable between the specified function
% spaces, as
%       POP = PMAT*[ZOP{1};...;ZOP{p}]
% where each Zop{i} is a tensor-PI operator, built by `monopmat`,
%   Zop{i}: L2^{n1}[Omega_in{1}] x ... x L2^{nq}[Omega_in{q}]
%               --> L2^{mi}[Omega_out{i}]
% for
%   Omega_in{j} = [a1,b1]^{d_in(j,1)} x ... x [aN,bN]^{d_in(j,N)}
%   Omega_out{i} = [a1,b1]^{d_out(i,1)} x ... x [aM,bM]^{d_out(i,M)}
%
% INPUTS
% - prog:   'struct' specifying a PIESOS program in which to declare the
%           operator variable;
% - dims:   struct with fields
%     in:   q x 1 array specifying the number of state variables (ni) on
%           which the operator acts;
%    out:   p x 1 array specifying the dimensions of the output to the
%           operator;
% - vars:   struct with fields
%     in:   1 x N cellstr specifying the names of the input spatial
%           variables to the operator;
%    out:   1 x M cellstr specifying the names of the output spatial
%           variables to the operator;
% - dom:    struct with fields
%    in:    N x 2 array specifying the intervals on which the input
%           variables are defined;
%   out:    M x 2 array specifying the intervals on which the output
%           variables are defined;
% - depmat: struct with fields
%       in: q x N integer array specifying the degree to which each of the
%           states depends on each of the input variables (d_in);
%      out: p x M logical (or 0/1) array specifying whether or not each of
%           the outputs depends on each of the output variables (d_out
%           for otype 1);
% - degs:   (optional) struct with fields 'in' and 'out', the degrees of
%           the monomials in each of the q blocks of each Zop{i}, passed
%           to `monopmat` as DEGS:
%     in:   1 x q cell, element degs.in{i} a di x N cell of the monomial
%           degrees (including 0) in each factor of each input variable
%           of block i, or a di x N array of maximal degrees, for di the
%           number of factors, max(depmat.in(i,:));
%    out:   1 x q cell, element degs.out{i} a 1 x M cell of the monomial
%           degrees in each output variable, or a 1 x M array of maximal
%           degrees.
%           If degs is a numeric array or a cell, the same degrees are
%           used for 'in' and 'out'. Defaults to 0;
% - opts:   (optional) struct with fields, passed to `monopmat` as OPTS:
% exclude:  1 x q cell, specifying for each of the operator blocks whether
%           any terms need to be excluded. Element opts.exclude{i} can be a
%           3^N3(i) logical array, where N3(i) is the number of 3-PI
%           variables (appearing both in the input and output) of the
%           operator block i, with element k along direction j set to true
%           to exclude the multiplier term (k=1), lower integral term
%           (k=2), or upper integral term (k=3), along spatial variable j
%           in the list of 3-PI variables (intersect(vars.in,vars.out)).
%           A single array, or a 3 x 1 array, may also be used for all
%           blocks. By default, no terms are excluded;
%   sep:    1 x q cell, specifying for each of the operator blocks whether
%           to enforce separability along each of the variables. Element
%           opts.sep{i} is a 1 x N3(i) logical array, specifying for each
%           of the 3-PI variables whether the operator is separable along
%           that variable, meaning that the kernel in the upper and lower
%           integrals is the same. A single array may also be used for
%           all blocks. By default, no separability is enforced;
%  otype:   1, 2 or 3, specifying whether to take the pointwise product
%           (1), the proper tensor product (2), or both (3) of the 3-PI
%           operators in constructing the tensor-PI operator basis. If
%           otype==1, then d_out = depmat.out. If otype==2, then d_out =
%           d*depmat.out, where d = max(depmat.in(i,:)). Defaults to 1;
%
% OUTPUTS
% - prog:   'struct' specifying the same PIESOS program as the input, but
%           with the indefinite operator variable added to the structure;
% - Pop:    'struct' representing the indefinite operator variable, with
%           fields P=Pmat and Z=Zop, so that the operator is Pop.P*Pop.Z;
% - Pmat:   p x r cell, for r the number of block rows of Zop, with
%           Pmat{i,j} a dims.out(i) x rj 'dpvar' object, constant in the
%           spatial variables, for rj the number of rows of block row j
%           of Zop;
% - Zop:    'tensopmat' object [Zop{1};...;Zop{p}] representing the
%           monomial basis operator in terms of which Pop is defined,
%               Pop = Pmat * Zop;
%
% NOTES
% - Unrecognized fields of opts give a warning (PIESOS:unknownOption).
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - piesos_lpivar
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% DJ, 11/24/2026: Initial coding;


% Check that the output dependence is boolean
if ~islogical(depmat.out)
    if any(any(~ismember(depmat.out,[0;1])))
        error("Distributed monomial degree of the output cannot be set explicitly.")
    end
    depmat.out = logical(depmat.out);
end
% Set default monomial degrees of 0 in the independent variables
if nargin<=5
    degs = 0;
end
if nargin<=6
    opts = struct();
end
check_opts(opts,'piesos_lpivar');

% Declare a basis of operators on this domain, one block row per output
Zop = [];
ndims_out = size(depmat.out,1);
for i=1:ndims_out
    depmat_i = depmat;
    depmat_i.out = repmat(depmat.out(i,:),size(depmat.in,1),1);
    Zop_i = monopmat(dims.in,vars,dom,depmat_i,degs,opts);
    Zop = [Zop; Zop_i];
end

% Declare an indefinite matrix acting on the operator basis
Z1c = num2cell(ones(1,numel(dims.out)));
Z2c = num2cell(ones(1,numel(Zop.dims.out)));
[prog,Pmat] = sosquadvar(prog,Z1c,Z2c,dims.out,Zop.dims.out);

% Declare the output operator as a struct
Pop = struct();
Pop.P = Pmat;
Pop.Z = Zop;

end