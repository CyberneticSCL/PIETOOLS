function [prog,Pop,Pmat,Zop] = piesos_lpivar(prog,dims,vars,dom,depmat,degs,opts)
% [PROG,POP,PMAT,ZOP] = PIESOS_LPIVAR(PROG,DIMS,VARS,DOM,DEPMAT,DEGS,OPTS)
% declares an indefinite operator variable between the specified function
% spaces, as
%       POP = PMAT*[ZOP{1};...;ZOP{q}]
% where each Zop{i} is a tensor-PI operator,
%   Zop{i}: L2^{n1}[Omega_in1] x ... x L2^{nq}[Omega_in{q}]
%               --> L2^{mi}[Omega_out{i}]
% for
%   Omega_in = [a1,b1]^{d_in(i,1)} x ... x [aN,bN]^{d_in(i,N)}
%   Omega_out = [a1,b1]^{d_out(i,1)} x ... x [aM,bM]^{d_out(i,M)}
%
% INPUTS
% - prog:   'struct' specifying a PIESOS program in which to declare the
%           operator variable (see `piesos_program`);
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
%      out: p x M logical array specifying whether or not each of the
%           outputs depends on each of the output variables (d_out);
% - degs:   (optional) struct with fields
%     in:   1 x q cell, specifying the degrees of the input variables in
%           each operator block. Each element degs.in{i} can be a di x N
%           cell, where di depends on depmat.in(i,:), specifying the
%           degrees of the monomial basis in each of the factors of each of
%           the input variables in block i. Element degs.in{i} can also be
%           a di x N array, specifying the maximal monomial degree in each
%           of the factors of each variable. If degs.in{i} is scalar, the
%           same degree is used for all variables;
%    out:   1 x q cell, specifying the degrees of the output variables in
%           each operator block. Each element degs.out{i} can be a 1 x M
%           cell, specifying the degrees of the monomial basis in each of
%           the output variables in block i, or a 1 x M array, specifying
%           the maximal monomial degree in each variable.
%           If degs is a numeric array or a cell, the same degrees are
%           used for 'in' and 'out'. Defaults to 0;
% - opts:   (optional) struct with fields
% exclude:  1 x q cell, specifying for each of the operator blocks whether
%           any terms need to be excluded. Element opts.exclude{i} can be a
%           3^N3(i) logical array, where N3(i) is the number of 3-PI
%           variables (appearing both in the input and output) of the
%           operator block i, with element k along direction j set to true
%           to exclude the multiplier term (k=1), lower integral term
%           (k=2), or upper integral term (k=3), along spatial variable j
%           in the list of 3-PI variables (intersect(vars.in,vars.out));
%   sep:    1 x q cell, specifying for each of the operator blocks whether
%           to enforce separability along each of the variables. Element
%           opts.sep{i} is a 1 x N3(i) logical array, specifying for each
%           of the 3-PI variables whether the operator is separable along
%           that variable, meaning that the kernel in the upper and lower
%           integrals is the same;
%  otype:   integer in {1,2} specifying whether to take the proper tensor
%           product (otype = 2) or the pointwise product (otype = 1) of the
%           3-PI operators in constructing the tensor-PI operator basis. If
%           otype==1, then d_out = depmat.out. If otype==2, then d_out =
%           d*depmat.out, where d = max(depmat.in(i,:));
%
% OUTPUTS
% - prog:   'struct' specifying the same PIESOS program as the input, but
%           with the indefinite operator variable added to the structure;
% - Pop:    'struct' representing the indefinite operator variable, with
%           fields P=Pmat and Z=Zop, so that the operator is Pop.P*Pop.Z;
% - Pmat:   m x nZ 'dpvar' object representing the matrix variable defining
%           the operator;
% - Zop:    nZ x n 'tensopmat' object representing the monomial basis
%           operator in terms of which Pop is defined,
%               Pop = Pmat * Zop;
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
% DJ, 10/15/2026: Initial coding;


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

% Declare a basis of operators on this domain
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