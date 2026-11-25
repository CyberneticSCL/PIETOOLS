function [prog,V,Pmat,Zop] = piesos_sosvar(prog,Z,opts)
% [PROG,V,PMAT,ZOP] = PIESOS_SOSVAR(PROG,Z,OPTS) computes a distributed SOS
% polynomial decision variable based on the monomial basis Z,
%   V = <Zop*Z, Pmat*Zop*Z>
% V is declared by `piesos_quadvar`, with ZR=ZL and opts.pos=true.
%
% INPUTS
% - prog:   'struct' specifying a PIESOS program in which to declare the
%           decision variable;
% - Z:      nZ x 1 'polyopvar' object representing a basis of distributed
%           monomials, or a vector of other distributed polynomials, e.g.
%           T*x for an operator T;
% - opts:   (optional) struct with fields:
%           - 'degs': struct with fields 'in' and 'out', or an array or
%               cell used for both, specifying the degrees of the
%               monomials in the dummy variables (kernels of Zop) and in
%               the primary variables (basis of the SOS matrix Pmat(s)),
%               respectively, as in `piesos_quadopvar`. Defaults to 0;
%           - 'exclude': 3 x 1 logical array indicating whether to exclude
%               the multiplier term (exclude(1)=1), lower integral
%               (exclude(2)=1), or upper integral (exclude(3)=1) from the
%               definition of basis operator Zop. For several spatial
%               variables, a term is excluded if its type is excluded
%               along any of the variables;
%           - 'psatz': 1 x N array of integers, or logical array,
%               specifying psatz multipliers g_k(s)>=0 on the domain,
%               with Pmat = sum_k g_k*Pmat_k and Pmat_k an SOS matrix, so
%               that Pmat(s) need only be positive semidefinite on the
%               domain. If psatz=0, no multiplier is used. If psatz(i)=qi,
%               the multipliers are (s_i-a_i)^j*(b_i-s_i)^(qi-j),
%               j=0,...,qi;
%           - 'sep': scalar boolean specifying whether to enforce
%               separability of the operator;
%           - 'otype': type of output of Zop: 1 for pointwise products, 2
%               for tensor products, 3 for both, see `piesos_quadopvar`.
%               Defaults to 1;
%           Other options of `piesos_quadvar` can be used, but the
%           operators on the left and right must be the same (opts.L,
%           opts.R, degs.L, degs.R);
%
% OUTPUTS
% - prog:   'struct' specifying the PIESOS program with the decision
%           variables defining Pmat and the condition Pmat>=0 added;
% - V:      scalar 'polyopvar' object representing a distributed SOS
%           polynomial functional;
% - Pmat:   r x r cell, for r the number of block rows of Zop, with
%           Pmat{i,j} an mi x mj 'dpvar' object in the primary variables,
%           representing the SOS matrix decision variable Pmat(s). With
%           regions (opts.intmat), a cell of such cells, as in
%           `piesos_quadvar`;
% - Zop:    m x n 'tensopmat' object representing the monomial basis
%           operator acting on the distributed monomial basis defined by Z;
%
% NOTES
% - Unrecognized fields of opts give a warning (PIESOS:unknownOption).
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - piesos_sosvar
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

% Verify that the monomials are appropriately specified
if ~isa(Z,'polyopvar')
    error("Basis of distributed monomials must be specified as 'polyopvar' object.")
end
if ~Z.chkval
    % Guarantee Z.Z/Z.C are mutually consistent
    Z = initialize(Z);
end

% Declare default monomial degree if not specified
if nargin<=2
    opts = struct();
end
if ~isa(opts,'struct')
    error("Options must be specified as a 'struct'.")
end
check_opts(opts,'piesos_sosvar');

% Declare the SOS functional using piesos_quadvar, by default with
% pointwise products only (otype 1)
opts.pos = true;
if ~isfield(opts,'otype')
    opts.otype = 1;
end
[prog,V,Pmat,ZopL,ZopR] = piesos_quadvar(prog,Z,[],opts);

% V is a sum of squares only if the operators on the left and right are the
% same, which opts.L, opts.R, opts.degs.L, opts.degs.R may violate.
if ~isequal(ZopL,ZopR)
    error("A distributed SOS polynomial requires the same operators on the left and right: opts.L and opts.R, and degs.L and degs.R, must be equal.")
end
Zop = ZopL;

end