function [prog,V,Pmat,Zop] = piesos_sosvar(prog,Z,opts)
% [PROG,V,PMAT,ZOP] = PIESOS_SOSVAR(PROG,Z,OPTS) computes a distributed SOS
% polynomial decision variable based on the monomial basis Z,
%   V = <Zop*Z, Pmat*Zop*Z>
% V is declared by `piesos_quadvar`, with ZR=ZL and opts.pos=true.
%
% INPUTS
% - prog:   'struct' specifying a PIESOS program in which to declare the
%           decision variable (see `piesos_program`);
% - Z:      n x 1 'polyopvar' object representing a basis of distributed
%           monomials;
% - opts:   (optional) struct with fields:
%           - 'degs': struct or array specifying the degrees of the
%               monomials to use in the parameterization of the operators
%               Zop. See `piesos_poslpivar` for more details;
%               The degrees in degs.in define the monomials in the dummy 
%               variables, in terms of which the kernels of Zop are   
%               defined, those in degs.out the monomials in the primary
%               variables, defining the basis of the SOS matrix Pmat(s).
%               If degs is an array or cell, the same degrees are used;
%           - 'exclude': 3x1 boolean array indicating whether to exclude
%               the multiplier term (exclude(1)=1), lower integral
%               (exclude(2)=1), or upper integral (exclude(3)=1) from the
%               definition of basis operator Zop. For several spatial
%               variables, a term is excluded if its type is excluded
%               along any of the variables (see `monopmat`);
%           - 'psatz': 1 x N array of integers, or logical array,   
%               specifying psatz multipliers g_k(s)>=0 on the domain,  
%               with Pmat = sum_k g_k*Pmat_k and Pmat_k an SOS matrix, so  
%               that V>=0 is enforced only for states on the interval over 
%               which the operator integrates. If psatz=0, no multiplier
%               is used. If psatz(i)=qi, the multipliers are (s_i-a_i)^j*
%               (b_i-s_i)^(qi-j), j=0,..,qi, see `piesos_poslpivar`;
%           - 'sep': scalar boolean specifying whether to enforce
%               separability of the operator;
%           Other options of `piesos_quadvar` can be used, but the
%           operators on the left and right must be the same (opts.L, 
%           opts.R, degs.L, degs.R);
%
% OUTPUTS
% - prog:   'struct' specifying the PIESOS program with the decision
%           variables defining Pmat and the condition Pmat>=0 added;
% - V:      scalar 'polyopvar' object representing a distributed polynomial
%           functional;
% - Pmat:   q x q cell, with Pmat{i,j} an mi x mj 'dpvar' object in the
%           primary variables, representing the SOS matrix decision
%           variable Pmat(s) (one block row for each of the q monomials
%           in Z, see `piesos_poslpivar`);
% - Zop:    m x n 'tensopmat' object representing the monomial basis
%           operator acting on the distributed monomial basis defined by Z;
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
% DJ, 10/15/2026: Initial coding;
% DJ, 10/21/2026: Describe 'exclude' for several spatial variables;

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

% Declare the SOS functional using piesos_quadvar
opts.pos = true;
[prog,V,Pmat,ZopL,ZopR] = piesos_quadvar(prog,Z,[],opts);

% V is a sum of squares only if the operators on the left and right are the
% same, which opts.L, opts.R, opts.degs.L, opts.degs.R may violate.
if ~isequal(ZopL,ZopR)
    error("A distributed SOS polynomial requires the same operators on the left and right: opts.L and opts.R, and degs.L and degs.R, must be equal.")
end
Zop = ZopL;

end