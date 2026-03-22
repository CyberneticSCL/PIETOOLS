function [prog,V,Pmat,Zop] = piesos_sosvar(prog,Z,opts)
% [PROG,V,PMAT,ZOP] = PIESOS_SOSVAR(PROG,Z,OPTS) computes a distributed SOS
% polynomial decision variable based on the monomial basis Z, 
%   V = <Zop*Z, Pmat*Zop*Z>
%
% INPUTS
% - prog:   LPI program structure representing the optimization program;
% - Z:      n x 1 'polyopvar' object representing a basis of distributed
%           monomials;
% - opts:   (optional) struct with fields:
%           - 'deg': cell or array specifying the degrees of the monomials
%               to use in the parameterization of the operators Zop.
%               See `piesos_poslpivar` for more details;
%           - 'exclude': 3x1 boolean array indicating whether to exclude 
%               the multiplier term (exclude(1)=1), lower integral 
%               (exclude(2)=1), or upper integral (exclude(3)=1) from the 
%               definition of basis operator Zop;
%           - 'psatz': 1xp array of integers specifying whether to use a
%               psatz multiplier to enforce positivity only on the interval
%               over which the operator integrates. If psatz=0, no such
%               multiplier is used. If psatz=1, the matrix Pmat will be
%               multiplied by (s-a)*(b-s), where [a,b]=Z.dom. If
%               psatz=[0,1], we have Pmat = Pmat1 + (s-a)*(b-s)*Pmat2;
%           - 'sep': scalar boolean specifying whether to enforce
%               separability of the operator;
%
% OUTPUTS:
% - prog:   LPI program structure with the decision variables defining Pmat
%           and the condition P>=0 added;
% - V:      scalar 'polyopvar' object representing a distributed polynomial
%           functional;
% - Pmat:   m x m 'dpvar' object representing the positive semidefinite
%           matrix decision variable;
% - Zop:    m x n 'tensopvar' object representing the monomial basis
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
% DJ, 03/05/2026: Initial coding


% Verify that the monomials are appropriately specified
if ~isa(Z,'polyopvar')
    error("Basis of distributed monomials must be specified as 'polyopvar' object.")
end
% Declare default monomial degree if not specified
if nargin<=2
    opts = struct();
end
if ~isa(opts,'struct')
    error("Options must be specified as a 'struct'.")
elseif ~isfield(opts,'deg')
    opts.deg = 0;
end

% Declare a positive semidefinite PI operator, Zop'*Pmat*Zop, acting on the
% desired monomial basis
pdegs = opts.deg;
[prog,Pmat,Zop] = piesos_poslpivar(prog,Z,pdegs,opts);

% Convert the distributed polynomial defined by <Zop*Z, Pmat*Zop*Z>
V = quad2lin(Pmat,Zop,Z); 

end