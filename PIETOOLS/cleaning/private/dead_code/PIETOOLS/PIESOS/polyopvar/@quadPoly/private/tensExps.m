function [E, Z3] = tensExps(Z1,Z2)
% Given two monomial column vectors Z1, Z2, this
% function calculates Z3 such that E*Z3 = Z1\otimes Z2.
%
% INPUTS
% - Z1, Z2: numeric column vectors;
%
% OUTPUTS
% - E:      numel(Z1)*numel(Z2) x numel(Z3) sparse 0/1 matrix with
%           E*Z3 = Z1\otimes Z2;
% - Z3:     column vector of the sorted unique entries of Z1\otimes Z2.
%
% NOTES
% - The entries of Z1\otimes Z2 are the products Z1(i)*Z2(j): applied to
%   vectors of exponents, they are not the exponents of the products of
%   the monomials.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - tensExps
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
% SS, 01/31/2026: Initial coding;

Zfull = kron(Z1,Z2);

[Z3,~, E] = unique(Zfull);

E = sparse(1:numel(Zfull), E, 1, numel(Zfull), numel(Z3));
end