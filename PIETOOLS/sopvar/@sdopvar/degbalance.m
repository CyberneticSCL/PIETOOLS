function [deg,report] = degbalance(P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [DEG,REPORT] = DEGBALANCE(P,OPTS) picks the degrees of the basis of a
% positive operator so that, declared with them, it carries monomials of
% similar degree to those of P: for every block, the maximal degree of P on
% every nonempty subset of variables, halved and rounded up (the
% Newton-polytope factor of a sum of squares in the Gram sense).
%
% INPUTS
% - P:      'sdopvar' object, the operator whose degrees are to be matched;
% - opts:   (optional) 'struct' with fields 'include' (nblk x n3 array, the
%           basis multi-indices to size; default every block), 'sep'
%           (1 x n3 logical, the separable directions) and 'weight'
%           ('uniform', default, halves every cap; 'slotwise' halves only
%           the integration slot);
%
% OUTPUTS
% - deg:    1 x nblk cell of degree structs (fields 'int', 'mult', 'joint',
%           'subset'), the degrees of the basis of a positive operator;
% - report: 'struct' with fields 'lines' (what was decided) and 'maxdeg'
%           (per-block per-subset maximal degrees read off P);
%
% See also DEGBALANCE_CORE, GET_EQ_OPTS, POSSOPVAR.
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
% Initial coding MMP, 09/21/2026

% Degrees read off P and halved, with or without options
if nargin<2
    [deg,report] = degbalance_core(P);
else
    [deg,report] = degbalance_core(P,opts);
end

end
