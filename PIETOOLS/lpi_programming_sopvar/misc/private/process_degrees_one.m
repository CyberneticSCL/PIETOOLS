function d = process_degrees_one(deg,n3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D = PROCESS_DEGREES_ONE(DEG,N3) completes the degree specification DEG of
% one monomial basis Z(theta,s) in the N3 integration variables theta_k and
% the N3 variables s_k: the caps int and mult as 1 x N3 arrays (int
% defaulting to 1, mult to int), the total cap joint (default
% sum(int)+sum(mult)) and the subset caps (default []).
%
% INPUT
% - deg:    'double' scalar, the cap int, or 'struct' with the optional
%           fields
%             int     double, a scalar or 1 x n3, the largest degree in
%                     each theta_k;
%             mult    double, a scalar or 1 x n3, the largest degree in
%                     each s_k;
%             joint   double scalar, the largest total degree in theta and
%                     s;
%             subset  double, 2^(2*n3) entries, the largest degree summed
%                     over each subset of the 2*n3 variables [theta_1 ..
%                     theta_n3, s_1 .. s_n3]: the subset with members b at
%                     linear index 1+sum(2.^(b-1)), entry 1 not used;
% - n3:     'double', the number of variables theta_k (and of s_k).
%
% OUTPUT
% - d:      'struct' with fields int, mult (1 x n3 double), joint (double)
%           and subset ([] or 1 x 2^(2*n3) double).
%
% NOTES
% DEG other than a numeric scalar or a 'struct', int or mult with other
% than 1 or n3 entries, a nonscalar joint and a subset with other than
% 2^(2*n3) entries are errors.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - process_degrees_one
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
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body; code unchanged. No help text deleted. The
%                file had no log: initial coding MP, 08/22/2026, as a local
%                function of 'possopvar' (the file since renamed
%                'sopquadvar'), moved to private/ by MMP, 09/21/2026; the
%                lines stamped MMP, 09/21/2026 add the subset caps.

% A scalar is the cap int
if isnumeric(deg) && isscalar(deg)
    deg = struct('int',deg);
elseif ~isa(deg,'struct')
    error("Monomial degrees should be specified as a scalar or 'struct' object.")
end

% Defaults: int 1, mult int
if ~isfield(deg,'int') || isempty(deg.int)
    deg.int = 1;
end
if ~isfield(deg,'mult') || isempty(deg.mult)
    deg.mult = deg.int;
end

% Per-variable caps as 1 x n3 rows
d = struct();
d.int = expand_caps(deg.int,n3,'int');
d.mult = expand_caps(deg.mult,n3,'mult');

% Total cap; default sum(int)+sum(mult)
if ~isfield(deg,'joint') || isempty(deg.joint)
    d.joint = sum(d.int)+sum(d.mult);
elseif ~isscalar(deg.joint)
    error("The joint degree should be specified as a scalar.")
else
    d.joint = deg.joint;
end

% Optional cap per variable SUBSET, over the order the basis is built in,    % MMP, 09/21/2026
% [theta_1..theta_n3, s_1..s_n3], indexed as 'build_exponent_grid'           % MMP, 09/21/2026
% documents. 'int', 'mult' and 'joint' are the singleton and full-set        % MMP, 09/21/2026
% special cases, so this subsumes them; empty leaves the basis determined    % MMP, 09/21/2026
% by those three alone, unchanged from before. Needed because a             % MMP, 09/21/2026
% 'poslpivar_2d' basis is cut out by 2^nvars subset caps, which per-variable % MMP, 09/21/2026
% plus total cannot express -- see 'settings2possopvar'.                     % MMP, 09/21/2026
if ~isfield(deg,'subset') || isempty(deg.subset)                            % MMP, 09/21/2026
    d.subset = [];                                                          % MMP, 09/21/2026
elseif numel(deg.subset)~=2^(2*n3)                                          % MMP, 09/21/2026
    error("A 'subset' degree array should have 2^(2*n3) entries for n3 "...
          +"spatial variables.")                                            % MMP, 09/21/2026
else                                                                        % MMP, 09/21/2026
    d.subset = reshape(deg.subset,1,[]);                                    % MMP, 09/21/2026
end                                                                         % MMP, 09/21/2026

end


