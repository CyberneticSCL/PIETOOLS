function [sp,dm] = opvar_space_list(P,side)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [SP,DM] = OPVAR_SPACE_LIST(P,SIDE) lists the output (SIDE = 'out') or
% input (SIDE = 'in') spaces of the container P: SP{k} the variable names
% of space k as a 1 x n cellstr in registry order (1 x 0 for R^q), DM the
% component counts.
%
% INPUT
% - P:      'opvar' or 'dopvar' object;
% - side:   'out' or 'in'.
%
% OUTPUT
% - sp:     1 x M cell of cellstr, M = number of spaces on that side;
% - dm:     M x 1 component counts.
%
% NOTES
% Cost: O(M * nv) logical indexing (nv registry variables); nothing in the
% number of decision variables.
%
% See also PARSE_OPVAR_SPACES, POSLPIVAR_LEGACY, LPIVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - opvar_space_list
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
% Initial coding MMP, 10/06/2026

% Space masks (one row per space over the registry) and component counts
% of the chosen side
switch side
    case 'out',     S = P.space_out;    dm = P.dim_out(:);
    case 'in',      S = P.space_in;     dm = P.dim_in(:);
    otherwise
        error('opvar_space_list:side',"SIDE must be 'out' or 'in'.")
end
% Variable names of each space, in registry order
sp = cell(1,size(S,1));
for k = 1:size(S,1),    sp{k} = reshape(P.vars(S(k,:)),1,[]);    end
end
