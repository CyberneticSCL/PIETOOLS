function [sp,dm] = copvar_space_list(P,side)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [SP,DM] = COPVAR_SPACE_LIST(P,SIDE) lists the output (SIDE = 'out') or
% input (SIDE = 'in') spaces of the container P in the form 'poscopvar',
% 'lpivar_cdopvar' and 'parse_copvar_spaces' take: SP{k} the variable names
% of space k as a 1 x n cellstr in registry order (1 x 0 for R^q), DM the
% component counts. It is the inverse of 'parse_copvar_spaces', and stands
% in for passing an opvar's .dim (e.g. poslpivar(prog,Top.dim,...)) to the
% legacy LPI routines.
%
% INPUT
% - P:      'copvar' or 'cdopvar' object;
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
% See also PARSE_COPVAR_SPACES, POSLPIVAR_SOP, LPIVAR_CDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - copvar_space_list
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
% Initial coding MMP, 10/06/2026: library version of the test-folder
%                'cx_space_list' (cx_exec, MMP 09/25/2026), same body, plus
%                a check on SIDE (cx_space_list read any other value as
%                'in'). cx_space_list stays, unchanged, for the frozen 2-D
%                helpers that call it.

switch side
    case 'out',     S = P.space_out;    dm = P.dim_out(:);
    case 'in',      S = P.space_in;     dm = P.dim_in(:);
    otherwise
        error('copvar_space_list:side',"SIDE must be 'out' or 'in'.")
end
sp = cell(1,size(S,1));
for k = 1:size(S,1),    sp{k} = reshape(P.vars(S(k,:)),1,[]);    end
end
