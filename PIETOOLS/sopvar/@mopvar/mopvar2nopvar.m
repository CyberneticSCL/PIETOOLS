function Pnop = mopvar2nopvar(Pmop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PNOP = MOPVAR2NOPVAR(PMOP) takes a 1 x 1 'mopvar' container whose block
% maps one L2 space to itself and returns the equivalent 'nopvar' object.
%
% INPUTS
% - Pmop:   1 x 1 'mopvar' whose single block has the same spatial
%           variables on both sides;
% OUTPUTS
% - Pnop:   'nopvar' object representing the same operator;
%
% NOTES
% Thin wrapper over 'sopvar2nopvar'. The grid must be 1 x 1 and the block's
% input and output spaces must agree, because 'nopvar' has no
% finite-dimensional component and no block structure to put one in; a
% container with an R space or more than one space per side has no 'nopvar'
% equivalent and is rejected rather than silently truncated.
%
% Converting INTO 'nopvar' is usually the profitable direction: it drops the
% container and the per-variable degree sets in favour of one shared degree
% vector, which is exactly what makes that class faster. See
% 'Bench_mopvar_vs_legacy' for the measured ratios.
%
% See also NOPVAR2MOPVAR, SOPVAR2NOPVAR, MOPVAR2OPVAR, MOPVAR2OPVAR2D.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mopvar2nopvar
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
% Initial coding MMP, 09/18/2026

if ~isa(Pmop,'mopvar')
    error('mopvar2nopvar:badInput','Input must be a mopvar object.')
end
[M,N] = size(Pmop);
if M~=1 || N~=1
    error('mopvar2nopvar:badGrid',...
        ['''nopvar'' has one space on each side and no block structure, so '...
         'only a 1x1 container converts; got %dx%d.'],M,N)
end
if isempty(Pmop.C{1,1})
    error('mopvar2nopvar:zeroBlock',...
        ['The single block is structurally zero, so there is no operator to '...
         'convert. Build an explicitly zero-valued sopvar block instead.'])
end
if ~isequal(Pmop.space_out(1,:),Pmop.space_in(1,:))
    error('mopvar2nopvar:spaceMismatch',...
        ['''nopvar'' maps a space to itself; this container maps {%s} to '...
         '{%s}.'],strjoin(Pmop.vars(Pmop.space_in(1,:)),','),...
        strjoin(Pmop.vars(Pmop.space_out(1,:)),','))
end
Pnop = sopvar2nopvar(Pmop.C{1,1});

end
