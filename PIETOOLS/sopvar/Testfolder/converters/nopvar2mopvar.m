function Pmop = nopvar2mopvar(Pnop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PMOP = NOPVAR2MOPVAR(PNOP) takes a 'nopvar' object and returns the
% equivalent 1 x 1 'mopvar' container.
%
% INPUTS
% - Pnop:   'nopvar' object, P: L2^n[s_1,...,s_N] -> L2^m[s_1,...,s_N];
% OUTPUTS
% - Pmop:   1 x 1 'mopvar' holding the single equivalent 'sopvar' block;
%
% NOTES
% A 'nopvar' has one space on each side and no finite-dimensional
% component, so its container is 1 x 1 and this routine is a thin wrapper
% over 'nopvar2sopvar'. It exists so that the mopvar converters cover every
% legacy class rather than most of them, and so a caller can move a
% 'nopvar' into the container world without knowing which block converter
% applies.
%
% Wrapping a 'nopvar' this way buys generality and costs speed, and that
% trade is the intended design of the two classes, not a defect in either:
% 'nopvar' is deliberately restricted - one space, no R component, and a
% single shared degree vector giving the full tensor basis of that degree on
% both sides - and those restrictions are what make its addition a plain
% matrix add where 'sopvar' must union per-variable degree sets. Measured
% ratios for the wrapped form against 'nopvar' itself are roughly 4.9x on
% plus and 1.0x on mtimes; see 'Bench_mopvar_vs_legacy'. Use 'nopvar' when
% the problem fits its restrictions.
%
% See also MOPVAR2NOPVAR, NOPVAR2SOPVAR, OPVAR2MOPVAR, BENCH_MOPVAR_VS_LEGACY.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - nopvar2mopvar
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

if ~isa(Pnop,'nopvar')
    error('nopvar2mopvar:badInput','Input must be a nopvar object.')
end
Pmop = mopvar({nopvar2sopvar(Pnop)});

end
