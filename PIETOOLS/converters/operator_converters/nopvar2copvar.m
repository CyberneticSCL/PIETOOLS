function Pmop = nopvar2copvar(Pnop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PMOP = NOPVAR2COPVAR(PNOP) takes a 'nopvar' object and returns the
% equivalent 1 x 1 'copvar' container.
%
% INPUTS
% - Pnop:   'nopvar' object, P: L2^n[s_1,...,s_N] -> L2^m[s_1,...,s_N];
% OUTPUTS
% - Pmop:   1 x 1 'copvar' holding the single equivalent 'sopvar' block;
%
% NOTES
% A 'nopvar' has one space on each side and no finite-dimensional
% component, so its container is 1 x 1 and this routine is a thin wrapper
% over 'nopvar2sopvar'. It exists so that the copvar converters cover every
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
% plus and 1.0x on mtimes; see 'Bench_copvar_vs_legacy'. Use 'nopvar' when
% the problem fits its restrictions.
%
% See also COPVAR2NOPVAR, NOPVAR2SOPVAR, OPVAR2COPVAR, BENCH_COPVAR_VS_LEGACY.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - nopvar2copvar
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
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change. Renamed here:
%                  Bench_mopvar_vs_legacy -> Bench_copvar_vs_legacy,
%                  nopvar2mopvar -> nopvar2copvar. File was 'nopvar2mopvar.m'.

if ~isa(Pnop,'nopvar')
    error('nopvar2copvar:badInput','Input must be a nopvar object.')
end
Pmop = copvar({nopvar2sopvar(Pnop)});

end
