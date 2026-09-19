function Pop = mopvar2opvar(Pmop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POP = MOPVAR2OPVAR(PMOP) takes a 'mopvar' container over the spaces
% (R^m, L2^n[s]) and returns the equivalent 'opvar' object.
%
% INPUTS
% - Pmop:   'mopvar' whose output and input spaces are the same list drawn
%           from {R^m, L2^n[s]} - that is, a 1 x 1 or 2 x 2 grid over one
%           spatial variable. Unlike 'sopvar2opvar', all four blocks may be
%           populated;
%
% OUTPUTS
% - Pop:    'opvar' object representing the same operator, with
%           P = C{1,1}, Q1 = C{1,2}, Q2 = C{2,1}, R = C{2,2} when both
%           spaces are present;
%
% NOTES
% The per-block conversion is delegated to 'sopvar2opvar' rather than
% reimplemented; each block comes back as a single-component 'opvar' and its
% one live component is copied into the corresponding slot. A structurally
% zero block ([]) contributes nothing and leaves that component empty.
%
% The grid is matched to the spaces by NAME, not by position: the container
% sorts its variable registry and a 1 x 1 grid may be either the R block or
% the L2 block, so which row is which is decided by whether that space's
% mask is empty. A row whose space has one variable is the L2 row; a row
% whose space is empty is the R row.
%
% Restrictions, each of which errors rather than producing a wrong operator:
% at most one spatial variable, since 'opvar' is the 1D class; the output
% and input space lists must agree, since 'opvar' has a single dim per
% space; and no more than two of each.
%
% See also OPVAR2MOPVAR, SOPVAR2OPVAR, MOPVAR2OPVAR2D, MOPVAR2NOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mopvar2opvar
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
    error('mopvar2opvar:badInput','Input must be a mopvar object.')
end
if numel(Pmop.vars)>1
    error('mopvar2opvar:tooManyVars',...
        ['''opvar'' is the 1D class; this container has %d spatial '...
         'variables (%s). Use ''mopvar2opvar2d'' or ''mopvar2nopvar''.'],...
        numel(Pmop.vars),strjoin(Pmop.vars,','))
end
[M,N] = size(Pmop);
% Identify each row and column as the R space (empty mask) or the L2 space.
row_is_L2 = any(Pmop.space_out,2);
col_is_L2 = any(Pmop.space_in ,2);
if M>2 || N>2 || numel(unique(row_is_L2))~=M || numel(unique(col_is_L2))~=N
    error('mopvar2opvar:badGrid',...
        ['An opvar has one R space and one L2 space, so the grid must have '...
         'at most one of each; got %dx%d.'],M,N)
end
if ~isequal(sort(row_is_L2(:)),sort(col_is_L2(:)))
    error('mopvar2opvar:asymmetric',...
        ['''opvar'' carries a single dimension per space, so the output and '...
         'input space lists must agree.'])
end

Pop = opvar();
Pop.I = Pmop.dom;
if isempty(Pop.I),      Pop.I = [0,1];      end
% opvar needs both a primary and a dummy variable name. The container stores
% only the primary, as 'sopvar' does, so the dummy follows the '<primary>_dum'
% convention that 'sopvar2opvar' also applies.
if ~isempty(Pmop.vars)
    Pop.var1 = polynomial(1,1,{Pmop.vars{1}},[1,1]);
    Pop.var2 = polynomial(1,1,{[Pmop.vars{1},'_dum']},[1,1]);
end
% dim = [out_R in_R; out_L2 in_L2], zero where that space is absent.
d = zeros(2,2);
for i = 1:M
    d(1+row_is_L2(i),1) = Pmop.dim_out(i);
end
for j = 1:N
    d(1+col_is_L2(j),2) = Pmop.dim_in(j);
end
Pop.dim = d;

nm = {'P' ,'Q1';
      'Q2','R'};
for i = 1:M
    for j = 1:N
        if isempty(Pmop.C{i,j}),    continue,   end
        Bk = sopvar2opvar(Pmop.C{i,j});
        slot = nm{1+row_is_L2(i),1+col_is_L2(j)};
        Pop.(slot) = Bk.(slot);
    end
end

end
