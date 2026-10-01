function Pop = copvar2opvar(Pc)                                             % MMP, 09/30/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POP = COPVAR2OPVAR(PC) takes a 'copvar' container over the spaces         % MMP, 09/30/2026
% (R^m, L2^n[s]) and returns the equivalent 'opvar' object.
%
% INPUTS
% - Pc:     'copvar' whose output and input spaces are the same list drawn  % MMP, 09/30/2026
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
% See also OPVAR2COPVAR, SOPVAR2OPVAR, COPVAR2OPVAR2D, COPVAR2NOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - copvar2opvar
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
%                  mopvar2nopvar -> copvar2nopvar,
%                  mopvar2opvar -> copvar2opvar,
%                  mopvar2opvar2d -> copvar2opvar2d. File was 'mopvar2opvar.m'.
% MMP, 09/30/2026: Renamed the input Pmop -> Pc, a name from before the
%                  09/25/2026 rename ('mopvar' is now the stub class in
%                  @mopvar). Mechanical, no functional change; each marked
%                  line differs from its old text only by that name.

if ~isa(Pc,'copvar')                                                        % MMP, 09/30/2026
    error('copvar2opvar:badInput','Input must be a copvar object.')
end
if numel(Pc.vars)>1                                                         % MMP, 09/30/2026
    error('copvar2opvar:tooManyVars',...
        ['''opvar'' is the 1D class; this container has %d spatial '...
         'variables (%s). Use ''copvar2opvar2d'' or ''copvar2nopvar''.'],...
        numel(Pc.vars),strjoin(Pc.vars,','))                                % MMP, 09/30/2026
end
[M,N] = size(Pc);                                                           % MMP, 09/30/2026
% Identify each row and column as the R space (empty mask) or the L2 space.
row_is_L2 = any(Pc.space_out,2);                                            % MMP, 09/30/2026
col_is_L2 = any(Pc.space_in ,2);                                            % MMP, 09/30/2026
if M>2 || N>2 || numel(unique(row_is_L2))~=M || numel(unique(col_is_L2))~=N
    error('copvar2opvar:badGrid',...
        ['An opvar has one R space and one L2 space, so the grid must have '...
         'at most one of each; got %dx%d.'],M,N)
end
if ~isequal(sort(row_is_L2(:)),sort(col_is_L2(:)))
    error('copvar2opvar:asymmetric',...
        ['''opvar'' carries a single dimension per space, so the output and '...
         'input space lists must agree.'])
end

Pop = opvar();
Pop.I = Pc.dom;                                                             % MMP, 09/30/2026
if isempty(Pop.I),      Pop.I = [0,1];      end
% opvar needs both a primary and a dummy variable name. The container stores
% only the primary, as 'sopvar' does, so the dummy follows the '<primary>_dum'
% convention that 'sopvar2opvar' also applies.
if ~isempty(Pc.vars)                                                        % MMP, 09/30/2026
    Pop.var1 = polynomial(1,1,{Pc.vars{1}},[1,1]);                          % MMP, 09/30/2026
    Pop.var2 = polynomial(1,1,{[Pc.vars{1},'_dum']},[1,1]);                 % MMP, 09/30/2026
end
% dim = [out_R in_R; out_L2 in_L2], zero where that space is absent.
d = zeros(2,2);
for i = 1:M
    d(1+row_is_L2(i),1) = Pc.dim_out(i);                                    % MMP, 09/30/2026
end
for j = 1:N
    d(1+col_is_L2(j),2) = Pc.dim_in(j);                                     % MMP, 09/30/2026
end
Pop.dim = d;

nm = {'P' ,'Q1';
      'Q2','R'};
for i = 1:M
    for j = 1:N
        if isempty(Pc.C{i,j}),    continue,   end                           % MMP, 09/30/2026
        Bk = sopvar2opvar(Pc.C{i,j});                                       % MMP, 09/30/2026
        slot = nm{1+row_is_L2(i),1+col_is_L2(j)};
        Pop.(slot) = Bk.(slot);
    end
end

end
