function Pop = copvar2opvar2d(Pc,vars_xy)                                   % MMP, 09/30/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POP = COPVAR2OPVAR2D(PC) takes a 'copvar' container over the spaces       % MMP, 09/30/2026
% R x L2[x] x L2[y] x L2[x,y] and returns the equivalent 'opvar2d'.
%
% POP = COPVAR2OPVAR2D(PC,VARS_XY) names which registry variable plays x    % MMP, 09/30/2026
% and which plays y, as a 1 x 2 cellstr.
%
% INPUTS
% - Pc:       'copvar' over at most two spatial variables, whose output and % MMP, 09/30/2026
%             input space lists agree. Unlike 'sopvar2opvar2d', all sixteen
%             blocks may be populated;
% - vars_xy:  optional 1 x 2 cellstr {xname,yname}. DEFAULTS TO THE SORTED
%             REGISTRY ORDER, Pc.vars, which is alphabetical - so a         % MMP, 09/30/2026
%             container over {'y','x'} would otherwise silently put 'x'
%             first. Pass this explicitly whenever the intended roles are
%             not alphabetical;
%
% OUTPUTS
% - Pop:      'opvar2d' object representing the same operator;
%
% NOTES
% Each block is delegated to 'sopvar2opvar2d' and its one live component is
% copied into the corresponding slot of the result. A structurally zero
% block contributes nothing and leaves that component empty.
% 'sopvar2opvar2d' picks x and y itself (Bk.var1); where that is the        % MMP, 09/28/2026
% reverse of vars_xy the component is taken from the mirror slot (x<->y)    % MMP, 09/28/2026
% and a cell-valued one transposed, cell axis d being direction d.          % MMP, 09/28/2026
%
% Blocks are matched to components BY SPACE, not by grid position: the
% container sorts its registry and may omit rows for spaces the operator
% does not use, so row 2 of a 3 x 3 grid is not necessarily L2[x]. Each
% row's space mask over (x,y) is read as one of {}, {x}, {y}, {x,y} and
% mapped to opvar2d's index 0, x, y, 2; the component name is then the
% (out,in) entry of
%
%     R00 R0x R0y R02
%     Rx0 Rxx Rxy Rx2
%     Ry0 Ryx Ryy Ry2
%     R20 R2x R2y R22
%
% See also OPVAR2D2COPVAR, SOPVAR2OPVAR2D, COPVAR2OPVAR, COPVAR2NOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - copvar2opvar2d
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
%                  mopvar2opvar2d -> copvar2opvar2d. File was
%                  'mopvar2opvar2d.m'.
% MMP, 09/28/2026: A block that 'sopvar2opvar2d' orients the reverse way to
%                  vars_xy (its rule: sorted, a lone 's2'/'y' is y) is read
%                  from the mirror slot, cells transposed. Before, the
%                  same-named slot was copied: 12 of the 16 components came
%                  out empty (block dropped with its dim zeroed, or rejected
%                  by opvar2d's setter), and R22 had its alpha
%                  axes swapped, or was rejected by opvar2d's setter when it
%                  depends on a dummy. Hit every non-sorted vars_xy, and the
%                  DEFAULT too wherever the lone-variable rule disagrees with
%                  sorted order (registry {a,b}: lone b; {y,z}: lone y, z).
% MMP, 09/30/2026: Renamed the input Pmop -> Pc, here and in 'dom_matrix',
%                  a name from before the 09/25/2026 rename ('mopvar' is now
%                  the stub class in @mopvar). Mechanical, no functional
%                  change; each marked line differs from its old text only
%                  by that name. The one line with a 09/28/2026 marker keeps
%                  its old text as '(was)'.

if ~isa(Pc,'copvar')                                                        % MMP, 09/30/2026
    error('copvar2opvar2d:badInput','Input must be a copvar object.')
end
if numel(Pc.vars)>2                                                         % MMP, 09/30/2026
    error('copvar2opvar2d:tooManyVars',...
        ['''opvar2d'' is the 2D class; this container has %d spatial '...
         'variables (%s). Use ''copvar2nopvar''.'],...
        numel(Pc.vars),strjoin(Pc.vars,','))                                % MMP, 09/30/2026
end
if nargin<2 || isempty(vars_xy)
    vars_xy = Pc.vars;                                                      % MMP, 09/30/2026
end
vars_xy = vars_xy(:)';
if ~all(ismember(Pc.vars,vars_xy))                                          % MMP, 09/30/2026
    error('copvar2opvar2d:varsMismatch',...
        'vars_xy {%s} does not cover the registry {%s}.',...
        strjoin(vars_xy,','),strjoin(Pc.vars,','))                          % MMP, 09/30/2026
end

nm = {'R00','R0x','R0y','R02';
      'Rx0','Rxx','Rxy','Rx2';
      'Ry0','Ryx','Ryy','Ry2';
      'R20','R2x','R2y','R22'};
[M,N] = size(Pc);                                                           % MMP, 09/30/2026
ri = space_index(Pc.space_out,Pc.vars,vars_xy,'output');                    % MMP, 09/30/2026
ci = space_index(Pc.space_in ,Pc.vars,vars_xy,'input');                     % MMP, 09/30/2026
if numel(unique(ri))~=M || numel(unique(ci))~=N
    error('copvar2opvar2d:repeatedSpace',...
        ['''opvar2d'' carries one dimension per space, so no two rows or '...
         'columns may be over the same space.'])
end

Pop = opvar2d();
Pop.I = dom_matrix(Pc,vars_xy);                                             % MMP, 09/30/2026
if ~isempty(vars_xy)
    v1 = polynomial(zeros(numel(vars_xy),1));
    v2 = polynomial(zeros(numel(vars_xy),1));
    for k = 1:numel(vars_xy)
        v1(k) = polynomial(1,1,{vars_xy{k}},[1,1]);
        v2(k) = polynomial(1,1,{[vars_xy{k},'_dum']},[1,1]);
    end
    Pop.var1 = v1;      Pop.var2 = v2;
end
d = zeros(4,2);
d(ri,1) = Pc.dim_out(:);                                                    % MMP, 09/30/2026
d(ci,2) = Pc.dim_in(:);                                                     % MMP, 09/30/2026
Pop.dim = d;

mir = [1 3 2 4];    % opvar2d space index with x and y exchanged            % MMP, 09/28/2026
for i = 1:M
    for j = 1:N
        if isempty(Pc.C{i,j}),    continue,   end                           % MMP, 09/30/2026
        Bk = sopvar2opvar2d(Pc.C{i,j});                                     % MMP, 09/30/2026
        slot = nm{ri(i),ci(j)};
%       Pop.(slot) = Bk.(slot);                                             % MMP, 09/28/2026 (was)
        % Bk's x,y follow sopvar2opvar2d's rule, not vars_xy. opvar2d has   % MMP, 09/28/2026
        % two directions, so one variable of the block decides: Bk agrees   % MMP, 09/28/2026
        % with vars_xy or is its reverse. Reversed: mirror slot, and cell   % MMP, 09/28/2026
        % axes swapped (R22{a,b} -> {b,a}; Ry2 1x3 -> Rx2 3x1, etc.).       % MMP, 09/28/2026
%       v = [Pmop.C{i,j}.vars.in, Pmop.C{i,j}.vars.out];                    % MMP, 09/30/2026 (was)
        v = [Pc.C{i,j}.vars.in, Pc.C{i,j}.vars.out];                        % MMP, 09/30/2026
        sw = ~isempty(v) && find(strcmp(pvar2varname(Bk.var1),v{1}),1)...
                         ~= find(strcmp(vars_xy,v{1}),1);                   % MMP, 09/28/2026
        if sw                                                               % MMP, 09/28/2026
            val = Bk.(nm{mir(ri(i)),mir(ci(j))});                           % MMP, 09/28/2026
            if iscell(val),     val = val.';    end                         % MMP, 09/28/2026
        else                                                                % MMP, 09/28/2026
            val = Bk.(slot);                                                % MMP, 09/28/2026
        end                                                                 % MMP, 09/28/2026
        Pop.(slot) = val;                                                   % MMP, 09/28/2026
    end
end

end

% ------------------------------------------------------------------------
function idx = space_index(mask,vars,vars_xy,side)
% opvar2d's space index for each row of a mask: 1 for R, 2 for L2[x], 3 for
% L2[y], 4 for L2[x,y], with x and y as named in vars_xy.
nsp = size(mask,1);
idx = zeros(nsp,1);
has_x = false(nsp,1);       has_y = false(nsp,1);
if numel(vars_xy)>=1
    k = find(strcmp(vars,vars_xy{1}),1);
    if ~isempty(k),     has_x = mask(:,k);      end
end
if numel(vars_xy)>=2
    k = find(strcmp(vars,vars_xy{2}),1);
    if ~isempty(k),     has_y = mask(:,k);      end
end
for i = 1:nsp
    if     ~has_x(i) && ~has_y(i),  idx(i) = 1;
    elseif  has_x(i) && ~has_y(i),  idx(i) = 2;
    elseif ~has_x(i) &&  has_y(i),  idx(i) = 3;
    else,                           idx(i) = 4;
    end
    if sum(mask(i,:))>2
        error('copvar2opvar2d:tooManyVarsInSpace',...
            '%s space %d has more than two variables.',side,i)
    end
end
end

function I = dom_matrix(Pc,vars_xy)                                         % MMP, 09/30/2026
% opvar2d's I is 2x2, row k the domain of vars_xy{k}. A direction the
% container never used has no recorded domain, so it defaults to [0,1].
I = [0 1;0 1];
for k = 1:min(2,numel(vars_xy))
    idx = find(strcmp(Pc.vars,vars_xy{k}),1);                               % MMP, 09/30/2026
    if ~isempty(idx),   I(k,:) = Pc.dom(idx,:);  end                        % MMP, 09/30/2026
end
end
