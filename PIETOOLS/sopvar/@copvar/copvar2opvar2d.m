function Pop = copvar2opvar2d(Pmop,vars_xy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POP = COPVAR2OPVAR2D(PMOP) takes a 'copvar' container over the spaces
% R x L2[x] x L2[y] x L2[x,y] and returns the equivalent 'opvar2d'.
%
% POP = COPVAR2OPVAR2D(PMOP,VARS_XY) names which registry variable plays x
% and which plays y, as a 1 x 2 cellstr.
%
% INPUTS
% - Pmop:     'copvar' over at most two spatial variables, whose output and
%             input space lists agree. Unlike 'sopvar2opvar2d', all sixteen
%             blocks may be populated;
% - vars_xy:  optional 1 x 2 cellstr {xname,yname}. DEFAULTS TO THE SORTED
%             REGISTRY ORDER, Pmop.vars, which is alphabetical - so a
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

if ~isa(Pmop,'copvar')
    error('copvar2opvar2d:badInput','Input must be a copvar object.')
end
if numel(Pmop.vars)>2
    error('copvar2opvar2d:tooManyVars',...
        ['''opvar2d'' is the 2D class; this container has %d spatial '...
         'variables (%s). Use ''copvar2nopvar''.'],...
        numel(Pmop.vars),strjoin(Pmop.vars,','))
end
if nargin<2 || isempty(vars_xy)
    vars_xy = Pmop.vars;
end
vars_xy = vars_xy(:)';
if ~all(ismember(Pmop.vars,vars_xy))
    error('copvar2opvar2d:varsMismatch',...
        'vars_xy {%s} does not cover the registry {%s}.',...
        strjoin(vars_xy,','),strjoin(Pmop.vars,','))
end

nm = {'R00','R0x','R0y','R02';
      'Rx0','Rxx','Rxy','Rx2';
      'Ry0','Ryx','Ryy','Ry2';
      'R20','R2x','R2y','R22'};
[M,N] = size(Pmop);
ri = space_index(Pmop.space_out,Pmop.vars,vars_xy,'output');
ci = space_index(Pmop.space_in ,Pmop.vars,vars_xy,'input');
if numel(unique(ri))~=M || numel(unique(ci))~=N
    error('copvar2opvar2d:repeatedSpace',...
        ['''opvar2d'' carries one dimension per space, so no two rows or '...
         'columns may be over the same space.'])
end

Pop = opvar2d();
Pop.I = dom_matrix(Pmop,vars_xy);
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
d(ri,1) = Pmop.dim_out(:);
d(ci,2) = Pmop.dim_in(:);
Pop.dim = d;

for i = 1:M
    for j = 1:N
        if isempty(Pmop.C{i,j}),    continue,   end
        Bk = sopvar2opvar2d(Pmop.C{i,j});
        slot = nm{ri(i),ci(j)};
        Pop.(slot) = Bk.(slot);
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

function I = dom_matrix(Pmop,vars_xy)
% opvar2d's I is 2x2, row k the domain of vars_xy{k}. A direction the
% container never used has no recorded domain, so it defaults to [0,1].
I = [0 1;0 1];
for k = 1:min(2,numel(vars_xy))
    idx = find(strcmp(Pmop.vars,vars_xy{k}),1);
    if ~isempty(idx),   I(k,:) = Pmop.dom(idx,:);  end
end
end
