function tl = exec_tools_sop(PIE,st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TL = EXEC_TOOLS_SOP(PIE,ST) returns the constructions that the settings
% ST define for the container executives on a 1-D or 2-D PIE: the positive
% operator and the stock slack of the settings, the margins, identity
% operators on the PIE's spaces, the degrees of a free operator, and the
% conversion of a solved container back to an 'opvar' or 'opvar2d' object.
%
% INPUT
% - PIE:  'pie_struct' object of a 1-D or 2-D PIE (initialized), with T a
%         'copvar' object;
% - st:   'struct' of LPI settings, completed by EXEC_SOP_SETTINGS:
%           1-D  eppos, eppos2, epneg, ddZ (each optional), and the degree
%                fields of POSLPIVAR_SETTINGS_SOP;
%           2-D  settings_2d, read by SETTINGS_2D_SOP and
%                POSLPIVAR_SETTINGS_2D_SOP.
%
% OUTPUT
% - tl:   'struct' with fields
%           nv             double, the spatial dimension, 1 or 2;
%           S              'struct', the 2-D settings of SETTINGS_2D_SOP
%                          ([] in 1-D);
%           dom            the domain argument of LPIVAR_CDOPVAR and
%                          POSCOPVAR: 1x2 double PIE.dom in 1-D; in 2-D a
%                          'struct' with fields vars (1x2 cellstr) and dom
%                          (2x2 double);
%           eppos          double column, one margin per space of PIE.T:
%                          1-D [eppos; eppos2] (default [1e-4; 1e-6]), 2-D
%                          the four entries of settings_2d.eppos;
%           epneg          double, the decay rate (default 0);
%           eppos_default  logical, true when the settings carry no eppos;
%           storage        function handle, [prog,P] = storage(prog,X,side):
%                          declares the positive operator of the settings
%                          on the SIDE ('out' or 'in') spaces of the
%                          'copvar' X; P is a 'cdopvar' object;
%           slack          function handle, [prog,N] = slack(prog,X):
%                          declares the stock slack of an inequality on the
%                          spaces of X;
%           eye            function handle, E = eye(side,e): the 'copvar'
%                          operator diag(e_k I) on the output (SIDE 1) or
%                          input (SIDE 2) spaces of PIE.T, e a scalar or one
%                          weight per space;
%           ident          function handle, I = ident(Bop,side): the
%                          identity 'copvar' on the output (SIDE 1) or input
%                          (SIDE 2) space of the 'copvar' Bop;
%           toop           function handle, Op = toop(C): the 'copvar' C as
%                          an 'opvar' (1-D) or 'opvar2d' (2-D) object;
%           zdeg           the degrees of a free operator for
%                          LPIVAR_CDOPVAR: 1-D the settings' ddZ (default 1);
%                          2-D a 'struct' with the largest entry of
%                          settings_2d.Zop_deg on every role.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - exec_tools_sop
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
% Initial coding MMP, 10/08/2026
% MMP, 10/28/2026: Moved from executives_sopvar/private/ to
%                utility_functions/ with the consolidation of the
%                executives; header to the header protocol, the licence
%                block added, step comments in the body; code unchanged.
% DJ (CC), 10/29/2026: Support the 'copvar' operators of the new pie_struct:
%                tl.ident and eye_state read the spaces of PIE.T and Bop
%                from the 'copvar' fields, as PIE.T.dim and Bop.dim no longer
%                exist; eye_state gives one block for each block of PIE.T,
%                and tl.ident one for each space of Bop;

% Dimension and spatial variables of the PIE
nv = PIE.dim;
vars = PIE.vars;
vnames = reshape(vars(:,1).varname,1,[]);
tl = struct('nv',nv,'S',[]);
if nv==1
    % 1-D: margins from eppos, eppos2, epneg; storage and slack from the
    % 1-D degree fields; conversion by copvar2opvar; Z degrees from ddZ
    tl.dom = PIE.dom;
    tl.eppos = [getf(st,'eppos',1e-4); getf(st,'eppos2',1e-6)];
    tl.epneg = getf(st,'epneg',0);
    tl.eppos_default = ~isfield(st,'eppos') || isempty(st.eppos);
    tl.storage = @(prog,X,side) poslpivar_settings_sop(prog,X,st,'lf',side,PIE.dom);
    tl.slack   = @(prog,X) poslpivar_settings_sop(prog,X,st,'slack','out',PIE.dom);
    tl.toop    = @(C) copvar2opvar(C);
    tl.zdeg    = getf(st,'ddZ',1);
elseif nv==2
    % 2-D: everything from settings_2d; conversion by copvar2opvar2d; one
    % degree cap for Z
    S = settings_2d_sop(st);
    tl.S = S;
    tl.dom = struct('vars',{vnames},'dom',PIE.dom);
    tl.eppos = S.eppos;     tl.epneg = S.epneg;     tl.eppos_default = S.eppos_default;
    tl.storage = @(prog,X,side) poslpivar_settings_2d_sop(prog,X,st,'lf',side,PIE);
    tl.slack   = @(prog,X) poslpivar_settings_2d_sop(prog,X,st,'slack','out',PIE);
    tl.toop    = @(C) copvar2opvar2d(C,vnames);
    tl.zdeg    = zcaps(S.Zop_deg);
else
    error('exec_tools_sop:dim','The stock executives exist for 1-D and 2-D PIEs only.');
end
% Identity operators: weighted on the state spaces, plain on an input or
% output space of an operator
tl.eye   = @(side,e) eye_state(PIE,side,e);
%   tl.ident = @(Bop,side) op2copvar_sop(mat2opvar(eye(size(Bop,side)),Bop.dim(:,side),vars,PIE.dom));   % DJ (CC), 10/29/2026 (was)
% The identity on the output (side 1) or input (side 2) spaces of a 'copvar'   % DJ (CC), 10/29/2026
% Bop, with one block for each of its spaces; an empty container if Bop has    % DJ (CC), 10/29/2026
% no space on that side (an empty channel), as the former opvar construction.  % DJ (CC), 10/29/2026
sd = {'out','in'};                                                          % DJ (CC), 10/29/2026
spc = @(B,s) arrayfun(@(k) B.vars(B.(['space_',sd{s}])(k,:)),1:numel(B.(['dim_',sd{s}])),'UniformOutput',false);   % DJ (CC), 10/29/2026
mkI = {@(B,s) copvar(), @(B,s) eye_copvar_sop(B.(['dim_',sd{s}]),spc(B,s),struct('vars',{B.vars},'dom',B.dom))};   % DJ (CC), 10/29/2026
tl.ident = @(Bop,side) mkI{1+~isempty(Bop.(['dim_',sd{side}]))}(Bop,side);  % DJ (CC), 10/29/2026
end


function Im = eye_state(PIE,side,e)
% IM = EYE_STATE(PIE,SIDE,E) returns diag(e_k I_{n_k}) on the state spaces
% of the PIE as a 'copvar' object.
% INPUT  PIE:   'pie_struct' object, 1-D or 2-D;
%        side:  1 (the output spaces of PIE.T, the PDE state) or 2 (its
%               input spaces, the fundamental state);
%        e:     double, one weight per space of that side (2 in 1-D, 4 in
%               2-D), or a scalar for every space.
% OUTPUT Im:    'copvar' object.
%   n = PIE.T.dim(:,side);                                                  % DJ (CC), 10/29/2026 (was)
% The type of space (1: R, 2: L2[s1], 3: L2[s2], 4: L2[s1,s2]) of each block   % DJ (CC), 10/29/2026
% of PIE.T, from the registry of its variables                              % DJ (CC), 10/29/2026
T = PIE.T;                                                                  % DJ (CC), 10/29/2026
if side==1                                                                  % DJ (CC), 10/29/2026
    dims = T.dim_out;  
    mask = T.space_out;  
else
    dims = T.dim_in;  
    mask = T.space_in;  
end
[~,loc] = ismember(T.vars,PIE.space.vars);                                  % DJ (CC), 10/29/2026
kind = ones(numel(dims),1);                                                 % DJ (CC), 10/29/2026
for k = 1:numel(loc),   kind = kind + loc(k)*mask(:,k);   end              % DJ (CC), 10/29/2026
%   if isscalar(e),     e = e*ones(numel(n),1);     end                     % DJ (CC), 10/29/2026 (was)
%   if numel(e)~=numel(n)                                                   % DJ (CC), 10/29/2026 (was)
if isscalar(e),     e = e*ones(2^PIE.dim,1);     end                        % DJ (CC), 10/29/2026
if numel(e)~=2^PIE.dim                                                      % DJ (CC), 10/29/2026
    error('exec_tools_sop:e','E should have one weight per space (%d) or be a scalar.',2^PIE.dim);   % DJ (CC), 10/29/2026
end
% The block-diagonal matrix of the weights, as a container with the blocks  % DJ (CC), 10/29/2026
% of PIE.T                                                                  % DJ (CC), 10/29/2026
blocks = cell(1,numel(dims));                                               % DJ (CC), 10/29/2026
spaces = cell(1,numel(dims));                                               % DJ (CC), 10/29/2026
for k = 1:numel(dims),  blocks{k} = e(kind(k))*eye(dims(k));  spaces{k} = T.vars(mask(k,:));   end   % DJ (CC), 10/29/2026
Im = mat2copvar_sop(blkdiag(blocks{:}),dims,spaces,struct('vars',{T.vars},'dom',T.dom));   % DJ (CC), 10/29/2026
end




function v = getf(s,f,d)
% V = GETF(S,F,D) returns the field F of the struct S, or D when S has no
% nonempty field F.
% INPUT  s: 'struct'; f: char, the field name; d: the default value.
% OUTPUT v: the field's value or D.
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)),  v = s.(f);  else,   v = d;  end
end


function d = zcaps(Z)
% D = ZCAPS(Z) returns one degree cap for every role of LPIVAR_CDOPVAR: the
% largest entry of the stock 2-D free-variable degrees Z.
% INPUT  Z: 'struct' with fields dx, dy, d2, each a nested cell of double
%           arrays (settings_2d.Zop_deg).
% OUTPUT d: 'struct' with fields mult (double), int (1x2 double), out and
%           in (double), each the cap.
m = max([cellmax(Z.dx), cellmax(Z.dy), cellmax(Z.d2), 0]);
d = struct('mult',m,'int',[m m],'out',m,'in',m);
end


function m = cellmax(c)
% M = CELLMAX(C) returns the largest entry of a nested cell of double
% arrays (0 when it holds none).
% INPUT  c: cell or double array.
% OUTPUT m: double.
m = 0;
if iscell(c)
    for k = 1:numel(c),     m = max(m,cellmax(c{k}));   end
elseif isnumeric(c) && ~isempty(c)
    m = max(c(:));
end
end
