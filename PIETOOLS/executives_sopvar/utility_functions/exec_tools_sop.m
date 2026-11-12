function tl = exec_tools_sop(PIE,st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TL = EXEC_TOOLS_SOP(PIE,ST) returns the constructions that the settings
% ST define for the container executives on a PIE in any number of spatial
% variables: the positive operator and the stock slack of the settings, the
% margins, identity operators on the PIE's spaces, the degrees of a free
% operator, the container options of the PIE's dimension and the
% per-direction bounds of the slack bases.
%
% INPUT
% - PIE:  'pie_struct' object as EXEC_PIE_SOP returns it: each operator a
%         'copvar' object or 0 (absent), in any number nv >= 1 of
%         variables;
% - st:   'struct' of LPI settings, completed by EXEC_SOP_SETTINGS:
%           1-D  eppos, eppos2, epneg (each optional), lpi.var.deg{1}, the
%                degrees of a free operator, and lpi.pos, the positive
%                operator of POSLPIVAR_SETTINGS_SOP;
%           2-D  settings_2d, read by SETTINGS_2D_SOP and
%                POSLPIVAR_SETTINGS_2D_SOP;
%           N-D  (nv >= 3) lpi.pos and lpi.var, entry 3: the positive
%                operator on a space in d >= 1 variables, the multiplier
%                basis operator (alpha = 1 in every direction) of degree
%                mult and the full-integral one (alpha = 4, separable in
%                every direction) with caps int in each variable and joint
%                in total, the rule of the stock light 2-D storage, the
%                identity basis on R^n; lpi.var.deg{3}, the degree of a free
%                operator; and settings_nd (optional), a 'struct' with the
%                optional fields
%                  LF_deg, LF_opts  instead, a POSCOPVAR degree
%                       specification and its options;
%                  LF_psatz  the term codes (default 0);
%                  eppos  the margin on every space, or 2^nv margins as
%                       for EYE below (default 1e-2);
%                  epneg  (default the field epneg of ST, else 0).
%
% OUTPUT
% - tl:   'struct' with fields
%           nv             double, the spatial dimension;
%           S              'struct', the 2-D settings of SETTINGS_2D_SOP,
%                          or in N-D the N-D settings above with field
%                          eppos_default ([] in 1-D);
%           dom            the domain argument of LPIVAR_CDOPVAR and
%                          POSCOPVAR: 1x2 double PIE.dom in 1-D; in 2-D
%                          and N-D a 'struct' with fields vars (1 x nv
%                          cellstr, in the order of PIE.vars) and dom
%                          (nv x 2 double);
%           eppos          double column, the margins of the spaces of
%                          PIE.T: 1-D [eppos; eppos2] (default [1e-4;
%                          1e-6]), 2-D the four entries of
%                          settings_2d.eppos, N-D settings_nd.eppos;
%           epneg          double, the decay rate (default 0);
%           eppos_default  logical, true when the settings carry no eppos;
%           storage        function handle, [prog,P] = storage(prog,X,side):
%                          declares the positive operator of the settings
%                          on the SIDE ('out' or 'in') spaces of the
%                          'copvar' X; P is a 'cdopvar' object;
%           slack          function handle, [prog,N] = slack(prog,X):
%                          declares the stock slack of an inequality on the
%                          spaces of X; in N-D an error (no stock settings
%                          exist there);
%           eye            function handle, E = eye(side,e): the 'copvar'
%                          operator diag(e_k I) on the output (SIDE 1) or
%                          input (SIDE 2) spaces of PIE.T; e a scalar, or
%                          2^nv weights, the weight of a space being entry
%                          1 + sum_k 2^(k-1) over the positions k in
%                          PIE.vars of its variables;
%           ident          function handle, I = ident(Bop,side): the
%                          identity 'copvar' on the output (SIDE 1) or input
%                          (SIDE 2) space of the 'copvar' or 'cdopvar' Bop;
%           zdeg           the degrees of a free operator for
%                          LPIVAR_CDOPVAR: 1-D lpi.var.deg{1}; 2-D a
%                          'struct' with the largest entry of
%                          settings_2d.Zop_deg on every role; N-D the same
%                          with lpi.var.deg{3};
%           sop            'struct', settings.sop with each per-dimension
%                          vector (dw_feas, dw_P, dw_Q, max_raise,
%                          max_raise_feas, cap_delta) reduced to its entry
%                          min(nv,end): the options the executives read;
%           capfor         function handle, cap = capfor(X): the opts.cap of
%                          LPI_INEQ_SOP for a slack of the target X, the
%                          bounds Jmax(T) + cap_delta of this dimension
%                          (Jmax of GET_LIFT_DEGS on PIE.T) on the variables
%                          of X; [] (no bound) when cap_delta is empty or
%                          NaN, X has no variable, or X has one T lacks.
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
% MMP, 10/29/2026: PIEs in any number of spatial variables and container
%                PIEs (T a 'copvar') in every dimension. N-D branch: no stock
%                settings exist in 3-D or higher, so the positive operator
%                is a POSCOPVAR of settings_nd (local SETTINGS_ND, POS_ND),
%                whose default has the structure of the stock light 2-D
%                storage, a multiplier and one separable full integral per
%                space (with all 3^nv basis operators at degree 1 the 3-D
%                heat stability SDP was past what MOSEK solves here); the
%                stock slack is an error, and a solved container is returned
%                as it is. eye and ident build on the spaces of a container
%                (local EYE_CONTAINER, IDENT_OP). Field capfor (local
%                CAP_FOR): the per-direction bounds of settings.sop.cap_delta
%                for opts.cap of LPI_INEQ_SOP, read from T (GET_LIFT_DEGS
%                info.Jmax) and laid over each target's variables, [] when
%                cap_delta is empty. Field sop (local SOP_OF_DIM): the
%                per-dimension vectors of the settings_sopvar tiers reduced
%                to this dimension's entry; a NaN cap_delta means no bound.
%                The 1-D and 2-D branches are unchanged.
% DJ (CC), 10/29/2026: Support the 'copvar' operators of the new pie_struct:
%                tl.ident and eye_state read the spaces of PIE.T and Bop
%                from the 'copvar' fields, as PIE.T.dim and Bop.dim no longer
%                exist; eye_state gives one block for each block of PIE.T,
%                and tl.ident one for each space of Bop;
% MMP, 10/30/2026: Merge of origin/ndopvar: the DJ tl.ident and EYE_STATE
%                bodies (origin commit 3df33dcc) and the container branch of
%                10/29/2026 both read the spaces from the 'copvar' fields;
%                the container branch is kept, the DJ lines commented out
%                where they stood. EYE_CONTAINER indexes the weights by
%                1 + sum_k 2^(k-1) over the variables of a space, equal to
%                the 1 + sum_k loc(k) of DJ in 1-D and 2-D (in 3 or more
%                variables that index gives {s3} and {s1,s2} the same
%                weight); IDENT_OP returns the 0 x 0 container of
%                EYE_COPVAR_SOP for an empty channel, where DJ returned
%                copvar().
% MMP, 11/04/2026: The opvar-PIE code is commented out, being unreachable
%                (EXEC_PIE_SOP returns a PIE whose operators are all 'copvar'
%                objects): EYE_STATE calls EYE_CONTAINER, IDENT_OP errors on
%                a class other than 'copvar' or 'cdopvar' (its arguments VARS
%                and PDOM dropped), capfor reads PIE.T without OP2COPVAR_SOP.
%                Supersedes the opvar-PIE parts of the entries of 10/29/2026
%                and 10/30/2026.
% MMP, 11/08/2026: Comment pass (maintainer, 11/08/2026): log condensed; the
%                help of toop and the 1-D, 2-D branch notes say that
%                copvar2opvar and copvar2opvar2d are retired (11/08/2026,
%                cleaning/private), so the handle errors when called (unused
%                since 11/04/2026); three commented-out notes removed; code
%                unchanged.
% MMP, 11/09/2026: Inlining pass (maintainer, 11/09/2026: locals that are
%                very short and called once go into the body): EYE_STATE
%                (one call of EYE_CONTAINER since 11/04/2026, with its
%                disabled opvar body of 10/29/2026) into the handle tl.eye,
%                SOP_OF_DIM into the construction of tl.sop, ZCAPS into the
%                2-D branch, each with its header as the comment of the
%                block. CAP_FOR and IDENT_OP stay: their handles need a
%                function body.
% MMP, 11/09/2026: tl.toop removed (maintainer, 11/09/2026: "why does
%                exec_tools need a converter?"): the handle converted a
%                container result back to an opvar for output; no executive
%                has called it since they return containers (11/04/2026),
%                and its targets copvar2opvar, copvar2opvar2d are retired.
% MMP, 11/12/2026: Native settings (maintainer, 11/12/2026: phase 2 of the
%                knob reduction): Z in 1-D from lpi.var.deg{1}, the former
%                ddZ; in N-D, LF_int, LF_mult and Zop_deg from lpi.pos and
%                lpi.var, entry 3, the former settings_nd fields, with the
%                total cap joint of lpi.pos in place of 2*LF_int (2*int in
%                every tier, so the same P). settings_nd keeps LF_deg,
%                LF_opts, LF_psatz, eppos and epneg.

% Dimension and spatial variables of the PIE
nv = PIE.dim;
vars = PIE.vars;
vnames = reshape(vars(:,1).varname,1,[]);
tl = struct('nv',nv,'S',[]);
if nv==1
%   % 1-D: margins from eppos, eppos2, epneg; storage and slack from the    % MMP, 11/12/2026 (was); MMP, 11/09/2026
%   % 1-D degree fields; Z degrees from ddZ                                 % MMP, 11/12/2026 (was); MMP, 11/09/2026
    % 1-D: margins from eppos, eppos2, epneg; storage from lpi.pos, stock   % MMP, 11/12/2026
    % slack from sop.stock (POSLPIVAR_SETTINGS_SOP); Z degrees from lpi.var % MMP, 11/12/2026
    tl.dom = PIE.dom;
    tl.eppos = [getf(st,'eppos',1e-4); getf(st,'eppos2',1e-6)];
    tl.epneg = getf(st,'epneg',0);
    tl.eppos_default = ~isfield(st,'eppos') || isempty(st.eppos);
    tl.storage = @(prog,X,side) poslpivar_settings_sop(prog,X,st,'lf',side,PIE.dom);
    tl.slack   = @(prog,X) poslpivar_settings_sop(prog,X,st,'slack','out',PIE.dom);
%   tl.toop    = @(C) copvar2opvar(C);                                      % MMP, 11/09/2026 (was)
%   tl.zdeg    = getf(st,'ddZ',1);                                          % MMP, 11/12/2026 (was)
    tl.zdeg    = st.lpi.var.deg{1};     % the former ddZ                    % MMP, 11/12/2026
elseif nv==2
    % 2-D: everything from settings_2d; one degree cap for Z                % MMP, 11/09/2026
    S = settings_2d_sop(st);
    tl.S = S;
    tl.dom = struct('vars',{vnames},'dom',PIE.dom);
    tl.eppos = S.eppos;     tl.epneg = S.epneg;     tl.eppos_default = S.eppos_default;
    tl.storage = @(prog,X,side) poslpivar_settings_2d_sop(prog,X,st,'lf',side,PIE);
    tl.slack   = @(prog,X) poslpivar_settings_2d_sop(prog,X,st,'slack','out',PIE);
%   tl.toop    = @(C) copvar2opvar2d(C,vnames);                             % MMP, 11/09/2026 (was)
%   tl.zdeg    = zcaps(S.Zop_deg);                                          % MMP, 11/09/2026 (was)
    % One degree cap for every role of LPIVAR_CDOPVAR: the largest entry of % MMP, 11/09/2026
    % the stock 2-D free-variable degrees (the former local ZCAPS)          % MMP, 11/09/2026
    m = max([cellmax(S.Zop_deg.dx), cellmax(S.Zop_deg.dy), cellmax(S.Zop_deg.d2), 0]); % MMP, 11/09/2026
    tl.zdeg    = struct('mult',m,'int',[m m],'out',m,'in',m);               % MMP, 11/09/2026
% else                                                                      % MMP, 10/29/2026 (was)
%     error('exec_tools_sop:dim','The stock executives exist for 1-D and 2-D PIEs only.'); % MMP, 10/29/2026 (was)
else                                                                        % MMP, 10/29/2026
    % N-D (a container PIE in 3 or more variables): the storage of          % MMP, 10/29/2026
    % settings_nd by POSCOPVAR; no stock slack; containers kept             % MMP, 10/29/2026
    S = settings_nd(st);                                                    % MMP, 10/29/2026
    tl.S = S;                                                               % MMP, 10/29/2026
    tl.dom = struct('vars',{vnames},'dom',PIE.dom);                         % MMP, 10/29/2026
    tl.eppos = S.eppos(:);  tl.epneg = S.epneg;     tl.eppos_default = S.eppos_default; % MMP, 10/29/2026
    tl.storage = @(prog,X,side) pos_nd(prog,X,side,S,tl.dom);               % MMP, 10/29/2026
    tl.slack   = @(prog,X) error('exec_tools_sop:slack',['The stock slack reads the stock 1-D ' ...
                     'and 2-D settings; for %d variables set settings.sop.slack = ''like''.'],nv); % MMP, 10/29/2026
%   tl.toop    = @(C) C;                                                    % MMP, 11/09/2026 (was); MMP, 10/29/2026
    z = S.Zop_deg;  tl.zdeg = struct('mult',z,'int',[z z],'out',z,'in',z);  % MMP, 10/29/2026
end
% Identity operators: weighted on the state spaces, plain on an input or
% output space of an operator
% tl.eye   = @(side,e) eye_state(PIE,side,e);                               % MMP, 10/29/2026 (was)
% tl.ident = @(Bop,side) op2copvar_sop(mat2opvar(eye(size(Bop,side)),Bop.dim(:,side),vars,PIE.dom)); % MMP, 10/29/2026 (was)
cdom = PIE.dom;     if nv>=2,   cdom = tl.dom;  end % the domain of a container % MMP, 10/29/2026
% tl.eye   = @(side,e) eye_state(PIE,side,e,cdom);                          % MMP, 11/09/2026 (was); MMP, 10/29/2026
% diag(e_k I) on the state spaces of the PIE, side 1 the PDE state, 2 the   % MMP, 11/09/2026
% fundamental state (EYE_CONTAINER directly: the former local EYE_STATE was % MMP, 11/09/2026
% that one call since 11/04/2026)                                           % MMP, 11/09/2026
tl.eye   = @(side,e) eye_container(PIE,side,e,cdom);                        % MMP, 11/09/2026
% tl.ident = @(Bop,side) ident_op(Bop,side,vars,PIE.dom,cdom);              % MMP, 11/04/2026 (was)
tl.ident = @(Bop,side) ident_op(Bop,side,cdom);    % copvar Bop only        % MMP, 11/04/2026
% The container options of this dimension, for the executives (the former   % MMP, 11/09/2026
% local SOP_OF_DIM): each per-dimension vector among dw_feas, dw_P, dw_Q,   % MMP, 11/09/2026
% max_raise, max_raise_feas and cap_delta reduced to its entry min(nv,end), % MMP, 11/09/2026
% the last entry serving every higher dimension; eppos_stab (one entry per  % MMP, 11/09/2026
% space) is kept                                                            % MMP, 11/09/2026
% tl.sop = [];    if isfield(st,'sop'),   tl.sop = sop_of_dim(st.sop,nv);     end % MMP, 11/09/2026 (was); MMP, 10/29/2026
tl.sop = [];                                                                % MMP, 11/09/2026
if isfield(st,'sop')                                                        % MMP, 11/09/2026
    tl.sop = st.sop;                                                        % MMP, 11/09/2026
    f = {'dw_feas','dw_P','dw_Q','max_raise','max_raise_feas','cap_delta'}; % MMP, 11/09/2026
    for k = 1:numel(f)                                                      % MMP, 11/09/2026
        if isfield(tl.sop,f{k}) && isnumeric(tl.sop.(f{k})) && numel(tl.sop.(f{k}))>1 % MMP, 11/09/2026
            v = tl.sop.(f{k});   tl.sop.(f{k}) = v(min(nv,end));            % MMP, 11/09/2026
        end                                                                 % MMP, 11/09/2026
    end                                                                     % MMP, 11/09/2026
end                                                                         % MMP, 11/09/2026
% Per-direction bounds on the slack bases (settings.sop.cap_delta): Jmax    % MMP, 10/29/2026
% of T plus the delta, laid over each target's variables; none by default   % MMP, 10/29/2026
% and where the delta of this dimension is NaN                              % MMP, 10/29/2026
tl.capfor = @(X) [];                                                        % MMP, 10/29/2026
cd = [];    if isstruct(tl.sop) && isfield(tl.sop,'cap_delta'),     cd = tl.sop.cap_delta;  end % MMP, 10/29/2026
if ~isempty(cd) && ~isnan(cd)                                               % MMP, 10/29/2026
%     Tc = op2copvar_sop(PIE.T);      [~,~,ti] = get_lift_degs(Tc,struct()); % MMP, 11/04/2026 (was)
    Tc = PIE.T;     [~,~,ti] = get_lift_degs(Tc,struct());  % a 'copvar'    % MMP, 11/04/2026
    ct = max(ti.Jmax,0) + cd;                                               % MMP, 10/29/2026
    tl.capfor = @(X) cap_for(X,Tc.vars,ct);                                 % MMP, 10/29/2026
end                                                                         % MMP, 10/29/2026
% The copvar tl.eye, tl.ident of DJ (CC), 10/29/2026, which the merge of    % MMP, 10/30/2026
% 10/30/2026 replaces by the container branch above (header entry):         % MMP, 10/30/2026
% tl.eye   = @(side,e) eye_state(PIE,side,e);                               % MMP, 10/30/2026 (was)
%   tl.ident = @(Bop,side) op2copvar_sop(mat2opvar(eye(size(Bop,side)),Bop.dim(:,side),vars,PIE.dom));   % DJ (CC), 10/29/2026 (was)
% % The identity on the output (side 1) or input (side 2) spaces of a 'copvar'   % DJ (CC), 10/29/2026
% % Bop, with one block for each of its spaces; an empty container if Bop has    % DJ (CC), 10/29/2026
% % no space on that side (an empty channel), as the former opvar construction.  % DJ (CC), 10/29/2026
% sd = {'out','in'};                                                        % MMP, 10/30/2026 (was)
% spc = @(B,s) arrayfun(@(k) B.vars(B.(['space_',sd{s}])(k,:)),1:numel(B.(['dim_',sd{s}])),'UniformOutput',false);   % MMP, 10/30/2026 (was)
% mkI = {@(B,s) copvar(), @(B,s) eye_copvar_sop(B.(['dim_',sd{s}]),spc(B,s),struct('vars',{B.vars},'dom',B.dom))};   % MMP, 10/30/2026 (was)
% tl.ident = @(Bop,side) mkI{1+~isempty(Bop.(['dim_',sd{side}]))}(Bop,side); % MMP, 10/30/2026 (was)
end


% function Im = eye_state(PIE,side,e)                                       % MMP, 10/29/2026 (was)
% % % BEGIN MMP, 10/29/2026: new local functions EYE_CONTAINER, IDENT_OP,
% % % SETTINGS_ND, POS_ND (container and N-D PIEs), CAP_FOR (the cap_delta
% % % bounds) and SOP_OF_DIM (the per-dimension options), header entries of
% % % 10/29/2026; everything from here to the END marker is new.
function Im = eye_container(PIE,side,e,dom)                                 % MMP, 10/29/2026
% IM = EYE_CONTAINER(PIE,SIDE,E,DOM) returns diag(e_k I_{n_k}) on the
% spaces of the container T of a container PIE as a 'copvar' object.
% INPUT  PIE:   'pie_struct' object, T a 'copvar' object;
%        side:  1 (the output spaces of PIE.T) or 2 (its input spaces);
%        e:     double, a scalar for every space, or 2^nv weights: the
%               weight of a space is entry 1 + sum_k 2^(k-1) over the
%               positions k in PIE.vars of its variables (R^n entry 1; in
%               2-D the rows R, L2[x], L2[y], L2[x,y] of an opvar2d dim);
%        dom:   the domain argument of MAT2COPVAR_SOP.
% OUTPUT Im:    'copvar' object, block k e_k I on space k.
% The spaces of that side, the variable names in the order of PIE.vars
sd = {'out','in'};
[sp,dm] = copvar_space_list(PIE.T,sd{side});
vn = reshape(PIE.vars(:,1).varname,1,[]);
% One weight per space: the scalar, or the entry of its variable subset
if isscalar(e)
    w = e*ones(numel(sp),1);
elseif numel(e)==2^numel(vn)
    w = zeros(numel(sp),1);
    for k = 1:numel(sp),    w(k) = e(1 + sum(2.^(find(ismember(vn,sp{k}))-1)));     end
else
    error('exec_tools_sop:e','E should be a scalar or have 2^nv = %d entries, one per subset of the variables.',2^numel(vn));
end
% The block-diagonal constant matrix as a container multiplier
blocks = cell(1,numel(sp));
for k = 1:numel(sp),    blocks{k} = w(k)*eye(dm(k));   end
Im = mat2copvar_sop(blkdiag(blocks{:}),dm,sp,dom);
end


% function I = ident_op(Bop,side,vars,pdom,cdom)                            % MMP, 11/04/2026 (was)
function I = ident_op(Bop,side,cdom)                                        % MMP, 11/04/2026
% I = IDENT_OP(BOP,SIDE,CDOM) returns the identity 'copvar' on the output
% (SIDE 1) or input (SIDE 2) spaces of BOP.
% INPUT  Bop:   'copvar' or 'cdopvar' object;
%        side:  1 or 2;
%        cdom:  the domain argument of EYE_COPVAR_SOP.
% OUTPUT I:     'copvar' object.
if isa(Bop,'copvar') || isa(Bop,'cdopvar')
    sd = {'out','in'};
    [sp,dm] = copvar_space_list(Bop,sd{side});
    I = eye_copvar_sop(dm,sp,cdom);
else
%     I = op2copvar_sop(mat2opvar(eye(size(Bop,side)),Bop.dim(:,side),vars,pdom)); % MMP, 11/04/2026 (was)
    % Every PIE operator is a 'copvar' (EXEC_PIE_SOP): anything else, such  % MMP, 11/04/2026
    % as the 0 of an absent channel, has no spaces to read                  % MMP, 11/04/2026
    error('exec_tools_sop:ident','BOP should be a ''copvar'' or ''cdopvar'' object.'); % MMP, 11/04/2026
end
end


function S = settings_nd(st)                                                % MMP, 10/29/2026
% S = SETTINGS_ND(ST) returns the N-D settings of the container executives:
% LF_int, LF_mult, LF_joint and Zop_deg from entry 3 of ST.lpi.pos and
% ST.lpi.var, the other fields from ST.settings_nd, each absent or empty one
% at its default.
% INPUT  st: 'struct' of LPI settings, with field lpi (pos, var) and the
%            optional field settings_nd ('struct': LF_deg, LF_opts,
%            LF_psatz, eppos, epneg).
% OUTPUT S:  'struct' with fields LF_int, LF_mult, LF_joint (lpi.pos int,
%            mult, joint), LF_deg ([]: the rule of LF_int, LF_mult), LF_opts
%            (struct()), LF_psatz (0), eppos (1e-2), epneg (ST.epneg, else
%            0), Zop_deg (lpi.var.deg{3}) and eppos_default (logical, true
%            when ST.settings_nd has no eppos).
% S = struct('LF_int',1,'LF_mult',0,'LF_deg',[],'LF_opts',struct(),'LF_psatz',0, ... % MMP, 11/12/2026 (was)
%            'eppos',1e-2,'epneg',getf(st,'epneg',0),'Zop_deg',1);          % MMP, 11/12/2026 (was)
ps = st.lpi.pos;    zd = st.lpi.var.deg;    k3 = @(v) v(min(3,end));        % MMP, 11/12/2026
S = struct('LF_int',k3(ps.int),'LF_mult',k3(ps.mult),'LF_joint',k3(ps.joint), ...
           'LF_deg',[],'LF_opts',struct(),'LF_psatz',0,'eppos',1e-2, ...
           'epneg',getf(st,'epneg',0),'Zop_deg',zd{min(3,end)});            % MMP, 11/12/2026
given = struct();
if isfield(st,'settings_nd') && isstruct(st.settings_nd),   given = st.settings_nd;     end
f = fieldnames(S);
for k = 1:numel(f)
    if isfield(given,f{k}) && ~isempty(given.(f{k})),   S.(f{k}) = given.(f{k});    end
end
S.eppos_default = ~(isfield(given,'eppos') && ~isempty(given.eppos));
end


function [prog,P] = pos_nd(prog,X,side,S,dom)                               % MMP, 10/29/2026
% [PROG,P] = POS_ND(PROG,X,SIDE,S,DOM) declares the positive operator of the
% N-D settings S on the SIDE spaces of X by POSCOPVAR: with S.LF_deg and
% S.LF_opts when S.LF_deg is given; else on each space in d >= 1 variables
% the multiplier basis operator of degree S.LF_mult and the separable
% full-integral one with caps S.LF_int per variable and 2*S.LF_int in total
% (the structure SETTINGS2POSSOPVAR gives the stock light 2-D storage), the
% identity basis on R^n; the term codes S.LF_psatz.
% INPUT  prog:  'struct', the LPI program;
%        X:     'copvar' or 'cdopvar' object;
%        side:  char, 'out' or 'in';
%        S:     'struct', as SETTINGS_ND returns;
%        dom:   'struct' with fields vars and dom, the domain of P.
% OUTPUT prog:  'struct', the program with the Gram matrices of P;
%        P:     'cdopvar' object.
[sp,dm] = copvar_space_list(X,side);
% An explicit POSCOPVAR specification
if ~isempty(S.LF_deg)
    co = S.LF_opts;     if ~isfield(co,'psatz'),    co.psatz = S.LF_psatz;  end
    [prog,P] = poscopvar(prog,dm,sp,dom,S.LF_deg,co);
    return
end
% The rule: per space, the identity basis on R^n (as POS2D in
% POSLPIVAR_SETTINGS_2D_SOP), else the multiplier (alpha 1...1) and the
% full-integral (alpha 4...4) basis operators, in the order of include
deg = cell(1,numel(sp));    incl = cell(1,numel(sp));
for k = 1:numel(sp)
    d = numel(sp{k});
    if d==0
        deg{k} = struct('int',zeros(1,numel(dom.vars)));   % one 0 per registry variable
    else
        deg{k} = {struct('int',S.LF_mult,'mult',0,'joint',S.LF_mult), ...
                  struct('int',S.LF_int,'mult',S.LF_int,'joint',S.LF_joint)}; % MMP, 11/12/2026
%                 struct('int',S.LF_int,'mult',S.LF_int,'joint',2*S.LF_int)}; % MMP, 11/12/2026 (was)
        incl{k} = [ones(1,d); 4*ones(1,d)];
    end
end
co = struct('psatz',S.LF_psatz,'sep',true);
co.include = incl;                  % assigned apart: a cell in struct() builds a struct array
[prog,P] = poscopvar(prog,dm,sp,dom,deg,co);
end


function c = cap_for(X,tv,tc)
% C = CAP_FOR(X,TV,TC) lays the bounds TC, one per variable TV of T, over
% the registry of X.
% INPUT  X:   'copvar' or 'cdopvar' object, the target of a slack;
%        tv:  1 x nv 'cellstr', the variables of T (sorted);
%        tc:  1 x nv double, the bounds.
% OUTPUT c:   1 x numel(X.vars) double, the bounds of the variables of X;
%             [] when X has no variable or one that TV lacks (no bound).
c = [];
xv = X.vars;
if isempty(xv),     return,     end
[ok,ix] = ismember(xv,tv);
if ~all(ok),    return,     end
c = reshape(tc(ix),1,[]);
end


% % % END MMP, 10/29/2026


function v = getf(s,f,d)
% V = GETF(S,F,D) returns the field F of the struct S, or D when S has no
% nonempty field F.
% INPUT  s: 'struct'; f: char, the field name; d: the default value.
% OUTPUT v: the field's value or D.
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)),  v = s.(f);  else,   v = d;  end
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
