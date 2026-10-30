function tl = exec_tools_sop(PIE,st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TL = EXEC_TOOLS_SOP(PIE,ST) returns the constructions that the settings
% ST define for the container executives on a PIE in any number of spatial
% variables: the positive operator and the stock slack of the settings, the
% margins, identity operators on the PIE's spaces, the degrees of a free
% operator, and the conversion of a solved container back to an 'opvar'
% or 'opvar2d' object (1-D, 2-D).
%
% INPUT
% - PIE:  'pie_struct' object as EXEC_PIE_SOP returns it: a 1-D or 2-D PIE
%         with T an 'opvar' (1-D) or 'opvar2d' (2-D) object, or a container
%         PIE (T a 'copvar' object) in any number nv >= 1 of variables;
% - st:   'struct' of LPI settings, completed by EXEC_SOP_SETTINGS:
%           1-D  eppos, eppos2, epneg, ddZ (each optional), and the degree
%                fields of POSLPIVAR_SETTINGS_SOP;
%           2-D  settings_2d, read by SETTINGS_2D_SOP and
%                POSLPIVAR_SETTINGS_2D_SOP;
%           N-D  (nv >= 3) settings_nd (optional), a 'struct' with the
%                optional fields
%                  LF_int, LF_mult  the positive operator on a space in
%                       d >= 1 variables: the multiplier basis operator
%                       (alpha = 1 in every direction) of degree LF_mult
%                       (default 0) and the full-integral one (alpha = 4,
%                       separable in every direction) with caps LF_int in
%                       each variable and 2*LF_int in total (default 1),
%                       the rule of the stock light 2-D storage; the
%                       identity basis on R^n;
%                  LF_deg, LF_opts  instead, a POSCOPVAR degree
%                       specification and its options;
%                  LF_psatz  the term codes (default 0);
%                  eppos  the margin on every space, or 2^nv margins as
%                       for EYE below (default 1e-2);
%                  epneg  (default the field epneg of ST, else 0);
%                  Zop_deg  the degree of a free operator (default 1).
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
%                          one weight per row of PIE.T.dim (opvar PIE), or
%                          2^nv weights (container PIE), the weight of a
%                          space being entry 1 + sum_k 2^(k-1) over the
%                          positions k in PIE.vars of its variables, which
%                          is the row order of PIE.T.dim in 1-D and 2-D;
%           ident          function handle, I = ident(Bop,side): the
%                          identity 'copvar' on the output (SIDE 1) or input
%                          (SIDE 2) space of the 'opvar', 'opvar2d' or
%                          'copvar' Bop;
%           toop           function handle, Op = toop(C): the 'copvar' C as
%                          an 'opvar' (1-D) or 'opvar2d' (2-D) object; C
%                          itself in N-D;
%           zdeg           the degrees of a free operator for
%                          LPIVAR_CDOPVAR: 1-D the settings' ddZ (default 1);
%                          2-D a 'struct' with the largest entry of
%                          settings_2d.Zop_deg on every role; N-D the same
%                          with settings_nd.Zop_deg;
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
% MMP, 10/29/2026: PIEs in any number of spatial variables (maintainer: the
%                executives must support n variables, at least 3-D), and
%                container PIEs (T a 'copvar') in every dimension.
%                - N-D branch: no stock settings exist in 3-D or higher, so
%                  the positive operator is a POSCOPVAR of settings_nd
%                  (local SETTINGS_ND, POS_ND), the stock slack is an error
%                  (the 'like' slack of LPI_INEQ_SOP is N-D), and a solved
%                  container is returned as it is.
%                - eye and ident build on the spaces of a container (local
%                  EYE_CONTAINER, IDENT_OP; MAT2COPVAR_SOP, EYE_COPVAR_SOP),
%                  since an N-D PIE has no opvar dims; an opvar PIE takes
%                  the former code.
%                The 1-D and 2-D branches are unchanged.
%                The default N-D storage has the structure of the stock
%                light 2-D storage (a multiplier and one separable full
%                integral per space): with all 3^nv basis operators at
%                degree 1, the 3-D heat stability SDP (heatNd, light) had m
%                145,123 (storage Gram 216, slack Gram 1728, build 439 s),
%                past what MOSEK solves on this workstation.
% MMP, 10/29/2026: Field capfor (local CAP_FOR): the per-direction bounds of
%                settings.sop.cap_delta for LPI_INEQ_SOP's opts.cap, read
%                from T (GET_LIFT_DEGS info.Jmax) and laid over each
%                target's variables; [] when cap_delta is empty, so every
%                program is unchanged by default.
% MMP, 10/29/2026: Field sop (local SOP_OF_DIM): the per-dimension vectors of
%                the settings_sopvar tiers ([1-D 2-D 3-D], since 10/29/2026)
%                reduced to this dimension's entry, so the executives keep
%                reading scalars; a NaN cap_delta means no bound.

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
    tl.toop    = @(C) C;                                                    % MMP, 10/29/2026
    z = S.Zop_deg;  tl.zdeg = struct('mult',z,'int',[z z],'out',z,'in',z);  % MMP, 10/29/2026
end
% Identity operators: weighted on the state spaces, plain on an input or
% output space of an operator
% tl.eye   = @(side,e) eye_state(PIE,side,e);                               % MMP, 10/29/2026 (was)
% tl.ident = @(Bop,side) op2copvar_sop(mat2opvar(eye(size(Bop,side)),Bop.dim(:,side),vars,PIE.dom)); % MMP, 10/29/2026 (was)
cdom = PIE.dom;     if nv>=2,   cdom = tl.dom;  end % the domain of a container % MMP, 10/29/2026
tl.eye   = @(side,e) eye_state(PIE,side,e,cdom);                            % MMP, 10/29/2026
tl.ident = @(Bop,side) ident_op(Bop,side,vars,PIE.dom,cdom);                % MMP, 10/29/2026
% The container options of this dimension: every per-dimension vector       % MMP, 10/29/2026
% reduced to its entry min(nv,end) (local SOP_OF_DIM), for the executives   % MMP, 10/29/2026
tl.sop = [];    if isfield(st,'sop'),   tl.sop = sop_of_dim(st.sop,nv);     end % MMP, 10/29/2026
% Per-direction bounds on the slack bases (settings.sop.cap_delta): Jmax    % MMP, 10/29/2026
% of T plus the delta, laid over each target's variables; none by default   % MMP, 10/29/2026
% and where the delta of this dimension is NaN                              % MMP, 10/29/2026
tl.capfor = @(X) [];                                                        % MMP, 10/29/2026
cd = [];    if isstruct(tl.sop) && isfield(tl.sop,'cap_delta'),     cd = tl.sop.cap_delta;  end % MMP, 10/29/2026
if ~isempty(cd) && ~isnan(cd)                                               % MMP, 10/29/2026
    Tc = op2copvar_sop(PIE.T);      [~,~,ti] = get_lift_degs(Tc,struct());  % MMP, 10/29/2026
    ct = max(ti.Jmax,0) + cd;                                               % MMP, 10/29/2026
    tl.capfor = @(X) cap_for(X,Tc.vars,ct);                                 % MMP, 10/29/2026
end                                                                         % MMP, 10/29/2026
end


% function Im = eye_state(PIE,side,e)                                       % MMP, 10/29/2026 (was)
function Im = eye_state(PIE,side,e,dom)                                     % MMP, 10/29/2026
% IM = EYE_STATE(PIE,SIDE,E,DOM) returns diag(e_k I_{n_k}) on the state
% spaces of the PIE as a 'copvar' object.
% INPUT  PIE:   'pie_struct' object: an opvar PIE (1-D or 2-D) or a
%               container PIE (T a 'copvar', any dimension);
%        side:  1 (the output spaces of PIE.T, the PDE state) or 2 (its
%               input spaces, the fundamental state);
%        e:     double, a scalar for every space, or for an opvar PIE one
%               weight per row of PIE.T.dim (2 in 1-D, 4 in 2-D), for a
%               container PIE 2^nv weights (EYE_CONTAINER);
%        dom:   the domain argument of MAT2COPVAR_SOP (container PIE).
% OUTPUT Im:    'copvar' object.
% A container PIE: on the spaces of its T                                   % MMP, 10/29/2026
if isa(PIE.T,'copvar'),     Im = eye_container(PIE,side,e,dom);     return,     end % MMP, 10/29/2026
n = PIE.T.dim(:,side);
if isscalar(e),     e = e*ones(numel(n),1);     end
if numel(e)~=numel(n)
    error('exec_tools_sop:e','E should have one weight per space (%d) or be a scalar.',numel(n));
end
% The block-diagonal matrix of the weights, as an opvar / opvar2d, then a
% container
blocks = cell(1,numel(n));
for k = 1:numel(n),     blocks{k} = e(k)*eye(n(k));    end
Imat = blkdiag(blocks{:});
if PIE.dim==1
    Pop = mat2opvar(Imat,n,PIE.vars,PIE.dom);
else
    Pop = opvar2d(Imat,[n,n],PIE.dom,PIE.vars);
end
Im = op2copvar_sop(Pop);
end


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


function I = ident_op(Bop,side,vars,pdom,cdom)                              % MMP, 10/29/2026
% I = IDENT_OP(BOP,SIDE,VARS,PDOM,CDOM) returns the identity 'copvar' on the
% output (SIDE 1) or input (SIDE 2) spaces of BOP.
% INPUT  Bop:   'opvar', 'opvar2d', 'copvar' or 'cdopvar' object;
%        side:  1 or 2;
%        vars:  PIE.vars, and pdom: PIE.dom (an opvar Bop);
%        cdom:  the domain argument of EYE_COPVAR_SOP (a container Bop).
% OUTPUT I:     'copvar' object.
if isa(Bop,'copvar') || isa(Bop,'cdopvar')
    sd = {'out','in'};
    [sp,dm] = copvar_space_list(Bop,sd{side});
    I = eye_copvar_sop(dm,sp,cdom);
else
    I = op2copvar_sop(mat2opvar(eye(size(Bop,side)),Bop.dim(:,side),vars,pdom));
end
end


function S = settings_nd(st)                                                % MMP, 10/29/2026
% S = SETTINGS_ND(ST) returns the N-D settings of the container executives:
% the fields of ST.settings_nd, each absent or empty one at its default.
% INPUT  st: 'struct' of LPI settings, with the optional field settings_nd
%            ('struct': LF_int, LF_mult, LF_deg, LF_opts, LF_psatz, eppos,
%            epneg, Zop_deg).
% OUTPUT S:  'struct' with fields LF_int (default 1), LF_mult (0), LF_deg
%            ([]: the rule of LF_int, LF_mult), LF_opts (struct()),
%            LF_psatz (0), eppos (1e-2), epneg (ST.epneg, else 0), Zop_deg
%            (1) and eppos_default (logical, true when ST.settings_nd has no
%            eppos).
S = struct('LF_int',1,'LF_mult',0,'LF_deg',[],'LF_opts',struct(),'LF_psatz',0, ...
           'eppos',1e-2,'epneg',getf(st,'epneg',0),'Zop_deg',1);
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
                  struct('int',S.LF_int,'mult',S.LF_int,'joint',2*S.LF_int)};
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


function s = sop_of_dim(s,nv)
% S = SOP_OF_DIM(S,NV) returns the container options S with each
% per-dimension field reduced to its entry for NV variables, the last
% entry serving every higher NV: dw_feas, dw_P, dw_Q, max_raise,
% max_raise_feas, cap_delta (eppos_stab, one entry per space, is kept).
% INPUT  s:  'struct', settings.sop (EXEC_SOP_SETTINGS);
%        nv: double, the spatial dimension.
% OUTPUT s:  'struct', S with those fields scalar.
f = {'dw_feas','dw_P','dw_Q','max_raise','max_raise_feas','cap_delta'};
for k = 1:numel(f)
    if isfield(s,f{k}) && isnumeric(s.(f{k})) && numel(s.(f{k}))>1
        v = s.(f{k});   s.(f{k}) = v(min(nv,end));
    end
end
end
% % % END MMP, 10/29/2026


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
