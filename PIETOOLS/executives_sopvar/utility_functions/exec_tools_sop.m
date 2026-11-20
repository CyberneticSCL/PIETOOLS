function tl = exec_tools_sop(PIE,st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TL = EXEC_TOOLS_SOP(PIE,ST) returns the constructions that the LPI
% settings ST define on a PIE in any number of spatial variables: the
% positive operator and the 'stock' slack of the settings, the margins,
% identity operators on the spaces of the PIE, the degrees of a free
% operator, the options of the executives for the dimension of the PIE and
% the per-direction bounds of the slack bases.
%
% INPUT
% - PIE:  'pie_struct' object, each operator a 'copvar' object or 0
%         (absent), in any number nv >= 1 of variables;
% - st:   'struct' of LPI settings:
%           1-D  eppos, eppos2, epneg (each optional), lpi.var.deg{1}, the
%                degrees of a free operator, and lpi.pos, the degrees of
%                the positive operator;
%           2-D  settings_2d, the 16-basis degrees, margins and 'stock'
%                slack of two variables;
%                without settings_2d, as N-D with lpi.var.deg{2} and no
%                'stock' slack;
%           N-D  (nv >= 3) lpi.pos and lpi.var, entry 3: the positive
%                operator on a space in d >= 1 variables, the multiplier
%                basis operator (alpha = 1 in every direction) of degree
%                mult and the full-integral one (alpha = 4, separable in
%                every direction) with caps int in each variable and joint
%                in total, the identity basis on R^n; lpi.var.deg{3}, the
%                degree of a free operator; and settings_nd (optional), a
%                'struct' with the optional fields
%                  LF_deg, LF_opts  instead, a degree specification of the
%                       positive operator and its options;
%                  LF_psatz  the term codes (default 0);
%                  eppos  the margin on every space, or 2^nv margins as
%                       for EYE below (default 1e-2);
%                  epneg  (default the field epneg of ST, else 0).
%
% OUTPUT
% - tl:   'struct' with fields
%           nv             double, the spatial dimension;
%           S              'struct', the completed 2-D settings, or in N-D
%                          the N-D settings above with field
%                          eppos_default ([] in 1-D);
%           dom            the domain argument of the operator
%                          declarations: 1x2 double PIE.dom in 1-D; in 2-D
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
%                          'copvar' X; P is a 'cdopvar' object; from the
%                          storage row in every dimension when lpi.pos of
%                          the program settings carries len_com, else as
%                          above;
%           slack          function handle, [prog,N] = slack(prog,X):
%                          declares the 'stock' slack of an inequality on
%                          the spaces of X; an error without 'stock'
%                          settings (N-D, and 2-D without settings_2d);
%           eye            function handle, E = eye(side,e): the 'copvar'
%                          operator diag(e_k I) on the output (SIDE 1) or
%                          input (SIDE 2) spaces of PIE.T; e a scalar, or
%                          2^nv weights, the weight of a space being entry
%                          1 + sum_k 2^(k-1) over the positions k in
%                          PIE.vars of its variables;
%           ident          function handle, I = ident(Bop,side): the
%                          identity 'copvar' on the output (SIDE 1) or input
%                          (SIDE 2) space of the 'copvar' or 'cdopvar' Bop;
%           zdeg           the degrees of a free operator: 1-D
%                          lpi.var.deg{1}; 2-D a 'struct' with the largest
%                          entry of settings_2d.Zop_deg on every role; N-D
%                          the same with lpi.var.deg{3};
%                          (2-D without settings_2d: lpi.var.deg{2});
%           sop            'struct', settings.sop with each per-dimension
%                          vector (dw_feas, dw_P, dw_Q, max_raise,
%                          max_raise_feas, cap_delta, max_raise_Q) reduced
%                          to its entry min(nv,end);
%           capfor         function handle, cap = capfor(X): the
%                          per-direction caps of the slack bases of an
%                          inequality on X, Jmax(T) + cap_delta of this
%                          dimension on the variables of X, with Jmax(T)
%                          the largest total degree per direction of a
%                          kernel monomial of PIE.T; [] (no bound) when
%                          cap_delta is empty or NaN, X has no variable,
%                          or X has one T lacks.
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

% Dimension and spatial variables of the PIE
nv = PIE.dim;
vars = PIE.vars;
vnames = reshape(vars(:,1).varname,1,[]);
tl = struct('nv',nv,'S',[]);
if nv==1
    % 1-D: margins from eppos, eppos2, epneg; storage from lpi.pos, 'stock'
    % slack from sop.stock; Z degrees from lpi.var
    tl.dom = PIE.dom;
    tl.eppos = [getf(st,'eppos',1e-4); getf(st,'eppos2',1e-6)];
    tl.epneg = getf(st,'epneg',0);
    tl.eppos_default = ~isfield(st,'eppos') || isempty(st.eppos);
    tl.storage = @(prog,X,side) poslpivar_settings_sop(prog,X,st,'lf',side,PIE.dom);
    tl.slack   = @(prog,X) poslpivar_settings_sop(prog,X,st,'slack','out',PIE.dom);
    tl.zdeg    = st.lpi.var.deg{1};     % the degrees of a free operator
elseif nv==2 && isfield(st,'settings_2d')
    % 2-D: everything from settings_2d; one degree cap for Z
    S = settings_2d_sop(st);
    tl.S = S;
    tl.dom = struct('vars',{vnames},'dom',PIE.dom);
    tl.eppos = S.eppos;     tl.epneg = S.epneg;     tl.eppos_default = S.eppos_default;
    tl.storage = @(prog,X,side) poslpivar_settings_2d_sop(prog,X,st,'lf',side,PIE);
    tl.slack   = @(prog,X) poslpivar_settings_2d_sop(prog,X,st,'slack','out',PIE);
    % One degree cap for every role of the free operator: the largest entry
    % of the 2-D free-operator degrees
    m = max([cellmax(S.Zop_deg.dx), cellmax(S.Zop_deg.dy), cellmax(S.Zop_deg.d2), 0]);
    tl.zdeg    = struct('mult',m,'int',[m m],'out',m,'in',m);
elseif nv==2
    % 2-D without settings_2d: margins and fallback storage of
    % SETTINGS_ND, POS_ND (the row: POS_ROW below); no 'stock' slack; Z
    % from lpi.var.deg{2}, built as N-D from entry 3
    S = settings_nd(st);    tl.S = S;
    tl.dom = struct('vars',{vnames},'dom',PIE.dom);
    tl.eppos = S.eppos(:);  tl.epneg = S.epneg;     tl.eppos_default = S.eppos_default;
    if isscalar(tl.eppos),  tl.eppos = tl.eppos*ones(4,1);  end  % one margin per space
    tl.storage = @(prog,X,side) pos_nd(prog,X,side,S,tl.dom);
    tl.slack   = @(prog,X) error('exec_tools_sop:slack',['No stock slack without ' ...
                     'settings_2d; set settings.sop.slack = ''like'', or ' ...
                     'settings.settings_2d = settings_2d_stock(level).']);
    z = st.lpi.var.deg{min(2,end)};
    if isnumeric(z) && isscalar(z),     tl.zdeg = struct('mult',z,'int',[z z],'out',z,'in',z);
    else,                               tl.zdeg = z;    end     % a struct or row: as given
else
    % N-D (a container PIE in 3 or more variables): the storage of
    % settings_nd; no 'stock' slack
    S = settings_nd(st);
    tl.S = S;
    tl.dom = struct('vars',{vnames},'dom',PIE.dom);
    tl.eppos = S.eppos(:);  tl.epneg = S.epneg;     tl.eppos_default = S.eppos_default;
    tl.storage = @(prog,X,side) pos_nd(prog,X,side,S,tl.dom);
    tl.slack   = @(prog,X) error('exec_tools_sop:slack',['The stock slack reads the stock 1-D ' ...
                     'and 2-D settings; for %d variables set settings.sop.slack = ''like''.'],nv);
    z = S.Zop_deg;  tl.zdeg = struct('mult',z,'int',[z z],'out',z,'in',z);
end
% The storage row (lpi.pos with len_com) in every dimension; otherwise the
% storage of the branch above (POS_ROW decides on the program settings).
% Not on a 2-D PIE with an explicit settings_2d: that branch is the
% 16-basis storage those settings define.
if ~(nv==2 && isfield(st,'settings_2d'))
    tl.storage = @(prog,X,side) pos_row(prog,X,side,tl.dom,tl.storage);
end
% Identity operators: weighted on the state spaces, plain on an input or
% output space of an operator
cdom = PIE.dom;     if nv>=2,   cdom = tl.dom;  end % the domain of a container
% diag(e_k I) on the state spaces of the PIE, side 1 the PDE state, 2 the
% fundamental state
tl.eye   = @(side,e) eye_container(PIE,side,e,cdom);
tl.ident = @(Bop,side) ident_op(Bop,side,cdom);    % copvar Bop only
% The options of this dimension: each per-dimension vector among dw_feas,
% dw_P, dw_Q, max_raise, max_raise_feas, cap_delta and max_raise_Q reduced
% to its entry min(nv,end), the last entry serving every higher
% dimension; eppos_stab (one entry per space) is kept
tl.sop = [];
if isfield(st,'sop')
    tl.sop = st.sop;
    f = {'dw_feas','dw_P','dw_Q','max_raise','max_raise_feas','cap_delta','max_raise_Q'};
    for k = 1:numel(f)
        if isfield(tl.sop,f{k}) && isnumeric(tl.sop.(f{k})) && numel(tl.sop.(f{k}))>1
            v = tl.sop.(f{k});   tl.sop.(f{k}) = v(min(nv,end));
        end
    end
end
% Per-direction bounds on the slack bases (settings.sop.cap_delta): Jmax
% of T plus the delta, laid over the variables of each target; none by
% default and where the delta of this dimension is NaN
tl.capfor = @(X) [];
cd = [];    if isstruct(tl.sop) && isfield(tl.sop,'cap_delta'),     cd = tl.sop.cap_delta;  end
if ~isempty(cd) && ~isnan(cd)
    Tc = PIE.T;     [~,~,ti] = get_lift_degs(Tc,struct());  % a 'copvar'
    ct = max(ti.Jmax,0) + cd;
    tl.capfor = @(X) cap_for(X,Tc.vars,ct);
end
end

function Im = eye_container(PIE,side,e,dom)
% IM = EYE_CONTAINER(PIE,SIDE,E,DOM) returns diag(e_k I_{n_k}) on the
% spaces of the container T of a container PIE as a 'copvar' object.
% INPUT  PIE:   'pie_struct' object, T a 'copvar' object;
%        side:  1 (the output spaces of PIE.T) or 2 (its input spaces);
%        e:     double, a scalar for every space, or 2^nv weights: the
%               weight of a space is entry 1 + sum_k 2^(k-1) over the
%               positions k in PIE.vars of its variables (R^n entry 1; in
%               2-D the spaces R, L2[x], L2[y], L2[x,y]);
%        dom:   the domain of the container (1 x 2 double, or a 'struct'
%               with fields vars and dom).
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

function I = ident_op(Bop,side,cdom)
% I = IDENT_OP(BOP,SIDE,CDOM) returns the identity 'copvar' on the output
% (SIDE 1) or input (SIDE 2) spaces of BOP.
% INPUT  Bop:   'copvar' or 'cdopvar' object;
%        side:  1 or 2;
%        cdom:  the domain of the container.
% OUTPUT I:     'copvar' object.
if isa(Bop,'copvar') || isa(Bop,'cdopvar')
    sd = {'out','in'};
    [sp,dm] = copvar_space_list(Bop,sd{side});
    I = eye_copvar_sop(dm,sp,cdom);
else
    % Every PIE operator is a 'copvar' or 0: anything else, such as the 0
    % of an absent channel, has no spaces to read
    error('exec_tools_sop:ident','BOP should be a ''copvar'' or ''cdopvar'' object.');
end
end

function S = settings_nd(st)
% S = SETTINGS_ND(ST) returns the N-D settings:
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
ps = st.lpi.pos;    zd = st.lpi.var.deg;    k3 = @(v) v(min(3,end));
% A storage row need not carry int, mult, joint (unused then): the
% defaults 1, 0, 2 in their place
S = struct('LF_int',k3(getf(ps,'int',1)),'LF_mult',k3(getf(ps,'mult',0)), ...
           'LF_joint',k3(getf(ps,'joint',2)), ...
           'LF_deg',[],'LF_opts',struct(),'LF_psatz',0,'eppos',1e-2, ...
           'epneg',getf(st,'epneg',0),'Zop_deg',zd{min(3,end)});
given = struct();
if isfield(st,'settings_nd') && isstruct(st.settings_nd),   given = st.settings_nd;     end
f = fieldnames(S);
for k = 1:numel(f)
    if isfield(given,f{k}) && ~isempty(given.(f{k})),   S.(f{k}) = given.(f{k});    end
end
S.eppos_default = ~(isfield(given,'eppos') && ~isempty(given.eppos));
end

function [prog,P] = pos_nd(prog,X,side,S,dom)
% [PROG,P] = POS_ND(PROG,X,SIDE,S,DOM) declares the positive operator of the
% N-D settings S on the SIDE spaces of X: with S.LF_deg and S.LF_opts when
% S.LF_deg is given; else on each space in d >= 1 variables the
% multiplier basis operator of degree S.LF_mult and the separable
% full-integral one with caps S.LF_int per variable and 2*S.LF_int in
% total, the identity basis on R^n; the term codes S.LF_psatz.
% INPUT  prog:  'struct', the LPI program;
%        X:     'copvar' or 'cdopvar' object;
%        side:  char, 'out' or 'in';
%        S:     'struct', as SETTINGS_ND returns;
%        dom:   'struct' with fields vars and dom, the domain of P.
% OUTPUT prog:  'struct', the program with the Gram matrices of P;
%        P:     'cdopvar' object.
[sp,dm] = copvar_space_list(X,side);
% An explicit degree specification
if ~isempty(S.LF_deg)
    co = S.LF_opts;     if ~isfield(co,'psatz'),    co.psatz = S.LF_psatz;  end
    [prog,P] = poscopvar(prog,dm,sp,dom,S.LF_deg,co);
    return
end
% The rule: per space, the identity basis on R^n, else the multiplier
% (alpha 1...1) and the full-integral (alpha 4...4) basis operators, in
% the order of include
deg = cell(1,numel(sp));    incl = cell(1,numel(sp));
for k = 1:numel(sp)
    d = numel(sp{k});
    if d==0
        deg{k} = struct('int',zeros(1,numel(dom.vars)));   % one 0 per registry variable
    else
        deg{k} = {struct('int',S.LF_mult,'mult',0,'joint',S.LF_mult), ...
                  struct('int',S.LF_int,'mult',S.LF_int,'joint',S.LF_joint)};
        incl{k} = [ones(1,d); 4*ones(1,d)];
    end
end
co = struct('psatz',S.LF_psatz,'sep',true);
co.include = incl;                  % assigned apart: a cell in struct() builds a struct array
[prog,P] = poscopvar(prog,dm,sp,dom,deg,co);
end

function [prog,P] = pos_row(prog,X,side,dom,stock)
% [PROG,P] = POS_ROW(PROG,X,SIDE,DOM,STOCK) declares the positive operator
% on the SIDE spaces of X: when prog.settings.lpi.pos has the field
% len_com (the storage row), the operator declared without degrees, from
% the row; otherwise STOCK(PROG,X,SIDE).
% INPUT  prog:  'struct', the LPI program, with settings;
%        X:     'copvar' or 'cdopvar' object;
%        side:  char, 'out' or 'in';
%        dom:   the domain of the operator (tl.dom);
%        stock: function handle, the storage of the branch of the dimension.
% OUTPUT prog:  'struct', the program with the Gram matrices of P;
%        P:     'cdopvar' object.
row = isfield(prog,'settings') && isfield(prog.settings,'lpi') && ...
      isfield(prog.settings.lpi,'pos') && isfield(prog.settings.lpi.pos,'len_com');
if row
    [sp,dm] = copvar_space_list(X,side);
    [prog,P] = poscopvar(prog,dm,sp,dom);
else
    [prog,P] = stock(prog,X,side);
end
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

function v = getf(s,f,d)
% V = GETF(S,F,D) returns the field F of the struct S, or D when S has no
% nonempty field F.
% INPUT  s: 'struct'; f: char, the field name; d: the default value.
% OUTPUT v: the value of the field, or D.
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
