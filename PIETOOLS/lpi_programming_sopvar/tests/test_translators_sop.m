function R = test_translators_sop()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = TEST_TRANSLATORS_SOP() checks poslpivar_sop, poslpivar_settings_sop,
% get_lpivar_degs_sop, trace_rn_sop and copvar_space_list in one spatial
% variable s1 on [0,1], against the stock routines or against their
% definitions, and errors at the first failed check:
%   (i)   on legacy inputs, field by field (same_val_sop): poslpivar_sop
%         returns what poslpivar returns at nargout 1 to 5 (degrees
%         {2,[1 2 3],[2 1 3]}, psatz 1, exclude [0 0 0 1]);
%         poslpivar_settings_sop(prog,dim,settings,'lf') returns the sum of
%         the poslpivar operators of dd1/options1 and dd12/options12
%         (lpisettings('custom') with override1 = override2 = 0);
%         get_lpivar_degs_sop returns what get_lpivar_degs returns; and
%         trace_rn_sop of a 'dopvar' returns the trace of its P part;
%   (ii)  each of 12 inputs built to fail raises an error, with the
%         identifier expected where one is listed (poslpivar_sop:side,
%         :dim, :excludeRn, :excludeL2, :nargout;
%         poslpivar_settings_sop:sosineq, :role; get_lpivar_degs_sop:dim;
%         trace_rn_sop:square; copvar_space_list:side), and with the
%         message of poslpivar for malformed degrees and options;
%   (iii) on the container R x L2[s1], poslpivar_settings_sop for the roles
%         'lf' and 'slack' returns the program and operator of two
%         poslpivar_sop calls summed;
%   (iv)  for the plants io1 and rd (frac 0.5) of cx_plant and the
%         lpisettings 'stripped', 'light', 'heavy', 'veryheavy', 'extreme'
%         and 'custom',
%         get_lpivar_degs_sop of the container operators R, R + eppos2 T'T
%         and R + eppos2 T T' returns what get_lpivar_degs returns for the
%         stock operators of the same settings; with the integral bases
%         excluded, both return -1 as third degree;
%   (v)   at random decision values, trace_rn_sop equals the sum of the
%         traces of the R^q diagonal blocks of the substituted operator, to
%         1e-12 (1 + |trace|), for a free operator on R^3, a free operator
%         W on R^2 x L2[s1], Mz*W*Mz' (Mz a random fixed operator), the
%         identity on R^3 x L2[s1] and the identity on L2^2[s1] (trace 0);
%   (vi)  for every pattern of zero dimensions of an 'opvar' with random
%         nonzero blocks (dimensions at most [2 3; 1 2]), copvar_space_list
%         of its 'copvar' (opvar2copvar) lists R before L2[s1], drops the
%         spaces of dimension 0, and returns the nonzero dimensions of each
%         side.
%
% INPUT
% none.
%
% OUTPUT
% - R:      'struct' with field checks ('double', the number of checks
%           passed); one printed line with that number.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - test_translators_sop
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
% Initial coding MMP, 10/06/2026
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body and headers of the local functions; code
%                unchanged. The deleted header text was a history note on
%                the copies these routines replace and a reference to
%                another test.

% No warnings for non-canonical multipliers; seed; s1 on [0,1]; check
% counter
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
rng(20261006);
pvar s1 s1_dum
dom = [0 1];
n = 0;

% % % (i) legacy dispatch
% poslpivar_sop against poslpivar at nargout 1 to 5
pl = lpiprogram(s1,s1_dum,dom);
dd = {2,[1 2 3],[2 1 3]};   op = struct('psatz',1,'exclude',[0 0 0 1],'sep',0);
for no = 1:5
    a = cell(1,no);     b = cell(1,no);
    [a{:}] = poslpivar_sop(pl,[1;1],dd,op);
    [b{:}] = poslpivar(pl,[1;1],dd,op);
    [tf,where] = same_val_sop(a,b,'out');
    assert(tf,'(i) poslpivar_sop nargout %d differs at %s',no,where);    n = n+1;
end
% poslpivar_settings_sop ('lf') against the two poslpivar terms
st = lpisettings('custom');
st.override1 = 0;   st.override2 = 0;           % both second terms declared
[pa,Pa] = poslpivar_settings_sop(pl,[1 1;1 1],st,'lf');
[pb,P1] = poslpivar(pl,[1 1;1 1],st.dd1,st.options1);
[pb,P2] = poslpivar(pb,[1 1;1 1],st.dd12,st.options12);
assert(same_val_sop({pa,Pa},{pb,P1+P2}),'(i) legacy lf pair differs');      n = n+1;
% get_lpivar_degs_sop and trace_rn_sop on legacy operators
PIE = cx_plant('rd',0.5);
[~,Rop] = poslpivar(lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom),PIE.T.dim,st.dd1,st.options1);
assert(isequal(get_lpivar_degs_sop(Rop,PIE.T),get_lpivar_degs(Rop,PIE.T)),'(i) legacy Qdeg');   n = n+1;
[~,Wop] = lpivar(pl,[2 2;0 0],1);
assert(same_val_sop(trace_rn_sop(Wop),trace(Wop.P)),'(i) legacy trace');   n = n+1;

% % % (ii) errors
% Containers: R^1 x L2[s1], L2[s1] -> R^2, one space over two variables
Xm = eye_copvar_sop([1;1],{{},{'s1'}},dom);         % R^1 x L2[s1]
Xr = rand_copvar(struct('out',{{{}}},'in',{{{'s1'}}}),struct('out',2,'in',1),dom,1,0.8);  % L2 -> R^2
X2 = eye_copvar_sop(1,{{'s1','s2'}},[0 1;0 1]);     % one space over two variables
pc = lpiprogram(s1,s1_dum,dom);
E = { @() poslpivar_sop(pc,Xr,1,[]),                                   'poslpivar_sop:side'
      @() poslpivar_sop(lpiprogram_sop({'s1','s2'},[0 1;0 1]),X2,1,[]), 'poslpivar_sop:dim'
      @() poslpivar_sop(pc,Xm,1,struct('exclude',[1 0 0 0])),         'poslpivar_sop:excludeRn'
      @() poslpivar_sop(pc,Xm,1,struct('exclude',[0 1 1 1])),         'poslpivar_sop:excludeL2'
      @() nout3(pc,Xm),                                                'poslpivar_sop:nargout'
      @() poslpivar_sop(pc,Xm,[1 2 3 4],[]),                           ''
      @() poslpivar_sop(pc,Xm,1,7),                                    ''
      @() poslpivar_settings_sop(pc,Xm,setfield(st,'sosineq_on',1),'slack'), 'poslpivar_settings_sop:sosineq'
      @() poslpivar_settings_sop(pc,Xm,st,'other'),                    'poslpivar_settings_sop:role'
      @() get_lpivar_degs_sop(X2),                                     'get_lpivar_degs_sop:dim'
      @() trace_rn_sop(rand_copvar(struct('out',{{{}}},'in',{{{}}}),struct('out',3,'in',2),dom,0,0.8)), 'trace_rn_sop:square'
      @() copvar_space_list(Xm,'both'),                                'copvar_space_list:side' };
% Each case raises an error, with the listed identifier where one is given
for k = 1:size(E,1)
    id = errid(E{k,1});
    assert(~isempty(id) || ~isempty(errmsg(E{k,1})),'(ii) case %d raised no error',k);
    if ~isempty(E{k,2})
        assert(strcmp(id,E{k,2}),'(ii) case %d: got ''%s'', expected ''%s''',k,id,E{k,2});
    end
    n = n+1;
end
% poslpivar's own messages for malformed degrees and options
assert(strcmp(errmsg(E{6,1}),errmsg(@() poslpivar(pl,[1;1],[1 2 3 4],struct()))),'(ii) degree message');
assert(strcmp(errmsg(E{7,1}),errmsg(@() poslpivar(pl,[1;1],1,7))),'(ii) options message');
n = n+2;

% % % (iii) poslpivar_settings_sop composition, container
% For each role, its two settings terms by poslpivar_sop, summed
for role = {'lf','slack'}
    if strcmp(role{1},'lf'),  f = {'dd1','options1','dd12','options12'};
    else,                     f = {'dd2','options2','dd3','options3'};
    end
    [pa,Pa] = poslpivar_settings_sop(pc,Xm,st,role{1},'out',dom);
    [pb,P1] = poslpivar_sop(pc,Xm,st.(f{1}),st.(f{2}),'out',dom);
    [pb,P2] = poslpivar_sop(pb,Xm,st.(f{3}),st.(f{4}),'out',dom);
    [tf,where] = same_val_sop({pa,Pa},{pb,P1+P2},role{1});
    assert(tf,'(iii) %s pair differs at %s',role{1},where);    n = n+1;
end

% % % (iv) Q degrees against stock get_lpivar_degs
% For each plant and preset: the stock R, and the container R on the
% output and on the input spaces of T; three forms each
for pk = {{'io1'},{'rd',0.5}}
    PIE = cx_plant(pk{1}{:});   T = PIE.T;  Tm = opvar2copvar(T);
    for ps = {'stripped','light','heavy','veryheavy','extreme','custom'}
        st6 = lpisettings(ps{1});   e2 = st6.eppos2;
        p0 = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
        [~,Rop] = poslpivar_settings_sop(p0,T.dim,st6,'lf');
        [~,Ro] = poslpivar_settings_sop(p0,Tm,st6,'lf','out',PIE.dom);
        [~,Ri] = poslpivar_settings_sop(p0,Tm,st6,'lf','in',PIE.dom);
        pairs = {get_lpivar_degs(Rop,T),          get_lpivar_degs_sop(Ro);
                 get_lpivar_degs(Rop+e2*T'*T,T),  get_lpivar_degs_sop(Ri+e2*Tm'*Tm);
                 get_lpivar_degs(Rop+e2*T*T',T),  get_lpivar_degs_sop(Ro+e2*Tm*Tm')};
        for r = 1:3
            assert(isequal(full(pairs{r,1}),pairs{r,2}),'(iv) %s/%s form %d: stock %s, container %s', ...
                   pk{1}{1},ps{1},r,mat2str(full(pairs{r,1})),mat2str(pairs{r,2}));
            n = n+1;
        end
    end
end
[~,R0] = poslpivar_sop(pc,Xm,1,struct('exclude',[0 0 1 1]));      % no R1, R2: deg2 = 0
[~,R0l] = poslpivar(pl,[1;1],1,struct('exclude',[0 0 1 1]));
Tl = opvar();  Tl.dim = [1 1;1 1];  Tl.I = dom;     % get_lpivar_degs needs an opvar T
q0 = get_lpivar_degs_sop(R0);
assert(isequal(q0,full(get_lpivar_degs(R0l,Tl))) && q0(3)==-1, ...
       '(iv) Z2, Z3 excluded: stock %s, container %s (expected deg2 = 0 -> -1)',mat2str(full(get_lpivar_degs(R0l,Tl))),mat2str(q0));
n = n+1;

% % % (v) trace against its definition
% The five operators: free on R^3, free W on R^2 x L2[s1], Mz*W*Mz', and
% two identities
[pw,W3] = lpivar_cdopvar(pc,struct('out',3,'in',3),struct('out',{{{}}},'in',{{{}}}),dom,[1 1 1]);
[pw,Wmx] = lpivar_cdopvar(pw,struct('out',[2;1],'in',[2;1]), ...
                          struct('out',{{{},{'s1'}}},'in',{{{},{'s1'}}}),dom,[1 1 1]);
Mz = rand_copvar(struct('out',{{{}}},'in',{{{},{'s1'}}}),struct('out',1,'in',[2;1]),dom,1,0.8);
Xs = {W3, Wmx, (Mz*Wmx)*Mz', eye_copvar_sop([3;1],{{},{'s1'}},dom), eye_copvar_sop(2,{{'s1'}},dom)};
for k = 1:numel(Xs)
    % trace_rn_sop, and the operator and its trace at random values
    X = Xs{k};
    tr = trace_rn_sop(X);
    if isa(X,'cdopvar')
        v = randn(numel(X.Zd),1);
        Xv = subs_dvar_sop(X,X.Zd,v);    trv = double(subs_dvar_sop(tr,X.Zd,v));
    else
        Xv = X;     trv = double(tr.C(1));
    end
    % The definition: traces of the diagonal blocks on spaces without a
    % spatial variable
    ref = 0;
    for kk = find(~any(Xv.space_out,2)).'
        if ~isempty(Xv.C{kk,kk}),   ref = ref + trace(full(Xv.C{kk,kk}.params{1}));   end
    end
    assert(abs(trv-ref) <= 1e-12*(1+abs(ref)),'(v) case %d: trace %g, definition %g',k,trv,ref);
    n = n+1;
end

% % % (vi) space lists against opvar dims
% Every pattern of zero dimensions with a nonzero input and output
for pat = 0:15
    dm = [bitand(pat,1),bitand(pat,2)/2; bitand(pat,4)/4,bitand(pat,8)/8].*[2 3;1 2];
    if ~any(dm(:,1)) || ~any(dm(:,2)),  continue,   end
    Pop = opvar();  Pop.dim = dm;   Pop.I = dom;    Pop.var1 = s1;  Pop.var2 = s1_dum;
    Pop.P = rand(dm(1,1),dm(1,2));      Pop.Q1 = rand(dm(1,1),dm(2,2));     % nonzero, so
    Pop.Q2 = rand(dm(2,1),dm(1,2));     Pop.R.R0 = rand(dm(2,1),dm(2,2));   % every block is kept
    Pc = opvar2copvar(Pop);
    % Each side: R before L2[s1], spaces of dimension 0 dropped
    for side = {'out','in'}
        c = 1 + strcmp(side{1},'in');
        keep = dm(:,c)>0;
        spx = {cell(1,0), {'s1'}};      spx = spx(keep);    dmx = dm(keep,c);
        [sp,dmv] = copvar_space_list(Pc,side{1});
        assert(isequal(sp,spx) && isequal(dmv,dmx),'(vi) dim %s side %s',mat2str(dm),side{1});
        n = n+1;
    end
end

% Summary and output
fprintf('test_translators_sop: %d checks passed.\n',n);
R = struct('checks',n);
end


function nout3(p,X)
% NOUT3(P,X) calls poslpivar_sop(P,X,1,[]) with three outputs.
% INPUT  p: 'struct', an LPI program; X: 'copvar' object.
% OUTPUT none.
[~,~,~] = poslpivar_sop(p,X,1,[]);
end

function id = errid(f)
% ID = ERRID(F) returns the identifier of the error that F() raises (''
% when it raises none).
% INPUT  f: 'function_handle', called with no argument.
% OUTPUT id: 'char'.
id = '';
try,        f();
catch ME,   id = ME.identifier;
end
end

function m = errmsg(f)
% M = ERRMSG(F) returns the message of the error that F() raises ('' when
% it raises none).
% INPUT  f: 'function_handle', called with no argument.
% OUTPUT m: 'char'.
m = '';
try,        f();
catch ME,   m = ME.message;
end
end
