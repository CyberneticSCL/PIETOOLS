function R = test_translators_sop()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = TEST_TRANSLATORS_SOP() tests the 1-D executive translators of
% lpi_programming_sopvar against the stock routines they mirror or against
% their definitions (never against the test-folder copies they replace):
%
% (i)   legacy dispatch: poslpivar_sop, poslpivar_settings_sop,
%       get_lpivar_degs_sop and trace_rn_sop on legacy operators return,
%       field by field, what poslpivar / get_lpivar_degs / trace(.P) return
%       (poslpivar at nargout 1..5);
% (ii)  errors: every documented identifier, from an input built to raise it;
% (iii) poslpivar_settings_sop on a container: the program and operator of
%       its explicit composition (two poslpivar_sop calls, P1 first, then
%       P1 + P2), for both roles, with override = 0;
% (iv)  get_lpivar_degs_sop on the container against stock get_lpivar_degs
%       on the stock operator declared from the same settings: plants io1
%       and rd, all six presets, the forms R, R + eppos2 T'T, R + eppos2 TT';
%       and deg2 = 0 gives -1, as the stock rule;
% (v)   trace_rn_sop against its definition: at random decision values, the
%       trace of the R^q diagonal blocks of the substituted operator, for a
%       free R^3 operator, a mixed R^2 x L2 one, Mz*W*Mz', a fixed
%       container and one without R^q (trace 0);
% (vi)  copvar_space_list against the opvar dims of opvar2copvar's input:
%       R first, zero rows dropped, for every pattern of zero dimensions.
% The span comparison of poslpivar_sop with poslpivar (degree forms, d2~=d3,
% psatz, exclude, sep) is the second part of test_poscopvar_vs_poslpivar.
%
% Initial coding MMP, 10/06/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
rng(20261006);
pvar s1 s1_dum
dom = [0 1];
n = 0;

% % % (i) legacy dispatch
pl = lpiprogram(s1,s1_dum,dom);
dd = {2,[1 2 3],[2 1 3]};   op = struct('psatz',1,'exclude',[0 0 0 1],'sep',0);
for no = 1:5
    a = cell(1,no);     b = cell(1,no);
    [a{:}] = poslpivar_sop(pl,[1;1],dd,op);
    [b{:}] = poslpivar(pl,[1;1],dd,op);
    [tf,where] = same_val_sop(a,b,'out');
    assert(tf,'(i) poslpivar_sop nargout %d differs at %s',no,where);    n = n+1;
end
st = lpisettings('custom');
st.override1 = 0;   st.override2 = 0;           % both second terms declared
[pa,Pa] = poslpivar_settings_sop(pl,[1 1;1 1],st,'lf');
[pb,P1] = poslpivar(pl,[1 1;1 1],st.dd1,st.options1);
[pb,P2] = poslpivar(pb,[1 1;1 1],st.dd12,st.options12);
assert(same_val_sop({pa,Pa},{pb,P1+P2}),'(i) legacy lf pair differs');      n = n+1;
PIE = cx_plant('rd',0.5);
[~,Rop] = poslpivar(lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom),PIE.T.dim,st.dd1,st.options1);
assert(isequal(get_lpivar_degs_sop(Rop,PIE.T),get_lpivar_degs(Rop,PIE.T)),'(i) legacy Qdeg');   n = n+1;
[~,Wop] = lpivar(pl,[2 2;0 0],1);
assert(same_val_sop(trace_rn_sop(Wop),trace(Wop.P)),'(i) legacy trace');   n = n+1;

% % % (ii) errors
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
[pw,W3] = lpivar_cdopvar(pc,struct('out',3,'in',3),struct('out',{{{}}},'in',{{{}}}),dom,[1 1 1]);
[pw,Wmx] = lpivar_cdopvar(pw,struct('out',[2;1],'in',[2;1]), ...
                          struct('out',{{{},{'s1'}}},'in',{{{},{'s1'}}}),dom,[1 1 1]);
Mz = rand_copvar(struct('out',{{{}}},'in',{{{},{'s1'}}}),struct('out',1,'in',[2;1]),dom,1,0.8);
Xs = {W3, Wmx, (Mz*Wmx)*Mz', eye_copvar_sop([3;1],{{},{'s1'}},dom), eye_copvar_sop(2,{{'s1'}},dom)};
for k = 1:numel(Xs)
    X = Xs{k};
    tr = trace_rn_sop(X);
    if isa(X,'cdopvar')
        v = randn(numel(X.Zd),1);
        Xv = subs_dvar_sop(X,X.Zd,v);    trv = double(subs_dvar_sop(tr,X.Zd,v));
    else
        Xv = X;     trv = double(tr.C(1));
    end
    ref = 0;
    for kk = find(~any(Xv.space_out,2)).'
        if ~isempty(Xv.C{kk,kk}),   ref = ref + trace(full(Xv.C{kk,kk}.params{1}));   end
    end
    assert(abs(trv-ref) <= 1e-12*(1+abs(ref)),'(v) case %d: trace %g, definition %g',k,trv,ref);
    n = n+1;
end

% % % (vi) space lists against opvar dims
for pat = 0:15
    dm = [bitand(pat,1),bitand(pat,2)/2; bitand(pat,4)/4,bitand(pat,8)/8].*[2 3;1 2];
    if ~any(dm(:,1)) || ~any(dm(:,2)),  continue,   end
    Pop = opvar();  Pop.dim = dm;   Pop.I = dom;    Pop.var1 = s1;  Pop.var2 = s1_dum;
    Pop.P = rand(dm(1,1),dm(1,2));      Pop.Q1 = rand(dm(1,1),dm(2,2));     % nonzero, so
    Pop.Q2 = rand(dm(2,1),dm(1,2));     Pop.R.R0 = rand(dm(2,1),dm(2,2));   % every block is kept
    Pc = opvar2copvar(Pop);
    for side = {'out','in'}
        c = 1 + strcmp(side{1},'in');
        keep = dm(:,c)>0;
        spx = {cell(1,0), {'s1'}};      spx = spx(keep);    dmx = dm(keep,c);
        [sp,dmv] = copvar_space_list(Pc,side{1});
        assert(isequal(sp,spx) && isequal(dmv,dmx),'(vi) dim %s side %s',mat2str(dm),side{1});
        n = n+1;
    end
end

fprintf('test_translators_sop: %d checks passed.\n',n);
R = struct('checks',n);
end


function nout3(p,X)
[~,~,~] = poslpivar_sop(p,X,1,[]);
end

function id = errid(f)
id = '';
try,        f();
catch ME,   id = ME.identifier;
end
end

function m = errmsg(f)
m = '';
try,        f();
catch ME,   m = ME.message;
end
end
