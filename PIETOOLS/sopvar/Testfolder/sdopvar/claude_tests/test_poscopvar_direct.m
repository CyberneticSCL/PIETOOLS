%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_POSCOPVAR_DIRECT checks the direct Gram-to-kernel map of
% 'poscopvar_direct' against 'poscopvar' (copquadvar), and compares its two
% paths and its cost.
%
% (1) Identity with poscopvar: on the same program the two declare the
%     same decision variables, in the same order, and every coefficient of
%     every block and cell agrees to rounding (keyed by block, cell,
%     entries, output and input exponents; a basis both routes do not
%     share is padding on one side and zero on the other). 1-D and 2-D,
%     R x L2 and single spaces, shifted domains, sep, include, joint and
%     subset caps, the product and face Positivstellensatz terms, and 3-D.
% (2) The 'fast' path (canonical weight as a Gram shift, g = 1 tables) and
%     the 'general' path (the weight's own tables) agree to rounding, on
%     domains with a ~= 0.
% (3) The cache: a second call with the same degrees reuses the position
%     map (cache_hit) and gives the same operator.
% (4) End to end: V'*V + L'*L certified at int 0, mult 0 with the lower and
%     upper components, the solved weight read back from the Gram.
% (5) Build time of poscopvar, poscopvar_direct cold and warm, and the
%     general path, printed with nnz(B).
% (6) Several terms in one call (psatz a row, psatz_offset per term) equal
%     the sum of the single-term calls on the same program, both through
%     poscopvar_direct and through poscopvar + plus_batch: same decision
%     variables, coefficients to rounding; Qcell one cell per term.
%
% Initial coding MMP, 10/07/2026
% MMP, 10/08/2026: part (6), the several-term call.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
for nm = {'poscopvar','copquadvar','poscopvar_direct','lift_copvar','lpi_eq_sop'}
    if numel(which(nm{1},'-all'))~=1
        error('test_poscopvar_direct: %s resolves to %d files; fix the path first.',...
              nm{1},numel(which(nm{1},'-all')));
    end
end
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
poscopvar_direct('clear');
npass = 0;
parts = getappdata(0,'direct_test_parts');
if isempty(parts),  parts = 1:6;    end
tol = 1e-12;

%% (1) identity with poscopvar
if any(parts==1)
sub2d = 3*ones(1,16);   sub2d(16) = 4;      % 4 basis variables [th1 th2 s1 s2]: pairs/triples <= 3, all <= 4
cases = {
  {'1D R x L2, 1/1',                 [1;1],  {{},{'s1'}},             [0 1],     {struct('int',0),struct('int',1,'mult',1)}, struct()}
  {'1D R^2 x L2^2, 2/1',             [2;2],  {{},{'s1'}},             [0 1],     {struct('int',0),struct('int',2,'mult',1)}, struct()}
  {'1D L2, 2/2',                     1,      {{'s1'}},                [0 1],     struct('int',2,'mult',2), struct()}
  {'1D L2^3, 1/2, dom [-1 2]',       3,      {{'s1'}},                [-1 2],    struct('int',1,'mult',2), struct()}
  {'1D R x L2, R weight 1',          [1;1],  {{},{'s1'}},             [0 1],     {struct('int',1),struct('int',1,'mult',1)}, struct()}
  {'1D R x L2, scalar deg 2',        [1;1],  {{},{'s1'}},             [0 1],     2, struct()}
  {'1D L2, product, dom [-1 2]',     1,      {{'s1'}},                [-1 2],    struct('int',1,'mult',1), struct('psatz',1)}
  {'1D L2, face 3, dom [1 3]',       1,      {{'s1'}},                [1 3],     struct('int',2,'mult',1), struct('psatz',3)}
  {'1D L2, face 4',                  1,      {{'s1'}},                [0 1],     struct('int',1,'mult',2), struct('psatz',4)}
  {'1D R x L2, sep',                 [1;1],  {{},{'s1'}},             [0 1],     struct('int',1,'mult',1), struct('sep',true)}
  {'1D R x L2, include {[],[1;2]}',  [1;1],  {{},{'s1'}},             [0 1],     struct('int',1,'mult',1), struct('include',{{[],[1;2]}})}
  {'1D L2, joint 3',                 1,      {{'s1'}},                [0 1],     struct('int',3,'mult',3,'joint',3), struct()}
  {'1D L2^2, per-component',         2,      {{'s1'}},                [0 1],     {{struct('int',1),struct('int',2,'mult',1),struct('int',1,'mult',2)}}, struct()}
  {'2D L2[s1,s2], 1/1',              1,      {{'s1','s2'}},           [0 1;0 2], struct('int',1,'mult',1), struct()}
  {'2D L2^2[s1,s2], 2/1',            2,      {{'s1','s2'}},           [0 1;0 1], struct('int',2,'mult',1), struct()}
  {'2D R x L2[s1] x L2[s1,s2]',      [1;1;1],{{},{'s1'},{'s1','s2'}}, [0 1;0 2], struct('int',1,'mult',1), struct()}
  {'2D L2[s1] x L2[s2], int [1 2]',  [1;1],  {{'s1'},{'s2'}},         [0 1;0 2], struct('int',[1 2],'mult',[1 1]), struct()}
  {'2D L2[s1,s2], product',          1,      {{'s1','s2'}},           [0 1;0 2], struct('int',1,'mult',1), struct('psatz',1)}
  {'2D L2[s1,s2], face 5, dom',      1,      {{'s1','s2'}},           [-1 1;0 2],struct('int',1,'mult',1), struct('psatz',5)}
  {'2D L2[s1,s2], sep [1 0]',        1,      {{'s1','s2'}},           [0 1;0 2], struct('int',1,'mult',1), struct('sep',[true,false])}
  {'2D L2[s1,s2], joint 3',          1,      {{'s1','s2'}},           [0 1;0 1], struct('int',2,'mult',2,'joint',3), struct()}
  {'2D L2[s1,s2], subset caps',      1,      {{'s1','s2'}},           [0 1;0 1], struct('int',2,'mult',2,'subset',sub2d), struct()}
  {'2D R x L2[s1,s2], 2/1',          [1;1],  {{},{'s1','s2'}},        [0 1;0 2], {struct('int',0),struct('int',1,'mult',2)}, struct()}
  {'3D L2[s1,s2,s3], 1/1',           1,      {{'s1','s2','s3'}},      [0 1;0 1;0 1], struct('int',1,'mult',1), struct()}
  };
for ic = 1:numel(cases)
    [lbl,dims,spaces,dom,deg,opts] = deal(cases{ic}{:});
    vars = reshape(unique([spaces{:}]),1,[]);
    prog = lpiprogram_sop(vars,dom);
    t0 = tic;   [~,Pc] = poscopvar(prog,dims,spaces,dom,deg,opts);              tc = toc(t0);
    prog = lpiprogram_sop(vars,dom);
    t0 = tic;   [~,Pd,~,id] = poscopvar_direct(prog,dims,spaces,dom,deg,opts);  td = toc(t0);
    if ~isequal(Pc.Zd(:),Pd.Zd(:))
        error('test_poscopvar_direct: case ''%s'': the decision variables differ.',lbl);
    end
    [mx,ref] = coef_diff(Pc,Pd,vars);
    if ~(mx<=tol*max(ref,1))
        error('test_poscopvar_direct: case ''%s'': coefficients differ, max |diff| %.3g.',lbl,mx);
    end
    v = verify(Pd);
    if ~v.true
        error('test_poscopvar_direct: case ''%s'' fails verify: %s',lbl,strjoin(v.flags,' | '));
    end
    fprintf('  passed: %-32s ndec %6d  max|diff| %.1e  poscopvar %.2fs  direct %.2fs\n',lbl,numel(Pd.Zd),mx,tc,td);
    npass = npass+1;
end
fprintf('(1) poscopvar_direct = poscopvar on %d cases.\n\n',numel(cases));
end

%% (2) fast against general
if any(parts==2)
cases = {
  {'1D L2, product, dom [-1 2]',     1,      {{'s1'}},         [-1 2],     struct('int',2,'mult',1), 1}
  {'1D R x L2, face 3, dom [1 3]',   [1;1],  {{},{'s1'}},      [1 3],      struct('int',1,'mult',1), 3}
  {'2D L2[s1,s2], product, dom',     1,      {{'s1','s2'}},    [-1 1;0 2], struct('int',1,'mult',1), 1}
  {'2D L2[s1,s2], face 6, dom',      1,      {{'s1','s2'}},    [-1 1;0 2], struct('int',1,'mult',1), 6}
  {'2D L2[s1,s2], joint 3, product', 1,      {{'s1','s2'}},    [0 1;0 1],  struct('int',2,'mult',1,'joint',3), 1}
  };
for ic = 1:numel(cases)
    [lbl,dims,spaces,dom,deg,code] = deal(cases{ic}{:});
    vars = reshape(unique([spaces{:}]),1,[]);
    prog = lpiprogram_sop(vars,dom);
    [~,Pf] = poscopvar_direct(prog,dims,spaces,dom,deg,struct('psatz',code,'path','fast'));
    prog = lpiprogram_sop(vars,dom);
    [~,Pg] = poscopvar_direct(prog,dims,spaces,dom,deg,struct('psatz',code,'path','general'));
    [mx,ref] = coef_diff(Pf,Pg,vars);
    if ~isequal(Pf.Zd(:),Pg.Zd(:)) || ~(mx<=tol*max(ref,1))
        error('test_poscopvar_direct: fast and general differ on ''%s'': max |diff| %.3g.',lbl,mx);
    end
    fprintf('  passed: %-32s fast = general, max|diff| %.1e\n',lbl,mx);
    npass = npass+1;
end
fprintf('(2) fast = general on %d cases.\n\n',numel(cases));
end

%% (3) the cache
if any(parts==3)
vars = {'s1','s2'};     dom = [0 1;0 2];    spaces = {vars};
deg = struct('int',2,'mult',2);
poscopvar_direct('clear');
prog = lpiprogram_sop(vars,dom);
t0 = tic;   [~,P1,~,i1] = poscopvar_direct(prog,1,spaces,dom,deg);    t1 = toc(t0);
prog = lpiprogram_sop(vars,dom);
t0 = tic;   [~,P2,~,i2] = poscopvar_direct(prog,1,spaces,dom,deg);    t2 = toc(t0);
if i1.cache_hit || ~i2.cache_hit
    error('test_poscopvar_direct: the cache should miss on the first call and hit on the second.');
end
[mx,ref] = coef_diff(P1,P2,vars);
if ~(mx==0)
    error('test_poscopvar_direct: a cache hit changed the operator (%.3g).',mx);
end
% the product term at int 1 lands on the plain term's index set at int 2
prog = lpiprogram_sop(vars,dom);
[~,~,~,i3] = poscopvar_direct(prog,1,spaces,dom,struct('int',1,'mult',2),struct('psatz',1));
if ~i3.cache_hit
    error('test_poscopvar_direct: the product term at int 1 should reuse the map of the plain term at int 2.');
end
fprintf('  passed: cache miss %.2fs (map %.2fs), hit %.2fs (map %.2fs), identical; the product term at int 1 hits the plain map at int 2\n',...
        t1,i1.t_map,t2,i2.t_map);
npass = npass+1;
fprintf('(3) cache passed.\n\n');
end

%% (4) end to end: V'*V + L'*L
if any(parts==4)
vars = {'s1'};  dom = [0 1];
ZL = {[0;1]};   ZR = {[0;1]};
C2 = sparse([2 0;-1 0]);    C3 = sparse([2 -1;0 0]);
params = reshape({sparse(2,2),C2,C3},[3 1 1]);
Pt = copvar({sopvar(params,struct('out',{vars},'in',{vars}),ZL,ZR,struct('out',dom,'in',dom),[1 1])});
prog = lpiprogram_sop(Pt);
[prog,Pop,~,i4] = poscopvar_direct(prog,1,{vars},dom,struct('int',0,'mult',0),struct('include',{{[2;3]}}));
prog = lpi_eq_sop(prog,Pop-Pt,'symmetric');
sopts = struct('solver','mosek','simplify',false);
if isempty(which('mosekopt')),  sopts.solver = 'sedumi';   end
prog = lpisolve(prog,sopts);
rel_b = cx_resid(prog);
% the Gram read back from its dpvar: Q = [q11 q12; q12 q22]
W = double(sosgetsol(prog,i4.Qdpvar));
fprintf('(4) V''*V + L''*L at int 0, mult 0, components {2,3}: %d decision variables, rel_b %.1e, W = [%.3f %.3f; %.3f %.3f], eig %.3f %.3f\n',...
        numel(Pop.Zd),rel_b,W(1,1),W(1,2),W(2,1),W(2,2),min(eig(W)),max(eig(W)));
if rel_b>1e-6 || min(eig(W))<-1e-6 || norm(W-[2 1;1 1])>1e-5
    error('test_poscopvar_direct: the end-to-end certificate is not reproduced.');
end
npass = npass+1;
fprintf('\n');
end

%% (5) build time
if any(parts==5)
tcases = {
  {'1D R^2 x L2^2, 2/2',          [2;2],  {{},{'s1'}},                       [0 1],         {struct('int',0),struct('int',2,'mult',2)}, 0}
  {'1D R^2 x L2^3, 3/3',          [2;3],  {{},{'s1'}},                       [0 1],         {struct('int',0),struct('int',3,'mult',3)}, 0}
  {'1D L2^16, 2/2',               16,     {{'s1'}},                          [0 1],         struct('int',2,'mult',2), 0}
  {'2D L2[s1,s2], 1/1',           1,      {{'s1','s2'}},                     [0 1;0 1],     struct('int',1,'mult',1), 0}
  {'2D L2[s1,s2], 2/2',           1,      {{'s1','s2'}},                     [0 1;0 1],     struct('int',2,'mult',2), 0}
  {'2D L2^4[s1,s2], 2/2',         4,      {{'s1','s2'}},                     [0 1;0 1],     struct('int',2,'mult',2), 0}
  {'2D four spaces, 1/1',         [1;1;1;1], {{},{'s1'},{'s2'},{'s1','s2'}}, [0 1;0 1],     {struct('int',0),struct('int',1,'mult',1),struct('int',1,'mult',1),struct('int',1,'mult',1)}, 0}
  {'2D L2[s1,s2], 2/2, product',  1,      {{'s1','s2'}},                     [0 1;0 1],     struct('int',2,'mult',2), 1}
  {'3D L2[s1,s2,s3], 1/1',        1,      {{'s1','s2','s3'}},                [0 1;0 1;0 1], struct('int',1,'mult',1), 0}
  {'3D L2[s1,s2,s3], 1/1, product',1,     {{'s1','s2','s3'}},                [0 1;0 1;0 1], struct('int',1,'mult',1), 1}
  };
tsel = getappdata(0,'direct_test_tcases');
if isempty(tsel),   tsel = 1:numel(tcases);     end
fprintf('%-32s %7s %9s %10s %10s %10s %9s\n','(5) build time','ndec','poscopvar','direct cold','direct warm','general','nnz(B)');
for ic = reshape(tsel,1,[])
    [lbl,dims,spaces,dom,deg,code] = deal(tcases{ic}{:});
    vars = reshape(unique([spaces{:}]),1,[]);
    opts = struct('psatz',code);
    poscopvar_direct('clear');
    prog = lpiprogram_sop(vars,dom);
    t0 = tic;   [~,Pc] = poscopvar(prog,dims,spaces,dom,deg,opts);              tc = toc(t0);
    prog = lpiprogram_sop(vars,dom);
    t0 = tic;   [~,Pd] = poscopvar_direct(prog,dims,spaces,dom,deg,opts);       td1 = toc(t0);
    prog = lpiprogram_sop(vars,dom);
    t0 = tic;   [~,Pd] = poscopvar_direct(prog,dims,spaces,dom,deg,opts);       td2 = toc(t0);
    optg = opts;    optg.path = 'general';  optg.cache = false;
    prog = lpiprogram_sop(vars,dom);
    t0 = tic;   [~,Pg] = poscopvar_direct(prog,dims,spaces,dom,deg,optg);       tg = toc(t0);
    if numel(Pc.Zd)~=numel(Pd.Zd) || numel(Pg.Zd)~=numel(Pd.Zd)
        error('test_poscopvar_direct: timing case ''%s'': ndec differ.',lbl);
    end
    nnzB = 0;
    for i = 1:numel(Pd.C)
        if ~isempty(Pd.C{i}),   nnzB = nnzB+sum(cellfun(@nnz,Pd.C{i}.params.B(:)));    end
    end
    fprintf('%-32s %7d %8.2fs %9.2fs %9.2fs %9.2fs %9d\n',lbl,numel(Pd.Zd),tc,td1,td2,tg,nnzB);
end
fprintf('\n');
end

%% (6) several terms in one call = the sum of the single-term calls
if any(parts==6)
cases = {
  {'1D R x L2, [0 1] at [0 1]',       [1;1],  {{},{'s1'}},     [0 1],     {struct('int',0),struct('int',2,'mult',1)}, [0 1],       [0 1]}
  {'1D L2, [0 3 4] at [0 1 1]',       1,      {{'s1'}},        [-1 2],    struct('int',2,'mult',2),                     [0 3 4],     [0 1 1]}
  {'2D L2, plain + 4 faces',          1,      {{'s1','s2'}},   [0 1;0 2], struct('int',1,'mult',1),                     [0 3 4 5 6], 0}
  {'2D L2, joint 3, [0 1] at [0 1]',  1,      {{'s1','s2'}},   [0 1;0 1], struct('int',2,'mult',1,'joint',3),           [0 1],       [0 1]}
  {'2D R x L2, [0 1 3] at [0 1 0]',   [1;1],  {{},{'s1','s2'}},[0 1;0 2], {struct('int',0),struct('int',1,'mult',2)},   [0 1 3],     [0 1 0]}
  };
for ic = 1:numel(cases)
    [lbl,dims,spaces,dom,deg,codes,offs] = deal(cases{ic}{:});
    vars = reshape(unique([spaces{:}]),1,[]);   nv = numel(vars);
    if isscalar(offs),  offs = repmat(offs,size(codes));    end
    % the one call
    prog = lpiprogram_sop(vars,dom);
    t0 = tic;   [~,Pm,Qm,im] = poscopvar_direct(prog,dims,spaces,dom,deg,struct('psatz',codes,'psatz_offset',offs));   tm = toc(t0);
    % the single-term calls on one program, each at the term's own degrees
    prog = lpiprogram_sop(vars,dom);    Ps = [];   Pcq = [];
    for c = 1:numel(codes)
        degc = reduced(deg,offs(c),codes(c),nv);
        [prog,Pc] = poscopvar_direct(prog,dims,spaces,dom,degc,struct('psatz',codes(c)));
        if isempty(Ps),     Ps = Pc;  else,   Ps = Ps+Pc;    end
    end
    % the same through poscopvar (copquadvar) and plus_batch
    prog = lpiprogram_sop(vars,dom);    terms = cell(1,numel(codes));
    for c = 1:numel(codes)
        degc = reduced(deg,offs(c),codes(c),nv);
        [prog,terms{c}] = poscopvar(prog,dims,spaces,dom,degc,struct('psatz',codes(c)));
    end
    if numel(terms)==1,     Pcq = terms{1};    else,   Pcq = plus_batch(terms{:});  end
    % The one call sorts all names at once; a sum keeps the merge order of
    % its terms, so the lists are compared as sets and the rows aligned by
    % name in coef_diff.
    if ~isequal(sort(Pm.Zd(:)),sort(Ps.Zd(:))) || ~isequal(sort(Pm.Zd(:)),sort(Pcq.Zd(:)))
        error('test_poscopvar_direct: case ''%s'': the decision variables of the one call and the sums differ.',lbl);
    end
    [m1,ref] = coef_diff(Pm,Ps,vars);   [m2,~] = coef_diff(Pm,Pcq,vars);
    if ~(m1<=tol*max(ref,1)) || ~(m2<=tol*max(ref,1))
        error('test_poscopvar_direct: case ''%s'': one call differs from the sums, max |diff| %.3g / %.3g.',lbl,m1,m2);
    end
    if numel(codes)>1 && ~(iscell(Qm) && numel(Qm)==numel(codes))
        error('test_poscopvar_direct: case ''%s'': Qcell should hold one cell per term.',lbl);
    end
    fprintf('  passed: %-32s %d terms, ndec %6d, one call %.2fs (routes %s), max|diff| %.1e / %.1e vs direct sum / poscopvar sum\n',...
            lbl,numel(codes),numel(Pm.Zd),tm,strjoin(im.route,','),m1,m2);
    npass = npass+1;
end
fprintf('(6) several terms in one call = the sum of single calls on %d cases.\n\n',numel(cases));
end

fprintf('\ntest_poscopvar_direct passed (%d checks, parts %s).\n',npass,mat2str(parts));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function degc = reduced(deg,off,code,nv)
% The degree specification of one term at its own degrees: 'int' lowered
% by OFF, 'mult' kept, a 'joint' cap lowered per direction the weight acts
% in (as poscopvar_direct does inside the one call).
if code==0,     ng = 0;     elseif code==1,     ng = nv;    else,   ng = 1;     end
if iscell(deg)
    degc = cellfun(@(s) reduced(s,off,code,nv),deg,'UniformOutput',false);
    return
end
if off==0,  degc = deg;     return,     end
if isnumeric(deg),  deg = struct('int',deg,'mult',deg);     end
w = 1;  if isfield(deg,'int') && ~isempty(deg.int),     w = deg.int;    end
if ~isfield(deg,'mult') || isempty(deg.mult),   deg.mult = w;   end
degc = deg;     degc.int = max(w-off,0);
if isfield(deg,'joint') && ~isempty(deg.joint),     degc.joint = max(deg.joint-off*ng,0);    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mx,ref] = coef_diff(P1,P2,vars)
% Largest difference between the coefficients of two decision containers
% on the same decision variables, keyed by (block, cell, entries, output
% and input exponents) so that padding on one side meets zero on the other.
[X1,k1] = coef_rows(P1,vars);
[X2,k2] = coef_rows(P2,vars);
% The same decision variables in any order: rows of X2 in the order of P1.
if ~isequal(P1.Zd(:),P2.Zd(:))
    [tf,loc2] = ismember(P1.Zd(:),P2.Zd(:));
    if ~all(tf) || numel(P1.Zd)~=numel(P2.Zd),  mx = inf;   ref = 0;   return,     end
    X2 = X2(loc2,:);
end
[keys,~,loc] = unique([k1;k2],'rows');
n1 = size(k1,1);    nk = size(keys,1);
X1 = X1*sparse(1:n1,loc(1:n1),1,n1,nk);
X2 = X2*sparse(1:size(k2,1),loc(n1+1:end),1,size(k2,1),nk);
ref = full(max(max(abs(X1))));      if isempty(ref),    ref = 0;    end
if ~isequal(size(X1),size(X2)),     mx = inf;   return,     end
mx = full(max(max(abs(X1-X2))));    if isempty(mx),     mx = 0;     end
end


function [X,keys] = coef_rows(Pop,vars)
nv = numel(vars);
ngen = numel(Pop.Zd);
X = sparse(ngen,0);     keys = zeros(0,5+2*nv);
for i = 1:size(Pop.C,1)
    for j = 1:size(Pop.C,2)
        B = Pop.C{i,j};
        if isempty(B) || ~isa(B,'sdopvar'),     continue,   end
        m = B.dims(1);  n = B.dims(2);
        EL = multiindex_grid(B.ZL,'first_slowest');  NL = size(EL,1);
        ER = multiindex_grid(B.ZR,'first_slowest');  NR = size(ER,1);
        [~,po] = ismember(B.vars.out,vars);     [~,pi_] = ismember(B.vars.in,vars);
        eo = zeros(NL,nv);  eo(:,po) = EL;
        ei = zeros(NR,nv);  ei(:,pi_) = ER;
        nc = m*NL*n*NR;
        [ridx,cidx] = ind2sub([m*NL,n*NR],(1:nc)');
        p = floor((ridx-1)/NL)+1;   l = ridx-(p-1)*NL;
        r = floor((cidx-1)/NR)+1;   t = cidx-(r-1)*NR;
        for q = 1:numel(B.params.B)
            Bq = B.params.B{q};
            if nnz(Bq)==0,  continue,   end
            kq = [repmat([i,j,q],nc,1),p,r,eo(l,:),ei(t,:)];
            X = [X, Bq];                                                    %#ok<AGROW>
            keys = [keys; kq];                                              %#ok<AGROW>
        end
    end
end
end
