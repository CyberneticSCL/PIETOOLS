function results = test_high_dim_pde(dims,names,nrep)
% TEST_HIGH_DIM_PDE Tests the PDE parser and converter on systems with 4 and
% 5 spatial variables (those of HIGH_DIM_SYSTEMS), against the PDE written
% by hand rather than against the routines that build the PIE.
%
% Each operator of the PIE is applied, block by block, to polynomials by
% its kernels (apply_sopvar), exactly, and compared with the PDE:
%  (A) for a PIE with constraints K, on states built to satisfy the BCs:
%      T xf = x, A xf = right-hand side, K xf = 0, the outputs, and the map
%      pde2pie.x.C, where xf = F D^alpha x for each block;
%  (B) for a PIE without K, on random fundamental states xf and inputs
%      compatible with the BCs: x = T xf + Tw wf + Tu uf satisfies the BCs,
%      D^alpha x = xf, A xf + Bw wf + Bu uf = right-hand side at x, and the
%      outputs. The fundamental inputs wf, uf are built from pde2pie.w/u.
% A system that must be rejected passes if the conversion errors with the
% expected message. Prints one line per system, with the time of the
% conversion and of the checks.
%
% Run as: results = test_high_dim_pde;  or  results = test_high_dim_pde(5);
%
% INPUT
% - dims:   (optional) array, the dimensions of the systems to run, from 4
%           and 5 (default 4: about 5 min). A 5D system with dynamics takes
%           1-3 min to convert, and about 3.5 GB of memory; it is skipped if
%           less than 4 GB is available;
% - names:  (optional) cellstr, the systems to run (overrides dims);
% - nrep:   (optional) positive integer, the number of random test
%           functions per check (default 1). Random fundamental states have
%           degree 2 in each variable in 4D, 1 in 5D;
%
% OUTPUT
% - results: 'struct' array with fields name, dim, status ('pass',
%            'rejected', 'skipped' or 'FAIL'), passed (logical; true for
%            pass and rejected), t_convert, t_check (seconds) and message.
%
% NOTES
% The folder must be on the path: testing_codes/private/ is off the
% default path, e.g.
%   addpath(fullfile(<PIETOOLS>,'testing_codes','private','converters','tests','high_dim'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - test_high_dim_pde
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding DJ - 11/25/2026

if nargin<1 || isempty(dims),  dims = 4;  end
if nargin<3 || isempty(nrep),  nrep = 1;  end
catl = high_dim_systems('list');
if nargin<2 || isempty(names)
    names = {catl(ismember([catl.dim],dims)).name};
end
warning('off','sopvar:noncanonicalMultiplier');
tol = 1e-9;
results = struct('name',{},'dim',{},'status',{},'passed',{},'t_convert',{},'t_check',{},'message',{});
for i=1:numel(names)
    S = high_dim_systems(names{i});
    r = struct('name',names{i},'dim',S.dim,'status','','passed',false,'t_convert',NaN,'t_check',NaN,'message','');
    % 5D with dynamics: skip if memory is short (memory is Windows only)
    if S.dim==5 && ~isempty(S.rhs)
        try [~,sv] = memory;  avail = sv.PhysicalMemory.Available;
        catch, avail = Inf;
        end
        if avail<4e9
            r.status = 'skipped';  r.message = sprintf('%.1f GB of memory available, 4 GB needed.',avail/1e9);
            results(end+1) = report(r);                                     %#ok<AGROW>
            continue
        end
    end
    % Convert, and classify a conversion error
    t0 = tic;
    try
        PIE = convert(S.PDE,'silent');
        r.t_convert = toc(t0);
    catch ME
        r.t_convert = toc(t0);
        if ~isempty(S.expect) && contains(ME.message,S.expect)
            r.status = 'rejected';  r.passed = true;
        else
            r.status = 'FAIL';
        end
        r.message = ME.message;
        results(end+1) = report(r);                                         %#ok<AGROW>
        continue
    end
    if ~isempty(S.expect)
        r.status = 'FAIL';  r.message = 'converted, but must be rejected';
        results(end+1) = report(r);                                         %#ok<AGROW>
        continue
    end
    blk = arrayfun(@(b) sprintf('%dD',nnz(PIE.T.space_in(b,:))),1:size(PIE.T.space_in,1),'UniformOutput',false);
    r.message = sprintf('xf blocks [%s]',strjoin(blk,' '));
    if ~S.check
        % Only the conversion: here a 3D state and its 4D history
        r.passed = isequal(sort(cellfun(@(c) sscanf(c,'%dD'),blk)),[3 4]);
        r.status = ternary(r.passed,'pass','FAIL');
        results(end+1) = report(r);                                         %#ok<AGROW>
        continue
    end
    % Semantic checks
    t0 = tic;
    try
        e = check_pie(S,PIE,6-S.dim,nrep);
        r.t_check = toc(t0);
        v = struct2cell(e);  v = [v{:}];
        r.passed = all(v(~isnan(v))<tol);
        r.status = ternary(r.passed,'pass','FAIL');
        f = fieldnames(e);  ok = ~isnan(v);
        r.message = [r.message,' | ',strjoin(cellfun(@(n,x) sprintf('%s %.1e',n,x),f(ok)',num2cell(v(ok)),'UniformOutput',false),', ')];
    catch ME
        r.t_check = toc(t0);
        r.status = 'FAIL';  r.message = sprintf('%s (%s:%d)',ME.message,ME.stack(1).name,ME.stack(1).line);
    end
    results(end+1) = report(r);                                             %#ok<AGROW>
end
st = {results.status};
fprintf('\n%d systems: %d pass, %d rejected as required, %d skipped, %d FAIL.\n',numel(results), ...
    nnz(strcmp(st,'pass')),nnz(strcmp(st,'rejected')),nnz(strcmp(st,'skipped')),nnz(strcmp(st,'FAIL')));
end



%% ------------------------------------------------------------------------
function e = check_pie(S,PIE,deg,nrep)
% E = CHECK_PIE(S,PIE,DEG,NREP) checks the PIE of system S (see the help of
% the file), with random polynomials of degree DEG in each variable.
%
% INPUT
% - S:      'struct' of HIGH_DIM_SYSTEMS;
% - PIE:    'pie_struct', the converted PIE;
% - deg:    nonnegative integer;
% - nrep:   positive integer, the number of random test functions;
%
% OUTPUT
% - e:      'struct' of relative errors (NaN if not checked): T, A, K, Z, Y,
%           C on states satisfying the BCs, and BC, D, A_rand, Z_rand,
%           Y_rand on random fundamental states.
T = PIE.T;  vreg = T.vars;  dom = T.dom;
pvars = PIE.space.vars;                    % column order of pde2pie.*.diff
idx = PIE.pde2pie.x.idx(:,1);  dmat = PIE.pde2pie.x.diff;  nb = numel(idx);
% The declared state of each block, by ID; an auxiliary block lacks some of
% the variables of its state, and holds its integral over them
myx = arrayfun(@(r) find(S.xid==PIE.space_pde.x.ID(r)),idx);
bvars = cell(nb,1);  drop = cell(nb,1);
for b=1:nb
    bvars{b} = vreg(T.space_in(b,:));
    drop{b} = setdiff(S.svars{myx(b)},bvars{b});
end
hasK = ~isempty(PIE.K) && ~isempty(PIE.K.C);  isBVP = isempty(S.rhs);
% The declared output of each PIE output, by ID
outs = {'z','Cz','Dzw','Dzu',S.zid,S.z; 'y','Cy','Dyw','Dyu',S.yid,S.y};
myo = cell(2,1);
for k=1:2
    if isfield(PIE.pde2pie,outs{k,1}) && ~isempty(PIE.pde2pie.(outs{k,1}).idx)
        myo{k} = arrayfun(@(r) find(outs{k,5}==PIE.space_pde.(outs{k,1}).ID(r)),PIE.pde2pie.(outs{k,1}).idx(:));
    end
end
e = struct('T',NaN,'A',NaN,'K',NaN,'Z',NaN,'Y',NaN,'C',NaN,'BC',NaN,'D',NaN,'A_rand',NaN,'Z_rand',NaN,'Y_rand',NaN);
if hasK
    % (A) states built to satisfy the BCs
    [e.T,e.A,e.K,e.C] = deal(0);
    for rep=1:nrep
        rng(100+rep);  c = randn(1,20);
        W = S.wgen(c);  U = S.ugen(c);  X = S.xs(c);
        bres = S.bc(X,W,U);
        for j=1:numel(bres)
            if pmax(bres{j})>1e-12*max(1,maxc(X)), error('The hand-built state violates BC %d (test error).',j); end
        end
        wf = fund_in(PIE,'w',W,S.wid);  uf = fund_in(PIE,'u',U,S.uid);
        if ~isBVP, F = S.rhs(X,W,U); end
        xf = cell(nb,1);
        for b=1:nb, xf{b} = Fop(Dop(X{myx(b)},pvars,dmat(b,:)),drop{b},vreg,dom); end
        Txf = lin(PIE,{'T','Tw','Tu'},xf,wf,uf,nb);
        Kxf = lin(PIE,{'K','Kw','Ku'},xf,wf,uf,sum(PIE.K.dim_out));
        if ~isBVP, Axf = lin(PIE,{'A','Bw','Bu'},xf,wf,uf,nb); end
        Cm = PIE.pde2pie.x.C;
        for b=1:nb
            k = myx(b);
            e.T = max(e.T,rel(Txf{b},Fop(X{k},drop{b},vreg,dom)));
            if ~isBVP, e.A = max(e.A,rel(Axf{b},Fop(F{k},drop{b},vreg,dom))); end
            % pde2pie.x.C maps D^alpha x to the block
            v = arrayfun(@(kk) polynomial(zeros(Cm.dim_in(kk),1)),1:size(Cm.C,2),'UniformOutput',false);
            v{b} = Dop(X{k},pvars,dmat(b,:));
            y = apply_cop(Cm,v);
            e.C = max(e.C,rel(y{b},xf{b}));
        end
        for j=1:numel(Kxf), e.K = max(e.K,pmax(Kxf{j})/max(1,maxc(xf))); end
        if isBVP, e.A = NaN; end
        for k=1:2
            if isempty(myo{k}), continue, end
            Zr = outs{k,6}(X,W,U);
            Zx = lin(PIE,outs(k,2:4),xf,wf,uf,numel(myo{k}));
            fn = upper(outs{k,1});
            if isnan(e.(fn)), e.(fn) = 0; end
            for j=1:numel(myo{k}), e.(fn) = max(e.(fn),rel(Zx{j},Zr{myo{k}(j)})); end
        end
    end
    return
end
% (B) random fundamental states and compatible inputs
[e.BC,e.D] = deal(0);
if ~isBVP, e.A_rand = 0; end
for rep=1:nrep
    rng(200+rep);  c = randn(1,20);
    xf = cell(nb,1);
    for b=1:nb, xf{b} = randpoly(bvars{b},deg,T.dim_in(b)); end
    W = S.wgen(c);  U = S.ugen(c);
    wf = fund_in(PIE,'w',W,S.wid);  uf = fund_in(PIE,'u',U,S.uid);
    x = lin(PIE,{'T','Tw','Tu'},xf,wf,uf,nb);
    X = cell(1,numel(S.svars));
    for b=1:nb, X{myx(b)} = x{b}; end
    bres = S.bc(X,W,U);
    for j=1:numel(bres), e.BC = max(e.BC,pmax(bres{j})/max(1,maxc(x))); end
    for b=1:nb, e.D = max(e.D,rel(Dop(X{myx(b)},pvars,dmat(b,:)),xf{b})); end
    if ~isBVP
        Axf = lin(PIE,{'A','Bw','Bu'},xf,wf,uf,nb);  F = S.rhs(X,W,U);
        for b=1:nb, e.A_rand = max(e.A_rand,rel(Axf{b},F{myx(b)})); end
    end
    for k=1:2
        if isempty(myo{k}), continue, end
        Zr = outs{k,6}(X,W,U);
        Zx = lin(PIE,outs(k,2:4),xf,wf,uf,numel(myo{k}));
        fn = [upper(outs{k,1}),'_rand'];
        if isnan(e.(fn)), e.(fn) = 0; end
        for j=1:numel(myo{k}), e.(fn) = max(e.(fn),rel(Zx{j},Zr{myo{k}(j)})); end
    end
end
end

function wf = fund_in(PIE,typ,W,myids)
% WF = FUND_IN(PIE,TYP,W,MYIDS) returns the fundamental inputs of type TYP
% ('w' or 'u') for the PDE inputs W: for each PIE input i, the PDE input
% with ID PIE.space_pde.(TYP).ID(idx(i)), differentiated diff(i,:) times,
% evaluated at loc(i,:) (-1 lower, 1 upper boundary), and mapped by C.
%
% INPUT
% - PIE:    'pie_struct';
% - typ:    'w' or 'u';
% - W:      1 x n cell of 'polynomial', the PDE inputs, in the order of
% - myids:  1 x n array, their IDs;
%
% OUTPUT
% - wf:     cell of 'polynomial', one per PIE input ({} if none).
wf = {};
if ~isfield(PIE.pde2pie,typ) || isempty(PIE.pde2pie.(typ).idx), return, end
p2p = PIE.pde2pie.(typ);  ids = PIE.space_pde.(typ).ID;
pvars = PIE.space.vars;  pdom = PIE.space.dom;
v = cell(numel(p2p.idx),1);
for i=1:numel(v)
    q = Dop(W{myids==ids(p2p.idx(i))},pvars,p2p.diff(i,:));
    for k=find(p2p.loc(i,:))
        q = subs(q,polynomial(pvars(k)),pdom(k,(p2p.loc(i,k)+3)/2));
    end
    v{i} = q;
end
if isfield(p2p,'C') && ~isempty(p2p.C) && ~isempty(p2p.C.C)
    wf = apply_cop(p2p.C,v);
else
    wf = v;
end
end

function y = lin(PIE,ops,xf,wf,uf,M)
% Y = LIN(PIE,OPS,XF,WF,UF,M) returns P*xf + Pw*wf + Pu*uf, with
% {P,Pw,Pu} the fields OPS of PIE, as a cell of M 'polynomial' blocks.
y = apply_cop(PIE.(ops{1}),xf,M);
if ~isempty(wf), y = cellfun(@plus,y,apply_cop(PIE.(ops{2}),wf,M),'UniformOutput',false); end
if ~isempty(uf), y = cellfun(@plus,y,apply_cop(PIE.(ops{3}),uf,M),'UniformOutput',false); end
end

function y = apply_cop(P,v,M)
% Y = APPLY_COP(P,V,M) applies the 'copvar' P to the cell V of 'polynomial'
% blocks, each block by its kernels: y{i} = sum_j P.C{i,j} v{j}. An absent
% operator gives M scalar zeros.
if isempty(P) || ~isa(P,'copvar') || isempty(P.C)
    if nargin<3, M = 0; end
    y = repmat({polynomial(0)},M,1);
    return
end
y = cell(size(P.C,1),1);
for i=1:size(P.C,1)
    acc = polynomial(zeros(P.dim_out(i),1));
    for j=1:size(P.C,2)
        if ~isempty(P.C{i,j}), acc = acc+apply_sopvar(P.C{i,j},polynomial(v{j})); end
    end
    y{i} = acc;
end
end

function y = Dop(p,pvars,d)
% Y = DOP(P,PVARS,D) differentiates P D(j) times in the variable PVARS{j}.
y = polynomial(p);
for j=1:numel(d)
    for i=1:d(j), y = diff(y,polynomial(pvars(j))); end
end
end

function y = Fop(p,drop,vreg,dom)
% Y = FOP(P,DROP,VREG,DOM) integrates P over the variables DROP, each on its
% interval (DOM rows ordered as VREG).
y = polynomial(p);
for j=1:numel(drop)
    kv = find(strcmp(vreg,drop{j}));
    y = int(y,polynomial(drop(j)),dom(kv,1),dom(kv,2));
end
end

function e = rel(a,b)
% E = REL(A,B) returns max|coef(A-B)| / max(1,max|coef(B)|).
e = pmax(polynomial(a)-polynomial(b))/max(1,pmax(b));
end

function m = pmax(p)
% M = PMAX(P) returns the largest absolute coefficient of the 'polynomial' P.
c = polynomial(p).coefficient;
if isempty(c), m = 0; else, m = full(max(abs(c(:)))); end
end

function m = maxc(X)
% M = MAXC(X) returns the largest PMAX over the cell X.
m = max([0;cellfun(@pmax,X(:))]);
end

function p = randpoly(vars,deg,m)
% P = RANDPOLY(VARS,DEG,M) returns an M x 1 'polynomial' in the variables
% VARS (cellstr), of degree at most DEG in each, with normal coefficients.
mons = {polynomial(1)};
for j=1:numel(vars)
    v = polynomial(vars(j));  new = cell(1,numel(mons)*(deg+1));
    for i=1:numel(mons)
        for q=0:deg, new{(i-1)*(deg+1)+q+1} = mons{i}*v^q; end
    end
    mons = new;
end
p = polynomial(zeros(m,1));
for i=1:numel(mons), p = p+randn(m,1)*mons{i}; end
end

function r = report(r)
% R = REPORT(R) prints the result R and returns it.
fprintf('%-14s %dD  %-8s  convert %7.2f s  check %7.2f s  %s\n',r.name,r.dim,r.status,r.t_convert,r.t_check,r.message);
end

function s = ternary(c,a,b)
% S = TERNARY(C,A,B) returns A if C, else B.
if c, s = a; else, s = b; end
end
