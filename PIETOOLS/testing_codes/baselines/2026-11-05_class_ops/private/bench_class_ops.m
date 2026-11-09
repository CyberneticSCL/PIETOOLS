function R = bench_class_ops(part,varargin)
% R = BENCH_CLASS_OPS(PART,...) times the standard class operations of
% 'opvar' against 'sopvar', and of 'dopvar' against 'sdopvar', on the same
% operators, and records the storage of each result.
% Operators act on L2^n[0,1] (3-PI: R0 multiplier, R1 lower, R2 upper kernel).
% Each operator is built once in the legacy class and converted once
% (OPVAR2SOPVAR, DOPVAR2SDOPVAR), so both classes hold the same operator.
% INPUT
% - part: 'check' semantic agreement of every operation at n = 2, d = 2;
%         'const' opvar vs sopvar over (n,d), varargin = {nlist, dlist};
%         'dvarB' dopvar vs sdopvar, each coefficient its own decision
%                 variable (the free-variable structure), over (n,d),
%                 varargin = {nlist, dlist};
%         'dvarQ' dopvar vs sdopvar at fixed n, d, each coefficient one
%                 decision variable drawn from q (nnz fixed as q grows),
%                 varargin = {qlist, n, d};
%         'dvarN' as 'dvarQ' with each coefficient a combination of
%                 ceil(4q/ncoef) variables (nnz about 4q),
%                 varargin = {qlist, n, d};
%         'dvarS' as 'dvarQ' with both operands on one shared list of
%                 decision variables, varargin = {qlist, n, d};
%         'inv'   inverse of I + 0.1 K, over (n,d), varargin = {nlist, dlist}.
% OUTPUT
% - R: struct array, one row per (case, operation, class), fields part, n,
%      d, q, nnz, op, cls, t (median seconds), nrep, bytes (whos of the
%      result), err (message if the operation failed), agree (relative
%      pointwise difference of the two results, 'check' and 'inv' only).
% Timing: a warm-up call, then the median of 3-25 calls filling about 0.3 s;
% an operation slower than 5 s is timed once.

% Initial coding MMP, 11/05/2026: baseline of the class operations before
% the opvar classes are retired (README.md in this folder).

pvar s s_dum
th = s_dum;     % the converters require the dummy variable '<primary>_dum'
switch part
    case 'check',   R = run_check(s,th);
    case 'const',   R = run_const(varargin{1},varargin{2},s,th);
    case 'dvarB',   R = run_dvar('dvarB',varargin{1},varargin{2},[],'own',s,th);
    case 'dvarQ',   R = run_dvar('dvarQ',varargin{2},varargin{3},varargin{1},'one',s,th);
    case 'dvarN',   R = run_dvar('dvarN',varargin{2},varargin{3},varargin{1},'many',s,th);
    case 'dvarS',   R = run_dvar('dvarS',varargin{2},varargin{3},varargin{1},'one',s,th);
    case 'inv',     R = run_inv(varargin{1},varargin{2},s,th);
    otherwise,      error('unknown part %s',part);
end
end

%% ---------------------------------------------------------------- parts
function R = run_const(nlist,dlist,s,th)
% opvar vs sopvar over (n,d): two independent random operators per case
R = [];
for n = nlist
    for d = dlist
        [A,S] = make_const(n,d,1,s,th);     [B,T] = make_const(n,d,2,s,th);
        ops = const_ops(n);
        for k = 1:size(ops,1)
            R = add(R,'const',n,d,0,0,ops{k,1},'opvar', @() ops{k,2}(A,B));
            R = add(R,'const',n,d,0,0,ops{k,1},'sopvar',@() ops{k,3}(S,T));
        end
        fprintf('##const n %d d %d done\n',n,d);
    end
end
end

function R = run_dvar(lbl,nlist,dlist,qlist,mode,s,th)
% 'dvarS': E on the same decision variable names as D (one shared list)
pfx = 'b';  if strcmp(lbl,'dvarS'),    pfx = 'a';  end
% dopvar vs sdopvar: D, E with disjoint decision variables, A/S constant
R = [];
if isempty(qlist),  qlist = NaN;  end
for n = nlist
    for d = dlist
        for q = qlist
            [D,SD,~,nnzD] = make_dvar(n,d,q,mode,1,s,th,'a');
            [E,SE] = make_dvar(n,d,q,mode,2,s,th,pfx);
            [A,S] = make_const(n,d,3,s,th);
            qq = numel(SD.Zd);
            ops = dvar_ops(n);
            for k = 1:size(ops,1)
                R = add(R,lbl,n,d,qq,nnzD,ops{k,1},'dopvar', @() ops{k,2}(D,E,A));
                R = add(R,lbl,n,d,qq,nnzD,ops{k,1},'sdopvar',@() ops{k,3}(SD,SE,S));
            end
            fprintf('##%s n %d d %d q %d nnz %d done\n',lbl,n,d,qq,nnzD);
        end
    end
end
end

function R = run_inv(nlist,dlist,s,th)
% inverse of I + 0.1 K; agreement = sup |P*Pinv*f - f| / sup|f| per class
R = [];
for n = nlist
    for d = dlist
        rng(4);
        A = rand_opvar([0 0;n n],d,s,th,[0 1]);
        A.R.R0 = 2*eye(n) + 0.1*A.R.R0;   A.R.R1 = 0.1*A.R.R1;   A.R.R2 = 0.1*A.R.R2;
        S = opvar2sopvar(A);
        R = add(R,'inv',n,d,0,0,'inv','opvar', @() inv(A));
        R(end).agree = inv_resid(@(f) op_apply(A,op_apply(R(end).out,f,s,th),s,th),n,s);
        R = add(R,'inv',n,d,0,0,'inv','sopvar',@() inv(S));
        R(end).agree = inv_resid(@(f) apply_sopvar(S,apply_sopvar(R(end).out,f)),n,s);
        R(end-1).out = [];  R(end).out = [];
        fprintf('##inv n %d d %d done\n',n,d);
    end
end
end

function R = run_check(s,th)
% Semantic agreement of every operation at n = 2, d = 2
n = 2;  d = 2;  R = [];
[A,S] = make_const(n,d,1,s,th);     [B,T] = make_const(n,d,2,s,th);
ops = const_ops(n);
for k = 1:size(ops,1)
    Po = ops{k,2}(A,B);     Ps = ops{k,3}(S,T);
    e = cmp_apply(Po,Ps,s,th);
    R = [R; struct('part','check','op',ops{k,1},'cls','opvar/sopvar','agree',e)]; %#ok<AGROW>
    fprintf('##check const %-12s %.2e\n',ops{k,1},e);
end
[D,SD] = make_dvar(n,d,NaN,'own',1,s,th,'a');   [E,SE] = make_dvar(n,d,NaN,'own',2,s,th,'b');
[A,S] = make_const(n,d,3,s,th);
nm = cellstr(string([SD.Zd(:); SE.Zd(:)]));     rng(7);  vals = randn(numel(nm),1);
% the evaluation itself, on the inputs
e0 = cmp_apply(dop_eval(D,nm,vals,s,th),sdop_eval(SD,nm,vals),s,th);
fprintf('##check dvar  %-12s %.2e\n','input',e0);
ops = dvar_ops(n);
for k = 1:size(ops,1)
    Po = ops{k,2}(D,E,A);   Ps = ops{k,3}(SD,SE,S);
    e = cmp_apply(dop_eval(Po,nm,vals,s,th),sdop_eval(Ps,nm,vals),s,th);
    R = [R; struct('part','check','op',ops{k,1},'cls','dopvar/sdopvar','agree',e)]; %#ok<AGROW>
    fprintf('##check dvar  %-12s %.2e\n',ops{k,1},e);
end
end

%% ---------------------------------------------------------------- operations
function ops = const_ops(n)
% {name, opvar handle, sopvar handle}
r = 1;  if n>1, r = 1:floor(n/2); end
ops = {'plus',      @(A,B) A+B,         @(S,T) S+T;
       'minus',     @(A,B) A-B,         @(S,T) S-T;
       'uminus',    @(A,B) -A,          @(S,T) -S;
       'scalar',    @(A,B) 2.5*A,       @(S,T) 2.5*S;
       'mtimes',    @(A,B) A*B,         @(S,T) S*T;
       'ctranspose',@(A,B) A',          @(S,T) S';
       'horzcat',   @(A,B) [A B],       @(S,T) [S T];
       'vertcat',   @(A,B) [A;B],       @(S,T) [S;T];
       'blkdiag',   @(A,B) blkdiag(A,B),@(S,T) blkdiag(S,T);
       'subsref',   @(A,B) A(r,:),      @(S,T) S(r,:)};
end

function ops = dvar_ops(n)
% {name, dopvar handle (D,E,A), sdopvar handle (SD,SE,S)}
r = 1;  if n>1, r = 1:floor(n/2); end
ops = {'plus',      @(D,E,A) D+E,          @(D,E,A) D+E;
       'plus_const',@(D,E,A) D+A,          @(D,E,A) D+A;
       'minus',     @(D,E,A) D-E,          @(D,E,A) D-E;
       'uminus',    @(D,E,A) -D,           @(D,E,A) -D;
       'scalar',    @(D,E,A) 2.5*D,        @(D,E,A) 2.5*D;
       'mtimes_DA', @(D,E,A) D*A,          @(D,E,A) D*A;
       'mtimes_AD', @(D,E,A) A*D,          @(D,E,A) A*D;
       'ctranspose',@(D,E,A) D',           @(D,E,A) D';
       'horzcat',   @(D,E,A) [D E],        @(D,E,A) [D E];
       'vertcat',   @(D,E,A) [D;E],        @(D,E,A) [D;E];
       'blkdiag',   @(D,E,A) blkdiag(D,E), @(D,E,A) blkdiag(D,E);
       'subsref',   @(D,E,A) D(r,:),       @(D,E,A) D(r,:)};
end

%% ---------------------------------------------------------------- builders
function [A,S] = make_const(n,d,seed,s,th)
% Random 3-PI opvar on L2^n[0,1] (RAND_OPVAR) and its sopvar
rng(seed);
A = rand_opvar([0 0;n n],d,s,th,[0 1]);
S = opvar2sopvar(A);
end

function [D,SD,names,nnzD] = make_dvar(n,d,q,mode,seed,s,th,prefix)
% Synthetic 3-PI dopvar on L2^n[0,1]: R0 with monomials s^0..s^d, R1 and R2
% with s^i th^j, i,j = 0..d. Each coefficient is a random constant plus a
% decision part: 'own' its own variable (q = number of coefficients), 'one'
% one variable drawn from q, 'many' ceil(4q/ncoef) variables drawn from q.
% Built in the dpvar layout: rows = (matrix row, [1; decision variables])
% with the decision variable inner, columns = (matrix column, monomial)
% with the monomial inner.
rng(seed);
dm0 = (0:d)';
[e1,e2] = ndgrid(0:d,0:d);  dm12 = [e1(:), e2(:)];
nZ = [numel(dm0), size(dm12,1), size(dm12,1)];
nc = sum(n*n*nZ);
if strcmp(mode,'own'),  q = nc;     end
sper = max(1,ceil(4*q/nc));
names = arrayfun(@(k) sprintf('%s%d',prefix,k),(1:q)','UniformOutput',false);
F = cell(1,3);  next = 0;   nnzD = 0;
for f = 1:3
    nz = nZ(f);
    [kk,ii,jj] = ndgrid(1:nz,1:n,1:n);
    ii = ii(:);     jj = jj(:);     kk = kk(:);     m = numel(ii);
    col = (jj-1)*nz + kk;
    rowc = (ii-1)*(q+1) + 1;
    switch mode
        case 'own',  lv = next + (1:m)';   next = next + m;    rr = rowc + lv;    cc = col;   vv = randn(m,1);
        case 'one',  lv = randi(q,m,1);    rr = rowc + lv;     cc = col;          vv = randn(m,1);
        case 'many', lv = randi(q,m,sper); rr = repmat(rowc,1,sper) + lv;  cc = repmat(col,1,sper);  vv = randn(m,sper);
    end
    C = sparse([rowc; rr(:)],[col; cc(:)],[rand(m,1); vv(:)],n*(q+1),n*nz);
    nnzD = nnzD + nnz(C);
    if f==1,    F{f} = dpvar(C,dm0,{'s'},names,[n n]);
    else,       F{f} = dpvar(C,dm12,{'s','s_dum'},names,[n n]);
    end
end
D = dopvar();
D.I = [0 1];    D.var1 = s;     D.var2 = th;
D.P = zeros(0,0);   D.Q1 = zeros(0,n);  D.Q2 = zeros(n,0);
D.R.R0 = F{1};  D.R.R1 = F{2};  D.R.R2 = F{3};
SD = dopvar2sdopvar(D);
end

%% ---------------------------------------------------------------- semantics
function g = op_apply(A,f,s,th)
% (A f)(s) = R0(s) f(s) + int_a^s R1(s,th) f(th) dth + int_s^b R2(s,th) f(th) dth
I = A.I;    ft = subs(f,s,th);
g = A.R.R0*f + int(A.R.R1*ft,th,I(1),s) + int(A.R.R2*ft,th,s,I(2));
end

function P = dop_eval(D,nm,vals,s,th)
% The dopvar D at the decision values VALS of the names NM: an opvar
P = opvar();    P.I = D.I;  P.var1 = D.var1;    P.var2 = D.var2;
[mo,ni] = size(D.R.R0);
P.P = zeros(0,0);   P.Q1 = zeros(0,ni);     P.Q2 = zeros(mo,0);
P.R.R0 = fld(D.R.R0,nm,vals);   P.R.R1 = fld(D.R.R1,nm,vals);   P.R.R2 = fld(D.R.R2,nm,vals);
end
function p = fld(F,nm,vals)
if isa(F,'dpvar'),  p = pi_dpvar_eval(F,nm,vals);   else,   p = F;  end
end

function P = sdop_eval(SD,nm,vals)
% The sdopvar SD at the decision values: C_gam = A_gam + B_gam' v (sdopvar
% definition), v ordered as SD.Zd
[tf,loc] = ismember(cellstr(string(SD.Zd(:))),nm);
v = zeros(numel(tf),1);     v(tf) = vals(loc(tf));
m = SD.dims(1)*prod([cellfun(@numel,SD.ZL),1]);
n = SD.dims(2)*prod([cellfun(@numel,SD.ZR),1]);
A = SD.params.A;    B = SD.params.B;    pC = cell(size(A));
for g = 1:numel(A)
    a = A{g};   b = [];     if g<=numel(B), b = B{g};   end
    if isempty(a) || (isscalar(a) && a==0),     a = sparse(m*n,1);  end
    w = a(:);
    if ~isempty(b) && nnz(b)>0,     w = w + (v.'*b).';  end
    pC{g} = reshape(w,m,n);
end
P = sopvar(pC,SD.vars,SD.ZL,SD.ZR,SD.dom,SD.dims);
end

function e = cmp_apply(Po,Ps,s,th)
% sup over 3 random polynomial inputs and 60 points of |Po f - Ps f| / sup|Ps f|
rng(11);    nin = Po.dim(2,2);  pts = rand(1,60);   num = 0;    den = 0;
for t = 1:3
    f = polynomial(zeros(nin,1));
    for i = 1:nin,  f(i) = randn + randn*s + randn*s^2;    end
    g1 = double(subs(op_apply(Po,f,s,th),s,pts));
    g2 = double(subs(apply_sopvar(Ps,f),s,pts));
    num = max(num,max(abs(g1(:)-g2(:))));     den = max(den,max(abs(g2(:))));
end
e = num/max(den,eps);
end

function e = inv_resid(fh,n,s)
% sup over 3 random polynomial inputs of |P Pinv f - f| / sup|f|
rng(12);    pts = rand(1,60);   num = 0;    den = 0;
for t = 1:3
    f = polynomial(zeros(n,1));
    for i = 1:n,    f(i) = randn + randn*s + randn*s^2;    end
    g = double(subs(fh(f),s,pts));     f0 = double(subs(f,s,pts));
    num = max(num,max(abs(g(:)-f0(:))));     den = max(den,max(abs(f0(:))));
end
e = num/den;
end

%% ---------------------------------------------------------------- timing
function R = add(R,part,n,d,q,nz,op,cls,fh)
% Time FH (warm-up, then the median of 3-25 calls filling ~0.3 s) and record it
row = struct('part',part,'n',n,'d',d,'q',q,'nnz',nz,'op',op,'cls',cls, ...
             't',NaN,'nrep',0,'bytes',NaN,'err','','agree',NaN,'out',[]);
try
    tic;    out = fh();     t1 = toc;
    if t1>5
        row.t = t1;     row.nrep = 1;
    else
        tic;    out = fh();     t2 = toc;
        nrep = max(3,min(25,ceil(0.3/max(t2,1e-6))));  ts = zeros(1,nrep);
        for k = 1:nrep,     tic;    out = fh();     ts(k) = toc;    end
        row.t = median(ts);     row.nrep = nrep;
    end
    w = whos('out');    row.bytes = w.bytes;    row.out = out;
catch ME
    row.err = ME.message;
end
if ~strcmp(part,'inv'),     row.out = [];   end
R = [R; row];
end
