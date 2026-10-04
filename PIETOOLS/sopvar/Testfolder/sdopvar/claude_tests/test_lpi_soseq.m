function test_lpi_soseq()
% TEST_LPI_SOSEQ checks 'lpi_soseq' against two independent references:
% (1) SOSTOOLS soseq on the same equations, soseq(prog,dpvar(C,zeros(1,0),
%     {},Zd,[1,n])): every field of the program must be isequal (compared
%     field by field: isequal on the whole struct is false even against
%     itself) and the new entry's At, b, Z of the same class and sparsity. soseq is SOSTOOLS'
%     own implementation, not an inverse of lpi_soseq, and the test also
%     catches a later SOSTOOLS change that lpi_soseq would not follow.
% (2) The equations themselves: for random values d of the decision
%     variables, x(first table row of Zd{i}) = d(i), sossolve's residual
%     At'*x - b must be -(C(1,l) + sum_i C(1+i,l) d(i)) on every column l
%     of C that is not zero, and At must vanish on every other row.
%
% Inputs: an lpiprogram with lpidecvar names and the Gram names of a
% poscopvar (sosquadvar lists a symmetric Gram entry twice, so the first-
% row rule matters), random sparse and full C with zero columns, constant-
% only columns and unused rows, q = 0 and 1, an all-zero C, a name list
% repeating a name (one copy cancelling a column), a name absent from the
% program with a zero row (accepted) and with a nonzero row (an error),
% NaN entries, complex and logical C, and pos given by the caller.
%
% Initial coding MMP, 10/02/2026

warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
rng(20261002);
prog = lpiprogram(polynomial({'s'}),[],[0 1]);
prog = lpidecvar(prog,[6,1]);
[prog,~] = poscopvar(prog,1,{{'s'}},[0 1],1);
dvt = prog.decvartable(:);
[nm,first] = unique(dvt,'first');                   % distinct names, first rows
rep = nm(accumarray(ic_of(dvt,nm),1)>1);            % names listed twice
assert(~isempty(rep),'the program should list some Gram name twice');
nchk = 0;

% (a) random equations, pos computed and pos given
for trial = 1:60
    q = randi([0,min(14,numel(nm))]);
    if trial<=4,    q = trial-1;    end
    pick = randperm(numel(nm),q);
    if q>0 && mod(trial,3)==0                       % make sure repeats are hit
        pick(1) = find(strcmp(nm,rep{randi(numel(rep))}),1);
        pick = unique(pick,'stable');   q = numel(pick);
    end
    Zd = nm(pick);
    n = randi([1,9]);
    C = sprand(q+1,n,0.45);
    C(:,randi(n)) = 0;                              % a 0 = 0 equation
    if q>0 && n>1,  C(2:end,randi(n)) = 0;  end     % a constant-only equation
    if q>1,         C(1+randi(q),:) = 0;    end     % an unused name
    if mod(trial,5)==0,     C = full(C);    end
    nchk = nchk + check_one(prog,C,Zd,first(pick),sprintf('trial %d',trial));
end

% (b) all-zero C: soseq's own path
nchk = nchk + check_one(prog,sparse(3,4),nm(1:2),first(1:2),'all-zero C');

% (c) a repeated name, one copy cancelling a column: soseq's own path
Zd = [nm(1); nm(2); nm(1)];
C = sparse([0 1 2; 1 0 0; 0 3 0; -1 0 1]);         % column 1: d1 - d1 = 0
p_ref = soseq(prog,dpvar(C,zeros(1,0),{},Zd,[1,3]));
p_new = lpi_soseq(prog,C,Zd);
assert(same_prog(p_new,p_ref),'repeated name: program differs from soseq');
assert(size(p_ref.expr.At{end},2)==2,'repeated name: soseq should drop the cancelled column');
nchk = nchk+2;

% (d) a name absent from the program: zero row accepted, nonzero row an error
Zd = [nm(1); {'not_a_decision_variable'}];
C = sparse([1 2; 3 0; 0 0]);
nchk = nchk + check_one(prog,C,Zd,[],'absent name, zero row');
C(3,1) = 1;
assert(~isempty(errmsg(@() soseq(prog,dpvar(C,zeros(1,0),{},Zd,[1,2])))),'soseq should refuse');
assert(~isempty(errmsg(@() lpi_soseq(prog,C,Zd))),'lpi_soseq should refuse an absent name');
nchk = nchk+2;

% (d2) inputs soseq treats its own way, handed to it: NaN entries (any()
% ignores NaN, so a NaN-only row or column is dropped), complex C (b is
% conjugated), a non-double C (refused by the dpvar constructor)
C = sparse([1 NaN 2 0; NaN 0 0 0; 0 3 NaN 1; 0 0 0 NaN]);   % NaN-only row 2, column 4 below
C(1,4) = 0;     C(3,4) = 0;
p_ref = soseq(prog,dpvar(C,zeros(1,0),{},nm(1:3),[1,4]));
assert(same_prog(lpi_soseq(prog,C,nm(1:3)),p_ref) && same_prog(lpi_soseq(prog,C,nm(1:3),first(1:3)),p_ref),...
       'NaN entries: program differs from soseq');
C = sparse([1i 2; 1 0; 0 1-2i]);
p_ref = soseq(prog,dpvar(C,zeros(1,0),{},nm(1:2),[1,2]));
assert(same_prog(lpi_soseq(prog,C,nm(1:2),first(1:2)),p_ref),'complex C: program differs from soseq');
C = sparse(logical([1 0; 1 1]));
assert(~isempty(errmsg(@() soseq(prog,dpvar(C,zeros(1,0),{},nm(1),[1,2])))) ...
       && ~isempty(errmsg(@() lpi_soseq(prog,C,nm(1),first(1)))),'logical C: both should refuse');
nchk = nchk+3;

% (e) appended after earlier entries: num and the cell positions
p0 = lpi_soseq(prog,sparse([1 0; 2 1]),nm(3));
p_ref = soseq(p0,dpvar(sparse([0 1; 1 1]),zeros(1,0),{},nm(4),[1,2]));
p_new = lpi_soseq(p0,sparse([0 1; 1 1]),nm(4));
assert(same_prog(p_new,p_ref) && p_new.expr.num==prog.expr.num+2,'second entry differs');
nchk = nchk+1;

fprintf('test_lpi_soseq passed (%d checks).\n',nchk);
end


function k = check_one(prog,C,Zd,pos,lbl)
% Compare lpi_soseq (pos computed; and pos given, if any) with soseq and
% with the equations. Returns the number of checks made.
n = size(C,2);
p_ref = soseq(prog,dpvar(C,zeros(1,0),{},Zd,[1,n]));
cand = {lpi_soseq(prog,C,Zd)};
if ~isempty(pos) || isempty(Zd),    cand{end+1} = lpi_soseq(prog,C,Zd,pos(:));    end
k = 0;
for t = 1:numel(cand)
    p = cand{t};
    assert(same_prog(p,p_ref),'%s (%d): program differs from soseq',lbl,t);
    e = p.expr.num;
    for f = {'At','b','Z'}
        a = p.expr.(f{1}){e};   r = p_ref.expr.(f{1}){e};
        assert(issparse(a)==issparse(r) && strcmp(class(a),class(r)),...
               '%s (%d): %s class/sparsity differs',lbl,t,f{1});
    end
    k = k+2;
end
% The equations, from the definition (independent of soseq).
if nnz(C)==0,   return,     end
e = p_ref.expr.num;     At = p_ref.expr.At{e};  b = p_ref.expr.b{e};
kept = full(any(C,1));
inT = ~cellfun(@isempty,cellfun(@(z) find(strcmp(prog.decvartable,z),1),Zd(:),'UniformOutput',false));
d = randn(numel(Zd),1);     d(~inT) = 0;           % absent names have zero rows
x = zeros(numel(prog.decvartable),1);
for i = find(inT(:))'
    x(find(strcmp(prog.decvartable,Zd{i}),1)) = d(i);
end
res = At'*x - b;
expect = -(full(C(1,kept))' + full(C(2:end,kept))'*d);
assert(numel(res)==nnz(kept) && norm(res-expect,inf)<=1e-12*(1+norm(expect,inf)),...
       '%s: residual differs from the equations',lbl);
rows = unique(cellfun(@(z) find(strcmp(prog.decvartable,z),1),Zd(inT)));
other = true(size(At,1),1);     other(rows) = false;
assert(nnz(At(other,:))==0,'%s: At has entries off the first rows of Zd',lbl);
k = k+2;
end

function tf = same_prog(p,r)
% Field by field: isequal on the whole program is false even for
% isequal(prog,prog), because the polynomial overload of isequal (vartable)
% returns one logical per entry; compare that field with all().
tf = isequal(fieldnames(p),fieldnames(r)) && isequal(fieldnames(p.expr),fieldnames(r.expr));
for f = setdiff(fieldnames(p)','expr')
    t = isequal(p.(f{1}),r.(f{1}));
    tf = tf && isequal(size(p.(f{1})),size(r.(f{1}))) && all(t(:));
end
for f = fieldnames(p.expr)'         % isequaln: NaN entries compare equal
    tf = tf && isequal(size(p.expr.(f{1})),size(r.expr.(f{1}))) && isequaln(p.expr.(f{1}),r.expr.(f{1}));
end
end

function ic = ic_of(dvt,nm)
[~,ic] = ismember(dvt,nm);
end

function m = errmsg(f)
m = '';
try,        f();
catch ME,   m = ME.message;
end
end
