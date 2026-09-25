%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_LPIVAR_CDOPVAR checks that 'lpivar_cdopvar' declares the operator
% family it claims: every admissible kernel monomial of every block and
% gamma cell carries exactly one free decision variable, and nothing else.
%
% Against the stock 'lpivar' (1-D, square and rectangular, several degrees):
%   COUNT    numel(Pop.Zd) = the variables lpivar declares = the closed form
%            n11*n12 + (d1+1)(n11*n22 + n21*n12 + n21*n22)
%                    + 2(d3+1)(d2+1) n21*n22;
%   SUPPORT  at random decision values, the kernel monomials of each slot
%            P, Q1, Q2, R0, R1, R2 coincide with the stock operator's.
% Both sides are evaluated independently of lpivar_cdopvar's own layout:
% the container through 'pi_blk_kernels' (each class's kernel definition),
% the stock operator through 'dpvar2poly' and polynomial 'subs'. Equal count
% and equal support, with one variable per monomial, give equal families.
%
% Beyond 1-D (2-D, 3-D, cross, rectangular, zero blocks), where there is no
% stock counterpart, against the specification itself:
%   every kernel monomial respects its role's cap (mult, int, out, in); a
%   multiplier direction carries no dummy variable (canonical form); the
%   support is FULL, so the count equals the number of admissible monomial
%   entries and the map from variables to operators is injective; and the
%   constructor had nothing to fold (no noncanonical warning).
%
% MMP, 09/25/2026: Initial coding

rng(20260926);
nchk = 0;

% ------------------------------------------------------------ vs lpivar, 1-D
DIMS = { [2 2;3 3], [1 1;1 1], [1 2;0 1], [0 1;2 0], [2 0;1 2], [0 0;2 3] };
DEGS = { [1 1 1], [2 1 3], [0 2 1], 2 };
for id = 1:numel(DIMS)
    for ig = 1:numel(DEGS)
        n = DIMS{id};   d = DEGS{ig};
        if isscalar(d),     d3 = [d d d];   else,   d3 = d;     end
        prog = lpiprogram(polynomial({'s1'}),[],[0,1]);
        nd0 = numel(prog.decvartable);
        [prog,Zop] = lpivar(prog,n,d);
        nleg = numel(prog.decvartable) - nd0;
        nform = n(1,1)*n(1,2) + (d3(1)+1)*(n(1,1)*n(2,2) + n(2,1)*n(1,2) + n(2,1)*n(2,2)) ...
                + 2*(d3(3)+1)*(d3(2)+1)*n(2,1)*n(2,2);

        [sp,dm,rmap,cmap] = spaces_1d(n);
        prog2 = lpiprogram(polynomial({'s1'}),[],[0,1]);
        [~,Pc] = lpivar_cdopvar(prog2,dm,sp,[0,1],d);
        assert(numel(Pc.Zd)==nleg && nleg==nform, ...
            'count: container %d, lpivar %d, formula %d for n=%s d=%s', ...
            numel(Pc.Zd),nleg,nform,mat2str(n),mat2str(d));

        % Support, slot by slot. rmap/cmap give the container row/column of
        % the R (1) and L2 (2) parts, 0 where that part is absent.
        leg = legacy_supports(Zop);
        con = container_supports(Pc,rmap,cmap);
        slots = fieldnames(leg);
        for k = 1:numel(slots)
            assert(isequal(leg.(slots{k}),con.(slots{k})), ...
                'support of %s differs for n=%s d=%s',slots{k},mat2str(n),mat2str(d));
        end
        [ok,msg] = unit_responses(Pc,25);
        assert(ok,'n=%s d=%s: %s',mat2str(n),mat2str(d),msg);
        nchk = nchk + 2 + numel(slots);
    end
end
fprintf('  passed: 1-D against lpivar, %d dims x %d degrees (%d checks)\n', ...
        numel(DIMS),numel(DEGS),nchk);

% ------------------------------------------------------------ n-D, by spec
CASES = {
  % label, spaces, dims, deg, occ
  '2D square R x L2[s] x L2[s,t]', ...
      {{},{'s'},{'s','t'}}, [1;2;1], struct('mult',2,'int',[1 2],'out',1,'in',2), []
  '3D square R x L2[s,t,u]', ...
      {{},{'s','t','u'}}, [1;1], 1, []
  'cross L2[s] <-> L2[t]', ...
      {{'s'},{'t'}}, [2;1], struct('mult',1,'int',[2 1],'out',2,'in',1), []
  'rectangular R^2 x L2[s] -> L2[s,t] x R', ...
      struct('out',{{ {'s','t'}, {} }},'in',{{ {}, {'s'} }}), ...
      struct('out',[1;2],'in',[2;1]), 2, logical([1 1; 0 1])
  'per-block degrees', ...
      {{'s'},{'s','t'}}, [1;1], {1, [2 1 3]; struct('int',[0 3]), 0}, []
};
for ic = 1:size(CASES,1)
    [lbl,sp,dm,dg,occ] = deal(CASES{ic,:});
    vars = unique([cat_spaces(sp,'out'), cat_spaces(sp,'in')]);
    prog = make_prog(vars,repmat([0 1],numel(vars),1));
    opts = struct();
    if ~isempty(occ),   opts.occ = occ;     end
    lastwarn('');
    [prog,Pc] = lpivar_cdopvar(prog,dm,sp,[0 1],dg,opts);
    [~,wid] = lastwarn;
    assert(~strcmp(wid,'sdopvar:noncanonicalMultiplier'),'%s: constructor folded a multiplier',lbl);
    assert(all(ismember(Pc.Zd,prog.decvartable)),'%s: a variable is not in the program',lbl);
    [ok,msg,nsupp] = check_spec(Pc,dg,occ);
    assert(ok,'%s: %s',lbl,msg);
    assert(nsupp==numel(Pc.Zd),'%s: support %d ~= variables %d (not injective)',lbl,nsupp,numel(Pc.Zd));
    [ok,msg] = unit_responses(Pc,40);
    assert(ok,'%s: %s',lbl,msg);
    nchk = nchk+5;
    fprintf('  passed: %-42s %5d variables\n',lbl,numel(Pc.Zd));
end

% ------------------------------------------------------------ arithmetic
% The result is an ordinary decision container: it composes with a fixed
% one, takes an adjoint, and its constraints reach the program.
prog = lpiprogram(polynomial({'s'}),[],[0,1]);
[prog,Q] = lpivar_cdopvar(prog,[2;3],{{},{'s'}},[0,1],[1 1 1]);
T = rand_copvar(struct('out',{{ {}, {'s'} }},'in',{{ {}, {'s'} }}), ...
                struct('out',[2;3],'in',[2;3]),[0,1],1,0.7);
R = T'*Q - Q'*T;                                            %#ok<NASGU>
prog = lpi_eq_cdopvar(prog,T'*Q);
assert(isa(Q','cdopvar') && isfield(prog,'expr'),'arithmetic / constraint on the result');
nchk = nchk+1;

% ------------------------------------------------------------ errors
prog = lpiprogram(polynomial({'s'}),[],[0,1]);
expect_error(@() lpivar_cdopvar(prog,[1;1],{{},{'s'}},[0,1],1,struct('occ',true(3))),'');
expect_error(@() lpivar_cdopvar(prog,[1;1],{{},{'s'}},[0,1],struct('bogus',1)),'');
expect_error(@() lpivar_cdopvar(prog,[1;1],{{},{'s_dum'}},[0,1],1),'');
expect_error(@() lpivar_cdopvar(prog,[1;1;1],{{},{'s'}},[0,1],1),'');
nchk = nchk+4;

fprintf('test_lpivar_cdopvar passed (%d checks).\n',nchk);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sp,dm,rmap,cmap] = spaces_1d(n)
% Container spaces for lpivar's n = [out_R in_R; out_L2 in_L2], dropping a
% zero-dimensional part; rmap(k)/cmap(k) is the container row/column of the
% R (k=1) or L2 (k=2) part, 0 if absent.
spo = {};   dmo = [];   rmap = [0 0];
spi = {};   dmi = [];   cmap = [0 0];
if n(1,1)>0,  spo{end+1} = {};      dmo(end+1,1) = n(1,1);  rmap(1) = numel(spo);  end
if n(2,1)>0,  spo{end+1} = {'s1'};  dmo(end+1,1) = n(2,1);  rmap(2) = numel(spo);  end
if n(1,2)>0,  spi{end+1} = {};      dmi(end+1,1) = n(1,2);  cmap(1) = numel(spi);  end
if n(2,2)>0,  spi{end+1} = {'s1'};  dmi(end+1,1) = n(2,2);  cmap(2) = numel(spi);  end
sp = struct('out',{spo},'in',{spi});
dm = struct('out',dmo,'in',dmi);
end


function S = legacy_supports(Zop)
% Kernel supports of a stock dopvar, by 'dpvar2poly' and 'subs' at random
% decision values; each slot as a sorted cellstr of 'i,j,deg_s,deg_th'.
flds = {'P','P'; 'Q1','Q1'; 'Q2','Q2'; 'R0','R0'; 'R1','R1'; 'R2','R2'};
S = struct();
for k = 1:size(flds,1)
    f = flds{k,1};
    if any(strcmp(f,{'R0','R1','R2'})),     X = Zop.R.(f);
    else,                                   X = Zop.(f);
    end
    S.(f) = cell(0,1);
    if isempty(X) || any(size(X)==0),  continue,   end
    Xp = dpvar2poly(X);
    dn = reshape(X.dvarname,[],1);
    if ~isempty(dn),    Xp = subs(Xp,dn,randn(numel(dn),1));  end
    % Q1 is integrated in s itself (lpivar's var1); name that the theta slot
    % so it lines up with the container, where it is the input's dummy.
    if strcmp(f,'Q1'),  S.(f) = supp(Xp,{'','s1'});
    else,               S.(f) = supp(Xp,{'s1','s1_dum'});
    end
end
end


function S = container_supports(Pc,rmap,cmap)
% The same slots read off the container through 'pi_blk_kernels'.
dval = randn(numel(Pc.Zd),1);
S = struct('P',{cell(0,1)},'Q1',{cell(0,1)},'Q2',{cell(0,1)}, ...
           'R0',{cell(0,1)},'R1',{cell(0,1)},'R2',{cell(0,1)});
slot = @(r,c) Pc.C{rmap(r),cmap(c)};
if rmap(1) && cmap(1),  K = pi_blk_kernels(slot(1,1),dval);  S.P  = supp(K{1},{'s1','s1_dum'});  end
if rmap(1) && cmap(2),  K = pi_blk_kernels(slot(1,2),dval);  S.Q1 = supp(K{1},{'','s1_dum'});    end
if rmap(2) && cmap(1),  K = pi_blk_kernels(slot(2,1),dval);  S.Q2 = supp(K{1},{'s1','s1_dum'});  end
if rmap(2) && cmap(2)
    K = pi_blk_kernels(slot(2,2),dval);
    S.R0 = supp(K{1},{'s1','s1_dum'});
    S.R1 = supp(K{2},{'s1','s1_dum'});
    S.R2 = supp(K{3},{'s1','s1_dum'});
end
end


function s = supp(X,vnames)
% Nonzero (entry, monomial) pairs of a matrix polynomial, as sorted
% 'i,j,d1,d2' strings with degrees in vnames (an empty name counts 0).
X = polynomial(X);
[m,n] = size(X);
cf = X.coefficient;     dg = X.degmat;      vn = X.varname;
s = cell(0,1);
for e = 1:m*n
    [i,j] = ind2sub([m,n],e);
    t = find(abs(cf(:,e))>1e-12);
    for r = t(:)'
        d = zeros(1,numel(vnames));
        for v = 1:numel(vnames)
            if isempty(vnames{v}),  continue,   end
            p = find(strcmp(vn,vnames{v}),1);
            if ~isempty(p),     d(v) = full(dg(r,p));   end
        end
        % any other variable with nonzero degree is a failure to report
        other = setdiff(find(full(dg(r,:))),find(ismember(vn,vnames)));
        if ~isempty(other)
            d(end+1) = -1;                                          %#ok<AGROW>
        end
        s{end+1,1} = sprintf('%d,%d,%s',i,j,mat2str(d));             %#ok<AGROW>
    end
end
s = sort(s);
end


function prog = make_prog(vars,dom)
% An LPI program for any number of spatial variables, as in
% 'test_poscopvar': 'lpiprogram' refuses more than two, and what it returns
% is 'sosprogram' plus the primary and dummy variables and a 'dom' field.
% lpivar_cdopvar only needs prog.decvartable, through 'lpidecvar'.
vars = reshape(vars,1,[]);
if numel(vars)<=2
    prog = lpiprogram(polynomial(vars(:)),[],dom);
    return
end
prog = sosprogram(polynomial([]));
dums = strcat(vars,'_dum');
prog.vartable = [prog.vartable; polynomial(vars(:)); polynomial(dums(:))];
prog.dom = dom;
end


function v = cat_spaces(sp,side)
if isstruct(sp),    L = sp.(side);  else,   L = sp;     end
v = cell(1,0);
for k = 1:numel(L)
    if ~isempty(L{k}),  v = [v, reshape(L{k},1,[])];    end        %#ok<AGROW>
end
end


function [ok,msg,nsupp] = check_spec(Pc,dg,occ)
% Every block and cell against the degree caps, the canonical form, and a
% full support; nsupp counts the admissible monomial entries found.
ok = true;  msg = '';   nsupp = 0;
[M,N] = size(Pc.C);
if isempty(occ),    occ = true(M,N);    end
dval = randn(numel(Pc.Zd),1);
for i = 1:M
    for j = 1:N
        B = Pc.C{i,j};
        if ~occ(i,j)
            if ~isempty(B),  ok = false;  msg = sprintf('block (%d,%d) should be []',i,j);  return,  end
            continue
        end
        d = block_deg(dg,i,j);
        S3 = intersect(B.vars.out,B.vars.in);
        S2 = setdiff(B.vars.out,B.vars.in);     S1 = setdiff(B.vars.in,B.vars.out);
        K = pi_blk_kernels(B,dval);
        for g = 1:numel(K)
            gam = gamma_of(g,numel(S3));
            X = polynomial(K{g});
            cf = X.coefficient;     D = full(X.degmat);     vn = X.varname;
            [m,n] = size(X);
            % expected count of monomial entries in this cell
            nl = 1;     nr = 1;
            for t = 1:numel(S3)
                if gam(t)==1,   nl = nl*(d.mult+1);
                else,           nl = nl*(d.int(1)+1);   nr = nr*(d.int(2)+1);
                end
            end
            nl = nl*(d.out+1)^numel(S2);    nr = nr*(d.in+1)^numel(S1);
            nfound = 0;
            for e = 1:m*n
                t = find(abs(cf(:,e))>1e-12);
                nfound = nfound + numel(t);
                for r = t(:)'
                    for p = 1:numel(vn)
                        v = vn{p};  deg = D(r,p);
                        if deg==0,  continue,   end
                        isdum = endsWith(v,'_dum');     base = erase(v,'_dum');
                        k3 = find(strcmp(S3,base),1);
                        if ~isempty(k3)
                            if isdum && gam(k3)==1
                                ok = false; msg = sprintf('(%d,%d) cell %d: dummy %s in a multiplier direction',i,j,g,v); return
                            end
                            if gam(k3)==1,  cap = d.mult;
                            elseif isdum,   cap = d.int(2);
                            else,           cap = d.int(1);
                            end
                        elseif any(strcmp(S2,base)) && ~isdum,  cap = d.out;
                        elseif any(strcmp(S1,base)) && isdum,   cap = d.in;
                        else
                            ok = false; msg = sprintf('(%d,%d) cell %d: unexpected variable %s',i,j,g,v); return
                        end
                        if deg>cap
                            ok = false; msg = sprintf('(%d,%d) cell %d: degree %d in %s exceeds %d',i,j,g,deg,v,cap); return
                        end
                    end
                end
            end
            if nfound ~= m*n*nl*nr
                ok = false; msg = sprintf('(%d,%d) cell %d: %d monomial entries, expected %d',i,j,g,nfound,m*n*nl*nr); return
            end
            nsupp = nsupp + nfound;
        end
    end
end
end


function [ok,msg] = unit_responses(Pc,nsample)
% Injectivity, directly: with decision variable k at 1 and the rest at 0,
% exactly ONE monomial entry of one kernel, in one block and one cell, is
% nonzero, and distinct variables land in distinct places. Count and full
% support cannot see a variable reused in two places (the count is what was
% declared, and random values still fill every position); this can.
ok = true;  msg = '';
q = numel(Pc.Zd);
ks = unique([1, q, randperm(q,min(nsample,q))]);
seen = containers.Map('KeyType','char','ValueType','double');
for k = ks
    dval = zeros(q,1);  dval(k) = 1;
    hits = {};
    for ii = 1:numel(Pc.C)
        if isempty(Pc.C{ii}),   continue,   end
        K = pi_blk_kernels(Pc.C{ii},dval);
        for g = 1:numel(K)
            X = polynomial(K{g});
            [r,e] = find(abs(X.coefficient)>1e-12);
            for t = 1:numel(r)
                hits{end+1} = sprintf('%d|%d|%d|%s',ii,g,e(t), ...
                    mono_key(X.varname,X.degmat(r(t),:)));                   %#ok<AGROW>
            end
        end
    end
    if numel(hits)~=1
        ok = false;  msg = sprintf('variable %d moves %d kernel coefficients, not 1',k,numel(hits));
        return
    end
    % The location names its variables: a bare degree vector is in the
    % polynomial's OWN variable order, which depends on its content, so s^1
    % and s_dum^1 would both print as [1].
    if isKey(seen,hits{1})
        ok = false;  msg = sprintf('variables %d and %d move the same coefficient',seen(hits{1}),k);
        return
    end
    seen(hits{1}) = k;
end
end


function k = mono_key(vn,deg)
% A monomial as 'name^deg' factors sorted by name, independent of the
% polynomial's variable order; '1' for the constant.
deg = full(deg);
f = {};
for p = find(deg)
    f{end+1} = sprintf('%s^%d',vn{p},deg(p));                              %#ok<AGROW>
end
if isempty(f),  k = '1';    else,   k = strjoin(sort(f),'*');   end
end


function d = block_deg(dg,i,j)
% The normalized degree specification of block (i,j), as lpivar_cdopvar
% reads it.
if iscell(dg),  x = dg{i,j};    else,   x = dg;     end
if isnumeric(x) && isscalar(x)
    d = struct('mult',x,'int',[x x],'out',x,'in',x);
elseif isnumeric(x)
    d = struct('mult',x(1),'int',[x(3) x(2)],'out',x(1),'in',x(1));
else
    d = struct('mult',1,'int',[1 1],'out',1,'in',1);
    f = fieldnames(x);
    for k = 1:numel(f),     d.(f{k}) = x.(f{k});    end
    if isscalar(d.int),     d.int = [d.int d.int];  end
end
end


function gam = gamma_of(g,n3)
if n3==0,   gam = zeros(1,0);   return,     end
c = cell(1,n3);     [c{:}] = ind2sub([3*ones(1,n3),1],g);    gam = cell2mat(c(1:n3));
end


function expect_error(f,id)
try
    f();
catch ME
    if ~isempty(id) && ~strcmp(ME.identifier,id)
        error('test_lpivar_cdopvar: expected ''%s'', got ''%s'': %s',id,ME.identifier,ME.message);
    end
    return
end
error('test_lpivar_cdopvar: expected an error, none was raised.');
end
