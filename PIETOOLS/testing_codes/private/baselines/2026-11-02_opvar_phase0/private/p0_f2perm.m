function [e, px, det] = p0_f2perm(A, B)
% [E,PX,DET] = P0_F2PERM(A,B) compares the operators of one F2 PIE in two
% passes up to a reordering of the scalar state components: the largest
% relative error over all operators present in both, minimized over one
% permutation PX of the state components, applied to the state side of T,
% A, Tw, Tu, Bw, Bu, Cz, Cy; the other side of each operator (w, u, z, y)
% is permuted freely (exhaustive for inputs, greedy for outputs). Only
% permutations that keep the variables of each component are tried.
%
% INPUT
% - A, B:   'struct', the contents of F2/pass<k>/<label>.mat of two passes.
%
% OUTPUT
% - e:      double, the largest relative error max|Ya-Yb| / max(1,max|Ya|)
%           over the operators, for the best PX (Inf: no admissible PX);
% - px:     1 x n double, the best permutation: state component c of
%           pass A is component px(c) of pass B;
% - det:    char, the error per operator for PX.
%
% NOTES
% - 2 test functions per input component, 4 points per variable set,
%   total degree <= 2; responses of B are computed once per (test function,
%   position) by linearity, so the search costs no further applications.
q = 2;  P = 4;
xin_ops = {'T','A','Cz','Cy'};      xout_ops = {'T','A','Tw','Tu','Bw','Bu'};
ops = intersect({A.fix.name}, {B.fix.name}, 'stable');
e = Inf;    px = [];    det = '';
iT = strcmp({A.fix.name},'T');  jT = strcmp({B.fix.name},'T');
if ~any(iT) || ~any(jT),    det = 'no T';   return,     end
sa = p0_comps(A.fix(iT).X);     sb = p0_comps(B.fix(jT).X);
if numel(sa)~=numel(sb),    det = 'state components differ in number';   return,     end
PX = p0_admissible(sa, sb);
S = cell(1,numel(ops));
for k = 1:numel(ops)
    S{k} = p0_prep(A.fix(strcmp({A.fix.name},ops{k})).X, B.fix(strcmp({B.fix.name},ops{k})).X, ...
        q, P, sum(double([A.label,ops{k}])));
end
for r = 1:size(PX,1)
    p = PX(r,:);    emax = 0;   dd = '';
    for k = 1:numel(ops)
        pin = [];   pout = [];
        if ismember(ops{k}, xin_ops),   pin = p;    end
        if ismember(ops{k}, xout_ops),  pout = p;   end
        ek = p0_operr(S{k}, pin, pout);
        dd = sprintf('%s %s:%.2g', dd, ops{k}, ek);
        emax = max(emax, ek);
        if emax >= e,   break,  end
    end
    if emax < e,    e = emax;   px = p;     det = strtrim(dd);  end
end
end

function S = p0_prep(Xa, Xb, q, P, seed)
% Responses of Xa to all test functions at their own positions, and of Xb
% to each test function at each input position with the same variables,
% on points shared by both.
[ia, oa] = p0_comps(Xa);    [ib, ob] = p0_comps(Xb);
S = struct('ia',{ia},'ib',{ib},'oa',{oa},'ob',{ob},'Ya',{{}},'R',{{}});
if numel(ia)~=numel(ib) || numel(oa)~=numel(ob),    return,     end
v = cellstr(string(Xa.vars(:)));
rng(seed,'twister');
t = cell(numel(ia),1);
for c = 1:numel(ia),    t{c} = p0_randpoly(1, q, ia{c});    end
keys = unique(cellfun(@(z) strjoin(z,','), [oa; ob], 'UniformOutput', false));
pts = containers.Map();
for j = 1:numel(keys)
    nm = strsplit(keys{j},',');     if isempty(keys{j}),    nm = {};    end
    di = Xa.dom(ismember(v, nm),:);
    pts(keys{j}) = di(:,1)' + rand(P, numel(nm)).*(di(:,2)-di(:,1))';
end
S.Ya = p0_resp(Xa, t, oa, pts);
S.R = cell(numel(ia), numel(ib));
for c = 1:numel(ia)
    for c2 = 1:numel(ib)
        if isequal(ia{c}, ib{c2})
            z = repmat({[]}, numel(ib), 1);     z{c2} = t{c};
            S.R{c,c2} = p0_resp(Xb, z, ob, pts);
        end
    end
end
end

function e = p0_operr(S, pin, pout)
% Error of one operator with the input (PIN) and output (POUT) permutation
% fixed when given, searched when empty.
e = Inf;
if isempty(S.Ya),   return,     end
if isempty(pin),    PI = p0_admissible(S.ia, S.ib);     else,   PI = pin;   end
sc = max(1, max(cellfun(@(y) max(abs(y(:))), S.Ya)));
for r = 1:size(PI,1)
    p = PI(r,:);
    if any(arrayfun(@(c) isempty(S.R{c,p(c)}), 1:numel(p))),   continue,   end
    Yb = S.R{1,p(1)};
    for c = 2:numel(p)
        Rc = S.R{c,p(c)};
        for o = 1:numel(Yb),    Yb{o} = Yb{o} + Rc{o};  end
    end
    ee = 0;
    if isempty(pout)                                % greedy output match
        used = false(numel(Yb),1);
        for o = 1:numel(S.Ya)
            bd = Inf;   bj = 0;
            for o2 = 1:numel(Yb)
                if used(o2) || ~isequal(S.oa{o}, S.ob{o2}),     continue,   end
                dd = max(abs(S.Ya{o}(:) - Yb{o2}(:)));
                if dd < bd,     bd = dd;    bj = o2;    end
            end
            if bj==0,   ee = Inf;   break,  end
            used(bj) = true;    ee = max(ee, bd/sc);
        end
    else
        for o = 1:numel(S.Ya)
            if ~isequal(S.oa{o}, S.ob{pout(o)}),    ee = Inf;   break,  end
            ee = max(ee, max(abs(S.Ya{o}(:) - Yb{pout(o)}(:))) / sc);
        end
    end
    e = min(e, ee);
end
end

function [ic, oc] = p0_comps(X)
% Variable names of each scalar input (IC) and output (OC) component.
v = cellstr(string(X.vars(:)));
ic = {};    oc = {};
for j = 1:numel(X.dim_in),  ic = [ic; repmat({v(logical(X.space_in(j,:)))'},X.dim_in(j),1)];      end %#ok<AGROW>
for i = 1:numel(X.dim_out), oc = [oc; repmat({v(logical(X.space_out(i,:)))'},X.dim_out(i),1)];    end %#ok<AGROW>
end

function Y = p0_resp(X, z, oc, pts)
% Output components of X at their points for the input components Z
% (empty = zero), block by block (apply_sopvar).
N = numel(X.dim_in);    M = numel(X.dim_out);
ej = [0; cumsum(X.dim_in(:))];  ei = [0; cumsum(X.dim_out(:))];
q = 0;
for c = 1:numel(z),     if ~isempty(z{c}),  q = size(z{c},2);   break,  end,    end
Y = cell(numel(oc),1);
for i = 1:M
    acc = polynomial(zeros(X.dim_out(i), q));
    for j = 1:N
        B = X.C{i,j};
        if isempty(B),  continue,   end
        idx = ej(j)+1:ej(j+1);
        if all(cellfun(@isempty, z(idx))),  continue,   end
        u = cell(numel(idx),1);
        for r = 1:numel(idx)
            if isempty(z{idx(r)}),  u{r} = polynomial(zeros(1,q));
            else,                   u{r} = polynomial(z{idx(r)});
            end
        end
        acc = acc + apply_sopvar(B, vertcat(u{:}));
    end
    for r = 1:X.dim_out(i)
        c = ei(i) + r;
        Y{c} = p0_evalpoly(acc(r,:), oc{c}, pts(strjoin(oc{c},',')));
    end
end
end

function P = p0_admissible(ia, ib)
% Rows: the permutations p with ib{p(c)} equal to ia{c} for every c, built
% group by group (components with equal variables), so their number is the
% product of the group factorials rather than numel(ia)!.
n = numel(ia);
if n == 0,  P = zeros(1,0);     return,     end
ka = cellfun(@(z) strjoin(z,','), ia, 'UniformOutput', false);
kb = cellfun(@(z) strjoin(z,','), ib, 'UniformOutput', false);
P = zeros(1,n);
for g = unique(ka)'
    pa = find(strcmp(ka, g{1}));    pb = find(strcmp(kb, g{1}));
    if numel(pa) ~= numel(pb),  P = zeros(0,n);     return,     end
    G = perms(pb(:)');
    if size(G,1)*size(P,1) > 1e5
        error('p0_f2perm:many','More than 1e5 admissible permutations.');
    end
    m0 = size(P,1);
    P = repelem(P, size(G,1), 1);
    P(:,pa) = repmat(G, m0, 1);
end
end
