function [prog,P,info] = pf_row_storage(prog,X,side,dom,row)
% [PROG,P,INFO] = PF_ROW_STORAGE(PROG,X,SIDE,DOM,ROW) declares the positive
% storage operator on the SIDE spaces of X from the row of the settings
% proposal (sopvar_settings_proposal.tex, Sec. "The settings"): for a
% space in e variables the entry e of ROW gives the lengths
%   len_inc  d1 (1 x e)        incomplete moments (pattern c)
%   len_com  d2 (1 x e)        complete moments (pattern b)
%   len_sos  d3 (3 x e)        multiplier lengths of the blocks a, b, c
%   patterns 'uniform' | 'separable' | 'all' | cellstr of pattern words
%            ('aab', ...) | numeric rows (1 a, 2 b, 3 c)
% and len_R the multiplier length of an R^q space (default 1, the identity).
% The reader (Sec. "The reader"): basis operators = the patterns; per
% operator and owned direction d, mult_d = length(alpha_d) - 1 (0, d2_d - 1,
% d1_d - 1), int_d = d3_{alpha_d,d} - 1; in a direction the space does not
% own, int_d = the evaluation length of the entry of the largest space
% owning d, less one (rule R-nonowned, this file); no total cap. A zero
% length drops the pattern (zero is absent). Codes of poscopvar_direct: a =
% 1; c = 2 (lower integral); b = 3 (upper) in a direction where some space
% carries c, else 4 (full, separable), by Prop. span of the proposal, which
% needs equal multiplier lengths for b and c. The storage carries no term.
% ROW: struct with fields len_R, len_inc, len_com, len_sos, patterns, each
% of the last four a 1 x 3 cell for entries 1, 2, 3 (the last entry serves
% every higher number of variables). SCRATCH, MMP, 11/16/2026.
[sp,dm] = copvar_space_list(X,side);
vars = dom.vars;    nv = numel(vars);    M = numel(sp);
nd = cellfun(@numel,sp);
own = false(M,nv);
for k = 1:M,    own(k,ismember(vars,sp{k})) = true;    end
e_of = @(k) min(max(nd(k),1),3);
if ~isfield(row,'len_R') || isempty(row.len_R),     row.len_R = 1;  end
% pattern rows per space (codes 1 a, 2 b, 3 c), lengths applied
pats = cell(1,M);
for k = 1:M
    if nd(k)==0,    continue;   end
    e = e_of(k);    n = nd(k);
    pt = row.patterns{e};
    if ischar(pt)
        switch pt
            case 'uniform',     A = [ones(1,n); 2*ones(1,n)];
            case 'separable',   A = multiindex_grid(repmat({[1 2]},1,n));
            case 'all',         A = multiindex_grid(repmat({[1 2 3]},1,n));
            otherwise,  error('pf_row_storage:patterns','unknown pattern set %s',pt);
        end
    elseif iscell(pt)
        A = zeros(numel(pt),n);
        for i = 1:numel(pt),    w = pt{i};    A(i,:) = (w=='b')*1 + (w=='c')*2 + 1;   end
    else
        A = pt;
    end
    d1 = fit(row.len_inc{e},n);   d2 = fit(row.len_com{e},n);   d3 = row.len_sos{e};
    if size(d3,2)~=n,   d3 = repmat(d3(:,1),1,n);    end
    keep = true(size(A,1),1);
    for i = 1:size(A,1)
        for d = 1:n
            a = A(i,d);
            if (a==2 && d2(d)==0) || (a==3 && d1(d)==0) || d3(a,d)==0,   keep(i) = false;  end
        end
    end
    pats{k} = A(keep,:);
end
% separable directions: no owning space carries a c pattern there
sepv = true(1,nv);
for k = 1:M
    if nd(k)==0,    continue;   end
    idx = find(own(k,:));
    for d = 1:nd(k)
        if any(pats{k}(:,d)==3),    sepv(idx(d)) = false;    end
    end
end
% the largest entry owning each direction, for rule R-nonowned
big = zeros(1,nv);
for d = 1:nv,   o = find(own(:,d));  if ~isempty(o),  big(d) = max(nd(o));  end,  end
% degrees and include rows
deg = cell(1,M);    incl = cell(1,M);
for k = 1:M
    if nd(k)==0,    deg{k} = struct('int',row.len_R-1);    continue;   end
    e = e_of(k);    n = nd(k);  idx = find(own(k,:));
    d1 = fit(row.len_inc{e},n);   d2 = fit(row.len_com{e},n);   d3 = row.len_sos{e};
    if size(d3,2)~=n,   d3 = repmat(d3(:,1),1,n);    end
    A = pats{k};    rows = zeros(size(A));   dk = cell(1,size(A,1));
    for i = 1:size(A,1)
        mu = zeros(1,nv);   in = zeros(1,nv);
        for d = 1:n
            a = A(i,d);     g = idx(d);
            switch a
                case 1,     rows(i,d) = 1;                      mu(g) = 0;
                case 2,     rows(i,d) = 4 - ~sepv(g);           mu(g) = d2(d)-1;   % 4 separable, 3 (upper) else
                case 3,     rows(i,d) = 2;                      mu(g) = d1(d)-1;
            end
            in(g) = d3(a,d)-1;
        end
        for g = find(~own(k,:))                   % rule R-nonowned (v2, 11/16 06:40): an operator carries its own largest multiplier degree into the directions its space does not own, as POSCOPVAR with DEG omitted does (int = the block degree in every theta); v1 (the evaluation length of the largest owner) lost S3
            in(g) = max(in(idx));
        end
        dk{i} = struct('int',in,'mult',mu);
        % per-direction caps (row field 'cap', one entry per direction of the
        % entry; Inf = none): owned directions from this entry, non-owned
        % from the entry of the largest owner (rule R-nonowned)
        if isfield(row,'cap') && ~isempty(row.cap)
            cv = Inf(1,nv);
            ce = fit(row.cap{e},n);
            for d = 1:n,    cv(idx(d)) = ce(d);     end
            for g = find(~own(k,:))
                eb = min(big(g),3);     cb = row.cap{eb};
                cv(g) = cb(1);
            end
            dk{i}.cap = cv;
        end
    end
    deg{k} = dk;    incl{k} = rows;
end
opts = struct('sep',sepv);  opts.include = incl;
[prog,P,~,bl] = poscopvar(prog,dm,sp,dom,deg,opts);
info = struct('basis_list',bl,'sep',sepv,'patterns',{pats},'deg',{deg});
end

function v = fit(v,n)
% V = FIT(V,N) a length row of n entries from a scalar or a row
v = reshape(v,1,[]);
if isscalar(v),     v = repmat(v,1,n);
elseif numel(v)~=n, error('pf_row_storage:len','a length row has %d entries for a space in %d variables',numel(v),n);
end
end
