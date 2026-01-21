function out = subsref(F,S)

if ~strcmp(S(1).type,'()')
    out = builtin('subsref',F,S);
    return;
end

m  = F.dim(1); n  = F.dim(2);
ds = size(F.Zs,1);
dt = size(F.Zt,1);

subs = S(1).subs;

if numel(subs) == 2
    r = expandIdx(subs{1}, m);
    c = expandIdx(subs{2}, n);

    rIdx = reshape(((r(:)-1)*ds + (1:ds)).', [], 1);
    cIdx = reshape(((c(:)-1)*dt + (1:dt)).', [], 1);

    Csub = F.C(rIdx, cIdx);
    out  = quadPoly(Csub, F.Zs, F.Zt, [numel(r) numel(c)], F.ns, F.nt);
elseif numel(subs) == 1
    idxRaw = subs{1};

    % F(:) should behave like MATLAB: column vector output
    if ischar(idxRaw) && strcmp(idxRaw,':')
        idx = (1:m*n).';
        outSz = [m*n 1];
    else
        idx = expandIdx(idxRaw, m*n);     % includes end already resolved
        outSz = size(idx);
        idx = idx(:);
    end

    [rb, cb] = ind2sub([m n], idx);
    K = numel(idx);
    [ro, co] = ind2sub(outSz, (1:K).');

    Icell = cell(K,1); Jcell = cell(K,1); Vcell = cell(K,1);

    for t = 1:K
        rSrc = (rb(t)-1)*ds + (1:ds);
        cSrc = (cb(t)-1)*dt + (1:dt);

        [ib,jb,vb] = find(F.C(rSrc, cSrc));
        if isempty(vb), continue; end

        Icell{t} = (ro(t)-1)*ds + ib;
        Jcell{t} = (co(t)-1)*dt + jb;
        Vcell{t} = vb;
    end

    I = vertcat(Icell{:});
    J = vertcat(Jcell{:});
    V = vertcat(Vcell{:});

    Cnew = sparse(I, J, V, outSz(1)*ds, outSz(2)*dt);
    out  = quadPoly(Cnew, F.Zs, F.Zt, outSz, F.ns, F.nt);

else
    error('quadPoly:subsref','Use F(r,c) or F(k).');
end

if numel(S) > 1
    out = subsref(out, S(2:end));
end

end

function idx = expandIdx(idx, N)
if ischar(idx) && strcmp(idx,':')
    idx = 1:N;
elseif isstring(idx) && isscalar(idx) && idx==":"
    idx = 1:N;
elseif islogical(idx)
    idx = find(idx);
end
idx = double(idx);
end
