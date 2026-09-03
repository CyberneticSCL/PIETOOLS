function p = leftPermuteVec(n, m)
% Return permutation vector p such that PL*x = x(p) for PL from leftPermute(n,m).
    n = n(:).'; m = m(:).';
    k = numel(n);
    Nm = prod(n.*m);

    dnm = n.*m;
    strideNM = kronStrides(dnm);

    % enumerate all alpha and beta tuples
    A = allSubs(n);
    B = allSubs(m);
    NA = size(A,1);
    NB = size(B,1);

    Aall = repelem(A, NB, 1);
    Ball = repmat(B,  NA, 1);

    Gamma = (Aall - 1).*m + Ball; % (Nm)×k

    % output index i for each input index j=1..Nm
    i = 1 + sum((Gamma - 1) .* strideNM, 2);
    j = (1:Nm).';

    % PL(i,j)=1 => PL*x has entries x(j) moved to i
    % So (PL*x)(i)=x(j). Therefore p(i)=j.
    p = zeros(Nm,1);
    p(i) = j;
end

function stride = kronStrides(d)
    k = numel(d);
    stride = ones(1,k);
    for ii = 1:k-1
        stride(ii) = prod(d(ii+1:end));
    end
    stride(k) = 1;
end

function S = allSubs(d)
    d = d(:).';
    k = numel(d);
    N = prod(d);
    S = zeros(N,k);
    idx = (0:N-1).';
    for ii = k:-1:1
        S(:,ii) = mod(idx, d(ii)) + 1;
        idx = floor(idx / d(ii));
    end
end
