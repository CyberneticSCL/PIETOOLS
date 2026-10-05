% Test_sopvar_clean_op.m
%

tol = 1e-8;

for trial = 1:12
    rng(7100+trial);

    N3 = mod(trial-1,3);           % number of shared (S3) variables
    N1 = randi(2);                 % number of input-only (S1) variables
    N2 = randi(2);                 % number of output-only (S2) variables

    S3vars = arrayfun(@(k)sprintf('z%d',k),1:N3,'UniformOutput',false);
    S1vars = arrayfun(@(k)sprintf('x%d',k),1:N1,'UniformOutput',false);
    S2vars = arrayfun(@(k)sprintf('y%d',k),1:N2,'UniformOutput',false);

    vars = struct('in',{[S3vars,S1vars]},'out',{[S2vars,S3vars]});
    dom = struct('in',repmat([0,1],N3+N1,1),'out',repmat([0,1],N2+N3,1));
    degs = struct('in',randi(2,[1,N3+N1]),'out',randi(2,[1,N2+N3]));
    matdim = [randi(2),randi(2)];
    dnsty = 0.8;

    Pop = rand_sopvar(matdim,vars,dom,degs,dnsty);
    x = rand_poly([matdim(2),1],polynomial(Pop.vars.in(:)),2,20);

    %%%% Case 1: a coefficient overwritten with sub-tolerance noise is
    %%%% zeroed; every other coefficient is left bit-identical.
    ztol = 1e-6;
    noise = ztol/1000;
    k = randi(numel(Pop.params));
    Pop_noisy = Pop;
    Pop_noisy.params{k}(1,1) = noise;
    [Pop_clean,is_zero1] = clean_op(Pop_noisy,ztol);
    if Pop_clean.params{k}(1,1)~=0
        error("Test_sopvar_clean_op failed: sub-tolerance noise not zeroed (trial=%d).",trial)
    end
    for kk=1:numel(Pop.params)
        diff_kk = Pop_clean.params{kk} - Pop_noisy.params{kk};
        if kk==k
            diff_kk(1,1) = 0;   % the one entry clean_op is expected to change
        end
        if nnz(diff_kk)~=0
            error("Test_sopvar_clean_op failed: cleaning changed an entry other than the noised one (trial=%d, cell=%d).",trial,kk)
        end
    end

    %%%% Case 2: entirely-zero operator -- is_zero true, acts as zero map.
    Pop_zero = Pop;
    for kk=1:numel(Pop_zero.params)
        Pop_zero.params{kk} = 0*Pop_zero.params{kk};
    end
    [Pop_zero_clean,is_zero2] = clean_op(Pop_zero);
    if ~is_zero2
        error("Test_sopvar_clean_op failed: entirely-zero operator not flagged is_zero (trial=%d).",trial)
    end
    Px_zero = apply_sopvar(Pop_zero_clean,x);
    if max(max(abs(full(double(cleanpoly(Px_zero,1e-10).C)))))>tol
        error("Test_sopvar_clean_op failed: entirely-zero operator did not act as the zero map (trial=%d).",trial)
    end

    %%%% Case 3: genuinely nonzero operator -- is_zero false, action
    %%%% exactly preserved (sprand coefficients are never below the
    %%%% default 1e-12 tolerance, so nothing should actually be zeroed).
    [Pop_clean3,is_zero3] = clean_op(Pop);
    if is_zero3
        error("Test_sopvar_clean_op failed: nonzero operator incorrectly flagged is_zero (trial=%d).",trial)
    end
    Px0 = apply_sopvar(Pop,x);
    Px3 = apply_sopvar(Pop_clean3,x);
    err3 = cleanpoly(Px0-Px3,1e-10);
    if max(max(abs(full(err3.C))))>tol
        error("Test_sopvar_clean_op failed: cleaning an already-clean operator changed its action (trial=%d).",trial)
    end
end

disp('Test_sopvar_clean_op passed.')
