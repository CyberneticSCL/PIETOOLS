% Test of clean_op for 'sdopvar': over 9 random operators (N3 = 0,1,2
% shared variables, 1-2 input-only and output-only variables), (1) entries
% of A and B set to 1e-14, ztol and 2*ztol: the first two are zeroed, the
% third kept, checked against the kernels (params2kernel), the decision
% variables and the canonical form (a rebuild through the constructor
% gives the same parameters); (2) an operator with all coefficients below
% ztol is flagged is_zero; (3) a nonzero operator is not flagged and is
% unchanged. Then (4) an 'sdopvar' from a 'dpvar' times a 'sopvar' (c*S)
% is cleaned, and 0*(c*S) flagged zero, and (5) the cost at q = 1e6
% decision variables is compared with uminus (both linear in the nonzeros
% of B). Errors on the first failure.
%
% DJ, 11/24/2026: Initial coding;
% DJ, 11/24/2026: Case 4 uses c*S, as c*x (a 'dpvar' times a state) is no
%   longer supported;

ztol = 1e-12;
nn = @(R) sum(cellfun(@nnz,R.params.A),'all') + sum(cellfun(@nnz,R.params.B),'all');

for trial = 1:9
    rng(8100+trial);
    N3 = mod(trial-1,3);    N1 = randi(2);  N2 = randi(2);
    S3vars = arrayfun(@(k)sprintf('z%d',k),1:N3,'UniformOutput',false);
    S1vars = arrayfun(@(k)sprintf('x%d',k),1:N1,'UniformOutput',false);
    S2vars = arrayfun(@(k)sprintf('y%d',k),1:N2,'UniformOutput',false);
    vars = struct('in',{[S3vars,S1vars]},'out',{[S2vars,S3vars]});
    P = rand_sdopvar([randi(2),randi(2)],vars,[0,1],1,5,0.6);

    %%%% Case 1: 1e-14 and ztol zeroed, 2*ztol kept, in A and in B
    Q = P;  ndrop = 0;
    for k=1:numel(Q.params.A)
        A = Q.params.A{k};  B = Q.params.B{k};
        [ia,~,va] = find(A);    [ib,jb,vb] = find(B);
        if numel(va)>=3, va(1:3) = [1e-14;-ztol;2*ztol];   ndrop = ndrop+2; end
        if numel(vb)>=3, vb(1:3) = [-1e-14;ztol;-2*ztol];  ndrop = ndrop+2; end
        Q.params.A{k} = sparse(ia,1,va,size(A,1),1);
        Q.params.B{k} = sparse(ib,jb,vb,size(B,1),size(B,2));
    end
    [Qc,is_zero1] = clean_op(Q,ztol);
    if nn(Q)-nn(Qc)~=ndrop || is_zero1
        error("Test_sdopvar_clean_op failed: %d coefficients removed, %d expected (trial=%d).",nn(Q)-nn(Qc),ndrop,trial)
    end
    % The kernels differ exactly by the dropped coefficients
    Kq = Q.kernels;     Kc = Qc.kernels;    ndiff = 0;
    for k=1:numel(Kq)
        D = Kc{k} - Kq{k};      D = D.C;
        if any(abs(nonzeros(D))>ztol)
            error("Test_sdopvar_clean_op failed: a kernel coefficient above ztol changed (trial=%d).",trial)
        end
        ndiff = ndiff + nnz(D);
    end
    if ndiff~=ndrop
        error("Test_sdopvar_clean_op failed: the kernels differ in %d coefficients, %d expected (trial=%d).",ndiff,ndrop,trial)
    end
    % Decision variables kept, canonical form kept
    R = sdopvar(Qc.params,Qc.vars,Qc.Zd,Qc.ZL,Qc.ZR,Qc.dom,Qc.dims);
    if ~isequal(Qc.Zd,Q.Zd) || ~isequal(R.params,Qc.params)
        error("Test_sdopvar_clean_op failed: decision variables or canonical form changed (trial=%d).",trial)
    end

    %%%% Case 2: all coefficients below ztol -- is_zero true
    Pz = P;
    for k=1:numel(Pz.params.A)
        Pz.params.A{k} = 1e-15*Pz.params.A{k};     Pz.params.B{k} = 1e-15*Pz.params.B{k};
    end
    [Pz_clean,is_zero2] = clean_op(Pz);
    if ~is_zero2 || nn(Pz_clean)~=0
        error("Test_sdopvar_clean_op failed: an operator below tolerance is not zero (trial=%d).",trial)
    end

    %%%% Case 3: nonzero operator -- not flagged, unchanged
    [P3,is_zero3] = clean_op(P);
    if is_zero3 || ~isequal(P3.params,P.params)
        error("Test_sdopvar_clean_op failed: a clean nonzero operator was changed or flagged (trial=%d).",trial)
    end
end

%%%% Case 4: an 'sdopvar' from a 'dpvar' times a 'sopvar', and 0 times it
dpvar c
S = mat2sopvar(polynomial(1),struct('in',{cell(1,0)},'out',{{'s'}}),struct('in',zeros(0,2),'out',[0,1]));
D = c*S;
[D1,is_zero4] = clean_op(D);
if is_zero4 || ~isa(D1,'sdopvar')
    error("Test_sdopvar_clean_op failed: c*S is not cleaned as a nonzero 'sdopvar'.")
end
[~,is_zero5] = clean_op(0*D);
if ~is_zero5
    error("Test_sdopvar_clean_op failed: 0*(c*S) is not flagged zero.")
end

%%%% Case 5: cost at q = 1e6, nothing dropped and 10% of B dropped, vs uminus
rng(8200);
q = 1e6;
P = rand_sdopvar([2,2],struct('in',{{'s1','s2'}},'out',{{'s1','s2'}}),[0,1],1,q,2/(64*9));
Ps = P;
for k=1:numel(Ps.params.B)
    [ib,jb,vb] = find(Ps.params.B{k});
    vb(rand(size(vb))<0.1) = 1e-15;
    Ps.params.B{k} = sparse(ib,jb,vb,size(Ps.params.B{k},1),size(Ps.params.B{k},2));
end
t0 = tic;   uminus(P);              tu = toc(t0);
t0 = tic;   clean_op(P,ztol);       tc = toc(t0);
t0 = tic;   clean_op(Ps,ztol);      ts = toc(t0);
fprintf('q = %g: uminus %.4f s, clean_op %.4f s, clean_op (10%% dropped) %.4f s\n',q,tu,tc,ts)
if tc>20*tu+0.2 || ts>20*tu+0.2
    error("Test_sdopvar_clean_op failed: clean_op is much slower than uminus at q = %g.",q)
end

disp('Test_sdopvar_clean_op passed.')
