% Test_intop_combine_terms_idempotent.m -- Tests that
% @intop/combine_terms.m (the form with arguments (Kop,degs,dims,depmat))
% puts the repeated factors of a state in a canonical order, so that
%   (1) combining the terms of a functional that was already combined
%       leaves it unchanged, exactly: the same intmat, monomials, and
%       parameters;
%   (2) the combined functional represents the same functional: it gives
%       the same value on x o ... o x as the original;
%   (3) in every term with nonzero parameters, the factors are ordered by
%       decreasing order of integration (the last variable integrated
%       first), ties broken by the order of integration along the next
%       spatial variable, and remaining ties keeping their order. This is
%       checked independently of the code, from the integration orders
%       (int_order) of the combined functional. Terms without parameters
%       are skipped: the terms are a product over the rows along each
%       variable, so they cannot all be in order, and combine_terms leaves
%       them as they are.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES
% - Part 1 expects the tied square [0,3,0,3] as the chain [0,3,0,1] with
%   kernel k(s1,s2)+k(s2,s1): combine_terms splits every integral into
%   cells (see Test_intop_combine_terms_canonical).
% - Part 3 draws tied orders of integration with rand_ties_intmat.m (seeds
%   50000+trial).

tol = 1e-8;

%% Part 1: N=1, degree 2: tied full square, and chains (ordered)
pvar s_dum1 s_dum2
intmat_arr = {[0,3,0,3],[0,3,0,1],[0,2,0,3]};
for j=1:numel(intmat_arr)
    Kop = intop({s_dum1*(1+s_dum2^2)},intmat_arr(j),{'s'},[0,1],false);
    K1 = combine_terms(Kop,2,1,true);
    K2 = combine_terms(K1,2,1,true);
    K3 = combine_terms(K2,2,1,true);
    if ~same_op(K1,K2) || ~same_op(K2,K3)
        error("Test_intop_combine_terms_idempotent failed: Part 1, intmat %s is not left unchanged by a second call.",mat2str(intmat_arr{j}))
    end
    % Expected kernel: k for [0,3,0,1], k(s2,s1) for [0,2,0,3], k(s1,s2)+k(s2,s1) for the square
    if ~isequal(K1.intmat{1},[0,3,0,1])
        error("Test_intop_combine_terms_idempotent failed: Part 1, intmat %s not given as [0,3,0,1].",mat2str(intmat_arr{j}))
    end
    ksw = s_dum2*(1+s_dum1^2);
    kexp = {s_dum1*(1+s_dum2^2)+ksw, s_dum1*(1+s_dum2^2), ksw};
    kk = kernels(K1);
    err = cleanpoly(kk{1}-kexp{j},1e-12);
    if ~isdouble(err) || abs(double(err))>tol
        error("Test_intop_combine_terms_idempotent failed: Part 1, kernel of intmat %s.",mat2str(intmat_arr{j}))
    end
end

%% Part 2: N=2, degree 2, hand-picked cases
vars = struct('in',{{'s1','s2'}},'out',{cell(1,0)});
dom = struct('in',[0,1;0,1],'out',zeros(0,2));
depmat = true(1,2);
pvar s1_dum1 s1_dum2 s2_dum1 s2_dum2
K1m = s1_dum1*s2_dum1^2 + 5*s1_dum2*s2_dum2^3 - 2*s1_dum2^2*s2_dum1;
K2m = 3*s1_dum1^2*s2_dum2 - s1_dum2*s2_dum1^2 + 7;
K3m = s1_dum1 + 2*s2_dum2 + s1_dum2*s2_dum1;
% intmat along s1, along s2, kernels (one per combination of terms)
cases = {{[0,2,0,3;0,3,0,1],[0,3,0,1],{K1m;K2m}}, ...              % 2 orders along s1, twisted along s2
         {[0,2,0,3;0,3,0,1;0,3,0,3],[0,3,0,1],{K1m;K2m;K3m}}, ...  % plus a tied term along s1
         {[0,3,0,3],[0,3,0,1;0,2,0,3],{K1m,K2m}}, ...              % tied along s1, orders along s2
         {[0,3,0,3],[0,3,0,3],{K3m}}, ...                          % tied along both
         {[0,2,0,3],[0,2,0,3;0,3,0,1],{K1m,K2m}}};                 % ordered along both, opposite along s2

pvar s t s1_1 s2_1 s1_2 s2_2
vx = s^2*t + 3*t^3 - s*t^2 + 1;
xval = subs(vx,[s;t],[s1_1;s2_1]) * subs(vx,[s;t],[s1_2;s2_2]);
for j=1:numel(cases)
    Kop = intop(cases{j}{3},cases{j}(1:2),vars,dom,false);
    K1 = combine_terms(Kop,2,1,depmat);
    K2 = combine_terms(K1,2,1,depmat);
    K3 = combine_terms(K2,2,1,depmat);
    if ~same_op(K1,K2) || ~same_op(K2,K3)
        error("Test_intop_combine_terms_idempotent failed: Part 2, case %d is not left unchanged by a second call.",j)
    end
    if ~is_canon(K1,2)
        error("Test_intop_combine_terms_idempotent failed: Part 2, case %d: factors not in canonical order.",j)
    end
    err = cleanpoly(apply_op(Kop,xval)-apply_op(K1,xval),1e-10);
    if ~isdouble(err) || abs(double(err))>tol
        error("Test_intop_combine_terms_idempotent failed: Part 2, case %d: the functional changed.",j)
    end
end

%% Part 3: random N in {1,2,3}, degree in {2,3,4}, with ties, scalar or vector state
for trial = 1:40
    rng(50000+trial);
    N = randi(3);
    d = randi([2,4]);
    if N==3,    d = min(d,3);   end
    if N==2,    d = min(d,3);   end
    dimval = randi(2);
    if d>=3,    dimval = 1;     end
    nterms_i = randi(3,1,N);

    vnames = cell(1,N);
    for i=1:N,  vnames{i} = ['s',num2str(i)];   end
    vars = struct('in',{vnames},'out',{cell(1,0)});
    dom = struct('in',repmat([0,1],N,1),'out',zeros(0,2));
    intmat = cell(1,N);
    dumvars = cell(1,N);
    for i=1:N
        intmat{i} = rand_ties_intmat(d,nterms_i(i),0.15);
        dv = cell(d,1);
        for j=1:d,  dv{j} = [vnames{i},'_dum',num2str(j)];  end
        dumvars{i} = polynomial(dv);
    end
    allpvars = vertcat(dumvars{:});
    params = cell([nterms_i,1]);
    for lin=1:numel(params)
        params{lin} = rand_poly([1,dimval^d],allpvars,1,6);
    end
    Kop = intop(params,intmat,vars,dom,false);
    K1 = combine_terms(Kop,d,dimval,true(1,N));
    K2 = combine_terms(K1,d,dimval,true(1,N));
    if ~same_op(K1,K2)
        error("Test_intop_combine_terms_idempotent failed: Part 3, trial=%d, N=%d, d=%d: not left unchanged by a second call.",trial,N,d)
    end
    if ~is_canon(K1,d)
        error("Test_intop_combine_terms_idempotent failed: Part 3, trial=%d, N=%d, d=%d: factors not in canonical order.",trial,N,d)
    end

    % Same-state input x^{od}: one base formula, renamed per copy
    tvars_p = polynomial(cellfun(@(a) ['t',num2str(a)],num2cell(1:N)','UniformOutput',false));
    vx = rand_poly([dimval,1],tvars_p,2,6);
    if all(isequal(vx,0)),  vx = ones(dimval,1);    end
    xval = 1;
    for k=1:d
        copyvars = cellfun(@(a) [a,'_',num2str(k)],vnames(:),'UniformOutput',false);
        xval = kron(xval,subs(vx,tvars_p,polynomial(copyvars)));
    end
    Y0 = apply_op(Kop,xval);
    Y1 = apply_op(K1,xval);
    err = cleanpoly(Y0-Y1,1e-10);
    nrmref = norm(double(cleanpoly(Y0,1e-10)),'fro');
    if ~isdouble(err) || norm(double(err),'fro')>tol*max(1,nrmref)
        error("Test_intop_combine_terms_idempotent failed: Part 3, trial=%d, N=%d, d=%d: the functional changed.",trial,N,d)
    end
end

%% Part 4: decision variables in the kernels (dpvar), degree 3, one variable, tied and ordered terms
pvar s_dum1 s_dum2 s_dum3
dpvar c1 c2
Kd = {c1*s_dum1+2*s_dum3; c2*s_dum2^2+s_dum1*s_dum3; c1*s_dum3+c2};
Kop = intop(Kd,{[0,4,0,4,0,4;0,2,0,3,0,4;0,4,0,1,0,2]},{'s'},[0,1],false);
K1 = combine_terms(Kop,3,1,true);
K2 = combine_terms(K1,3,1,true);
if ~same_op(K1,K2) || ~is_canon(K1,3)
    error("Test_intop_combine_terms_idempotent failed: Part 4, kernels with decision variables.")
end
pvar s s_1 s_2 s_3
vb = 1+s+2*s^2;
xval = subs(vb,s,s_1)*subs(vb,s,s_2)*subs(vb,s,s_3);
errd = apply_op(Kop,xval)-apply_op(K1,xval);
if ~isempty(errd.C) && max(abs(errd.C(:)))>tol
    error("Test_intop_combine_terms_idempotent failed: Part 4, the functional with decision variables changed.")
end

disp('Test_intop_combine_terms_idempotent passed.')


%% ---- local functions ----------------------------------------------------
function tf = same_op(A,B)
% A and B have the same intmat, monomials, decision variables, and parameters
tf = isequal(A.intmat,B.intmat) && isequal(A.Z.in,B.Z.in) && isequal(A.Zd,B.Zd) ...
    && isequal(size(A.params),size(B.params)) && all(cellfun(@isequal,A.params(:),B.params(:)));
end

function tf = is_canon(Kop,m)
% In each term with nonzero parameters (a combination of one intmat row per
% spatial variable), the m factors are in decreasing integration order, ties
% broken by the order along the next spatial variable
io = Kop.int_order;
N = numel(io);
nr = cellfun(@(a) size(a,1),io);
is_nz = cellfun(@(a) nnz(a)>0,Kop.params);
tf = true;
for c=1:prod(nr)
    if ~is_nz(c)
        continue
    end
    sub = cell(1,N);
    [sub{:}] = ind2sub([nr,1],c);
    key = zeros(m,N);
    for t=1:N
        % integration order of the variable each copy equals (itself, or the one a multiplier is evaluated at)
        row = Kop.intmat{t}(sub{t},:);
        is_dlt = row(1:2:end)==row(2:2:end);
        tg = 1:m;
        tg(is_dlt) = row(2*find(is_dlt)-1);
        key(:,t) = io{t}(sub{t},tg)';
    end
    for p=1:m-1
        df = key(p,:)-key(p+1,:);
        nz = find(df~=0,1);
        if ~isempty(nz) && df(nz)<0
            tf = false;
            return
        end
    end
end
end
