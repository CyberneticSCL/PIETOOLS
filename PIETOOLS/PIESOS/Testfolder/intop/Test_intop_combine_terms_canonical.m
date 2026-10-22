% Test_intop_combine_terms_canonical.m -- Tests that @intop/combine_terms.m
% (the form with arguments (Kop,degs,dims,depmat)) gives the same intmat
% and the same kernels for equal functionals, so that the coefficients of
% equal functionals can be compared cell by cell (as piesos_eq does).
%
% Every integral is written as a sum of cells, one for each ordering of the
% dummy variables consistent with its limits, in canonical form (each
% variable integrated from a up to the next larger one, the largest over
% the full domain, multipliers [t,t] at a variable that is not a
% multiplier). For copies of the same state, the cells are then relabeled
% to the same canonical order, adding the kernels of cells that coincide.
% This is checked
%   Part 1: for the square [0,3,0,3] with kernel k, which is the chain
%           [0,3,0,1] with k(s1,s2)+k(s2,s1), and for other forms of the
%           same region;
%   Part 2: for two different states, for which the cells cannot be
%           relabeled;
%   Part 3: for the cube [0,4,0,4,0,4], which becomes one chain, with the
%           kernel the sum over the 6 orderings;
%   Part 4: for multipliers, whose kernels are written in terms of the
%           variable at which they are evaluated: the kernel of [0,3,1,1]
%           is k(s1,s1);
%   Part 5: for random functionals with ties in N=1,2,3 spatial variables:
%           relabeling the copies of the functional (reorder_states) does
%           not change the result of combine_terms, the result is made of
%           cells, the value is unchanged, and calling combine_terms again
%           does nothing;
%   Part 6: the same, with multipliers, for kernels that depend on all
%           variables, also those of the multipliers, which are replaced
%           by the variable at which the multiplier is evaluated (the
%           lowest of the equal variables);
%   Part 7: for rows that are not ordered, e.g. s1 integrated from s3 to s2,
%           which the 'intop' constructor rejects: after splitting them with
%           canonicalize_intmat (the indicator of s3<=s1<=s2; an empty region
%           gives an error), the value is unchanged and the result is made of
%           cells;
%   Part 8: for groups of equal variables (multipliers), also for
%           different states;
%   Part 9: for vector-valued states: copies that are equal in all spatial
%           variables are interchangeable, x_i(s)*x_j(s) = x_j(s)*x_i(s),
%           so that only the sum of the kernels on the columns (i,j) and
%           (j,i) is determined, which is put on the column (i,j), i<=j;
%   Part 10: for kernels with decision variables, as the unknown
%           functionals of the SOS programs.
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/22/2026: Part 7: rows that are not ordered are split with
%   canonicalize_intmat before declaring the functional, as the 'intop'
%   constructor rejects them (before, they were integrated with signed limits);

% NOTES
% - Oracle: expected intmat and kernels are written out by hand per case;
%   the value check compares apply_op before vs. after combine_terms.
% - Random parts (5, 6, 9b) are seeded (rng(60000+..), rng(61000+..)).

tol = 1e-8;
pvar s_dum1 s_dum2 s_dum3
sw = @(k,a,b) subs(k,[a;b],[b;a]);

%% Part 1: d=2, one variable: the square and its chains
k = s_dum1^2+3*s_dum1*s_dum2^3;                 % not symmetric
ksw = sw(k,s_dum1,s_dum2);
mk = @(K,im) intop({K},{im},{'s'},[0,1],false);
% intmat, kernel, expected kernel on [0,3,0,1]
cases = {[0,3,0,3],k,k+ksw; [0,3,0,3],ksw,k+ksw; [0,3,0,1],k+ksw,k+ksw; [0,2,0,3],k+ksw,k+ksw; ...
         [2,3,0,3],k,k; [0,3,0,1],k,k; [0,2,0,3],k,ksw};
for j=1:size(cases,1)
    K1 = combine_terms(mk(cases{j,2},cases{j,1}),2,1,true);
    if ~isequal(K1.intmat{1},[0,3,0,1])
        error("Test_intop_combine_terms_canonical failed: Part 1, case %d: intmat %s, expected [0,3,0,1].",j,mat2str(K1.intmat{1}))
    end
    kk = kernels(K1);
    err = cleanpoly(kk{1}-cases{j,3},1e-12);
    if ~isdouble(err) || abs(double(err))>tol
        error("Test_intop_combine_terms_canonical failed: Part 1, case %d: wrong kernel.",j)
    end
    if ~same_op(combine_terms(K1,2,1,true),K1)
        error("Test_intop_combine_terms_canonical failed: Part 1, case %d: not idempotent.",j)
    end
end

%% Part 2: two different states: [0,3,0,3] is the sum of two cells, with the same kernel
K1 = combine_terms(mk(k,[0,3,0,3]),[1,1],[1,1],true(2,1));
kk = kernels(K1);
if ~isequal(K1.intmat{1},[0,2,0,3;0,3,0,1]) && ~isequal(K1.intmat{1},[0,3,0,1;0,2,0,3])
    error("Test_intop_combine_terms_canonical failed: Part 2, intmat %s.",mat2str(K1.intmat{1}))
end
for r=1:2
    err = cleanpoly(kk{r}-k,1e-12);
    if ~isdouble(err) || abs(double(err))>tol
        error("Test_intop_combine_terms_canonical failed: Part 2, kernel %d.",r)
    end
end
K1 = combine_terms(mk(k,[2,3,0,3]),[1,1],[1,1],true(2,1));
if ~isequal(K1.intmat{1},[0,3,0,1])
    error("Test_intop_combine_terms_canonical failed: Part 2, [2,3,0,3] is not [0,3,0,1].")
end

%% Part 3: the cube is one chain, with the kernel summed over all orderings
k3 = s_dum1+2*s_dum2^2+3*s_dum3^3;
K1 = combine_terms(intop({k3},{[0,4,0,4,0,4]},{'s'},[0,1],false),3,1,true);
if ~isequal(K1.intmat{1},[0,4,0,1,0,2])
    error("Test_intop_combine_terms_canonical failed: Part 3, intmat %s.",mat2str(K1.intmat{1}))
end
ksum = 0;
P = perms(1:3);
vv = [s_dum1;s_dum2;s_dum3];
for r=1:6
    ksum = ksum+subs(k3,vv,vv(P(r,:)));
end
kk = kernels(K1);
err = cleanpoly(kk{1}-ksum,1e-12);
if ~isdouble(err) || abs(double(err))>tol
    error("Test_intop_combine_terms_canonical failed: Part 3, kernel is not the sum over the orderings.")
end

%% Part 4: multipliers: [2,2,0,3] is [0,3,1,1] with the variables switched
K1 = combine_terms(mk(k,[0,3,1,1]),2,1,true);
K2 = combine_terms(mk(k,[2,2,0,3]),2,1,true);
if ~isequal(K1.intmat{1},[0,3,1,1]) || ~isequal(K2.intmat{1},[0,3,1,1])
    error("Test_intop_combine_terms_canonical failed: Part 4, intmat %s and %s, expected [0,3,1,1].",mat2str(K1.intmat{1}),mat2str(K2.intmat{1}))
end
% Both are evaluated at s2=s1: the kernel is k(s1,s1), without s_dum2
kd = subs(k,s_dum2,s_dum1);
kk1 = kernels(K1);      kk2 = kernels(K2);
err1 = cleanpoly(kk1{1}-kd,1e-12);     err2 = cleanpoly(kk2{1}-kd,1e-12);
if ~isdouble(err1) || abs(double(err1))>tol || ~isdouble(err2) || abs(double(err2))>tol
    error("Test_intop_combine_terms_canonical failed: Part 4, kernels.")
end

%% Part 5 and 6: random functionals, N in {1,2,3}, with ties; without and with multipliers
for variant = 1:2
    pdelta = 0.3*(variant==2);
    for trial = 1:25
        rng(60000+100*variant+trial);
        N = randi(3);
        d = randi([2,4]);
        if N>=2,    d = min(d,3);   end
        nterms_i = randi(3,1,N);
        vnames = arrayfun(@(a) ['s',num2str(a)],1:N,'UniformOutput',false);
        vars = struct('in',{vnames},'out',{cell(1,0)});
        dom = struct('in',repmat([0,1],N,1),'out',zeros(0,2));
        intmat = cell(1,N);
        dumvars = cell(1,N);
        for i=1:N
            intmat{i} = rand_ties_intmat(d,nterms_i(i),pdelta);
            dumvars{i} = polynomial(arrayfun(@(a) [vnames{i},'_dum',num2str(a)],(1:d)','UniformOutput',false));
        end
        % Kernels in all dummy variables, also those of multipliers
        allnames = cell(1,0);
        for i=1:N
            allnames = [allnames,dumvars{i}.varname(:)'];   %#ok<AGROW>
        end
        params = cell([nterms_i,1]);
        for lin=1:numel(params)
            params{lin} = rand_poly([1,1],allnames,1,6);
        end
        Kop = intop(params,intmat,vars,dom,false);
        K1 = combine_terms(Kop,d,1,true(1,N));
        % (a) idempotent
        if ~same_op(combine_terms(K1,d,1,true(1,N)),K1)
            error("Test_intop_combine_terms_canonical failed: Part %d, trial %d (N=%d,d=%d): not idempotent.",4+variant,trial,N,d)
        end
        % (b) all integrals are cells
        for i=1:N
            for r=1:size(K1.intmat{i},1)
                if ~is_cell_row(K1.intmat{i}(r,:))
                    error("Test_intop_combine_terms_canonical failed: Part %d, trial %d: row %s is not a cell.",4+variant,trial,mat2str(K1.intmat{i}(r,:)))
                end
            end
        end
        % (b2) the kernels do not depend on the variables of multipliers
        kk = kernels(K1);
        for lin=1:numel(kk)
            kc = cleanpoly(kk{lin},1e-12);      % a zero kernel keeps the names of its variables
            sb = cell(1,N);
            [sb{:}] = ind2sub([size(K1.params,1:N),1],lin);
            for i=1:N
                row = K1.intmat{i}(sb{i},:);
                is_dlt = row(1:2:end)==row(2:2:end);
                if any(ismember(kc.varname,dumvars{i}.varname(is_dlt(:))))
                    error("Test_intop_combine_terms_canonical failed: Part %d, trial %d (N=%d,d=%d): a kernel depends on a multiplier variable.",4+variant,trial,N,d)
                end
            end
        end
        % (c) the value is unchanged
        tvars = polynomial(arrayfun(@(a) ['t',num2str(a)],(1:N)','UniformOutput',false));
        vx = rand_poly([1,1],tvars,2,6);
        if all(isequal(vx,0)),  vx = polynomial(1);     end
        xval = 1;
        for kk=1:d
            copyvars = arrayfun(@(a) [vnames{a},'_',num2str(kk)],(1:N)','UniformOutput',false);
            xval = kron(xval,subs(vx,tvars,polynomial(copyvars)));
        end
        Y0 = apply_op(Kop,xval);
        err = cleanpoly(Y0-apply_op(K1,xval),1e-10);
        if ~isdouble(err) || norm(double(err),'fro')>tol*max(1,norm(double(cleanpoly(Y0,1e-10)),'fro'))
            error("Test_intop_combine_terms_canonical failed: Part %d, trial %d (N=%d,d=%d): the functional changed.",4+variant,trial,N,d)
        end
        % (d) relabeling the copies does not change the result
        perm = randperm(d);
        K2 = combine_terms(reorder_states(Kop,perm,ones(1,d),true(d,N)),d,1,true(1,N));
        if ~same_kernels(K1,K2)
            error("Test_intop_combine_terms_canonical failed: Part %d, trial %d (N=%d,d=%d): relabeled copies give another form.",4+variant,trial,N,d)
        end
    end
end

%% Part 7: rows that are not ordered: [3,2,0,4,0,4] is s1 from s3 to s2, nonzero only for s3<=s1<=s2
pvar s s_1 s_2 s_3
vb = {1+s+2*s^2, 2-s+s^3};
rows7 = {[3,2,0,3,0,4],[2,3,0,4,0,4],[3,2,0,4,0,4],[0,4,3,1,0,4]};
degs7 = {3,3,3,[2,1]};
for j=1:numel(rows7)
    for r=1:numel(degs7)
        dg = degs7{r};
        if numel(dg)==1,    state_idcs = [1,1,1];   else,   state_idcs = [1,1,2];  end
        k7 = rand_poly([1,1],{'s_dum1','s_dum2','s_dum3'},1,5);
        % The constructor rejects rows that are not ordered
        try
            [im7,p7] = intop.canonicalize_intmat(rows7(j),{k7});
        catch
            % [3,2,0,3,0,4]: s2<=s3, so that s3<=s1<=s2 only on a null set: the functional is zero
            if j~=1,    error('Test_intop_combine_terms_canonical failed: Part 7, row %s: empty region.',mat2str(rows7{j})),   end
            continue
        end
        Kop = intop(p7,im7,{'s'},[0,1],false);
        K1 = combine_terms(Kop,dg,ones(1,numel(dg)),true(numel(dg),1));
        if ~same_op(combine_terms(K1,dg,ones(1,numel(dg)),true(numel(dg),1)),K1)
            error("Test_intop_combine_terms_canonical failed: Part 7, row %s, degs %s: not idempotent.",mat2str(rows7{j}),mat2str(dg))
        end
        for q=1:size(K1.intmat{1},1)
            if ~is_cell_row(K1.intmat{1}(q,:))
                error("Test_intop_combine_terms_canonical failed: Part 7, row %s: %s is not a cell.",mat2str(rows7{j}),mat2str(K1.intmat{1}(q,:)))
            end
        end
        sv = [s_1,s_2,s_3];
        xval = 1;
        for c=1:3
            xval = kron(xval,subs(vb{state_idcs(c)},s,sv(c)));
        end
        Y0 = apply_op(Kop,xval);
        err = cleanpoly(Y0-apply_op(K1,xval),1e-10);
        if ~isdouble(err) || abs(double(err))>tol*max(1,abs(double(cleanpoly(Y0,1e-10))))
            error("Test_intop_combine_terms_canonical failed: Part 7, row %s, degs %s: the functional changed.",mat2str(rows7{j}),mat2str(dg))
        end
    end
end
%% Part 8: groups of equal variables: the lowest variable of the group is free, the others are evaluated at it
% [0,4,1,1,1,1]: s2=s3=s1; [3,3,3,3,0,4]: s1=s3, s2=s3, s3 free (the same group, other representative)
k3 = s_dum1+2*s_dum2^2+3*s_dum3^3;
kd3 = subs(k3,[s_dum2;s_dum3],[s_dum1;s_dum1]);             % s1+2*s1^2+3*s1^3
forms8 = {[0,4,1,1,1,1],[3,3,3,3,0,4],[0,4,1,1,2,2],[0,4,3,3,1,1]};
for j=1:numel(forms8)
    for r=1:2
        if r==1,    dg = 3;         dm = 1;         dp = true;
        else,       dg = [1,1,1];   dm = [1,1,1];   dp = true(3,1);    end
        K1 = combine_terms(intop({k3},{forms8{j}},{'s'},[0,1],false),dg,dm,dp);
        if ~isequal(K1.intmat{1},[0,4,1,1,1,1])
            error("Test_intop_combine_terms_canonical failed: Part 8, %s with degs %s: intmat %s, expected [0,4,1,1,1,1].",mat2str(forms8{j}),mat2str(dg),mat2str(K1.intmat{1}))
        end
        kk = kernels(K1);
        err = cleanpoly(kk{1}-kd3,1e-12);
        if ~isdouble(err) || abs(double(err))>tol
            error("Test_intop_combine_terms_canonical failed: Part 8, %s with degs %s: kernel is not k(s1,s1,s1).",mat2str(forms8{j}),mat2str(dg))
        end
        if ~same_op(combine_terms(K1,dg,dm,dp),K1)
            error("Test_intop_combine_terms_canonical failed: Part 8, %s: not idempotent.",mat2str(forms8{j}))
        end
    end
end
% two groups, d=4: [2,2,0,4,4,4,0,5] is s1=s2 (below the group s3=s4, s4 free): free s1 and s3
pvar s_dum4
k4 = s_dum1+2*s_dum2^2+3*s_dum3^3+s_dum4*s_dum2;
kd4 = subs(k4,[s_dum2;s_dum4],[s_dum1;s_dum3]);
K1 = combine_terms(intop({k4},{[2,2,0,4,4,4,0,5]},{'s'},[0,1],false),[1,1,1,1],[1,1,1,1],true(4,1));
kk = kernels(K1);
err = cleanpoly(kk{1}-kd4,1e-12);
if ~isequal(K1.intmat{1},[0,3,1,1,0,5,3,3]) || ~isdouble(err) || abs(double(err))>tol
    error("Test_intop_combine_terms_canonical failed: Part 8, two groups: intmat %s.",mat2str(K1.intmat{1}))
end
%% Part 9: vector-valued state (2 elements), copies equal at a point: the columns (1,2) and (2,1) are added
% kernel columns in the order (1,1),(1,2),(2,1),(2,2)
kc = {s_dum1+1, 2*s_dum1^2, -3*s_dum1, 5};
K9 = intop({[kc{1},kc{2},kc{3},kc{4}]},{[0,3,1,1]},{'s'},[0,1],false);
K1 = combine_terms(K9,2,2,true);
kk = kernels(K1);
kexp = [kc{1},kc{2}+kc{3},0,kc{4}];
err = cleanpoly(kk{1}-kexp,1e-12);
if ~isequal(K1.intmat{1},[0,3,1,1]) || ~isdouble(err) || any(abs(double(err))>tol)
    error("Test_intop_combine_terms_canonical failed: Part 9, the columns of equal copies are not combined.")
end
if ~same_op(combine_terms(K1,2,2,true),K1)
    error("Test_intop_combine_terms_canonical failed: Part 9, not idempotent.")
end
% the value does not change, and the symmetric and asymmetric kernels of the same functional agree
pvar s s_1 s_2
vx = [1+s; 2-s^2];
xval = kron(subs(vx,s,s_1),subs(vx,s,s_2));
% evaluated at s_2=s_1 through the multiplier: apply_op does it
Y0 = apply_op(K9,xval);     Y1 = apply_op(K1,xval);
kc_sym = {kc{1},kc{2}+kc{3},0,kc{4}};
K9s = intop({[kc_sym{1},kc_sym{2},kc_sym{3},kc_sym{4}]},{[0,3,1,1]},{'s'},[0,1],false);
Y2 = apply_op(K9s,xval);
if ~isdouble(cleanpoly(Y0-Y1,1e-10)) || ~isdouble(cleanpoly(Y0-Y2,1e-10)) || any(abs(double(cleanpoly(Y0-Y1,1e-10)))>tol) || any(abs(double(cleanpoly(Y0-Y2,1e-10)))>tol)
    error("Test_intop_combine_terms_canonical failed: Part 9, the functional changed.")
end
K2 = combine_terms(K9s,2,2,true);
if ~same_op(K1,K2)
    error("Test_intop_combine_terms_canonical failed: Part 9, equal functionals have different kernels.")
end

%% Part 9b: random vector-valued functionals with multipliers: idempotent, the value is unchanged, relabeling the copies does not change the result
for trial = 1:15
    rng(61000+trial);
    N = randi(2);
    d = randi([2,3]);
    nterms_i = randi(3,1,N);
    vnames = arrayfun(@(a) ['s',num2str(a)],1:N,'UniformOutput',false);
    vars = struct('in',{vnames},'out',{cell(1,0)});
    dom = struct('in',repmat([0,1],N,1),'out',zeros(0,2));
    intmat = cell(1,N);
    allnames = cell(1,0);
    for i=1:N
        intmat{i} = rand_ties_intmat(d,nterms_i(i),0.4);
        allnames = [allnames,arrayfun(@(a) [vnames{i},'_dum',num2str(a)],1:d,'UniformOutput',false)];  %#ok<AGROW>
    end
    params = cell([nterms_i,1]);
    for lin=1:numel(params)
        params{lin} = rand_poly([1,2^d],allnames,1,5);
    end
    Kop = intop(params,intmat,vars,dom,false);
    K1 = combine_terms(Kop,d,2,true(1,N));
    if ~same_op(combine_terms(K1,d,2,true(1,N)),K1)
        error("Test_intop_combine_terms_canonical failed: Part 9b, trial %d (N=%d,d=%d): not idempotent.",trial,N,d)
    end
    tvars = polynomial(arrayfun(@(a) ['t',num2str(a)],(1:N)','UniformOutput',false));
    vx = rand_poly([2,1],tvars,2,5);
    xval = 1;
    for kcopy=1:d
        copyvars = arrayfun(@(a) [vnames{a},'_',num2str(kcopy)],(1:N)','UniformOutput',false);
        xval = kron(xval,subs(vx,tvars,polynomial(copyvars)));
    end
    Y0 = apply_op(Kop,xval);
    err = cleanpoly(Y0-apply_op(K1,xval),1e-10);
    if ~isdouble(err) || norm(double(err),'fro')>tol*max(1,norm(double(cleanpoly(Y0,1e-10)),'fro'))
        error("Test_intop_combine_terms_canonical failed: Part 9b, trial %d (N=%d,d=%d): the functional changed.",trial,N,d)
    end
    perm = randperm(d);
    K2 = combine_terms(reorder_states(Kop,perm,2*ones(1,d),true(d,N)),d,2,true(1,N));
    if ~same_kernels(K1,K2)
        error("Test_intop_combine_terms_canonical failed: Part 9b, trial %d (N=%d,d=%d): relabeled copies give another form.",trial,N,d)
    end
end
%% Part 10: decision variables in the kernel of a multiplier
dpvar c1 c2
Kd = intop({c1*s_dum1+c2*s_dum2^2+3},{[0,3,1,1]},{'s'},[0,1],false);
Kc = combine_terms(Kd,2,1,true);
kk = kernels(Kc);
ed = kk{1}-(c1*s_dum1+c2*s_dum1^2+3);
if ~isequal(Kc.intmat{1},[0,3,1,1]) || (~isempty(ed.C) && max(abs(ed.C(:)))>tol)
    error("Test_intop_combine_terms_canonical failed: Part 10, kernel with decision variables on a multiplier.")
end
if ~same_op(combine_terms(Kc,2,1,true),Kc)
    error("Test_intop_combine_terms_canonical failed: Part 10, not idempotent.")
end
disp('Test_intop_combine_terms_canonical passed.')


%% ---- local functions ----------------------------------------------------
function tf = same_op(A,B)
% Same intmat, monomials, decision variables, and parameters
tf = isequal(A.intmat,B.intmat) && isequal(A.Z.in,B.Z.in) && isequal(A.Zd,B.Zd) ...
    && isequal(size(A.params),size(B.params)) && all(cellfun(@isequal,A.params(:),B.params(:)));
end

function tf = same_kernels(A,B)
% Same intmat and the same kernels (the order of the monomials Z.in is irrelevant)
tf = isequal(A.intmat,B.intmat) && isequal(size(A.params),size(B.params));
if ~tf
    return
end
ka = kernels(A);
kb = kernels(B);
for k=1:numel(ka)
    e = cleanpoly(ka{k}-kb{k},1e-10);
    if ~isdouble(e) || any(any(abs(double(e))>1e-9))
        tf = false;
        return
    end
end
end

function tf = is_cell_row(row)
% The variables that are not multipliers form a chain, each integrated from a up
% to the next, the largest over [a,b], and each multiplier is evaluated at one of them
d = numel(row)/2;
L = row(1:2:end);       U = row(2:2:end);
is_dlt = L==U;
fre = find(~is_dlt);
up = U(~is_dlt);
tf = all(L(~is_dlt)==0) && all(ismember(L(is_dlt),fre)) && nnz(up==d+1)==1 ...
    && numel(unique(up))==numel(up) && all(ismember(up(up<=d),fre));
end
