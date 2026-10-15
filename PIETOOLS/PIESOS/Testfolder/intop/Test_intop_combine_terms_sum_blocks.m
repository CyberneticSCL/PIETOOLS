% Test_intop_combine_terms_sum_blocks.m -- Tests that
% @intop/combine_terms.m, which sums the blocks of terms with the same
% order of the copies of a state in one pass, gives a functional that
%   (1) represents the same functional as the input: the same value on
%       x o ... o x, for kernels that are polynomials or depend on decision
%       variables (dpvar);
%   (2) keeps the decision variables Zd of the input (the blocks share
%       them, they are not compared or merged);
%   (3) has unique, sorted integrals, and is left unchanged by a second
%       call;
%   (4) does not depend on the order of the terms of the input: permuting
%       the rows of intmat (and the parameters accordingly) gives the same
%       functional (the same value, and the same intmat after combining).
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES
% - Does not compare the one-pass sum of the blocks with a sum of the
%   blocks by plus; it checks only the properties (1)-(4).

tol = 1e-8;

%% Part 1: N=2, degree 2, twisted terms with decision variables
vars = struct('in',{{'s1','s2'}},'out',{cell(1,0)});
dom = struct('in',[0,1;0,1],'out',zeros(0,2));
pvar s1_dum1 s1_dum2 s2_dum1 s2_dum2
dpvar c1 c2
K1m = c1*s1_dum1*s2_dum1^2 + 5*s1_dum2*s2_dum2^3 - 2*c2*s1_dum2^2*s2_dum1;
K2m = 3*s1_dum1^2*s2_dum2 - c2*s1_dum2*s2_dum1^2 + 7*c1;
K3m = s1_dum1 + c2*s2_dum2 + s1_dum2*s2_dum1;
cases = {{[0,2,0,3;0,3,0,1],[0,3,0,1],{K1m;K2m}}, ...
         {[0,2,0,3;0,3,0,1;0,3,0,3],[0,3,0,1],{K1m;K2m;K3m}}, ...
         {[0,3,0,3],[0,3,0,1;0,2,0,3],{K1m,K2m}}, ...
         {[0,2,0,3],[0,2,0,3;0,3,0,1],{K1m,K2m}}};
pvar s t s1_1 s2_1 s1_2 s2_2
vx = s^2*t + 3*t^3 - s*t^2 + 1;
xval = subs(vx,[s;t],[s1_1;s2_1]) * subs(vx,[s;t],[s1_2;s2_2]);
for j=1:numel(cases)
    Kop = intop(cases{j}{3},cases{j}(1:2),vars,dom,false);
    K1 = combine_terms(Kop,2,1,true(1,2));
    check_result(Kop,K1,xval,tol,sprintf("Part 1, case %d",j));
end

%% Part 2: random N in {1,2,3}, several rows per variable, dpvar and polynomial kernels
dvn = {'c1';'c2';'c3'};
for trial = 1:40
    rng(80000+trial);
    N = randi(3);
    d = randi([2,4]);
    if N>=2,    d = min(d,3);   end
    nterms_i = randi(5,1,N);
    if N==2,    nterms_i = min(nterms_i,4);     end
    if N==3,    nterms_i = min(nterms_i,3);     end
    vnames = cell(1,N);     intmat = cell(1,N);     dumvars = cell(1,N);
    for i=1:N
        vnames{i} = ['s',num2str(i)];
        intmat{i} = rand_ties_intmat(d,nterms_i(i),0.15);
        dv = cell(d,1);
        for j=1:d,  dv{j} = [vnames{i},'_dum',num2str(j)];  end
        dumvars{i} = polynomial(dv);
    end
    vars = struct('in',{vnames},'out',{cell(1,0)});
    dom = struct('in',repmat([0,1],N,1),'out',zeros(0,2));
    allpvars = vertcat(dumvars{:});
    nt = cellfun(@(a) size(a,1),intmat);
    params = cell([nt,1]);
    use_dp = mod(trial,2)==0;
    for lin=1:numel(params)
        if use_dp
            params{lin} = rand_dpvar([1,1],dvn,allpvars,1,6);
        else
            params{lin} = rand_poly([1,1],allpvars,1,6);
        end
    end
    Kop = intop(params,intmat,vars,dom,false);
    K1 = combine_terms(Kop,d,1,true(1,N));

    tvars_p = polynomial(cellfun(@(a) ['t',num2str(a)],num2cell(1:N)','UniformOutput',false));
    vx = rand_poly([1,1],tvars_p,2,6);
    if all(isequal(vx,0)),  vx = 1+tvars_p(1);  end
    xval = 1;
    for k=1:d
        copyvars = cellfun(@(a) [a,'_',num2str(k)],vnames(:),'UniformOutput',false);
        xval = xval*subs(vx,tvars_p,polynomial(copyvars));
    end
    msg = sprintf("Part 2, trial=%d, N=%d, d=%d",trial,N,d);
    check_result(Kop,K1,xval,tol,msg);

    % (4): the order of the terms does not matter
    pm = arrayfun(@randperm,nt,'UniformOutput',false);
    if N==1,    ps = {pm{1},1};     else,   ps = pm;    end
    imat_p = cellfun(@(a,p) a(p,:),intmat,pm,'UniformOutput',false);
    Kperm = intop(params(ps{:}),imat_p,vars,dom,false);
    K1p = combine_terms(Kperm,d,1,true(1,N));
    if ~isequal(K1p.intmat,K1.intmat) || ~isequal(K1p.Z.in,K1.Z.in)
        error("Test_intop_combine_terms_sum_blocks failed: %s: result depends on the order of the terms.",msg)
    end
    for k=1:numel(K1.params)
        if norm(K1p.params{k}-K1.params{k},'fro')>1e-10*max(1,norm(K1.params{k},'fro'))
            error("Test_intop_combine_terms_sum_blocks failed: %s: parameters depend on the order of the terms.",msg)
        end
    end
end

disp('Test_intop_combine_terms_sum_blocks passed.')


%% ---- local functions ----------------------------------------------------
function check_result(Kop,K1,xval,tol,msg)
% (1) same value on xval, (2) same decision variables, (3) unique sorted
% integrals and idempotent
Y0 = apply_op(Kop,xval);
Y1 = apply_op(K1,xval);
err = Y0-Y1;
if isa(err,'dpvar')
    if ~isempty(err.C) && max(abs(err.C(:)))>tol*max(1,max(abs(Y0.C(:))))
        error("Test_intop_combine_terms_sum_blocks failed: %s: the functional changed.",msg)
    end
else
    err = cleanpoly(err,1e-10);
    nrmref = norm(double(cleanpoly(Y0,1e-10)),'fro');
    if ~isdouble(err) || norm(double(err),'fro')>tol*max(1,nrmref)
        error("Test_intop_combine_terms_sum_blocks failed: %s: the functional changed.",msg)
    end
end
if ~isequal(K1.Zd,Kop.Zd)
    error("Test_intop_combine_terms_sum_blocks failed: %s: decision variables changed.",msg)
end
for i=1:numel(K1.intmat)
    [~,imat_u] = uniquerows_integerTable(K1.intmat{i});
    if size(imat_u,1)~=size(K1.intmat{i},1)
        error("Test_intop_combine_terms_sum_blocks failed: %s: integrals are not unique.",msg)
    end
end
N = numel(K1.intmat);
K2 = combine_terms(K1,max(K1.depmat.in),1,true(1,N));
if ~isequal(K2.intmat,K1.intmat) || ~isequal(K2.Z.in,K1.Z.in) || ~isequal(size(K2.params),size(K1.params))
    error("Test_intop_combine_terms_sum_blocks failed: %s: not unchanged by a second call.",msg)
end
for k=1:numel(K1.params)
    if ~isequal(K2.params{k},K1.params{k})
        error("Test_intop_combine_terms_sum_blocks failed: %s: parameters changed by a second call.",msg)
    end
end
end
