% Test_intop_combine_terms_ND_random.m -- Randomized test of
% @intop/combine_terms.m for a repeated state spanning N>=2 spatial
% directions, complementing the hand-picked scenarios in
% Test_intop_combine_terms_2D_twisted.m.
%
% Each trial draws N in {2,3} spatial directions, degree d in {2,3} (number
% of copies of the single state), state size dims in {1,2}, and 1-3 terms
% per direction, each term's nested integration order drawn independently
% at random (rand_valid_intmat.m). Verifies that
% apply_op(Kop,xval) == apply_op(combine_terms(Kop,...),xval) for a
% same-state input (one base formula, renamed per copy), which is the input
% family combine_terms is used on.
%
% Oracle: combine_terms only permutes copies of the same state
%   (reorder_states.m, Kop_out*(z_{so(1)} o ... o z_{so(d)}) =
%   Kop*(z_1 o ... o z_d) for any z_i and permutation 'so', see
%   Test_intop_reorder_states.m). With z_1=...=z_d=x the input is invariant
%   under permutation of its factors, so the value on x^{od} is unchanged
%   however the terms are grouped or reordered.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES
% - Value check only: it does not check the canonical form (see
%   Test_intop_combine_terms_canonical).

tol = 1e-8;

for trial = 1:20
    rng(30000+trial);
    N = randi([2,3]);
    d = randi([2,3]);
    dimval = randi([1,2]);
    if d==3 && N==3
        dimval = 1;   % keep input dimension n=dims^d manageable
    end
    nterms_i = randi([1,3],1,N);

    vnames = cell(1,N);
    for i=1:N, vnames{i} = ['s',num2str(i)]; end
    vars = struct('in',{vnames},'out',{cell(1,0)});
    dom = struct('in',repmat([0,1],N,1),'out',zeros(0,2));

    intmat = cell(1,N);
    dumvars = cell(1,N);
    for i=1:N
        intmat{i} = rand_valid_intmat(d,nterms_i(i));
        dv = cell(d,1);
        for j=1:d, dv{j} = [vnames{i},'_dum',num2str(j)]; end
        dumvars{i} = polynomial(dv);
    end
    allpvars = vertcat(dumvars{:});

    n_in = dimval^d;
    params = cell(nterms_i);
    for lin=1:numel(params)
        params{lin} = rand_poly([1,n_in],allpvars,2,8);
    end

    Kop = intop(params,intmat,vars,dom,false);
    Kop2 = combine_terms(Kop,d,dimval,true(1,N));

    % Same-state input x^{od}: one base formula, renamed per copy
    tvars = cell(1,N);
    for i=1:N, tvars{i} = ['t',num2str(i)]; end
    tvars_p = polynomial(tvars(:));
    vx = rand_poly([dimval,1],tvars_p,2,6);
    if all(isequal(vx,0)), vx = ones(dimval,1); end

    xval = 1;
    for k=1:d
        copyvars = cell(1,N);
        for i=1:N, copyvars{i} = [vnames{i},'_',num2str(k)]; end
        vcopy = subs(vx,tvars_p,polynomial(copyvars(:)));
        xval = kron(xval,vcopy);
    end

    Y0 = apply_op(Kop,xval);
    Y1 = apply_op(Kop2,xval);
    while isa(Y0,'cell'), Y0 = Y0{1}; end
    while isa(Y1,'cell'), Y1 = Y1{1}; end
    err = cleanpoly(Y0-Y1,1e-10);
    nrmref = norm(double(cleanpoly(Y0,1e-10)),'fro');
    if ~isdouble(err) || norm(double(err),'fro') > tol*max(1,nrmref)
        error("Test_intop_combine_terms_ND_random failed: trial=%d, N=%d, d=%d, dims=%d, nterms=%s.", ...
            trial,N,d,dimval,mat2str(nterms_i))
    end
end

disp('Test_intop_combine_terms_ND_random passed.')
