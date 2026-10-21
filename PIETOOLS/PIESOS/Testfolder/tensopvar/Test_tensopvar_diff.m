% TEST_TENSOPVAR_DIFF
% Tests @tensopvar/diff, for operators with 'sopvar' factors and for nested
% operators (a 'tensopvar' as factor of another). Oracle: for a tensor-PI
% operator Cop and inputs X, one for each of the (leaf) factors,
% apply_op(diff(Cop,vars,d),X) must equal the derivative of
% apply_op(Cop,X), taken by 'polynomial' calculus, stacked over the rows
% of d. apply_op applies the leaf factors with apply_sopvar, which shares
% no code with diff for sopvar or with the product rule.
% Cases, flat operators: pointwise products (otype(1)=1) of 1-3 factors
% with 1-3 terms, first and higher order, mixed partials in several
% variables, stacked rows of d, factors independent of the variable,
% errors, and single-factor operators of both otypes.
% Cases, nested operators (function test_nested): (1) random trees of
% pointwise products of up to three (top level) or two (nested) factors,
% nested to depth two, with one or two terms, factors depending on
% different subsets of the variables (1-3 variables), as the root of a
% pointwise product and as a sum of nested single factors (the form of a
% sum of monomials); first and higher order, mixed partials, stacked rows
% of d, absent variables; (2) structure of the derivative; (3) errors;
% (4) nested blocks in a 'tensopmat' and in a 'polyopvar'.
%
% Run pietools_path_update first.

% DJ (CC), 10/20/2026: Initial coding;

rng(11);
tol = 1e-8;
npass = 0;      nfail = 0;      n_nontriv = 0;      n_err = 0;
dom = [0.5,2];

%% Random pointwise products
% variables s1,s2 (shared); factors depend on random subsets of them
for trial=1:40
    J = randi(3);                       % number of factors
    K = randi(3);                       % number of terms
    m = randi(2);
    nj = randi(2,J,1);
    dims = [m*ones(J,1),nj];
    vn = {'s1','s2'};
    vars = struct('in',{vn},'out',{vn});
    % Each factor maps from/to a random nonempty subset of the variables
    depmat.in = false(J,2);     depmat.out = false(J,2);
    for j=1:J
        depmat.in(j,:) = [true,true];
        depmat.out(j,:) = [true,true];
        if rand<0.5, depmat.out(j,randi(2)) = false; end    % nonempty: >=1 left
    end
    depmat.in(:,:) = true;
    Cop = rand_tensopvar(dims,vars,dom,depmat,K,[1,2],2,0.3);
    % Remove multipliers in the differentiated variables, which are not PI derivatives
    Cop = zero_mult_all(Cop,{'s1','s2'});
    % Distributed monomial
    X = cell(1,J);
    for j=1:J
        X{j} = rand_poly([nj(j),1],polynomial(Cop.ops{1,j}.vars.in(:)),2,6);
    end
    Y = apply_op(Cop,X);
    s1 = polynomial({'s1'});    s2 = polynomial({'s2'});
    % Orders: single derivatives, higher orders, mixed partials, stacked rows
    cases = {{{'s1'},1,       @(Y) diff(Y,s1)}, ...
             {{'s2'},1,       @(Y) diff(Y,s2)}, ...
             {{'s1'},2,       @(Y) diff(diff(Y,s1),s1)}, ...
             {{'s1','s2'},[1,1], @(Y) diff(diff(Y,s1),s2)}, ...
             {{'s1','s2'},[2,1], @(Y) diff(diff(diff(Y,s1),s1),s2)}, ...
             {{'s1','s2'},[1,0;0,1;1,1], @(Y) [diff(Y,s1);diff(Y,s2);diff(diff(Y,s1),s2)]}, ...
             {{'s2','s1'},[1,0;0,1;1,1], @(Y) [diff(Y,s2);diff(Y,s1);diff(diff(Y,s1),s2)]}, ...      % unsorted names
             {{'s1'},[1;1], @(Y) [diff(Y,s1);diff(Y,s1)]}};                                          % duplicate rows: splits shared
    for c=1:numel(cases)
        vc = cases{c}{1};   dc = cases{c}{2};   ref_fun = cases{c}{3};
        % Derivatives of order >1 in a variable require the boundary terms to vanish
        % or the multiplier parts of the factors to be absent; these may throw
        try
            dCop = diff(Cop,vc,dc);
        catch ME
            % Only second-order derivatives (cases 3 and 5) may fail, when a factor
            % has a nonvanishing boundary term; first order must always work
            if contains(ME.message,"cannot be expressed as a PI operator") && any(c==[3,5])
                n_err = n_err+1;
                continue
            else
                rethrow(ME)
            end
        end
        Yref = ref_fun(Y);
        dY = apply_op(dCop,X);
        n_nontriv = n_nontriv + (norm_poly(Yref)>tol);
        [npass,nfail] = tally(norm_poly(dY-Yref)<tol && isequal(size(dY),size(Yref)),npass,nfail, ...
            sprintf('trial %d (J=%d,K=%d), case %d',trial,J,K,c));
    end
end

%% First-order derivatives never throw for multiplier-free factors; second-order do in general
% Use J=2, one term, factors with kernels vanishing on the diagonal: d^2 is a PI operator
vars1 = struct('in',{{'s1'}},'out',{{'s1'}});
dm1 = struct('in',[0,1],'out',[0,1]);
C2 = sparse([2,1],[1,2],[1,-1],2,2);            % kernel s-t
Pk = sopvar({sparse(2,2);C2;sparse(2,2)},vars1,{[0;1]},{[0;1]},dm1,[1,1]);
Cop = sopvar2tensopvar({Pk,Pk});                % (P x)(P y), pointwise
X = {rand_poly([1,1],polynomial({'s1'}),2,3),rand_poly([1,1],polynomial({'s1'}),2,3)};
Y = apply_op(Cop,X);    s1 = polynomial({'s1'});
dY = apply_op(diff(Cop,{'s1'},2),X);
[npass,nfail] = tally(norm_poly(dY-diff(diff(Y,s1),s1))<tol,npass,nfail,'(s-t) kernels, d^2 of product');
% third derivative: error from the factor (d^3 of a single factor)
[npass,nfail] = tally(throws(@() diff(Cop,{'s1'},3),"cannot be expressed as a PI operator"),npass,nfail,'third derivative error');

%% Errors: multiplier part in a differentiated variable, bad input
Pm = rand_sopvar([1,1],vars1,[0,1],1,0.8);
Pm = add_mult(Pm);
Cm = sopvar2tensopvar({Pm,Pk});
[npass,nfail] = tally(throws(@() diff(Cm,{'s1'},1),"cannot be expressed as a PI operator"),npass,nfail,'multiplier factor error');
[npass,nfail] = tally(throws(@() diff(Cop,{'s1'},-1),"nonnegative"),npass,nfail,'negative order');
[npass,nfail] = tally(throws(@() diff(Cop,{'s1'},[1,1]),"each of the spatial"),npass,nfail,'order size');
% nested operators are differentiated recursively, see test_nested below
Cnest = tensopvar();    Cnest.ops = {Cop};      Cnest.otype = [1,2];    Cnest = initialize(Cnest);
dY = apply_op(diff(Cnest,{'s1'},1),X);
[npass,nfail] = tally(norm_poly(dY-diff(Y,s1))<tol && norm_poly(dY)>0,npass,nfail,'nested operator');
Ct2 = sopvar2tensopvar({Pk,Pk},[2,2]);
[npass,nfail] = tally(throws(@() diff(Ct2,{'s1'},1),"type-2"),npass,nfail,'type-2 multi-factor error');

%% Variable that is no output of any factor: zero operator
Cz = diff(Cop,{'s7'},1);
dY = apply_op(Cz,X);
[npass,nfail] = tally(norm_poly(dY)<tol && isequal(size(Cz),[1,1]),npass,nfail,'non-output variable');
Cz = diff(Cop,{'s7','s1'},[1,0;0,1]);
dY = apply_op(Cz,X);
[npass,nfail] = tally(norm_poly(dY-[0*Y;diff(Y,s1)])<tol,npass,nfail,'zero and nonzero rows');
% zero order returns the operator
[npass,nfail] = tally(norm_poly(apply_op(diff(Cop,{'s1'},0),X)-Y)<tol,npass,nfail,'zero order');
% default order, scalar order, char name
Cd1 = diff(Cop,{'s1'});    Cd2 = diff(Cop,'s1',1);
[npass,nfail] = tally(norm_poly(apply_op(Cd1,X)-apply_op(Cd2,X))<tol,npass,nfail,'default order, char name');

%% Single factor, both otypes: reduces to diff for sopvar
for ot=[1,2]
    P1 = zero_mult(rand_sopvar([2,1],vars1,dom,2,0.5));
    Ci = sopvar2tensopvar(P1,[ot,2]);
    x = rand_poly([1,1],polynomial({'s1'}),2,3);
    dC = diff(Ci,{'s1'},1);
    % apply_op renames the output variable s1 -> s1_1 for a type-2 product, so compare
    % through the pointwise version of the same operator
    otype_kept = isequal(dC.otype,[ot,2]);
    dC.otype = [1,2];
    [npass,nfail] = tally(otype_kept && norm_poly(apply_op(dC,x)-diff(apply_sopvar(P1,x),s1))<tol,npass,nfail,sprintf('single factor otype(1)=%d',ot));
end
% sum of terms
P2 = zero_mult(rand_sopvar([2,1],vars1,dom,2,0.5));
Cs = sopvar2tensopvar({P1;P2});
x = rand_poly([1,1],polynomial({'s1'}),2,3);
[npass,nfail] = tally(norm_poly(apply_op(diff(Cs,{'s1'},1),x)-diff(apply_sopvar(P1,x)+apply_sopvar(P2,x),s1))<tol,npass,nfail,'two terms, one factor');

% Nested operators
[npass,nfail,n_nontriv,n_err] = test_nested(npass,nfail,n_nontriv,n_err);

fprintf('\nTest_tensopvar_diff: %d passed, %d failed (%d nonzero comparisons, %d expected second-order errors).\n',npass,nfail,n_nontriv,n_err);
if nfail>0
    error('Test_tensopvar_diff failed.');
end


% TEST_NESTED tests @tensopvar/diff for nested operators, see the header
function [npass,nfail,n_nontriv,n_err] = test_nested(npass,nfail,n_nontriv,n_err)
    rng(41);
    tol = 1e-8;
    vn_all = {'s1','s2','s3'};

    %% (1) Random nested operators
    n_nested = 0;
    for trial=1:30
        nv = randi(3);
        vn = vn_all(1:nv);
        m = randi(2);
        kind = mod(trial-1,3)+1;
        switch kind
            case 1      % pointwise product of 1-3 factors, nested to depth 2, 1-2 terms
                spec = rand_spec(2,m,vn,true);
                Rop = inst(spec,[1,2]);
            case 2      % sum of two nested pointwise products, one factor per term
                spec = rand_spec(2,m,vn,true);
                Rop = tensopvar();      Rop.otype = [1,2];
                Rop.ops = {inst(spec,[1,2]);inst(spec,[1,2])};
                Rop = initialize(Rop);
            case 3      % the same with a Kronecker product (otype(1)=2) of a single factor
                spec = rand_spec(2,m,vn,true);
                Rop = tensopvar();      Rop.otype = [2,2];
                Rop.ops = {inst(spec,[1,2]);inst(spec,[1,2])};
                Rop = initialize(Rop);
        end
        is_nested = any(cellfun('isclass',Rop.ops(:),'tensopvar'));
        n_nested = n_nested+is_nested;
        % Inputs for each leaf factor, and the operator for evaluation: a single factor with
        % otype(1)=2 names the output variables s_1,..., so evaluate the pointwise version
        leaves = collect_leaves(Rop);
        X = cell(1,numel(leaves));
        for l=1:numel(leaves)
            X{l} = rand_poly([leaves{l}.dims(2),1],polynomial(leaves{l}.vars.in(:)),2,5);
        end
        Y = apply_op(as_pointwise(Rop),X);
        cases = {{vn(1),1,       @(Y) diff(Y,polynomial(vn(1)))}, ...
                 {vn(nv),1,      @(Y) diff(Y,polynomial(vn(nv)))}, ...
                 {vn(1),[1;1],   @(Y) [diff(Y,polynomial(vn(1)));diff(Y,polynomial(vn(1)))]}, ...
                 {[vn(nv),{'s9'}],[1,0;0,0;1,1], @(Y) [diff(Y,polynomial(vn(nv)));Y;0*Y]}, ...
                 {vn(1),2,       @(Y) diff(diff(Y,polynomial(vn(1))),polynomial(vn(1)))}};
        if nv>1
            sA = polynomial(vn(1));     sB = polynomial(vn(nv));
            cases = [cases, {{vn([nv,1]),[1,1],  @(Y) diff(diff(Y,sA),sB)}, ...
                             {vn([nv,1]),[1,0;0,1;1,1], @(Y) [diff(Y,sB);diff(Y,sA);diff(diff(Y,sA),sB)]}}];
        end
        for c_i=1:numel(cases)
            vc = cases{c_i}{1};     dc = cases{c_i}{2};     ref_fun = cases{c_i}{3};
            is_second = max(dc(:))>=2;
            try
                dRop = diff(Rop,vc,dc);
            catch ME
                % second-order derivatives may fail for a leaf with nonvanishing boundary terms
                if is_second && contains(ME.message,"cannot be expressed as a PI operator")
                    n_err = n_err+1;
                    continue
                else
                    rethrow(ME)
                end
            end
            Yref = ref_fun(Y);
            dY = apply_op(as_pointwise(dRop),X);
            n_nontriv = n_nontriv + (norm_poly(Yref)>tol);
            ok = norm_poly(dY-Yref)<tol && isequal(size(dY),size(Yref)) && isequal(size(dRop),[size(dc,1)*m,Rop.matdim(2)]) ...
                && isequal(dRop.otype,Rop.otype) && num_factors(dRop)==num_factors(Rop);
            [npass,nfail] = tally(ok,npass,nfail,sprintf('trial %d (kind %d, nv=%d), case %d',trial,kind,nv,c_i));
        end
    end
    fprintf('(%d of 30 random operators nested, %d nonzero comparisons)\n',n_nested,n_nontriv);

    %% (2) Structure: the nested factor of the derivative has one term per nonzero split
    vars1 = struct('in',{{'s1','s2'}},'out',{{'s1','s2'}});
    rng(5);
    A = zero_mult(rand_sopvar([1,1],vars1,[0,1],2,0.5));
    B = zero_mult(rand_sopvar([1,1],vars1,[0,1],2,0.5));
    C = zero_mult(rand_sopvar([1,1],vars1,[0,1],2,0.5));
    N = pointwise(A,B);                 % nested factor A.*B
    R = pointwise(N,C);                 % (A.*B).*C
    dR = diff(R,{'s1'},1);              % (A'.*B + A.*B').*C + (A.*B).*C'
    % two terms: N unchanged with C', and N' (two terms) with C
    nt = [size(dR.ops{1,1}.ops,1),size(dR.ops{2,1}.ops,1)];
    [npass,nfail] = tally(size(dR.ops,1)==2 && isa(dR.ops{1,1},'tensopvar') && isequal(sort(nt),[1,2]) && num_factors(dR)==3,npass,nfail,'terms of the derivative of (A.*B).*C');
    % a factor that does not depend on the differentiated variable is dropped
    varsA = struct('in',{{'s1'}},'out',{{'s1'}});       varsB = struct('in',{{'s2'}},'out',{{'s2'}});
    A1 = zero_mult(rand_sopvar([1,1],varsA,[0,1],2,0.8));
    B2 = zero_mult(rand_sopvar([1,1],varsB,[0,1],2,0.8));
    R1 = pointwise(pointwise(A1,B2),A1);
    dR1 = diff(R1,{'s2'},1);            % only B2 depends on s2: a single term A1.*B2'.*A1
    [npass,nfail] = tally(size(dR1.ops,1)==1 && size(dR1.ops{1,1}.ops,1)==1,npass,nfail,'factors independent of the variable are dropped');
    % stacked rows, some of which do not apply to a term (invalid rows of a term are zero), with a
    % zero row, for a nested factor in the first and in the last position
    X0 = {rand_poly([1,1],polynomial({'s1'}),2,4),rand_poly([1,1],polynomial({'s2'}),2,4),rand_poly([1,1],polynomial({'s1'}),2,4)};
    sA = polynomial({'s1'});    sB = polynomial({'s2'});
    % (variant 3: nested to depth three, with a nested first factor in the nested factor)
    for variant=1:3
        Xs = X0;
        if variant==1,      Rs = pointwise(pointwise(A1,B2),A1);
        elseif variant==2,  Rs = pointwise(A1,pointwise(B2,A1));
        else,               Rs = pointwise(pointwise(pointwise(A1,B2),A1),A1);
                            Xs = [X0,{rand_poly([1,1],sA,2,4)}];    end
        Ys = apply_op(Rs,Xs);
        dRs = diff(Rs,{'s1','s2'},[1,0;0,1;1,1;0,0;1,0]);
        Yref = [diff(Ys,sA);diff(Ys,sB);diff(diff(Ys,sA),sB);Ys;diff(Ys,sA)];
        [npass,nfail] = tally(norm_poly(apply_op(dRs,Xs)-Yref)<tol && norm_poly(Yref)>0 && isequal(size(dRs),[5,1]),npass,nfail,sprintf('stacked rows with invalid rows and a zero row, variant %d',variant));
    end
    % stacking does not add terms: for each of the three valid splits, a single row of d applies, and
    % the nested factor of the term is a single product (the other rows are zero)
    dR1s = diff(R1,{'s1','s2'},[1,0;0,1]);
    [npass,nfail] = tally(size(dR1s.ops,1)==3 && all(cellfun(@(F) isa(F,'tensopvar') && size(F.ops,1)==1,dR1s.ops(:,1))),npass,nfail,'terms of stacked derivative of a nested factor');
    % a nested factor that does not depend on the variable gives no term
    R2 = pointwise(pointwise(A1,A1),B2);
    dR2 = diff(R2,{'s2'},1);
    [npass,nfail] = tally(size(dR2.ops,1)==1,npass,nfail,'zero derivative of a nested factor gives no term');
    % second order: splits over a nested factor have multinomial coefficient 2. Kernels (s-t) vanish
    % on the diagonal, so that second derivatives are PI operators
    vars_s = struct('in',{{'s1'}},'out',{{'s1'}});      dom_s = struct('in',[0,1],'out',[0,1]);
    Pk = sopvar({sparse(2,2);sparse([2,1],[1,2],[1,-1],2,2);sparse(2,2)},vars_s,{[0;1]},{[0;1]},dom_s,[1,1]);
    Xk = {rand_poly([1,1],sA,2,3),rand_poly([1,1],sA,2,3),rand_poly([1,1],sA,2,3)};
    for variant=1:2
        if variant==1,  Rk = pointwise(pointwise(Pk,Pk),Pk);
        else,           Rk = pointwise(Pk,pointwise(Pk,Pk));    end
        Yk = apply_op(Rk,Xk);
        dRk = diff(Rk,{'s1'},[2;1]);
        Yref = [diff(diff(Yk,sA),sA);diff(Yk,sA)];
        [npass,nfail] = tally(norm_poly(apply_op(dRk,Xk)-Yref)<tol && norm_poly(Yref)>0,npass,nfail,sprintf('second order over a nested factor, variant %d',variant));
    end
    % no factor depends on the variable: zero operator of the same size and structure
    dR0 = diff(R1,{'s9'},1);
    [npass,nfail] = tally(isequal(size(dR0),size(R1)) && norm_poly(apply_op(dR0,X0))==0 && num_factors(dR0)==num_factors(R1),npass,nfail,'absent variable: zero operator');
    % zero order: the operator itself, nested or not
    dRz = diff(R1,{'s1'},0);
    [npass,nfail] = tally(norm_poly(apply_op(dRz,X0)-apply_op(R1,X0))<tol && norm_poly(apply_op(dRz,X0))>0,npass,nfail,'zero order');
    dRz2 = diff(R1,{'s1','s2'},[0,0;1,0]);
    Yz = apply_op(R1,X0);
    [npass,nfail] = tally(norm_poly(apply_op(dRz2,X0)-[Yz;diff(Yz,polynomial({'s1'}))])<tol,npass,nfail,'zero row stacked with nonzero row');

    %% (3) Errors
    % type-2 (Kronecker) product of the outputs of several factors, at the top or nested
    T2 = tensopvar();   T2.otype = [2,2];   T2.ops = {A1,B2};   T2 = initialize(T2);
    [npass,nfail] = tally(throws(@() diff(T2,{'s1'},1),"type-2"),npass,nfail,'type-2 product at the top level');
    Rn2 = tensopvar();  Rn2.otype = [1,2];  Rn2.ops = {T2,A1};  Rn2 = initialize(Rn2);
    [npass,nfail] = tally(throws(@() diff(Rn2,{'s1'},1),"type-2"),npass,nfail,'type-2 product nested');
    % a multiplier part in a differentiated variable of a nested factor is not a PI derivative
    Am = add_mult(A1);
    Rm = pointwise(pointwise(Am,B2),A1);
    [npass,nfail] = tally(throws(@() diff(Rm,{'s1'},1),"cannot be expressed as a PI operator"),npass,nfail,'multiplier part in nested factor');
    % ...but is fine for a variable the nested factor does not depend on
    dRm = diff(Rm,{'s2'},1);
    [npass,nfail] = tally(norm_poly(apply_op(dRm,X0)-diff(apply_op(Rm,X0),polynomial({'s2'})))<tol,npass,nfail,'multiplier in other variable');
    % input checks of the top level apply to nested operators
    [npass,nfail] = tally(throws(@() diff(R1,{'s1'},-1),"nonnegative"),npass,nfail,'negative order');
    [npass,nfail] = tally(throws(@() diff(R1,{'s1'},[1,1]),"each of the spatial"),npass,nfail,'order size');

    %% (4) Nested blocks of a tensopmat, and nested operators in a polyopvar
    rng(77);
    for nv=1:3
        vn = vn_all(1:nv);      sv = polynomial(vn);
        for trial=1:3
            m1 = randi(2);      n2 = randi(2);
            spec = rand_spec(2,m1,vn,true);
            T = tensopvar();    T.otype = [1,2];    T.ops = {inst(spec,[1,2]);inst(spec,[1,2])};     T = initialize(T);
            % the variables of the block on R must be output variables of T
            M = rand_poly([m1,n2],polynomial(T.vars.out(:)),2,4);
            Cop = tensopmat({T,M});
            leaves = collect_leaves(T);
            X1 = cell(1,numel(leaves));
            for l=1:numel(leaves)
                X1{l} = rand_poly([leaves{l}.dims(2),1],polynomial(leaves{l}.vars.in(:)),2,5);
            end
            x2 = randn(n2,1);
            YT = apply_op(tensopmat({T}),{X1});
            sA = polynomial(vn(1));     sB = polynomial(vn(nv));
            vc = vn([nv,1]);            dc = [1,0;0,1;1,1];
            if nv==1,  vc = vn(1);      dc = [1;1];     end
            dCop = diff(Cop,vc,dc);
            dY = apply_op(dCop,{X1,x2});
            if nv==1
                Yref = [diff(YT,sA)+diff(M,sA)*x2; diff(YT,sA)+diff(M,sA)*x2];
            else
                Yref = [diff(YT,sB)+diff(M,sB)*x2; diff(YT,sA)+diff(M,sA)*x2; diff(diff(YT,sA),sB)+diff(diff(M,sA),sB)*x2];
            end
            n_nontriv = n_nontriv + (norm_poly(Yref)>tol);
            [npass,nfail] = tally(norm_poly(dY-Yref)<tol && isequal(size(dY),size(Yref)),npass,nfail,sprintf('nested block in tensopmat, nv=%d, trial %d',nv,trial));
        end
    end
    for trial=1:12
        N = randi(3);
        vn = vn_all(1:N);       sv = polynomial(vn);
        ma = randi(2);          mb = randi(2);
        f = build_f_nested(vn,ma,mb);
        xa = rand_poly([ma,1],sv,2,5);      xb = rand_poly([mb,1],sv,2,5);
        Y = eval_polyopvar(f,{xa;xb});
        s_first = polynomial(vn(1));        s_last = polynomial(vn(N));
        cases = {{vn(1),1,       @(Y) diff(Y,s_first)}, ...
                 {vn(N),1,       @(Y) diff(Y,s_last)}, ...
                 {vn(1),[1;1],   @(Y) [diff(Y,s_first);diff(Y,s_first)]}};
        for c_i=1:numel(cases)
            vc = cases{c_i}{1};     dc = cases{c_i}{2};     ref_fun = cases{c_i}{3};
            df = diff(f,vc,dc);
            Yref = ref_fun(Y);
            dY = eval_polyopvar(df,{xa;xb});
            n_nontriv = n_nontriv + (norm_poly(Yref)>tol);
            [npass,nfail] = tally(norm_poly(dY-Yref)<tol && isequal(size(dY),size(Yref)) && df.chkval,npass,nfail,sprintf('nested polyopvar, trial %d (N=%d), case %d',trial,N,c_i));
        end
    end
end


%% Helpers
function e = norm_poly(p)
if isa(p,'double')
    e = max([0;abs(p(:))]);
    return
end
p = cleanpoly(polynomial(p),1e-12);
if isempty(p.C) || nnz(p.C)==0
    e = 0;
else
    e = max(abs(p.C(:)));
end
end

function [npass,nfail] = tally(ok,npass,nfail,name)
if ok
    npass = npass+1;
else
    nfail = nfail+1;    fprintf('FAIL: %s\n',name);
end
end

function tf = throws(f,msg)
tf = false;
try
    f();
catch ME
    tf = contains(ME.message,msg);
end
end

function P = zero_mult(P)
% Remove the parameters with a multiplier in some variable: derivatives of PI operators
vS3 = intersect(P.vars.in,P.vars.out);
params = P.params;
for k=1:numel(params)
    g = gamma_of_cell(k,numel(vS3));
    if any(g==1)
        params{k} = sparse(size(params{k},1),size(params{k},2));
    end
end
P = sopvar(params,P.vars,P.ZL,P.ZR,P.dom,P.dims);
end

function P = add_mult(P)
% Add a multiplier part to the 1D sopvar P
params = P.params;
params{1}(1,1) = 1;
P = sopvar(params,P.vars,P.ZL,P.ZR,P.dom,P.dims);
end

function Cop = zero_mult_all(Cop,vn)
% Zero all parameters of all factors with a multiplier in one of the variables vn
for l=1:numel(Cop.ops)
    P = Cop.ops{l};
    vS3 = intersect(P.vars.in,P.vars.out);
    params = P.params;
    for k=1:numel(params)
        g = gamma_of_cell(k,numel(vS3));
        if any(g(ismember(vS3,vn))==1)
            params{k} = sparse(size(params{k},1),size(params{k},2));
        end
    end
    Cop.ops{l} = sopvar(params,P.vars,P.ZL,P.ZR,P.dom,P.dims);
end
Cop = initialize(Cop);
end

function spec = rand_spec(depth,m,vn,is_root)
% Random tree: a leaf (m x n 'sopvar' with random input and output variables) or a pointwise
% product of 1-3 factors, summed over 1-2 terms, each factor again a tree. The root is a node
if ~is_root && (depth==0 || rand<0.35)
    nv = numel(vn);
    sin = vn(rand(1,nv)<0.6);       if isempty(sin),   sin = vn(randi(nv));    end
    if rand<0.5,    sout = sin;
    else,           sout = vn(rand(1,nv)<0.6);     if isempty(sout),  sout = vn(randi(nv));   end
    end
    spec = struct('leaf',true,'m',m,'n',randi(2),'in',{sin},'out',{sout},'kids',{{}},'K',1);
else
    J = randi(2+is_root);
    kids = cell(1,J);
    for j=1:J
        kids{j} = rand_spec(depth-1,m,vn,false);
    end
    spec = struct('leaf',false,'m',m,'n',0,'in',{{}},'out',{{}},'kids',{kids},'K',randi(2));
end
end

function R = inst(spec,otype)
% Random operator for the tree spec: all terms have the same structure
if spec.leaf
    vars = struct('in',{spec.in},'out',{spec.out});
    R = zero_mult(rand_sopvar([spec.m,spec.n],vars,[0,1],2,0.5));
    return
end
J = numel(spec.kids);
ops = cell(spec.K,J);
for k=1:spec.K
    for j=1:J
        ops{k,j} = inst(spec.kids{j},otype);
    end
end
R = tensopvar();    R.otype = otype;    R.ops = ops;    R = initialize(R);
end

function leaves = collect_leaves(R)
% The leaf 'sopvar' operators of the first term of R, in the order of the inputs of apply_op
leaves = {};
for j=1:size(R.ops,2)
    F = R.ops{1,j};
    if isa(F,'tensopvar')
        leaves = [leaves,collect_leaves(F)];
    else
        leaves = [leaves,{F}];
    end
end
end

function R = as_pointwise(R)
% apply_op names the output variables of a Kronecker product of the outputs of a single factor
% s_1,..., so evaluate the pointwise version of the same operator
R.otype = [1,R.otype(2)];
end

function N = pointwise(A,B)
% Pointwise product of the outputs of the operators A and B
N = tensopvar();        N.otype = [1,2];        N.ops = {A,B};      N = initialize(N);
end

function f = build_f_nested(vn,ma,mb)
% f = C*Z(x), Z = [1,0; 1,2], monomials a and a o b, as in Test_polyopvar_diff, but with the
% operators acting on each monomial a sum of two nested operators of the same structure
N = numel(vn);
vars = struct('in',{vn},'out',{vn});
dep1 = struct('in',true(1,N),'out',true(1,N));
dep2 = struct('in',true(2,N),'out',true(2,N));
T1a = zero_mult_all(rand_tensopvar([2,ma],vars,[0,1],dep1,2,[1,2],2,0.3),vn);
T1b = zero_mult_all(rand_tensopvar([2,ma],vars,[0,1],dep1,2,[1,2],2,0.3),vn);
T2a = zero_mult_all(rand_tensopvar([2,ma;2,mb],vars,[0,1],dep2,2,[1,2],2,0.3),vn);
T2b = zero_mult_all(rand_tensopvar([2,ma;2,mb],vars,[0,1],dep2,2,[1,2],2,0.3),vn);
T1 = tensopvar();   T1.otype = [1,2];   T1.ops = {T1a;T1b};     T1 = initialize(T1);
T2 = tensopvar();   T2.otype = [1,2];   T2.ops = {T2a;T2b};     T2 = initialize(T2);
f = polyopvar();
f.Z = [1,0;1,2];            f.varname = {'a';'b'};
f.dims.in = [ma;mb];
f.vars.in = vn;             f.dom.in = repmat([0,1],N,1);
f.depmat.in = true(2,N);
f.C = tensopmat({T1,T2});
f = initialize(f);
end
