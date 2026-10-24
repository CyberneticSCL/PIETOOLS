function Test_intop_plus()
% Test_intop_plus.m -- one local function for each case:
%   case_basic              plus for 'intop'
%   case_factor_order       plus on two 'intop' objects whose dummy copy #k
%                           represents different factors of the monomial
%                           they act on
%   case_N0                 the sum of 'intop' objects without spatial
%                           variables (N==0: a plain constant/polynomial
%                           multiplier, no integration)
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/23/2026: Merge Test_intop_plus_factor_order,
%   Test_intop_plus_N0 into this test, as cases;

cases = {@case_basic, @case_factor_order, @case_N0};
nfail = 0;
for i=1:numel(cases)
    tstart = tic;
    try
        cases{i}();
        fprintf('PASS  %s (%.1f s)\n',func2str(cases{i}),toc(tstart));
    catch err
        nfail = nfail+1;
        fprintf('FAIL  %s\n      %s\n',func2str(cases{i}),err.message);
        for k=1:min(3,numel(err.stack))
            fprintf('        in %s, line %d\n',err.stack(k).name,err.stack(k).line);
        end
    end
end
if nfail>0
    error('Test_intop_plus: %d of %d cases failed.',nfail,numel(cases));
end
fprintf('Test_intop_plus: all %d cases passed.\n',numel(cases));
end

%% --------------------------------------------------------------------------
function case_basic()
% case_basic -- Tests plus for 'intop'.
%
% For random functionals Aop, Bop on N variables, (Aop+Bop)(v) is compared
% to Aop(v)+Bop(v) on a random polynomial v, both evaluated with apply_op.
%

    m = 2;  n = 3;      % row and column dimensions of functional
    N = 3;              % number of variables
    dom = randi([-2,2],[N,1]) +[zeros(N,1),2*rand([N,1])];
    degs = randi(4,[1,N]);
    pdegs = 3;

    p = 2;              % column dimension of function
    pdegv = 2;

    Aop = rand_intop([m,n],N,dom,degs,pdegs);
    Bop = rand_intop([m,n],N,dom,degs,pdegs);

    Cop = Aop+Bop;

    pvars = Aop.pvars.in;
    pvars_full = polynomial(zeros(0,1));
    for i=1:N
        pvars_full = [pvars_full; pvars{i}(:)];
    end
    v = rand_poly([n,p],pvars_full,pdegv,50);

    Av = apply_op(Aop,v);
    Bv = apply_op(Bop,v);
    Cv = apply_op(Cop,v);

    % Residual must be constant and, relative to |Av+Bv|, below tol
    tol = 1e-8;
    ApBv = Av+Bv;
    err_C = cleanpoly(Cv - ApBv,1e-10);
    nrmref = norm(full(ApBv.C),'fro');
    if ~isdouble(err_C)
        error("At least one of the outputs is not constant.")
    elseif norm(double(err_C),'fro')>tol*max(1,nrmref)
        error("Outputs do not match: residual exceeds tolerance.")
    else
        disp('err_C = ')
        disp(double(err_C))
    end
end

%% --------------------------------------------------------------------------
function case_factor_order()
% case_factor_order -- Tests plus on two 'intop' objects whose
% dummy copy #k represents different factors of the monomial they act on.
%
% @intop/plus.m only checks that the operands share variable names and
% dummy counts (isequal(A.vars,B.vars) && isequal(A.depmat,B.depmat));
% nothing records which factor each dummy copy stands for. Callers that sum
% 'intop' objects (piesos_sosvar, @polyopvar/combine_terms) must
% canonicalize the operands to a common per-factor order first, e.g. with
% reorder_vars. Two operands with different per-factor splits of the same
% 3-dummy functional are built directly with the constructor, and:
%   (1) the naive sum A+B gives a wrong value;
%   (2) reorder_vars(B,{tau^{-1}}) before summing gives the right value,
%       with tau(j) the target factor that native dummy j represents.
%
% Oracle: the values are compared to a symbolic ground truth from int() of
% the kernels times the target functions, not to another intop routine.
%
%
% NOTES
% - The test fails if the naive sum A+B matches the ground truth, and with
%   "did not correctly canonicalize" if reorder_vars does not give the
%   right value for B.
% - var_order = tau^{-1}: B's native dummy j represents target
%   tau=[2,3,1], so the order passed to reorder_vars is [3,1,2].

    tol = 1e-8;

    pvar s1
    pvar s1_dum1 s1_dum2 s1_dum3
    d1 = s1_dum1; d2 = s1_dum2; d3 = s1_dum3;

    % Two distinct 3-dummy kernels. A's native dummy order matches the target factor
    % order [1,2,3]; B's native dummy j represents target tau(j)=[2,3,1] (a different
    % split of the same 3 targets).
    KA = (s1_dum1 + 2*s1_dum2^2)*(3*s1_dum3 + 1);
    KB = (s1_dum1^2 + 1)*(s1_dum2*s1_dum3 + 2);

    intmat_full = [0,4,0,4,0,4];    % d=3, each dummy integrated over full [a,b]=[0,1]
    A = intop({KA}, {intmat_full}, {'s1'}, [0,1]);
    B = intop({KB}, {intmat_full}, {'s1'}, [0,1]);

    % A+B passes the vars/depmat check of @intop/plus.m, so the hazard is silent
    if ~isequal(A.vars,B.vars) || ~isequal(A.depmat,B.depmat)
        error("Test setup failed: A and B should pass intop/plus.m's own vars/depmat check.")
    end

    % Three distinct target state functions (ground truth by symbolic int(), independent of intop)
    f1 = 2*s1^2 - s1 + 1;
    f2 = s1^2 + 3*s1 + 0.5;
    f3 = -s1^2 + 2*s1 + 2;

    % trueA: A's native order already is target order [1,2,3].
    integrand_A = KA * subs(f1,s1,d1) * subs(f2,s1,d2) * subs(f3,s1,d3);
    trueA = double(int(int(int(integrand_A,d3,0,1),d2,0,1),d1,0,1));

    % trueB: B's native dummy1,2,3 mean target 2,3,1 respectively (tau=[2,3,1]).
    integrand_B = KB * subs(f2,s1,d1) * subs(f3,s1,d2) * subs(f1,s1,d3);
    trueB = double(int(int(int(integrand_B,d3,0,1),d2,0,1),d1,0,1));
    trueSum = trueA + trueB;

    % Setup check: apply_op reproduces each operand's own ground truth when fed its
    % own (native-order) substitution
    vterm_canon = subs(f1,s1,d1)*subs(f2,s1,d2)*subs(f3,s1,d3);
    valA = double(apply_op(A,vterm_canon));
    vterm_B_own = subs(f2,s1,d1)*subs(f3,s1,d2)*subs(f1,s1,d3);
    valB_own = double(apply_op(B,vterm_B_own));
    if abs(valA-trueA)>tol*max(1,abs(trueA))
        error("Test setup failed: apply_op(A,...) does not match its own ground truth.")
    elseif abs(valB_own-trueB)>tol*max(1,abs(trueB))
        error("Test setup failed: apply_op(B,...) does not match its own ground truth.")
    end

    % ---- (1) NAIVE sum is WRONG --------------------------------------------
    % By linearity the value is A(vterm_canon)+B(vterm_canon), and vterm_canon is
    % not B's own substitution
    C_naive = A + B;
    val_naive = double(apply_op(C_naive,vterm_canon));
    if abs(val_naive-trueSum) < tol*max(1,abs(trueSum))
        error(['case_factor_order: expected the NAIVE sum to be ' ...
            'WRONG (this is exactly the hazard this test exists to catch), ' ...
            'but it matched ground truth -- the repro no longer demonstrates ' ...
            'the hazard; the kernels/targets may need to be revisited.'])
    end

    % ---- (2) Canonicalize via reorder_vars, THEN sum: CORRECT --------------
    B_canon = reorder_vars(B,{[3,1,2]});
    valB_canon = double(apply_op(B_canon,vterm_canon));
    if abs(valB_canon-trueB)>tol*max(1,abs(trueB))
        error("reorder_vars(B,{[3,1,2]}) did not correctly canonicalize B to target order.")
    end
    C_fixed = A + B_canon;
    val_fixed = double(apply_op(C_fixed,vterm_canon));
    if abs(val_fixed-trueSum) > tol*max(1,abs(trueSum))
        error(['case_factor_order failed: canonicalized sum %.10f ' ...
            'does not match ground truth %.10f'], val_fixed, trueSum)
    end

    fprintf('case_factor_order: naive sum = %.10f (WRONG, as expected; true = %.10f)\n', val_naive, trueSum);
    fprintf('case_factor_order: canonicalized sum = %.10f (matches true = %.10f)\n', val_fixed, trueSum);
    disp('case_factor_order passed.')
end

%% --------------------------------------------------------------------------
function case_N0()
% case_N0 -- Tests the sum of 'intop' objects without spatial
% variables (N==0: a plain constant/polynomial multiplier, no integration);
% see Test_apply_intop_N0.m for the single-operand case. Sections:
%  (1) sum and difference of such operators (common_intmat with N=0 empty
%      index lists);
%  (2) the 'intop' returned by int(...) of an operator without variables
%      (as in the constant term of piesos_sosvar) must have vars.in and
%      intmat of size 1x0, as isvalid requires of intop(c), so that it can
%      be added to intop(c);
%  (3) the same 1x0 shapes for an 'intop' declared with a vars struct with
%      empty fields {}.
%
% Oracle: the sums are compared with the direct matrix products (K1+K2)*v,
% not with another 'intop' routine.
%
%
% NOTES
% - Multipliers without spatial variables must be constant (they can still
%   act on a function of the free variable t); kernels are compared as
%   returned by 'kernels'.

    tol = 1e-8;
    p = 2;
    pvar t

    %%%% (1) sum and difference of constants
    a = intop(2);       b = intop(3);
    c = a+b;
    if ~isequal(c.matdim,[1,1]) || ~isscalar(c.params) || double(c.params{1})~=5
        error("intop(2)+intop(3) must be the constant 5.")
    end
    d = a-b;
    if double(d.params{1})~=-1
        error("intop(2)-intop(3) must be the constant -1.")
    end
    v1 = randn(1,p);
    Cv1 = apply_op(c,v1);
    if ~isdouble(Cv1) || norm(double(Cv1)-5*v1,'fro')>tol*max(1,norm(5*v1,'fro'))
        error("(intop(2)+intop(3)) applied to v must equal 5*v.")
    end
    % Sum with 0 and with itself
    a0 = a+0;       aa = a+a;
    if ~isequal(a0.params{1},a.params{1}) || double(aa.params{1})~=4
        error("Sum with 0 or with itself is incorrect.")
    end

    %%%% (1) matrix-valued multipliers depending on the function variable
    K1 = polynomial(randn(2,2));
    K2 = polynomial(randn(2,2));
    v2 = rand_poly([2,p],t,2,10);
    for sgn = [1,-1]
        if sgn==1
            Kop = intop(K1)+intop(K2);      ref = (K1+K2)*v2;
        else
            Kop = intop(K1)-intop(K2);      ref = (K1-K2)*v2;
        end
        Kv = apply_op(Kop,v2);
        err = cleanpoly(Kv-ref,1e-10);
        nrmref = norm(full(ref.C),'fro');
        if ~isdouble(err) || norm(double(err),'fro')>tol*max(1,nrmref)
            error("Sum (sign %d) of matrix multipliers does not match (K1+/-K2)*v.",sgn)
        end
    end

    %%%% (1) kernels with decision variables
    dpvar c1 c2
    Kd = intop(c1)+intop(2*c2-1);
    Kd_kernels = kernels(Kd);
    dif = Kd_kernels{1}-(c1+2*c2-1);
    if nnz(dif.C)~=0
        error("Sum of 'dpvar' kernels is incorrect.")
    end

    %%%% (2) operator without variables obtained by integrating a tensor-PI operator
    % The constant function 2 on [0,1], an operator from R to L2[0,1] (as the
    % constant blocks of the basis in piesos_sosvar), integrates to the constant 2
    pvars = struct('in',{cell(1,0)},'out',{{'s'}});
    pdom = struct('in',zeros(0,2),'out',[0,1]);
    Kint = int(mat2tensopvar(polynomial(2),pvars,pdom));
    if double(Kint.params{1})~=2
        error("int of the constant function 2 on [0,1] must be the constant 2.")
    elseif ~isequal(size(Kint.vars.in),[1,0]) || ~isequal(size(Kint.intmat),[1,0])
        error("int of an operator without variables must have 1x0 vars.in and intmat.")
    end
    Kc = Kint+a;
    Cv = apply_op(Kc,v1);
    if ~isdouble(Cv) || norm(double(Cv)-4*v1,'fro')>tol*max(1,norm(4*v1,'fro'))
        error("int(2)+intop(2) applied to v must equal 4*v.")
    end

    %%%% (3) declaration with empty-cell variables {}
    K0 = intop({polynomial(7)},{},struct('in',{{}},'out',{cell(1,0)}),zeros(0,2),true);
    if ~isequal(size(K0.vars.in),[1,0]) || ~isequal(size(K0.vars.out),[1,0]) || ~isequal(size(K0.intmat),[1,0])
        error("The empty vars and intmat must be stored as 1x0 cells.")
    end
    K9 = K0+a;
    if double(K9.params{1})~=9
        error("intop(7)+intop(2) must be the constant 9.")
    end
    if ~isequal(K0.vars,a.vars) || ~isequal(K0.intmat,a.intmat)
        error("Constants must have the same vars and intmat.")
    end

    %%%% Different function spaces must still be rejected
    pvar s
    Kvar = intop({polynomial(1)},{[0,2]},{'s'},[0,1]);     % int_0^1 v(s) ds
    try
        Kvar+a;
        error("case_N0:noerror","Sum of operators on different function spaces must error.")
    catch ME
        if strcmp(ME.identifier,"case_N0:noerror")
            rethrow(ME)
        end
    end

    disp('case_N0 passed: sums of operators without variables are correct.')
end
