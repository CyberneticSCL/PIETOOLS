% Test_intop_plus_factor_order.m -- Tests plus on two 'intop' objects whose
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
% DJ (CC), 10/15/2026: Initial coding;
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
    error(['Test_intop_plus_factor_order: expected the NAIVE sum to be ' ...
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
    error(['Test_intop_plus_factor_order failed: canonicalized sum %.10f ' ...
        'does not match ground truth %.10f'], val_fixed, trueSum)
end

fprintf('Test_intop_plus_factor_order: naive sum = %.10f (WRONG, as expected; true = %.10f)\n', val_naive, trueSum);
fprintf('Test_intop_plus_factor_order: canonicalized sum = %.10f (matches true = %.10f)\n', val_fixed, trueSum);
disp('Test_intop_plus_factor_order passed.')
