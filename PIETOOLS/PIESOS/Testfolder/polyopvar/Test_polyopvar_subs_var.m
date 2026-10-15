% Test_polyopvar_subs_var.m -- exercise the subs_var (variable-name
% substitution) path of @polyopvar/subs.m, as distinct from
% Test_polyopvar_subs_fctr.m, which only exercises subs_factor
% (position-index substitution).
%
% Oracle: fnew(y) is compared with Aop applied to the same y (apply_op),
% and f's own output space must be unchanged by the substitution.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES:
% - subs_var must align state names with common_states, not common_vars:
%   common_vars also unifies the output spaces, which would overwrite the
%   (empty) vars.out of the functional f with that of val (check (1)).
% - Evaluation uses eval_polyopvar (test-only, does not use subs), not
%   apply_poly, which evaluates by subs and would make the check circular.

pvar s1
svars = {'s1'};
dom = [0,1];

% f(x) = int_0^1 x(s1) ds1 : scalar functional, linear in state x
rng(8001);
Aop = rand_intop([1,1],svars,dom,1,3);
f = polyopvar();
f.Z = 1;
f.C = tensopmat(Aop);
f.varname = {'x'};
f.dims.in = 1;          f.dims.out = 1;
f.vars.in = svars;      f.vars.out = cell(1,0);
f.dom.in = dom;         f.dom.out = zeros(0,2);
f.depmat.in = true;     f.depmat.out = zeros(1,0);
f.chkval = true;

y = polyopvar('y',svars,dom);

% fnew(y) = int_0^1 y(s1) ds1
fnew = subs(f,{'x'},y);

% (1) f's own output space must be preserved exactly (f is still a
% scalar functional, regardless of val's own output space)
if ~isequal(fnew.vars.out,cell(1,0)) || ~isequal(size(fnew.dom.out),[0,2]) ...
        || ~isequal(size(fnew.depmat.out),[1,0])
    error("Test_polyopvar_subs_var failed: substitution corrupted f's own output space (vars.out/dom.out/depmat.out).")
end

% (2) the substituted-out variable 'x' must no longer appear in Z
x_idx = find(strcmp(fnew.varname,'x'));
if ~isempty(x_idx) && any(fnew.Z(:)==x_idx)
    error("Test_polyopvar_subs_var failed: fnew.Z still references the substituted-out variable x.")
end

% (3) the result must independently validate
fnew.chkval = false;
fnew = initialize(fnew);

% (4) numerically: fnew(y) must equal Aop*y for random y
yval = randn*s1 + randn;
xvals = cell(numel(fnew.varname),1);
for kk = 1:numel(fnew.varname)
    if strcmp(fnew.varname{kk},'y')
        xvals{kk} = yval;
    else
        xvals{kk} = polynomial(0);   % unused slot (x), value irrelevant
    end
end
Yfnew = eval_polyopvar(fnew,xvals);
Yref = apply_op(Aop,yval);
err = double(Yfnew-Yref);
if norm(err) > 1e-8
    error("Test_polyopvar_subs_var failed: fnew(y) does not match Aop*y numerically (residual %.3e).",norm(err))
end

disp('Test_polyopvar_subs_var passed.')
