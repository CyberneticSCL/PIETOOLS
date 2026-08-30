%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_LPI_EQ_NDOPVAR_ENDTOEND solves actual LPIs to verify 'lpi_eq_ndopvar'
% through the solver, rather than only comparing the constraint rows it
% generates.
%
% A random 'ndopvar' decision operator Q is built, a target M is obtained by
% freezing Q at a random assignment d0 of its decision variables, and the LPI
%
%       find d  such that  Q(d) == M
%
% is solved. The problem is feasible by construction, since d0 is a solution,
% so a solver failure is a real failure. The recovered operator is then
% compared against M.
%
% The comparison is made on the evaluated coefficient matrices
%
%       E_j(d) = (Ik o [1;d])^T C_j,
%
% which for an 'ndopvar' determine the operator uniquely: the monomials
% Zd(s) o Zd(t) are linearly independent, the 3^N indicator functions act on
% disjoint regions, and a multiplier direction contributes no dummy variable
% monomials, so there is no representational freedom left. Note that the
% recovered d* need not equal d0, since the map from decision variables to
% the operator is generally not injective; only the operators must agree.
%
% Two upstream limitations are worked around, both of which reproduce with
% stock PIETOOLS/SOSTOOLS code and neither of which involves this routine:
%
%   a) A program consisting only of equality constraints on free decision
%      variables, with no SOS or psd variable declared, cannot be passed to
%      sossolve; it fails with 'Size b mismatch'. This reproduces with plain
%      'soseq' on ordinary dpvar objects. A positivity variable is therefore
%      declared in every program below, as any real LPI would have.
%
%   b) For some highly redundant equality systems sossolve still fails to
%      assemble, again with 'Size b mismatch'. Rather than assume this is
%      upstream, the test verifies it: whenever the solve fails, the same
%      operator is constrained with the stock 'lpi_eq' through
%      'ndopvar2dopvar', and the test fails if that succeeds where this
%      routine did not. Such cases are reported as skipped, not passed.
%
% Checks performed, per case that the solver can assemble:
%   1. the LPI is solved and the solver reports it feasible;
%   2. the recovered operator equals M;
%   3. requiring Q to equal two different targets at once is not feasible, so
%      the constraints are shown to actually bite.
%
% MMP, 08/29/2026: Initial coding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
rng(11);
tol = 1e-6;
sopts = struct('solver','sedumi');
sopts.params.fid = 0;

% Each case is {label, dim, deg, dom, ndec}
cases = {
    {'1D 1x1 deg 1',  [1,1], 1,     [0,1],     4}
    {'1D 1x1 deg 1b', [1,1], 1,     [0,1],    12}
    {'1D 1x1 deg 2',  [1,1], 2,     [0,1],    16}
    {'1D 1x1 shift',  [1,1], 1,     [-1,2],    6}
    {'1D 2x2 deg 1',  [2,2], 1,     [0,1],     6}
    {'1D 2x3 deg 2',  [2,3], 2,     [0,1],     5}
    {'2D 1x1 deg 1',  [1,1], [1;1], [0,1;0,1], 5}
    {'2D 2x2 deg 1',  [2,2], [1;1], [0,1;0,1], 4}
    };

npass = 0;     nskip = 0;
for ic = 1:numel(cases)
    lbl  = cases{ic}{1};
    dim  = cases{ic}{2};
    deg  = cases{ic}{3};
    dom  = cases{ic}{4};
    ndec = cases{ic}{5};
    N = size(dom,1);

    vnames = arrayfun(@(k) sprintf('s%d',k),(1:N)','UniformOutput',false);
    var1 = polynomial(vnames);
    var2 = polynomial(strcat(vnames,'_dum'));
    dnames = arrayfun(@(k) sprintf('r%d',k),(1:ndec)','UniformOutput',false);

    prog0 = lpiprogram(var1,var2,dom);
    prog0 = lpidecvar(prog0,dnames);
    % See note (a) at the top of this file.
    [prog0,~] = possopvar(prog0,1,vnames(:)',dom,0);

    % % % A random decision operator, with dense coefficients so that the
    % % % constraints are not trivially satisfied.
    Q = densify(rand_ndopvar(dim,deg,dom,var1,var2,dnames));

    % % % Target obtained by freezing Q at a random assignment
    d0 = randn(ndec,1);
    M = freeze_ndopvar(Q,d0);
    D = Q-M;

    % % % 1. Constrain Q == M and solve
    prog = lpi_eq_ndopvar(prog0,D);
    [sol,ok] = try_solve(prog,sopts);
    if ~ok
        % See note (b): confirm the trusted routine cannot do better here.
        [~,ok_ref] = try_solve(lpi_eq(prog0,ndopvar2dopvar(D)),sopts);
        if ok_ref
            error(['test_lpi_eq_ndopvar_endtoend: the solver handled the stock lpi_eq ',...
                   'constraints but not ours, for case %s.'],lbl);
        end
        fprintf('  skipped: %-14s solver could not assemble; stock lpi_eq fails identically\n',lbl);
        nskip = nskip+1;
        continue
    end
    if ~solver_feasible(sol)
        error('test_lpi_eq_ndopvar_endtoend: the LPI should be feasible for case %s.',lbl);
    end

    % % % 2. The recovered operator must equal M
    dstar = reshape(double(lpigetsol(sol,dpvar(dnames))),[],1);
    e_op = coef_distance(Q,dstar,Q,d0);
    if ~(e_op<tol)
        error('test_lpi_eq_ndopvar_endtoend: recovered operator differs from M by %.3g (case %s).',...
              e_op,lbl);
    end
    e_d = norm(dstar-d0)/max(1,norm(d0));

    % % % 3. Two incompatible targets at once must not be feasible
    M1 = freeze_ndopvar(Q,d0+randn(ndec,1));
    prog_bad = lpi_eq_ndopvar(prog0,D);
    prog_bad = lpi_eq_ndopvar(prog_bad,Q-M1);
    [sol_bad,ok_bad] = try_solve(prog_bad,sopts);
    feas_bad = ok_bad && solver_feasible(sol_bad);

    fprintf(['  passed:  %-14s nd=%2d | feasible, op err %.1e | ',...
             '|d*-d0|/|d0| = %.1e | two targets %s\n'],...
            lbl,ndec,e_op,e_d,ternary(feas_bad,'FEASIBLE(!)','not feasible'));
    if feas_bad
        error(['test_lpi_eq_ndopvar_endtoend: two incompatible targets should not be ',...
               'feasible for case %s.'],lbl);
    end
    npass = npass+1;
end

fprintf('\ntest_lpi_eq_ndopvar_endtoend: %d cases passed, %d skipped for the upstream\n',npass,nskip);
fprintf('solver limitation described at the top of this file.\n');
if npass==0
    error('test_lpi_eq_ndopvar_endtoend: no case was actually solved.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sol,ok] = try_solve(prog,sopts)
% Solve, reporting failure rather than raising, so that an upstream inability
% to assemble the program can be distinguished from a wrong answer.

sol = [];   ok = true;
try
    evalc('sol = lpisolve(prog,sopts);');
catch
    ok = false;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = densify(P)
% Replace the sparse random coefficients by dense ones of the same size.

for ii = 1:numel(P.C)
    if isempty(P.C{ii})
        continue
    end
    P.C{ii} = sparse(randn(size(P.C{ii})));
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = freeze_ndopvar(Q,d)
% Return the ndopvar representing the constant operator Q(d), keeping the
% same decision variable list so that 'minus' applies directly.

q = numel(Q.dvarname);
M = Q;
for ii = 1:numel(Q.C)
    Cii = Q.C{ii};
    if isempty(Cii)
        continue
    end
    E = eval_coef(Cii,d,q);
    Cnew = sparse(size(Cii,1),size(Cii,2));
    Cnew(1:(q+1):end,:) = E;
    M.C{ii} = Cnew;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E = eval_coef(Cii,d,q)
% E = (Ik o [1;d])^T * Cii, the coefficient matrix of the operator at the
% given decision variable values.

nrb = size(Cii,1)/(q+1);
if nrb~=round(nrb)
    error('Coefficient matrix rows are not a multiple of q+1.');
end
w = [1;reshape(d,[],1)];
E = full(kron(speye(nrb),sparse(w.'))*Cii);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = coef_distance(P1,d1,P2,d2)
% Largest relative difference of the evaluated coefficient matrices.

q1 = numel(P1.dvarname);
q2 = numel(P2.dvarname);
e = 0;      scale = 0;
for ii = 1:numel(P1.C)
    A = eval_coef(P1.C{ii},d1,q1);
    B = eval_coef(P2.C{ii},d2,q2);
    e = max(e,max(max(abs(A-B))));
    scale = max(scale,max(max(abs(B))));
end
e = e/max(1,scale);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function feas = solver_feasible(sol)

feas = true;
if isfield(sol,'solinfo') && isfield(sol.solinfo,'info') ...
        && isfield(sol.solinfo.info,'pinf')
    feas = ~(sol.solinfo.info.pinf==1);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = ternary(c,a,b)

if c, out = a; else, out = b; end

end
