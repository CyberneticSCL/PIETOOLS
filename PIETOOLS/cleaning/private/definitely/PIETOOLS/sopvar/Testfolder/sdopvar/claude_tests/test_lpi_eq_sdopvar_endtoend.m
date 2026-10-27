%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_LPI_EQ_SDOPVAR_ENDTOEND solves actual LPIs to verify the whole chain:
% possopvar declares a positive operator Q, lpi_eq_sdopvar constrains Q to
% equal a known target M, the program is solved, and the recovered operator
% is compared against M.
%
% The target is built as M = R*R, for R a random operator of the form
%
%       (R x)(theta) = W * (Z x)(theta),
%
% where Z is the stack of basis operators that possopvar itself uses,
%
%       (Z_alpha x)(theta) = int_s I_alpha(theta-s) Z^alpha(theta,s) x(s) ds,
%
% and W is a random matrix. Then
%
%       M = R*R = Z* (W'W) Z,
%
% which is positive semidefinite and, being of exactly the form possopvar
% parameterizes, is attainable. That makes the LPI feasible, so a solver
% failure is a real failure rather than a badly posed problem. The identity
% <x,Mx> = ||R x||^2 is verified symbolically before each solve, so M really
% is R*R and not merely something asserted to be.
%
% Note that R is NOT built with '@sopvar/mtimes'. That method currently
% errors ('Size inputs must be integers') on a random 1D sopvar from
% 'rand_sopvar', so composing two sopvar objects is not available; the
% construction above reaches the same operator without it.
%
% Checks performed, per case:
%   1. <x,Mx> = ||R x||^2, confirming M = R*R;
%   2. the LPI is solved and the solver reports it feasible;
%   3. the recovered operator equals M, compared as operators via <y,.x> for
%      random polynomial test functions, not by comparing coefficients;
%   4. the recovered Gram matrix Q is positive semidefinite;
%   5. the same solve with the 'symmetric' option recovers the same operator,
%      which exercises the adjoint pairing (Q-M is self-adjoint here);
%   6. constraining Q = -M instead is detected as infeasible, so that the
%      constraints are shown to actually bite.
%
% Two spatial variables are included, so the multi-index contraction of the
% multiplier directions is exercised through a real solve and not only in
% isolation.
%
% MMP, 08/29/2026: Initial coding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
rng(5);
tol = 1e-6;
npair = 8;
sopts = struct('solver','sedumi');
sopts.params.fid = 0;

% Each case is {label, vars, m, deg, rank of W, domain}
cases = {
    {'1D m=1 deg 0 rk1', {'s1'},      1, 0, 1, [0,1]}
    {'1D m=1 deg 0 rk3', {'s1'},      1, 0, 3, [0,1]}
    {'1D m=1 deg 1 rk2', {'s1'},      1, 1, 2, [0,1]}
    {'1D m=1 deg 1 rk5', {'s1'},      1, 1, 5, [0,1]}
    {'1D m=2 deg 0 rk3', {'s1'},      2, 0, 3, [0,1]}
    {'1D m=1 deg 1 shft',{'s1'},      1, 1, 3, [-1,2]}
    {'2D m=1 deg 0 rk2', {'s1','s2'}, 1, 0, 2, [0,1;0,1]}
    {'2D m=1 deg 0 rk5', {'s1','s2'}, 1, 0, 5, [0,1;0,1]}
    {'2D m=2 deg 0 rk4', {'s1','s2'}, 2, 0, 4, [0,1;0,1]}
    };

npass = 0;
for ic = 1:numel(cases)
    lbl  = cases{ic}{1};
    vars = cases{ic}{2};
    m    = cases{ic}{3};
    deg  = cases{ic}{4};
    kr   = cases{ic}{5};
    dom  = cases{ic}{6};
    N = numel(vars);

    % % % Declare the positive operator variable
    v1 = polynomial(vars(:));
    v2 = polynomial(strcat(vars(:),'_dum'));
    prog0 = lpiprogram(v1,v2,dom);
    [prog0,Q,Qcell,alpha_list] = possopvar(prog0,m,vars,dom,deg);
    dnames = cellstr(string(Q.Zd(:)));

    % % % Build a random target M = R*R with R = W*Z
    [names,nQ] = q_name_matrix(Qcell);
    W = randn(kr,nQ);
    Q0 = W.'*W;
    d0 = assign_from_matrix(names,Q0,dnames);
    M = const_sdopvar(Q,d0);

    % % % 1. Confirm that M really is R*R
    xf = rand_poly(vars,1,m);
    lhs = pair_sdopvar(M,xf,xf,vars,dom);
    rhs = norm_sq_Rx(W,alpha_list,deg,m,xf,vars,dom);
    e_RR = abs(lhs-rhs)/max(1,abs(rhs));
    if ~(e_RR<tol)
        error('test_lpi_eq_sdopvar_endtoend: <x,Mx> = %.10g but ||Rx||^2 = %.10g (case %s).',...
              lhs,rhs,lbl);
    end

    % % % 2. Constrain Q == M and solve
    D = Q + (-1)*M;
    prog = lpi_eq_sdopvar(prog0,D);
    evalc('sol = lpisolve(prog,sopts);');
    if ~solver_feasible(sol)
        error('test_lpi_eq_sdopvar_endtoend: the LPI should be feasible for case %s.',lbl);
    end

    % % % 3. The recovered operator must equal M
    dstar = recover_dvals(sol,dnames);
    e_op = op_distance(const_sdopvar(Q,dstar),M,vars,dom,m,npair);
    if ~(e_op<tol)
        error('test_lpi_eq_sdopvar_endtoend: recovered operator differs from M by %.3g (case %s).',...
              e_op,lbl);
    end

    % % % 4. The recovered Gram matrix must be positive semidefinite
    Qsol = matrix_from_assign(names,dstar,dnames);
    lmin = min(eig((Qsol+Qsol.')/2));
    if lmin < -1e-6*max(1,norm(Qsol))
        error('test_lpi_eq_sdopvar_endtoend: recovered Q is not psd (min eig %.3g, case %s).',...
              lmin,lbl);
    end

    % % % 5. The same solve using the symmetric option. Q and M are both
    % % % self-adjoint, so their difference is too and the option applies.
    prog_s = lpi_eq_sdopvar(prog0,D,'symmetric');
    evalc('sol_s = lpisolve(prog_s,sopts);');
    if ~solver_feasible(sol_s)
        error('test_lpi_eq_sdopvar_endtoend: the symmetric variant should be feasible (case %s).',lbl);
    end
    ds = recover_dvals(sol_s,dnames);
    e_sym = op_distance(const_sdopvar(Q,ds),M,vars,dom,m,npair);
    if ~(e_sym<tol)
        error(['test_lpi_eq_sdopvar_endtoend: the symmetric variant recovered an operator ',...
               'differing from M by %.3g (case %s).'],e_sym,lbl);
    end

    % % % 6. The infeasible direction must be detected
    prog_bad = lpi_eq_sdopvar(prog0,Q + M);
    evalc('sol_bad = lpisolve(prog_bad,sopts);');
    feas_bad = solver_feasible(sol_bad);

    fprintf(['  passed: %-18s nd=%5d nQ=%3d | M=R*R %.0e | op err %.0e | sym %.0e | ',...
             'min eig %+.1e | Q=-M %s\n'],...
            lbl,numel(dnames),nQ,e_RR,e_op,e_sym,lmin,ternary(feas_bad,'FEASIBLE(!)','infeasible'));
    if feas_bad
        error(['test_lpi_eq_sdopvar_endtoend: constraining Q = -M should be infeasible ',...
               'for case %s, but the solver reported it feasible.'],lbl);
    end
    npass = npass+1;
end

fprintf('\ntest_lpi_eq_sdopvar_endtoend passed (%d cases).\n',npass);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = op_distance(P,M,vars,dom,m,npair)
% Largest relative difference of <y,Px> and <y,Mx> over random test pairs.

e = 0;
for r = 1:npair
    xr = rand_poly(vars,1,m);
    yr = rand_poly(vars,1,m);
    a = pair_sdopvar(P,xr,yr,vars,dom);
    b = pair_sdopvar(M,xr,yr,vars,dom);
    e = max(e,abs(a-b)/max(1,abs(b)));
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = pair_sdopvar(P,xf,yf,vars,dom)
% <y,P x> in closed form from the stored kernels of an sdopvar whose decision
% variables have already been frozen.

N = numel(vars);
vars_dum = strcat(vars,'_dum');
m = P.dims(1);      n = P.dims(2);
NL = prod(cellfun(@numel,P.ZL));
NR = prod(cellfun(@numel,P.ZR));
ZLp = monom_vec(P.ZL,vars);
ZRp = monom_vec(P.ZR,vars_dum);

xt = xf;
for k = 1:N
    xt = subs(xt,polynomial(vars(k)),polynomial(vars_dum(k)));
end

Px = polynomial(zeros(m,1));
for k = 1:numel(P.params.A)
    Ck = reshape(full(P.params.A{k}),m*NL,n*NR);
    if nnz(Ck)==0
        continue
    end
    Kk = polynomial(zeros(m,n));
    for p = 1:m
        for r = 1:n
            Kk(p,r) = ZLp.'*(Ck((p-1)*NL+(1:NL),(r-1)*NR+(1:NR))*ZRp);
        end
    end
    gam = gamma_index(k,N);
    expr = Kk*xt;
    for i = 1:N
        sp = polynomial(vars_dum(i));
        sv = polynomial(vars(i));
        switch gam(i)
            case 1
                expr = subs(expr,sp,sv);
            case 2
                expr = int(expr,sp,dom(i,1),sv);
            case 3
                expr = int(expr,sp,sv,dom(i,2));
        end
    end
    Px = Px + expr;
end

expr = yf.'*Px;
for i = 1:N
    expr = int(expr,polynomial(vars(i)),dom(i,1),dom(i,2));
end
val = double(expr);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = norm_sq_Rx(W,alpha_list,deg,m,xf,vars,dom)
% ||R x||^2 for R = W*Z, with Z the stack of possopvar basis operators.

N = numel(vars);
Zx = [];
for i = 1:size(alpha_list,1)
    Zx = [Zx; basis_apply(alpha_list(i,:),deg,m,xf,vars,dom)];              %#ok<AGROW>
end
Rx = W*Zx;
expr = Rx.'*Rx;
for i = 1:N
    expr = int(expr,polynomial({[vars{i},'_int']}),dom(i,1),dom(i,2));
end
val = double(expr);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Zw = basis_apply(alpha,deg,m,w,vars,dom)
% (Z_alpha w)(theta) = int_s I_alpha(theta-s) (I_m o Z^alpha(theta,s)) w(s) ds,
% rebuilding the monomial vector that possopvar uses for this multi-index.

N = numel(vars);
caps_int = repmat(deg,1,N);
caps_mult = repmat(deg,1,N);
caps_mult(alpha==1) = 0;

E = exponent_grid([caps_int,caps_mult]);
T = size(E,1);
vnames = [strcat(vars(:),'_int'); vars(:)];
Zi = polynomial(speye(T),E,vnames,[T,1]);

Zw = polynomial(zeros(m*T,1));
for p = 1:m
    for u = 1:T
        expr = Zi(u)*w(p);
        for k = 1:N
            sv = polynomial(vars(k));
            th = polynomial({[vars{k},'_int']});
            switch alpha(k)
                case 1
                    expr = subs(expr,sv,th);
                case 2
                    expr = int(expr,sv,dom(k,1),th);
                case 3
                    expr = int(expr,sv,th,dom(k,2));
            end
        end
        Zw((p-1)*T+u) = expr;
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E = exponent_grid(caps)

n = numel(caps);
grids = cell(1,n);
for i = 1:n
    grids{i} = (0:caps(i))';
end
subsc = cell(1,n);
[subsc{:}] = ndgrid(grids{:});
E = zeros(numel(subsc{1}),n);
for i = 1:n
    E(:,i) = subsc{i}(:);
end
E = sortrows(E);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function idx = gamma_index(k,N)
% Multi-index of a parameter cell, first variable the fastest index.

idx = ones(1,N);
k = k-1;
for i = 1:N
    idx(i) = mod(k,3)+1;
    k = floor(k/3);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Z = monom_vec(Zc,vnames)
% kron(Z{1}(x_1),...,Z{N}(x_N)) as a polynomial column vector.

degmat = zeros(1,0);
for i = 1:numel(Zc)
    Zi = reshape(Zc{i},[],1);
    degmat = [kron(degmat,ones(numel(Zi),1)), kron(ones(size(degmat,1),1),Zi)];  %#ok<AGROW>
end
T = size(degmat,1);
Z = polynomial(speye(T),degmat,vnames(:),[T,1]);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [names,nQ] = q_name_matrix(Qcell)
% Assemble the full matrix of decision variable names of the Gram matrix.

rows = cell(size(Qcell,1),1);
for i = 1:size(Qcell,1)
    row = {};
    for j = 1:size(Qcell,2)
        row = [row, cellstr(string(Qcell{i,j}))];                           %#ok<AGROW>
    end
    rows{i} = row;
end
names = vertcat(rows{:});
nQ = size(names,1);
if size(names,2)~=nQ
    error('Qcell does not describe a square matrix.');
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = assign_from_matrix(names,Q0,dnames)
% Assign decision variable values realizing the symmetric matrix Q0.

[tf,loc] = ismember(names,dnames);
if ~all(tf(:))
    error('A name in Qcell is not among the operator decision variables.');
end
d = zeros(numel(dnames),1);
d(loc(:)) = Q0(:);
if norm(reshape(d(loc),size(Q0))-Q0,'fro') > 1e-12*max(1,norm(Q0,'fro'))
    error('The naming of Q is not symmetric; cannot realize a symmetric matrix.');
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Q0 = matrix_from_assign(names,d,dnames)

[~,loc] = ismember(names,dnames);
Q0 = reshape(d(loc),size(names));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = const_sdopvar(Q,d)
% Freeze an sdopvar at the given decision variable values, returning an
% sdopvar with no decision variables representing the same operator.

pA = cell(size(Q.params.A));
pB = cell(size(Q.params.B));
for k = 1:numel(pA)
    pA{k} = Q.params.A{k} + Q.params.B{k}.'*d;
    pB{k} = sparse(0,numel(pA{k}));
end
M = sdopvar(struct('A',{pA},'B',{pB}),Q.vars,{},Q.ZL,Q.ZR,Q.dom,Q.dims);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = recover_dvals(sol,dnames)
% Read the solved values of the named decision variables.

dv = dpvar(dnames);
d = double(lpigetsol(sol,dv));
d = reshape(d,[],1);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function feas = solver_feasible(sol)
% Read the primal feasibility flag the solver reports.

feas = true;
if isfield(sol,'solinfo') && isfield(sol.solinfo,'info') ...
        && isfield(sol.solinfo.info,'pinf')
    feas = ~(sol.solinfo.info.pinf==1);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = rand_poly(vars,deg,m)
% Random m x 1 polynomial vector of degree at most deg in each variable.

E = exponent_grid(repmat(deg,1,numel(vars)));
T = size(E,1);
w = polynomial(randn(T,m),E,vars(:),[m,1]);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = ternary(c,a,b)

if c, out = a; else, out = b; end

end
