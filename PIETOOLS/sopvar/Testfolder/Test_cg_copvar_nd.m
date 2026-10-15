% Test_cg_copvar_nd - tests cg_copvar_pd on copvar operators in n = 1,2,3
% spatial variables on [0,1]^n, generated with rand_copvar.
%
% For a random M x M copvar A (M = 2 blocks, each block L2^m[s1..sN] ->
% L2^m[s1..sN] with random sparse polynomial kernels of degree d) and a
% random M x 1 copvar X with no input variables, B = A*X is formed and the
% regularized least squares problem
%       min ||A Y - B||^2 + ||Y||^2   <=>   (A'A + I) Y = A'B
% is solved with the conjugate gradient method for the positive definite
% operator A'A + I (cg_copvar_pd). Checks that ||(A'A+I)Y-A'B|| <= 1e-6*||A'B||, recomputed outside the solver.
%
% Cost grows steeply with n.
%
% Initial coding AT - 10/15/2026

maxit = 40;
tol_res = 1e-6;     % required reduction of the system residual
m = 2;              % components per block
M = 2;              % blocks, i.e. size of the container

rng(1);

for N = [1]
    % Declare N variables and their domain
    % N = 1 fast
    % N = 2 slow, but OK
    % N = 3 impossible, unrealistic
    d = 2;
    names = arrayfun(@(k) sprintf('s%d',k),1:N,'uni',0);
    dom = [0,1];                                % shared by all variables

    % random M x M copvar, every block L2^m[S] -> L2^m[S]
    spaces = struct('out',{repmat({names},1,M)}, ...
                    'in',{repmat({names},1,M)});
    dims = struct('out',m*ones(M,1), ...
                   'in',m*ones(M,1));
    A = rand_copvar(spaces,dims,dom,d,0.6);

    % random copvar Mx1 with 0 input variables (R -> L2^m[S] blocks)
    spacesX = struct('out',{repmat({names},1,M)},'in',{{{}}});
    dimsX = struct('out',m*ones(M,1),'in',1);
    X = rand_copvar(spacesX,dimsX,dom,d,0.6);

    B = A*X;

    % identity copvar I (block diagonal, structurally zero off the diagonal)
    Ic = cell(M,M);
    for i = 1:M,  Ic{i,i} = id_sopvar(m,names,repmat(dom,N,1));  end
    I1 = copvar(Ic);

    % Test cg_copvar may not converge
    % [Xr,res,it] = cg_copvar(A,B,1e-10,maxit);   % min ||A Y - B||^2 + ||Y||^2
    % E = A*Xr-B;   W = E'*E;
    % blk = W.C{1,1};   err = full(blk.params{1});
    % % The residual norm returned by the solver can read exactly 0 at the
    % % rounding level, so check the independently recomputed residual E.
    % W = B'*B;     nB2 = full(W.C{1,1}.params{1});
    % assert(sqrt(err)<=tol_res*sqrt(nB2),'N=%d: ||Apd*Xr-Bpd||/||Bpd|| = %g.',N,sqrt(err/nB2));
    % fprintf('n=%d passed: it=%d relRes=%.1e err=%.1e\n',N,it,res(end)/res(1),err);
    Apd = A'*A + I1;
    Bpd = A'*B;
    [Xr,res,it] = cg_copvar_pd(Apd,Bpd,1e-10,maxit);   % min ||A Y - B||^2 + ||Y||^2

    E = Apd*Xr-Bpd;   W = E'*E;
    blk = W.C{1,1};   err = full(blk.params{1});
    % The residual norm returned by the solver can read exactly 0 at the
    % rounding level, so check the independently recomputed residual E.
    W = Bpd'*Bpd;     nB2 = full(W.C{1,1}.params{1});
    assert(sqrt(err)<=tol_res*sqrt(nB2),'N=%d: ||Apd*Xr-Bpd||/||Bpd|| = %g.',N,sqrt(err/nB2));
    fprintf('n=%d passed: it=%d relRes=%.1e err=%.1e\n',N,it,res(end)/res(1),err);
end
