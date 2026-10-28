% Test_qp_copvar_nd - tests qp_copvar on copvar operators in n = 1,2 spatial
% variables on [0,1]^n, generated with rand_copvar, in the way
% Test_cg_copvar_nd tests cg_copvar_pd.
%
% For a random M x M copvar A0 (M = 2 blocks, each block L2^m[s1..sN] ->
% L2^m[s1..sN], m = 2 components, with random sparse polynomial kernels of
% degree d = 2) and a random M x 1 copvar X0 with no input variables (blocks
% R -> L2^m[s1..sN]), the least squares problem
%       min_X ||A X - B||^2     over X of degree <= deg in every variable,
% deg the degree of X0, is a QP in the coefficients of X (qp_copvar), solved
% in the Legendre basis with 'backslash'. The operator A and the right-hand
% side B are set by the variable caseid:
%   caseid = 1: A = A0'A0 + I, strictly positive definite, B = X0 is a random operator;
%   caseid = 2: A = A0, indefinite, B = A*X0, in the range of A;
%   caseid = 3: A = A0, indefinite, B = X0, in general not in the range of A.
% In cases 1 and 2 the minimum is 0 and is attained at X = X0 (case 2: X0 is
% the unique minimizer when A is injective on the polynomials of degree <=
% deg, as it is for the random A0 used), in case 3 the objective is not 0.
% A warning, not an error (asserts before 10/28/2026), is issued if
% ||A Xr - B||^2 > tol_res*||B||^2 (the objective is not 0) or
% ||Xr - X0||^2 > tol_err*||X0||^2, both recomputed outside the solver; for
% caseid = 3 both warnings are expected. For each n one line is printed: the
% assembly, deg, the number of unknowns, the time of qp_copvar, cond(H) and
% the relative objective and error.
%
% Assembly: 'container' for n = 1 and 'apply' for n >= 2, where the container
% products need tens of GB (measured for case 1, 2-D, 2 blocks: container
% 41.6 s and 13.6 GB at deg 2, and 96 s and 30.6 GB at deg 8; 'apply' about 16 s
% and 1.1 GB at deg 2). To use the container for n = 2, remove the switch to
% 'apply' in the line setting assembly.
%
% Cost grows steeply with n (n = 2: about 10 to 17 s at deg 2 with 'apply').
%
% Initial coding AT - 10/27/2026
% AT, 10/28/2026: eye_copvar_sop for the identity, warnings in place of the
%                asserts, caseid 1 to 3; header describes them.

tol_err = 1e-8;     % tolerance on ||X-X0||^2/||X0||^2
tol_res = 1e-8;     % tolerance on ||A X - B||^2/||B||^2
m = 2;              % components per block
M = 2;              % blocks, i.e. size of the container

rng(1);


% caseid = 1; % A is strictly PD, B random operator
% caseid = 3; % A is indefinite, B = X0: not in the range of A, the objective is not 0
caseid = 3; % A is indefinite, B = A*X0 is in the range of A

for N = [1 2]
    % Declare N variables and their domain
    % N = 1 fast
    % N = 2 slow, but OK
    d = 2;
    names = arrayfun(@(k) sprintf('s%d',k),1:N,'uni',0);
    dom = [0,1];                                % shared by all variables

    % random M x M copvar, every block L2^m[S] -> L2^m[S]
    spaces = struct('out',{repmat({names},1,M)}, ...
                    'in',{repmat({names},1,M)});
    dims = struct('out',m*ones(M,1), ...
                   'in',m*ones(M,1));
    A0 = rand_copvar(spaces,dims,dom,d,0.6);

    % random copvar Mx1 with 0 input variables (R -> L2^m[S] blocks), the solution
    spacesX = struct('out',{repmat({names},1,M)},'in',{{{}}});
    dimsX = struct('out',m*ones(M,1),'in',1);
    X0 = rand_copvar(spacesX,dimsX,dom,d,0.6);
    deg = 0;                                    % X0 has degree <= deg in every variable
    for i = 1:M,  deg = max([deg,cellfun(@max,X0.C{i,1}.ZL(:)')]);  end

    % identity copvar I (block diagonal, structurally zero off the diagonal) 
    I1 = eye_copvar_sop(dims,spaces,dom);                                   

    switch caseid
        case 1 % the objective is 0
            A = A0'*A0 + I1; % strictly positive definite
            B = X0;
        case 2 % the objective is 0
            A = A0;          % indefinite
            B = A*X0;       % B is in range of A
        case 3 % the objective is not 0
            A = A0;
            B = X0;
    end
    assembly = 'container';   %if N>1,  assembly = 'apply';  end
    tic;
    [Xr,info] = qp_copvar(A,B, ...
        struct('deg',6,'basis','legendre','assembly',assembly,'solver','backslash'));   % min ||A X - B||^2
    t = toc;

    E = A*Xr-B;     W = E'*E;     err_res = full(W.C{1,1}.params{1});
    D = Xr-X0;      W = D'*D;     err_X   = full(W.C{1,1}.params{1});
    W = B'*B;       nB2 = full(W.C{1,1}.params{1});
    W = X0'*X0;     nX2 = full(W.C{1,1}.params{1});
%
    % a warning if the objective is not 0 for caseid = 3 by design   % AT, 10/28/2026
    if err_res>tol_res*nB2                                                 
        warning('Test_qp_copvar_nd:objective','N=%d: the objective ||A*Xr-B||^2/||B||^2 = %g is not 0.',N,err_res/nB2);   % AT, 10/28/2026
    end                                                                     
    
    if caseid == 2
        if err_X>tol_err*nX2                                                   
            warning('Test_qp_copvar_nd:solution','N=%d: ||Xr-X0||^2/||X0||^2 = %g, Xr is not X0.',N,err_X/nX2);   % AT, 10/28/2026
        end                                                                     
    end
    fprintf('n=%d done: %s assembly, deg %d, %d unknowns, %.1f s, cond(H)=%.1e, res=%.1e, err=%.1e\n', ...
        N,assembly,deg,info.n,t,info.cond,err_res/nB2,err_X/nX2);
end
