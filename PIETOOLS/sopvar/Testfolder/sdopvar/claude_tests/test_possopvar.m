%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_POSSOPVAR verifies the operator returned by 'possopvar'.
%
% For each test case, the script checks the defining identity of Sec. 9,
%
%   <y, Pop x> = sum_{i,j} int_theta g(theta) (Z_{alpha_i} y)(theta)' ...
%                                      Q_{ij} (Z_{alpha_j} x)(theta) dtheta
%
% by evaluating both sides exactly, in closed form, for random polynomial
% test functions x, y and a random assignment of the decision variables:
%
%   * the left-hand side is computed from the kernel form actually stored in
%     the returned sdopvar object, i.e. from params.A, params.B, ZL and ZR;
%   * the right-hand side is computed directly from the definition of the
%     basis operators Z_alpha and the matrix Q named by Qcell.
%
% Since the identity is linear in Q, the decision variables are assigned
% random (not necessarily positive semidefinite) values; this tests the
% construction more sharply than a feasible assignment would.
%
% The script also checks that the returned operator is self-adjoint. Note
% that this is tested at the level of the operator, i.e. by verifying
% <y,Pop x> = <x,Pop y>, rather than by comparing Pop with Pop' using 'eq'.
% The reason is that the representation of a multiplier kernel is not unique:
% the factor delta(s_k-s'_k) identifies s_k with s'_k, so degrees may be
% moved freely between ZL and ZR in that coordinate. 'possopvar' returns the
% canonical form, in which a multiplier cell carries no ZR degree in its
% delta coordinates, but @sdopvar/ctranspose transposes the coefficient
% matrix without restoring that normalization, so eq(Pop,Pop') returns false
% for any operator with a nonzero multiplier part. The integral cells are
% compared against the transpose directly, since those are unambiguous.
%
% MP, 08/22/2026: Initial coding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
rng(1);
tol = 1e-8;

% Each case is {label, m, vars, dom, deg, options}
cases = {
    {'1D, m=1, deg 1',                  1, {'s1'},      [0,1],       1, struct()}
    {'1D, m=1, deg 2',                  1, {'s1'},      [0,1],       2, struct()}
    {'1D, m=2, deg 1',                  2, {'s1'},      [0,1],       1, struct()}
    {'1D, m=1, deg 1, psatz',           1, {'s1'},      [0,1],       1, struct('psatz',1)}
    {'1D, m=1, shifted domain',         1, {'s1'},      [-1,2],      1, struct()}
    {'1D, m=1, int/mult degrees',       1, {'s1'},      [0,1],       struct('int',2,'mult',1), struct()}
    {'1D, m=1, multiplier only',        1, {'s1'},      [0,1],       1, struct('include',1)}
    {'1D, m=1, integrals only',         1, {'s1'},      [0,1],       1, struct('include',[2;3])}
    {'1D, m=2, lower integral only',    2, {'s1'},      [0,1],       1, struct('include',2)}
    {'2D, m=1, deg 1',                  1, {'s1','s2'}, [0,1;0,1],   1, struct()}
    {'2D, m=1, deg 1, psatz',           1, {'s1','s2'}, [0,1;0,1],   1, struct('psatz',1)}
    {'2D, m=1, mixed domains',          1, {'s1','s2'}, [0,1;-1,3],  1, struct()}
    {'2D, m=2, joint degree cap',       2, {'s1','s2'}, [0,1;0,1],   struct('int',1,'mult',1,'joint',2), struct()}
    {'2D, m=1, mult and lower only',    1, {'s1','s2'}, [0,1;0,1],   1, struct('include',[1 1;2 2])}
    };

npass = 0;
for ic = 1:numel(cases)
    [lbl,m,vars,dom,deg,options] = deal(cases{ic}{:});
    n3 = numel(vars);
    if size(dom,1)==1
        dom = repmat(dom,n3,1);
    end

    % % % Declare the operator
    prog = lpiprogram(polynomial(vars(:)),[],dom);
    [prog,Pop,Qcell,alpha_list] = possopvar(prog,m,vars,dom,deg,options);

    % % % Assign random values to the decision variables
    dvars = cellstr(string(Pop.Zd(:)));
    dval = randn(numel(dvars),1);

    % % % Random polynomial test functions
    xfun = rand_polyvec(vars,m);
    yfun = rand_polyvec(vars,m);

    % % % Check the defining identity
    lhs = pair_kernel(Pop,xfun,yfun,m,vars,dom,dval);
    rhs = pair_factored(Qcell,alpha_list,xfun,yfun,m,vars,dom,deg,options,dval,dvars);
    err = abs(lhs-rhs)/max(1,abs(rhs));
    if ~(err<tol)
        error('test_possopvar: identity failed for case ''%s'': lhs=%.12g, rhs=%.12g, rel err=%.3g.',...
              lbl,lhs,rhs,err);
    end

    % % % Check self-adjointness of the operator
    swp = pair_kernel(Pop,yfun,xfun,m,vars,dom,dval);
    err_adj = abs(lhs-swp)/max(1,abs(lhs));
    if ~(err_adj<tol)
        error('test_possopvar: operator is not self-adjoint for case ''%s'': <y,Px>=%.12g, <x,Py>=%.12g.',...
              lbl,lhs,swp);
    end

    % % % Check that the transpose reproduces every non-multiplier cell
    err_ct = compare_integral_cells(Pop,Pop');
    if ~(err_ct<tol)
        error('test_possopvar: transpose differs on an integral cell for case ''%s'' (max err %.3g).',...
              lbl,err_ct);
    end

    % % % Check that a positive semidefinite Q gives a positive operator.
    % This is not implied by the identity check, which uses the same
    % convention on both sides; it would however fail if the blocks Q_{ij}
    % were paired with the wrong basis operators.
    dpos = assign_psd(Qcell,dvars);
    quad = pair_kernel(Pop,xfun,xfun,m,vars,dom,dpos);
    if quad < -tol*max(1,abs(quad))
        error('test_possopvar: <x,Px> = %.12g < 0 for a positive semidefinite Q, case ''%s''.',...
              quad,lbl);
    end

    fprintf('  passed: %-34s (n_dec = %5d, identity %.1e, adjoint %.1e, <x,Px> = %.3e)\n',...
            lbl,numel(dvars),err,err_adj,quad);
    npass = npass+1;
end

fprintf('possopvar test passed (%d of %d cases).\n',npass,numel(cases));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = pair_kernel(Pop,xfun,yfun,m,vars,dom,dval)                   %#ok<INUSD>
% Evaluate <y,Pop x> from the kernel form stored in the sdopvar object.

val = pi_kernel_pairing(pi_sdopvar_kernels(Pop,dval),xfun,yfun,vars,dom);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = pair_factored(Qcell,alpha_list,xfun,yfun,m,vars,dom,deg,options,dval,dvars)
% Evaluate the factored form directly from the definition,
%   sum_{i,j} int_theta g(theta) (Z_alpha_i y)(theta)' Q_ij (Z_alpha_j x)(theta) dtheta.

n3 = numel(vars);
vars_int = strcat(vars,'_int');
nblk = size(alpha_list,1);

Zy = cell(1,nblk);
Zx = cell(1,nblk);
for i = 1:nblk
    Zi = basis_monomials(alpha_list(i,:),deg,i,vars_int,vars);
    Zy{i} = apply_basis_op(Zi,yfun,alpha_list(i,:),vars,vars_int,dom,m);
    Zx{i} = apply_basis_op(Zi,xfun,alpha_list(i,:),vars,vars_int,dom,m);
end

if isfield(options,'psatz') && ~isempty(options.psatz) && options.psatz
    gfun = polynomial(1);
    for k = 1:n3
        thk = polynomial(vars_int(k));
        gfun = gfun*(thk-dom(k,1))*(dom(k,2)-thk);
    end
else
    gfun = polynomial(1);
end

expr = polynomial(0);
for i = 1:nblk
    for j = 1:nblk
        Qij = q_block(Qcell{i,j},dval,dvars);
        expr = expr + Zy{i}.'*Qij*Zx{j};
    end
end
val = full_integral(gfun*expr,vars_int,dom);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err = compare_integral_cells(P,Pt)
% Compare P and its transpose on every cell that has no delta factor. Those
% cells have a unique representation, so a self-adjoint operator must
% reproduce them exactly. The remaining (multiplier) cells are skipped, see
% the note at the top of this file.

if ~isequal(cellstr(string(P.Zd(:))),cellstr(string(Pt.Zd(:)))) ...
        || ~isequal(P.ZL,Pt.ZL) || ~isequal(P.ZR,Pt.ZR)
    error('test_possopvar: transpose returned a different basis.');
end

n3 = numel(P.vars.in);
err = 0;
for k = 1:numel(P.params.A)
    if any(pi_gamma_index(k,n3)==1)
        continue
    end
    dA = P.params.A{k}-Pt.params.A{k};
    dB = P.params.B{k}-Pt.params.B{k};
    err = max([err;abs(dA(:));abs(dB(:))]);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Zi = basis_monomials(alpha,deg,iblk,vars_int,vars_out)
% Rebuild the monomial vector Z^alpha(theta,s) used by 'possopvar'.

n3 = numel(vars_out);
if iscell(deg)
    deg = deg{iblk};
end
if isnumeric(deg) && isscalar(deg)
    deg = struct('int',deg);
end
if ~isfield(deg,'int') || isempty(deg.int)
    deg.int = 1;
end
if ~isfield(deg,'mult') || isempty(deg.mult)
    deg.mult = deg.int;
end
caps_int = expand(deg.int,n3);
caps_mult = expand(deg.mult,n3);
caps_mult(alpha==1) = 0;
if ~isfield(deg,'joint') || isempty(deg.joint)
    jointcap = sum(caps_int)+sum(caps_mult);
else
    jointcap = deg.joint;
end

E = grid_exponents([caps_int,caps_mult],jointcap);
T = size(E,1);
Zi = polynomial(speye(T),E,[vars_int(:);vars_out(:)],[T,1]);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Zw = apply_basis_op(Zi,w,alpha,vars,vars_int,dom,m)
% Evaluate (Z_alpha w)(theta) = int_s I_alpha(theta-s) (I_m kron Z^alpha(theta,s)) w(s) ds
% in closed form. The result is an (m*T) x 1 polynomial in theta.

T = size(Zi,1);
Zw = polynomial(zeros(m*T,1));
for p = 1:m
    for u = 1:T
        expr = Zi(u)*w(p);
        for k = 1:numel(vars)
            sk = polynomial(vars(k));
            thk = polynomial(vars_int(k));
            switch alpha(k)
                case 1      % delta(theta_k-s_k)
                    expr = subs(expr,sk,thk);
                case 2      % I(theta_k-s_k), i.e. s_k <= theta_k
                    expr = int(expr,sk,dom(k,1),thk);
                case 3      % I(s_k-theta_k), i.e. s_k >= theta_k
                    expr = int(expr,sk,thk,dom(k,2));
            end
        end
        Zw((p-1)*T+u) = expr;
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = full_integral(expr,vartab,dom)
% Integrate a scalar polynomial over the full domain and return a double.

for k = 1:numel(vartab)
    expr = int(expr,polynomial(vartab(k)),dom(k,1),dom(k,2));
end
val = double(expr);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dval = assign_psd(Qcell,dvars)
% Assign values to the decision variables so that the matrix Q named by
% Qcell is positive semidefinite.

% Assemble the full matrix of decision variable names.
names = cell(size(Qcell,1),1);
for i = 1:size(Qcell,1)
    row = {};
    for j = 1:size(Qcell,2)
        row = [row, cellstr(string(Qcell{i,j}))];                            %#ok<AGROW>
    end
    names{i} = row;
end
names = vertcat(names{:});

nQ = size(names,1);
if size(names,2)~=nQ
    error('test_possopvar: Qcell does not describe a square matrix.');
end

M = randn(nQ,nQ);
S = (M'*M)/nQ;

[tf,loc] = ismember(names,dvars);
if ~all(tf(:))
    error('test_possopvar: decision variable in Qcell not found in Pop.Zd.');
end
dval = zeros(numel(dvars),1);
dval(loc(:)) = S(:);

% A name occurring at both (p,q) and (q,p) is written twice, so confirm that
% the assignment really does reproduce the intended symmetric matrix.
if norm(reshape(dval(loc),nQ,nQ)-S,'fro') > 1e-12*max(1,norm(S,'fro'))
    error('test_possopvar: the naming of Q is not symmetric, cannot assign a psd value.');
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Q = q_block(names,dval,dvars)
% Look up the numeric values of the decision variables naming a block of Q.

names = cellstr(string(names));
[tf,loc] = ismember(names,dvars);
if ~all(tf(:))
    error('test_possopvar: decision variable in Qcell not found in Pop.Zd.');
end
Q = reshape(dval(loc),size(names));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = rand_polyvec(vars,m)
% Random m x 1 polynomial vector of degree at most 1 in each variable.

n3 = numel(vars);
E = grid_exponents(ones(1,n3),[]);
T = size(E,1);
if n3==0
    w = polynomial(randn(m,1));
else
    w = polynomial(randn(T,m),E,vars(:),[m,1]);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E = grid_exponents(caps,jointcap)
% Local copy of 'build_exponent_grid', which lives in a private folder.

caps = reshape(caps,1,[]);
n = numel(caps);
if isempty(jointcap)
    jointcap = sum(caps);
end
if n==0
    E = zeros(1,0);
    return
end
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
E = E(sum(E,2)<=jointcap,:);
E = sortrows(E);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function caps = expand(caps,n3)

caps = reshape(caps,1,[]);
if isscalar(caps)
    caps = repmat(caps,1,n3);
end

end
