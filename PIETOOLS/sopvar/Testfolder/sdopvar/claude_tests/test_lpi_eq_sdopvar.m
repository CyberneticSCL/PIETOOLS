%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_LPI_EQ_SDOPVAR verifies 'lpi_eq_sdopvar'.
%
% The constraints generated for an operator P must cut out exactly the set of
% decision variable values for which P is the zero operator, no more and no
% less. That reference set is computed here without reusing any of the
% machinery in 'lpi_eq_sdopvar': the pairing <y,P x> is evaluated in closed
% form from the stored kernels for a batch of random polynomial test
% functions, once per unit decision vector. Since the pairing is affine in
% the decision variables, those evaluations assemble a matrix PSI with
%
%       <y_r, P(d) x_r> = PSI(r,:) * [1;d],
%
% and the operator is zero exactly when PSI*[1;d] = 0, for enough random test
% pairs. The test then checks that the affine set defined by the generated
% constraints coincides with that one, by comparing row spaces.
%
% The decisive case is P-P' for a self-adjoint P. That is the zero operator,
% so it must impose no restriction at all.
%
% MMP, 08/29/2026: This case used to be decisive because 'ctranspose' moved
% monomial degrees from ZL to ZR, which delta(s-s') makes equivalent, so
% P-P' had a multiplier parameter outside the canonical form and
% constraining the stored coefficients would wrongly have restricted the
% decision variables. 'lpi_eq_sdopvar' contracted the multiplier directions
% to avoid that. The canonical form is now a class invariant enforced by the
% constructor, and 'ctranspose' preserves it, so a non-canonical parameter
% can no longer be built here and the cases below marked 'non-canonical' are
% canonicalized on construction. The check that the generated constraints
% cut out exactly the operator-zero set is unaffected and still the point of
% this file; the invariant itself is verified in
% 'test_canonical_multiplier'.
%
% MMP, 08/28/2026: Initial coding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
rng(2);

vars = {'s1'};
dom = [0,1];
npair = 40;

fprintf('=== generated constraints vs the true operator-zero set ===\n');
npass = 0;
cases = {
    {'canonical multiplier',      true,  {[0;1]},   {[0;1]},   3}
    {'non-canonical multiplier',  false, {[0;1]},   {[0;1]},   3}
    {'non-canonical, deg 2',      false, {[0;1;2]}, {[0;1;2]}, 4}
    {'canonical, deg 2',          true,  {[0;1;2]}, {[0;1;2]}, 4}
    {'asymmetric bases',          false, {[0;1;2]}, {[0;1]},   3}
    };

for ic = 1:numel(cases)
    lbl = cases{ic}{1};
    canonical = cases{ic}{2};
    ZL = cases{ic}{3};
    ZR = cases{ic}{4};
    nd = cases{ic}{5};

    dvarname = arrayfun(@(k) sprintf('r%d',k),(1:nd)','UniformOutput',false);
    prog = lpiprogram(polynomial(vars(:)),polynomial(strcat(vars(:),'_dum')),dom);
    prog = lpidecvar(prog,dvarname);

    P = make_random_sdopvar(ZL,ZR,dvarname,vars,dom,canonical);

    % % % Constraints under test
    prog_new = lpi_eq_sdopvar(prog,P);
    [A1,b1] = collect_constraints(prog_new,prog);

    % % % Reference: the true operator-zero set
    PSI = build_psi(P,vars,dom,npair);
    [A2,b2] = psi_system(PSI,P,prog);

    r1 = rank_aug(A1,b1);
    r2 = rank_aug(A2,b2);
    rj = rank_aug([A1;A2],[b1;b2]);
    if ~(r1==r2 && rj==r1)
        error(['test_lpi_eq_sdopvar: constraints do not match the operator-zero set ',...
               'for case %s: rank(new)=%d, rank(true)=%d, rank(joint)=%d.'],lbl,r1,r2,rj);
    end

    fprintf('  passed: %-26s %3d rows, rank %2d == true rank %2d\n',...
            lbl,size(A1,1),r1,r2);
    npass = npass+1;
end

% % % The decisive case: P-P' for a self-adjoint P is the zero operator, so
% % % it must impose no restriction, even though its multiplier parameter is
% % % not in the canonical form.
fprintf('\n=== P-P'''' for a self-adjoint P must impose no restriction ===\n');
for dg = [0,1]
    prog = lpiprogram(polynomial({'s1'}),polynomial({'s1_dum'}),[0,1]);
    [prog,P] = possopvar(prog,1,{'s1'},[0,1],dg);
    D = P + (-1)*P';   % @sdopvar defines no 'minus'

    prog_c = lpi_eq_sdopvar(prog,D);
    [Ac,bc] = collect_constraints(prog_c,prog);
    rc = rank_aug(Ac,bc);

    % How restrictive would it be without the multiplier contraction?
    [An,bn] = naive_constraints(D);
    rn = rank_aug(An,bn);

    if rc~=0
        error(['test_lpi_eq_sdopvar: P-P'''' is the zero operator but the generated ',...
               'constraints have rank %d (deg %d).'],rc,dg);
    end
    fprintf(['  passed: deg %d   contracted rank %d (no restriction), ',...
             'naive C==0 would have rank %d over %d decision variables\n'],...
            dg,rc,rn,numel(P.Zd));
    if rn==0
        fprintf('          note: naive would also be vacuous here, so this case is not discriminating\n');
    end
    npass = npass+1;
end

% % % Constraining a possopvar operator to zero must force it to zero
fprintf('\n=== possopvar operator constrained to zero ===\n');
prog = lpiprogram(polynomial({'s1'}),polynomial({'s1_dum'}),[0,1]);
[prog,P] = possopvar(prog,1,{'s1'},[0,1],0);
prog_z = lpi_eq_sdopvar(prog,P);
[Az,bz] = collect_constraints(prog_z,prog);
PSI = build_psi(P,{'s1'},[0,1],npair);
rz = rank_aug(Az,bz);
[Ap,bp] = psi_system(PSI,P,prog);
rt = rank_aug(Ap,bp);
rjz = rank_aug([Az;Ap],[bz;bp]);
if ~(rz==rt && rjz==rz)
    error(['test_lpi_eq_sdopvar: constraints on the possopvar operator do not match ',...
           'the operator-zero set: %d vs %d (joint %d).'],rz,rt,rjz);
end
fprintf('  passed: %d decision variables, rank %d == true rank %d\n',numel(P.Zd),rz,rt);
npass = npass+1;

% % % Input validation
fprintf('\n=== input validation ===\n');
probes = { {'polynomial input', @() lpi_eq_sdopvar(prog,polynomial(1))}
           {'double input',     @() lpi_eq_sdopvar(prog,1)}
           {'ndopvar input',    @() lpi_eq_sdopvar(prog,ndopvar())}
           {'bad option',       @() lpi_eq_sdopvar(prog,P,'skew')} };
for k = 1:numel(probes)
    raised = false;
    try
        probes{k}{2}();
    catch ME
        raised = true;
        fprintf('  errors   : %-18s %s\n',probes{k}{1},strtrim(ME.message));
    end
    if ~raised
        error('test_lpi_eq_sdopvar: expected an error for %s.',probes{k}{1});
    end
end

fprintf('\nlpi_eq_sdopvar test passed (%d checks).\n',npass);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = make_random_sdopvar(ZL,ZR,dvarname,vars,dom,canonical)
% Build a scalar sdopvar on one spatial variable with random parameters. If
% canonical is true, the multiplier parameter carries no ZR degree.

nd = numel(dvarname);
NL = prod(cellfun(@numel,ZL));
NR = prod(cellfun(@numel,ZR));
nC = NL*NR;

A = cell(3,1);      B = cell(3,1);
for k = 1:3
    Ak = randn(NL,NR);
    Bk = randn(nd,nC);
    if k==1 && canonical
        % Keep only the first ZR degree, i.e. no dependence on the dummy
        % variable in the multiplier direction.
        mask = false(NL,NR);    mask(:,1) = true;
        Ak = Ak.*mask;
        Bk = Bk.*repmat(reshape(mask,1,[]),nd,1);
    end
    A{k} = sparse(Ak(:));
    B{k} = sparse(Bk);
end

vv = struct();      vv.in = vars;       vv.out = vars;
dd = struct();      dd.in = dom;        dd.out = dom;
P = sdopvar(struct('A',{A},'B',{B}),vv,dvarname,ZL,ZR,dd,[1,1]);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PSI = build_psi(P,vars,dom,npair)
% Assemble PSI with <y_r,P(d) x_r> = PSI(r,:)*[1;d], by evaluating the
% pairing once per unit decision vector.

nd = numel(P.Zd);
PSI = zeros(npair,nd+1);
X = cell(1,npair);      Y = cell(1,npair);
for r = 1:npair
    X{r} = rand_poly(vars,2);
    Y{r} = rand_poly(vars,2);
end

base = zeros(npair,1);
for r = 1:npair
    base(r) = pair_val(P,zeros(nd,1),X{r},Y{r},vars,dom);
end
PSI(:,1) = base;
for k = 1:nd
    ek = zeros(nd,1);   ek(k) = 1;
    for r = 1:npair
        PSI(r,k+1) = pair_val(P,ek,X{r},Y{r},vars,dom) - base(r);
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = pair_val(P,d,xf,yf,vars,dom)
% Evaluate <y,P(d) x> in closed form from the stored kernels.

s = polynomial(vars(1));
t = polynomial({[vars{1},'_dum']});
NL = prod(cellfun(@numel,P.ZL));
NR = prod(cellfun(@numel,P.ZR));

ZLp = monom_vec(P.ZL{1},vars{1});
ZRp = monom_vec(P.ZR{1},[vars{1},'_dum']);
xt = subs(xf,s,t);

Px = polynomial(0);
for k = 1:3
    Ck = reshape(full(P.params.A{k} + P.params.B{k}.'*d),NL,NR);
    Kk = ZLp.'*(Ck*ZRp);
    switch k
        case 1      % multiplier: delta(s-t)
            Px = Px + subs(Kk*xt,t,s);
        case 2      % lower integral: t <= s
            Px = Px + int(Kk*xt,t,dom(1,1),s);
        case 3      % upper integral: t >= s
            Px = Px + int(Kk*xt,t,s,dom(1,2));
    end
end
val = double(int(yf*Px,s,dom(1,1),dom(1,2)));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Z = monom_vec(degs,vname)

degs = degs(:);
T = numel(degs);
Z = polynomial(speye(T),degs,{vname},[T,1]);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = rand_poly(vars,deg)

E = (0:deg)';
w = polynomial(randn(deg+1,1),E,vars(1),[1,1]);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,b] = naive_constraints(P)
% The constraint system that setting every stored coefficient to zero would
% give, for comparison. Not used to constrain anything.

A = [];     b = [];
for k = 1:numel(P.params.A)
    Ak = full(P.params.A{k});
    Bk = full(P.params.B{k});
    A = [A; -Bk.'];                                                         %#ok<AGROW>
    b = [b; Ak(:)];                                                         %#ok<AGROW>
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,b] = collect_constraints(prog,prog0)

A = [];     b = [];
for k = prog0.expr.num+1:prog.expr.num
    A = [A; full(prog.expr.At{k}).'];                                       %#ok<AGROW>
    b = [b; full(prog.expr.b{k})];                                          %#ok<AGROW>
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = rank_aug(A,b)

M = [A,b];
if isempty(M)
    r = 0;
    return
end
nrm = vecnorm(M,2,2);
nrm(nrm<eps) = 1;
M = M./nrm;
M(abs(M)<1e-12) = 0;
r = rank(M,1e-8*max(size(M)));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,b] = psi_system(PSI,P,prog)
% Express the operator-zero conditions over the full list of program decision
% variables, so that they can be compared against the constraint rows that
% SOSTOOLS stores, which are indexed by prog.decvartable.

names = cellstr(string(P.Zd(:)));
known = cellstr(string(prog.decvartable(:)));
[tf,loc] = ismember(names,known);
if ~all(tf)
    error('Operator refers to a decision variable not in the program.');
end
A = zeros(size(PSI,1),numel(known));
A(:,loc) = PSI(:,2:end);
b = -PSI(:,1);

end
