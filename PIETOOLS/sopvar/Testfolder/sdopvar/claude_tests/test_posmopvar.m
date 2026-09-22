%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_POSMOPVAR checks that the container 'posmopvar' returns really is the
% quadratic form it claims to be. It is the container counterpart of
% 'test_possopvar', which does the same for the single-space case.
%
% The defining identity, Sec. 8.4:
%
%   <y,Pop x> = sum_{c,e} <Z_{alpha_c} y_{k(c)}, Q_{ce} Z_{alpha_e} x_{k(e)}>
%
% where the inner product on the right is over the auxiliary variable theta,
% and Z_alpha maps space k into L_2[theta] by
%
%   (Z_alpha w)(theta) = int_{s in dom(s^k)} I_alpha(theta-s) ...
%                           (I_m kron Z^alpha(theta,s)) w(s) ds.
%
% The two sides are computed by disjoint routes:
%
%   LEFT  from the returned operator's stored kernels, applied through the
%         sdopvar definition -- shared variables carry an indicator, input-only
%         variables are integrated over their whole domain, output-only ones
%         are free multipliers;
%   RIGHT from Qcell and basis_list ONLY, by integrating each basis operator
%         against the test function in closed form and pairing over theta.
%
% Nothing in the right-hand route reads the operator, so the two agree only if
% every Q block is paired with the correct basis operators and each block of
% the container assembles the correct subset of the sum.
%
% WHY THIS IS NEEDED ON TOP OF THE SPAN TESTS. Both
% 'test_posmopvar_vs_poslpivar' and 'test_posmopvar_vs_poslpivar_2d' compare
% the linear SPAN of the generator kernels, which is invariant under every
% invertible relabeling of the decision variables -- including relabelings
% that are not symmetries of the psd cone. A defect that sends the wrong Q
% entry to the wrong generator product leaves the span identical while
% destroying both self-adjointness and positivity. This test is what pins the
% pairing, and it needs no external oracle and no solver.
%
% Three further properties are checked on the same object, since the machinery
% is already in place:
%
%   * self-adjointness, as <y,Px> == <x,Py> at RANDOM decision values -- the
%     Gram is symmetric by naming, so the identity must hold for every
%     admissible d, not only for positive semidefinite ones;
%   * positivity, as <z,Pz> >= 0 once a psd Gram is substituted, judged
%     RELATIVE to the size of the block contributions being summed. That last
%     check is not implied by the identity, but it would fail if the blocks
%     Q_{ce} were paired with the wrong basis operators;
%   * that the Gram is a SINGLE matrix spanning every (space, multi-index)
%     pair, which is what makes the container positive as a whole rather than
%     blockwise, and is the file's headline claim.
%
% The exponent grid of each basis operator is rebuilt here from the degree
% specification. That duplicates a convention -- 'build_exponent_grid' lives
% in 'lpis_sopvar/private' and is not reachable from this folder -- but only a
% convention: the monomial ORDER has to agree for the Qcell blocks to line up,
% while the identity being checked is the semantic content.
%
% MMP, 09/21/2026: Initial coding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
rng(17);

% A second complete copy of the toolbox can sit under .claude/worktrees, and
% 'pietools_path_update' is a single addpath(genpath(...)) which does not skip
% dot-directories, so a result could otherwise come from the wrong copy.
for nm = {'posmopvar','mopquadvar','sopquadvar','lpi_eq_mdopvar'}
    if numel(which(nm{1},'-all'))~=1
        error('test_posmopvar: %s resolves to %d files; fix the path first.',...
              nm{1},numel(which(nm{1},'-all')));
    end
end

w = warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');

D_R  = struct('int',0);
D_L1 = struct('int',1,'mult',1,'joint',2);

% {label, dims, spaces, dom, deg, options}
cases = {
  {'R x L2[s1]',                [1;1],  {{},{'s1'}},              [0,1],      1,          struct()}
  {'R^2 x L2^2[s1]',            [2;2],  {{},{'s1'}},              [0,1],      1,          struct()}
  {'R x L2[s1], shifted',       [1;1],  {{},{'s1'}},              [-1,2],     1,          struct()}
  {'R x L2[s1], deg 2',         [1;1],  {{},{'s1'}},              [0,1],      2,          struct()}
  {'R x L2[s1], per-space deg', [1;1],  {{},{'s1'}},              [0,1],      {D_R,D_L1}, struct()}
  {'L2[s1] x L2[s2]',           [1;1],  {{'s1'},{'s2'}},          [0,1;0,2],  1,          struct()}
  {'L2[s1] x L2[s1,s2]',        [1;1],  {{'s1'},{'s1','s2'}},     [0,1;0,2],  1,          struct()}
  {'L2[s1,s2] x L2[s2]',        [1;1],  {{'s1','s2'},{'s2'}},     [0,1;0,2],  1,          struct()}
  {'R x L2[s1] x L2[s1,s2]',    [1;1;1],{{},{'s1'},{'s1','s2'}},  [0,1;0,2],  1,          struct()}
  {'R x L2[s1], psatz',         [1;1],  {{},{'s1'}},              [0,1],      1,          struct('psatz',1)}
  {'R x L2[s1], separable',     [1;1],  {{},{'s1'}},              [0,1],      1,          struct('sep',true)}
  {'L2[s1] x L2[s2], sep s1',   [1;1],  {{'s1'},{'s2'}},          [0,1;0,2],  1,          struct('sep',[true,false])}
  {'R x L2[s1], mult+lower',    [1;1],  {{},{'s1'}},              [0,1],      1,          struct('include',{{[],[1;2]}})}
  {'R x L2[s1], type sym',      [1;1],  {{},{'s1'}},              [0,1],      1,          struct('type','sym')}
  % Above two variables. These are not more of the same: THREE registry
  % variables are the fewest at which a single block pair can carry an
  % output-only, an input-only AND a shared direction at once, and the fewest
  % at which a side can carry more than one non-shared direction. Below that,
  % 'int_onesided's grouped index decode has at most one entry per group, so
  % an iL/iR transposition or a wrong radix order is unobservable. The last
  % case puts three non-shared directions on one side.
  % Degree 0 with a single basis per space, on a domain NOT starting at 0.
  % This is the minimal reproducer for a defect found on 09/21/2026: the
  % block's coefficient matrix then has exactly ONE nonzero, which is the
  % only case in which 'int_onesided's index arithmetic implicitly expanded
  % to a matrix and scaled the block by the one-sided fan-out. Degree 1 or
  % more hides it (several nonzeros), a multiplier direction hides it
  % (fan-out 1), and a domain starting at 0 hides half of it (the alpha=2
  % constant weight is then zero and drops out). Keep all three conditions.
  {'deg0 {s1,s2} x {s2}, a>0',  [1;1],  {{'s1','s2'},{'s2'}},     [1,2;0,2], 0,          struct('include',{{[2,1],[1]}})}
  {'deg0 {s1,s2} x {s2}, all',  [1;1],  {{'s1','s2'},{'s2'}},     [1,2;0,2], 0,          struct()}
  % Above two variables. THREE registry variables are the fewest at which one
  % block pair can carry an output-only, an input-only AND a shared direction
  % at once, and the fewest at which one side can carry more than one
  % non-shared direction; below that 'int_onesided's grouped index decode has
  % at most one entry per group, so a transposition of iL and iR or a wrong
  % radix order is unobservable. The second case puts three non-shared
  % directions on one side.
  {'nv=3  {s1,s2} x {s2,s3}',   [1;1],  {{'s1','s2'},{'s2','s3'}},[1,2;0,2;-1,1],0,      struct()}
  {'nv=3  R x L2[s1,s2,s3]',    [1;1],  {{},{'s1','s2','s3'}},    [1,2;0,2;-1,1],0,      struct()}
  {'nv=4  {s1,s2} x {s3,s4}',   [1;1],  {{'s1','s2'},{'s3','s4'}},[1,2;0,2;-1,1;0,3],0,  struct()}
  };

npass = 0;
for ic = 1:size(cases,1)
    [lbl,dims,spaces,dom,deg,options] = deal(cases{ic,1}{:});

    vars = reshape(unique([spaces{:}]),1,[]);
    nv = numel(vars);
    prog = make_prog(vars,dom);
    if strcmp(get_opt(options,'type',''),'sym')
        [prog,Pop,Qcell,basis_list] = mopquadvar(prog,dims,spaces,dom,deg,options);
    else
        [prog,Pop,Qcell,basis_list] = posmopvar(prog,dims,spaces,dom,deg,options);
    end
    dvars = cellstr(string(Pop.Zd(:)));
    ndec = numel(dvars);

    % % % The Gram is ONE matrix over every (space, multi-index) pair
    N = size(basis_list,1);
    if ~isequal(size(Qcell),[N,N])
        error(['test_posmopvar: Qcell is %s for %d basis operators; the Gram '...
               'is not a single matrix spanning every pair (case %s).'],...
               mat2str(size(Qcell)),N,lbl);
    end
    if size(unique(basis_list,'rows'),1)~=N
        error('test_posmopvar: basis_list repeats a row (case %s).',lbl);
    end

    % % % Test functions, one per space
    x = cell(1,numel(spaces));      y = cell(1,numel(spaces));
    for k = 1:numel(spaces)
        x{k} = rand_testfun(spaces{k},dims(k));
        y{k} = rand_testfun(spaces{k},dims(k));
    end

    % % % The identity, at random decision values
    dval = randn(ndec,1);
    a_ker = pair_container(Pop,y,x,dval,spaces,dom,vars);
    a_fac = pair_factored(Qcell,basis_list,y,x,dval,dvars,spaces,dims,dom,vars,deg,options);
    scl = max([abs(a_ker),abs(a_fac),1e-30]);
    e_id = abs(a_ker-a_fac)/scl;

    % % % Self-adjointness, from the same two routes' left one
    b_ker = pair_container(Pop,x,y,dval,spaces,dom,vars);
    e_adj = abs(a_ker-b_ker)/max([abs(a_ker),abs(b_ker),1e-30]);

    % % % Positivity with a psd Gram, judged relative to the terms summed
    if strcmp(get_opt(options,'type',''),'sym')
        q_rel = NaN;                % a 'sym' variable is not claimed positive
    else
        dpsd = assign_psd(Qcell,dvars);
        q_rel = inf;
        for r = 1:3
            z = cell(1,numel(spaces));
            for k = 1:numel(spaces),    z{k} = rand_testfun(spaces{k},dims(k));     end
            qz = pair_container(Pop,z,z,dpsd,spaces,dom,vars);
            mag = pair_scale(Pop,z,dpsd,spaces,dom,vars);
            q_rel = min(q_rel,qz/max(mag,1e-30));
        end
    end

    if e_id > 1e-8
        error(['test_posmopvar: the defining identity fails for case ''%s'': '...
               '<y,Px> = %.12g from the kernels but %.12g from Qcell '...
               '(relative %.3g).'],lbl,a_ker,a_fac,e_id);
    end
    if e_adj > 1e-8
        error(['test_posmopvar: the container is not self-adjoint for case '...
               '''%s'': <y,Px> = %.12g but <x,Py> = %.12g (relative %.3g).'],...
               lbl,a_ker,b_ker,e_adj);
    end
    if ~isnan(q_rel) && q_rel < -1e-9
        error(['test_posmopvar: <z,Pz> = %.3g relative to the summed block '...
               'magnitudes with a psd Gram for case ''%s''.'],q_rel,lbl);
    end

    fprintf(['  passed: %-28s N=%-3d ndec=%-6d identity %8.2e  adjoint %8.2e'...
             '  <z,Pz>/scale %s\n'],lbl,N,ndec,e_id,e_adj,relstr(q_rel));
    npass = npass+1;
end

warning(w);
fprintf('test_posmopvar passed (%d of %d cases).\n',npass,size(cases,1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prog = make_prog(vars,dom)
% An LPI program for ANY number of spatial variables.
%
% 'lpiprogram' refuses more than two (lpiprogram.m:117-118), but that check is
% the only obstacle: what it returns is 'sosprogram' plus a vartable holding
% the primary and dummy variables and a 'dom' field, and both are assembled
% here for any nv. Measured: the direct route builds a usable program at
% nv = 1..4, and 'mopquadvar' never reads prog.dom or prog.vartable anyway --
% it only forwards prog to 'sosquadvar'. Note that the '_int' integration
% variables are not in the vartable even in the two-variable case, so
% 'sosquadvar' does not require them to be.

vars = reshape(vars,1,[]);
if numel(vars)<=2
    prog = lpiprogram(polynomial(vars(:)),[],dom);
    return
end
prog = sosprogram(polynomial([]));
dums = strcat(vars,'_dum');
prog.vartable = [prog.vartable; polynomial(vars(:)); polynomial(dums(:))];
prog.dom = dom;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = pair_container(Pop,y,x,dval,spaces,dom,vars)
% <y,Px> = sum_{k,l} <y_k, P_{kl} x_l>, from the STORED kernels.

val = 0;
for k = 1:numel(spaces)
    for l = 1:numel(spaces)
        B = Pop.C{k,l};
        if isempty(B),  continue,   end
        val = val + pair_block(B,y{k},x{l},dval,dom,vars);
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mag = pair_scale(Pop,z,dval,spaces,dom,vars)
% The size of the block contributions summed to form <z,Pz>, with the signs
% removed. This is the denominator that makes the positivity check relative:
% <z,Pz> is a sum of terms of both signs, and an absolute threshold on a
% near-cancelling sum says nothing.

mag = 0;
for k = 1:numel(spaces)
    for l = 1:numel(spaces)
        B = Pop.C{k,l};
        if isempty(B),  continue,   end
        mag = mag + abs(pair_block(B,z{k},z{l},dval,dom,vars));
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = pair_block(B,yk,xl,dval,dom,vars)
% <y, B x> for one block, from the 'sdopvar' kernel definition:
%
%   (Bx)(S2,S3) = sum_gam int_{S1} int_{S3dum} I_gam(S3-S3dum)
%                           K_gam(S2,S3; S3dum,S1) x(S3dum,S1) dS3dum dS1
%
% with index 1 a multiplier delta(s-s'), 2 the lower integral s'<=s and 3 the
% upper -- the convention 'pi_kernel_pairing/apply_indicator' uses. The output
% variables are integrated last. 'pi_kernel_pairing' itself assumes
% vars.in == vars.out and so does not cover the off-diagonal blocks.

O = B.vars.out;     I = B.vars.in;
S3 = intersect(O,I);                % sorted, as the class indexes the cell
S1 = setdiff(I,S3);
Kc = pi_blk_kernels(B,dval);
if numel(Kc)~=3^numel(S3)
    error('test_posmopvar: %d parameter cells for %d shared variables.',...
          numel(Kc),numel(S3));
end

xdum = xl;
for t = 1:numel(I)
    xdum = subs(xdum,polynomial(I(t)),polynomial({[I{t} '_dum']}));
end

Bx = polynomial(zeros(B.dims(1),1));
for q = 1:numel(Kc)
    gam = pi_gamma_index(q,numel(S3));
    e = Kc{q}*xdum;
    for t = 1:numel(S3)
        d = S3{t};      lim = dom(strcmp(vars,d),:);
        sp = polynomial({[d '_dum']});      s = polynomial({d});
        switch gam(t)
            case 1,     e = subs(e,sp,s);
            case 2,     e = int(e,sp,lim(1),s);
            case 3,     e = int(e,sp,s,lim(2));
        end
    end
    for t = 1:numel(S1)
        d = S1{t};      lim = dom(strcmp(vars,d),:);
        e = int(e,polynomial({[d '_dum']}),lim(1),lim(2));
    end
    Bx = Bx + e;
end

expr = yk.'*Bx;
for t = 1:numel(O)
    d = O{t};   lim = dom(strcmp(vars,d),:);
    expr = int(expr,polynomial({d}),lim(1),lim(2));
end
val = double(expr);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = pair_factored(Qcell,basis_list,y,x,dval,dvars,spaces,dims,dom,vars,deg,options)
% <y,Px> computed from Qcell and basis_list alone, by applying each basis
% operator to the test function in closed form and pairing over theta:
%
%   sum_{c,e} int_theta g(theta) (Z_c y)(theta)' Q_{ce} (Z_e x)(theta) dtheta.
%
% This route never reads Pop, so agreement with 'pair_container' is what pins
% the assignment of Q blocks to basis-operator products.

nv = numel(vars);
vars_int = strcat(vars,'_int');
N = size(basis_list,1);

% The psatz weight enters at the integration variable, as in 'poslpivar'.
gfun = polynomial(1);
if get_opt(options,'psatz',0)
    for d = 1:nv
        th = polynomial(vars_int(d));
        gfun = gfun*(th-dom(d,1))*(dom(d,2)-th);
    end
end

U = cell(1,N);      V = cell(1,N);
for c = 1:N
    k = basis_list(c,1);
    U{c} = apply_basis(basis_list(c,:),x{k},spaces{k},dims(k),dom,vars,vars_int,deg,c,basis_list);
    V{c} = apply_basis(basis_list(c,:),y{k},spaces{k},dims(k),dom,vars,vars_int,deg,c,basis_list);
end

expr = polynomial(0);
for c = 1:N
    for e = 1:N
        Q = q_block(Qcell{c,e},dval,dvars);
        if nnz(Q)==0,   continue,   end
        if size(Q,1)~=numel(V{c}) || size(Q,2)~=numel(U{e})
            error(['test_posmopvar: Q block (%d,%d) is %s but the basis '...
                   'operators give %d and %d rows.'],c,e,mat2str(size(Q)),...
                   numel(V{c}),numel(U{e}));
        end
        expr = expr + V{c}.'*(Q*U{e});
    end
end
expr = gfun*expr;
for d = 1:nv
    expr = int(expr,polynomial(vars_int(d)),dom(d,1),dom(d,2));
end
val = double(expr);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = apply_basis(brow,w,own_vars,m,dom,vars,vars_int,deg,c,basis_list)
% (Z_alpha w)(theta) in closed form, an (m*T) x 1 polynomial in theta:
%
%   int_{s} I_alpha(theta-s) kron(w(s), Z^alpha(theta,s)) ds
%
%   alpha=1  delta(theta-s)      -> substitute s := theta
%   alpha=2  I_1(theta-s),  s<=theta -> int from a to theta
%   alpha=3  I_-1(theta-s), s>=theta -> int from theta to b
%   alpha=4  full domain         -> int from a to b
%
% Note kron(w,Z), not kron(Z,w): (I_m kron Z) has the matrix row as the OUTER
% index, which is the layout 'sosquadvar' gives the rows of a Gram block.

nv = numel(vars);
own = find(ismember(vars,own_vars));    nk = numel(own);
alpha = brow(1+own);

E = basis_grid(brow,own,nv,deg,c,basis_list);
T = size(E,1);
names = [reshape(vars_int,[],1); reshape(vars(own),[],1)];
Z = polynomial(speye(T),E,names,[T,1]);

u = kron(w,Z);
for t = 1:nk
    d = vars{own(t)};       lim = dom(own(t),:);
    s = polynomial({d});    th = polynomial(vars_int(own(t)));
    switch alpha(t)
        case 1,     u = subs(u,s,th);
        case 2,     u = int(u,s,lim(1),th);
        case 3,     u = int(u,s,th,lim(2));
        case 4,     u = int(u,s,lim(1),lim(2));
        otherwise
            error('test_posmopvar: multi-index entry %d is not in {1,2,3,4}.',alpha(t));
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E = basis_grid(brow,own,nv,deg,c,basis_list)
% The exponent grid of basis operator c over [theta_1..theta_nv, s^k].
%
% This reimplements the degree expansion rather than calling
% 'build_exponent_grid', which is private to 'lpis_sopvar'. Only the ORDER is
% a shared convention -- rows sorted ascending -- and it has to be shared for
% the Gram blocks to line up. 'subset' caps are not covered here; those are
% exercised by the 'poslpivar_2d' comparison.

k = brow(1);
spec = deg;
if iscell(spec)
    spec = spec{k};
    if iscell(spec)
        spec = spec{sum(basis_list(1:c,1)==k)};
    end
end
if isnumeric(spec) && isscalar(spec)
    spec = struct('int',spec);
end
if isfield(spec,'subset') && ~isempty(spec.subset)
    error('test_posmopvar: this test does not cover subset caps.');
end
cap_int  = expand_cap(getfld(spec,'int',1),nv);
cap_mult = expand_cap(getfld(spec,'mult',getfld(spec,'int',1)),nv);
cap_mult = cap_mult(own);
% A multiplier direction identifies s_d with theta_d, so s_d dependence of
% the basis would be redundant and is excluded.
cap_mult(brow(1+own)==1) = 0;
caps = [cap_int,cap_mult];
jc = getfld(spec,'joint',sum(cap_int)+sum(cap_mult));

E = zeros(1,0);
for i = 1:numel(caps)
    col = (0:caps(i)).';
    E = [repelem(E,numel(col),1), repmat(col,size(E,1),1)];                 %#ok<AGROW>
    E = E(sum(E,2)<=jc,:);
end
E = sortrows(E);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = expand_cap(c,nv)

c = reshape(c,1,[]);
if isscalar(c),     c = repmat(c,1,nv);     end
if numel(c)~=nv
    error('test_posmopvar: a degree cap should be scalar or have nv entries.');
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Q = q_block(names,dval,dvars)
% Numeric values of the decision variables naming one block of the Gram.

names = cellstr(string(names));
[tf,loc] = ismember(names,dvars);
if ~all(tf(:))
    error('test_posmopvar: a decision variable named in Qcell is not in Pop.Zd.');
end
Q = reshape(dval(loc),size(names));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dval = assign_psd(Qcell,dvars)
% Decision variable values making the Gram named by Qcell positive
% semidefinite, with the same symmetry check 'test_possopvar/assign_psd'
% makes: without it a naming defect would present as a positivity defect.

names = cell(size(Qcell,1),1);
for i = 1:size(Qcell,1)
    row = {};
    for j = 1:size(Qcell,2)
        row = [row, cellstr(string(Qcell{i,j}))];                           %#ok<AGROW>
    end
    names{i} = row;
end
names = vertcat(names{:});
nQ = size(names,1);
if size(names,2)~=nQ
    error('test_posmopvar: Qcell does not describe a square matrix.');
end
M = randn(nQ,nQ);
S = (M'*M)/nQ;
[tf,loc] = ismember(names,dvars);
if ~all(tf(:))
    error('test_posmopvar: a decision variable named in Qcell is not in Pop.Zd.');
end
dval = zeros(numel(dvars),1);
dval(loc(:)) = S(:);
if norm(reshape(dval(loc),nQ,nQ)-S,'fro') > 1e-12*max(1,norm(S,'fro'))
    error('test_posmopvar: the naming of Q is not symmetric.');
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = rand_testfun(v,m)
% An m x 1 test function on L_2[v], of degree 2 in each variable, or a
% constant vector when the space is R^m.

f = polynomial(randn(m,1));
for t = 1:numel(v)
    s = polynomial(v(t));
    f = f + randn(m,1)*s + randn(m,1)*s^2;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = get_opt(options,name,dflt)

v = dflt;
if isa(options,'struct') && isfield(options,name) && ~isempty(options.(name))
    v = options.(name);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = getfld(s,name,dflt)

v = dflt;
if isfield(s,name) && ~isempty(s.(name))
    v = s.(name);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = relstr(q)

if isnan(q),    s = '   n/a  (sym)';
else,           s = sprintf('%+8.2e',q);
end

end
