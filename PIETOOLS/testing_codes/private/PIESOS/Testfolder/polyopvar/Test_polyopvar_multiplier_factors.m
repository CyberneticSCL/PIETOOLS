% Test_polyopvar_multiplier_factors.m -- distributed monomials with known
% factors.
%
% A monomial Z(i,:) containing zeros, e.g. Z(i,:) = [1,0,2] for x*m(s)*y,
% is acted on by a tensor-PI operator with more factors than its degree
% nnz(Z(i,:)): factor j acts on state Z(i,j) or, if Z(i,j)=0, is a
% multiplier m(s) from R to a function space (a known factor).
% 'initialize' must not move the zeros of Z, which would destroy the
% correspondence between the factors and the positions in Z.
%
% Oracle: the polynomial P(x,y,z) is evaluated (eval_polyopvar) at
% polynomial states and compared to the product x(s)*m(s)*y(s) of
% 'polynomial' objects.
% Parts:
%   (1) construction/initialize: Z kept, num_factors, degree, matdim,
%       value;
%   (2) invalid combinations of Z and operator error; zeros in monomials
%       with as many factors as their degree are moved to the right;
%   (3) times: with states, with known factors and with constants;
%   (4) plus/minus/vertcat (common_basis, combine_terms);
%   (5) clean_poly (all-zero columns of Z);
%   (6) subs of states by name/position, by linear functions and by
%       constants;
%   (7) innerprod, with known factors in the first, second or both
%       vectors;
%   (8) a constant times a state without spatial variables (no
%       multiplier left);
%   (9) subs of a known function m(s), a degree-0 'polyopvar' (a
%       multiplier from R), for some of the factors of a monomial.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES:
% - Evaluation uses eval_polyopvar (Testfolder/polyopvar), not
%   apply_poly: apply_poly evaluates by subs, which this test checks.

tol = 1e-8;
pvar s
m1 = 1-s^2;     m2 = 2+s;
vals = struct('x',1+2*s,'y',3-s+s^2,'z',-2+s+0.5*s^3,'w',0.5-s+2*s^2,'v',1-3*s);
ncases = 0;

% Identity operators on the states
polyopvar x(s) y(s) z(s) w(s) v(s)
Id = struct('x',x.C.ops{1,1},'y',y.C.ops{1,1},'z',z.C.ops{1,1},'w',w.C.ops{1,1},'v',v.C.ops{1,1});

%% (1) Construction and initialize keep the zeros in place
specs = {   [1,0,2],        {'x',m1,'y'}
            [1,2,0],        {'x','y',m1}
            [0,1,2],        {m1,'x','y'}
            [1,0,0,2],      {'x',m1,m2,'y'}
            [0,1],          {m1,'x'}
            [1,0],          {'x',m1}
            [0,0],          {m1,m2}
            [3,0,1,0,2],    {'z',m2,'x',m1,'y'}
            [1,0,1],        {'x',m1,'x'}    };
for k=1:size(specs,1)
    Zrow = specs{k,1};      facs = specs{k,2};
    f = build_f(Zrow,facs,Id);
    if ~isequal(full(f.Z),Zrow)
        error("Test_polyopvar_multiplier_factors failed (1.%d): Z was changed to [%s].",k,num2str(full(f.Z)))
    end
    [nfac,nops] = num_factors(f);
    if nfac~=numel(facs) || nops~=numel(facs)
        error("Test_polyopvar_multiplier_factors failed (1.%d): num_factors returned %d,%d, expected %d.",k,nfac,nops,numel(facs))
    elseif ~isequal(f.matdim,[1,1]) || full(sum(sum(f.degmat)))~=nnz(Zrow)
        error("Test_polyopvar_multiplier_factors failed (1.%d): matdim or degmat.",k)
    end
    assert_close(eval_f(f,vals),prod_facs(Zrow,facs,vals),tol,sprintf('(1.%d)',k));
    ncases = ncases+1;
end

%% (2) Invalid combinations
% A zero in Z at the position of a factor that is not a multiplier from R
expect_error(@() build_f([1,2,0],{'x',m1,'y'},Id),'multiplier');
expect_error(@() build_f([1,0,2],{'x','y',m1},Id),'multiplier');
% State beyond the last factor
expect_error(@() build_f([1,0,0,2],{'x',m1,'y'},Id),'factors');
% Padding in monomials with as many factors as their degree is still moved to the right
f = build_f([0,1,0,2],{'x','y'},Id);
if ~isequal(full(f.Z),[1,2,0,0])
    error("Test_polyopvar_multiplier_factors failed (2): padding not moved to the right.")
end
assert_close(eval_f(f,vals),vals.x*vals.y,tol,'(2)');
ncases = ncases+5;

%% (3) Products
fm = build_f([1,0,2],{'x',m1,'y'},Id);
fm2 = build_f([1,2,0],{'x','y',m2},Id);
fx = build_f(1,{'x'},Id);
fz = build_f(3,{'z'},Id);
pm = prod_facs([1,0,2],{'x',m1,'y'},vals);
pm2 = prod_facs([1,2,0],{'x','y',m2},vals);
% (3a) known factors in either factor of the product: factors of the second
% factor start after ALL factors of the first, not after its states
tests = {   fm.*fz,     [1,0,2,3],      pm*vals.z
            fz.*fm,     [3,1,0,2],      vals.z*pm
            fm.*fm,     [1,0,2,1,0,2], pm*pm
            fm.*fm2,    [1,0,2,1,2,0], pm*pm2
            fm.*fx,     [1,0,2,1],      pm*vals.x
            fx.*fm,     [1,1,0,2],      vals.x*pm   };
for k=1:size(tests,1)
    if ~isequal(full(tests{k,1}.Z),tests{k,2})
        error("Test_polyopvar_multiplier_factors failed (3a.%d): Z = [%s], expected [%s].",k,num2str(full(tests{k,1}.Z)),num2str(tests{k,2}))
    end
    assert_close(eval_f(tests{k,1},vals),tests{k,3},tol,sprintf('(3a.%d)',k));
    ncases = ncases+1;
end
% (3b) the product of a function with a constant is a function with a known factor
c = constant2polyopvar(2);
tests = {   c.*fx,      [0,1],          2*vals.x
            fx.*c,      [1,0],          2*vals.x
            c.*fm,      [0,1,0,2],      2*pm        };
for k=1:size(tests,1)
    % Z may omit trailing zeros (factors beyond its width are multipliers from R)
    Zk = full(tests{k,1}.Z);    Zk(end+1:numel(tests{k,2})) = 0;
    if ~isequal(Zk,tests{k,2}) || num_factors(tests{k,1})~=numel(tests{k,2})
        error("Test_polyopvar_multiplier_factors failed (3b.%d): Z = [%s] with %d factors, expected [%s].",k,num2str(full(tests{k,1}.Z)),num_factors(tests{k,1}),num2str(tests{k,2}))
    end
    assert_close(eval_f(tests{k,1},vals),tests{k,3},tol,sprintf('(3b.%d)',k));
    ncases = ncases+1;
end

% (3c) a monomial without operators (a zero block of a product) has no factors to tell it apart
% from a monomial with the same Z that has operators, and is merged with it (factor_key), so no
% such empty column is left in the basis of a product of vectors of monomials, with a constant
vt = [x;y;z];
Zb2 = dmonomials(vt,0:3);
P2 = Zb2.*Zb2;
[~,~,is_empty] = num_factors(P2);
Zfull = full(P2.Z);
if any(is_empty) && any(ismember(Zfull(is_empty,:),Zfull(~is_empty,:),'rows'))
    error("Test_polyopvar_multiplier_factors failed (3c): an empty monomial is not merged with one with the same Z.")
end
ncases = ncases+1;

%% (4) Sums
% (4a) equal monomials are merged, different ones (also in the same states) are not
g = fm+fm;
if size(g.Z,1)~=1 || ~isequal(full(g.Z),[1,0,2])
    error("Test_polyopvar_multiplier_factors failed (4a): monomials not merged.")
end
assert_close(eval_f(g,vals),2*pm,tol,'(4a) fm+fm');
g = fm+fm2;
if size(g.Z,1)~=2
    error("Test_polyopvar_multiplier_factors failed (4a): different monomials merged.")
end
assert_close(eval_f(g,vals),pm+pm2,tol,'(4a) fm+fm2');
g = fm+fx.*fz+fx;
assert_close(eval_f(g,vals),pm+vals.x*vals.z+vals.x,tol,'(4a) fm+x*z+x');
g = fm-fm+fx;
if size(g.Z,1)~=2
    error("Test_polyopvar_multiplier_factors failed (4a): fm-fm+x has %d monomials, expected 2.",size(g.Z,1))
end
assert_close(eval_f(g,vals),vals.x,tol,'(4a) fm-fm+x');
% (4a2) the factor count is part of the identity of a monomial: x o c (Z=[1,0], two factors) and
% x (Z=[1], one factor) are the same row of Z once padded, but are not merged
fxc = build_f([1,0],{'x',m2},Id);
g = fxc+fx;
if size(g.Z,1)~=2
    error("Test_polyopvar_multiplier_factors failed (4a2): x o c and x are merged.")
end
assert_close(eval_f(g,vals),vals.x*m2+vals.x,tol,'(4a2) x o c + x');
g = fx+fxc;
assert_close(eval_f(g,vals),vals.x+vals.x*m2,tol,'(4a2) x + x o c');
g = fxc+fxc+fx;
if size(g.Z,1)~=2
    error("Test_polyopvar_multiplier_factors failed (4a2): %d monomials, expected 2.",size(g.Z,1))
end
assert_close(eval_f(g,vals),2*vals.x*m2+vals.x,tol,'(4a2) 2 x o c + x');
ncases = ncases+3;
% (4b) vertical concatenation
vcat = [fm; fx.*fz];
Y = eval_f(vcat,vals);
assert_close(Y(1),pm,tol,'(4b) [fm;x*z] row 1');
assert_close(Y(2),vals.x*vals.z,tol,'(4b) [fm;x*z] row 2');
ncases = ncases+7;

%% (5) clean_poly only removes trailing all-zero columns of Z
g = fm;     g.Z = [g.Z,sparse(1,2)];
g = clean_poly(g);
if ~isequal(full(g.Z),[1,0,2])
    error("Test_polyopvar_multiplier_factors failed (5): trailing zero columns not removed, Z = [%s].",num2str(full(g.Z)))
end
assert_close(eval_f(g,vals),pm,tol,'(5) trailing');
% Interior column of Z that is zero in all monomials: [1,0,2;1,0,1]
g = clean_poly(fm+build_f([1,0,1],{'x',m1,'x'},Id));
if size(g.Z,2)~=3 || any(g.Z(:,2))
    error("Test_polyopvar_multiplier_factors failed (5): interior zero column removed, Z = [%s].",num2str(full(g.Z(:)')))
end
assert_close(eval_f(g,vals),pm+vals.x*m1*vals.x,tol,'(5) interior');
ncases = ncases+2;

%% (6) Substitution
% New functions of the states w,v, with their values. The operator on a monomial of
% more than one factor can only be composed with operators of a single factor
% (see Test_polyopvar_subs/case_fctr_function), so the new functions are linear in w,v
valw = {    2*w,                            2*vals.w
            w+2*v,                          vals.w+2*vals.v     };
for k=[1,2,3,4,5,8,9]
    Zrow = specs{k,1};      facs = specs{k,2};
    f = build_f(Zrow,facs,Id);
    pos = find(Zrow);
    for kv=1:size(valw,1)
        names = {'x','y','z'};
        for st = unique(names(Zrow(pos)),'stable')
            % (6a) Substitute the state in all factors, by name
            v = vals;           v.(st{1}) = valw{kv,2};
            Yref = prod_facs(Zrow,facs,v);
            fnew = subs(f,st,valw{kv,1});
            assert_close(eval_f(fnew,vals),Yref,tol,sprintf('(6a) spec %d, state %s, val %d',k,st{1},kv));
            ncases = ncases+1;
            % (6b) Substitute a single factor, by position
            ip = find(strcmp(names(Zrow(pos)),st{1}));
            if numel(ip)==1
                fnew = subs(f,pos(ip),valw{kv,1});
                assert_close(eval_f(fnew,vals),Yref,tol,sprintf('(6b) spec %d, factor %d, val %d',k,pos(ip),kv));
                ncases = ncases+1;
            end
        end
    end
end
% Substituting a factor that is a multiplier from R, or beyond the last factor, is not possible
expect_error(@() subs(fm,2,valw{1,1}),'multiplier');
expect_error(@() subs(fm,4,valw{1,1}),'degree');
ncases = ncases+2;

% (6c) Several states at once
fnew = subs(fm,{'x';'y'},[valw{1,1};valw{2,1}]);
v = vals;   v.x = valw{1,2};    v.y = valw{2,2};
assert_close(eval_f(fnew,vals),prod_facs([1,0,2],{'x',m1,'y'},v),tol,'(6c)');
ncases = ncases+1;

% (6d) Constants: all states replaced, the known factors remain
cx = 1+s;   cy = 2-s^2;
v = vals;   v.x = cx;       v.y = cy;
Yref = prod_facs([1,0,2],{'x',m1,'y'},v);
fnew = subs(fm,{'x';'y'},[cx;cy]);
if any(any(fnew.Z))
    error("Test_polyopvar_multiplier_factors failed (6d): result should be a constant (degree-0) term.")
end
assert_close(eval_f(fnew,vals),Yref,tol,'(6d) by name');
fnew = subs(fm,[1,3],[cx;cy]);
assert_close(eval_f(fnew,vals),Yref,tol,'(6d) by position');
ncases = ncases+2;

%% (7) Inner product: the monomials of the second vector follow all factors of the first
% (sosprogramming/innerprod.m), compared with the integral computed directly
tests = {   innerprod(fm,fz),   int(pm*vals.z,s,0,1)
            innerprod(fz,fm),   int(vals.z*pm,s,0,1)
            innerprod(fm,fm2),  int(pm*pm2,s,0,1)       };
for k=1:size(tests,1)
    assert_close(eval_V(tests{k,1},vals),tests{k,2},tol,sprintf('(7.%d)',k));
    ncases = ncases+1;
end

%% (8) A constant times a state that does not depend on any spatial variable: otimes
% merges the two scalar factors, so there is no multiplier from R left in the operator,
% and Z must be [1] (checked through the value, as 'initialize' moves padding)
polyopvar a
c = constant2polyopvar(2);
pc = c.*a;
if ~isequal(full(pc.Z),[1,0]) && ~isequal(full(pc.Z),1)
    error("Test_polyopvar_multiplier_factors failed (8): Z = [%s], expected [1].",num2str(full(pc.Z)))
end
assert_close(eval_polyopvar(pc,{polynomial(3)}),polynomial(6),tol,'(8) c.*a');
ncases = ncases+1;

%% (9) Substitution of a known function m(s), a degree-0 polyopvar (a multiplier from R), for
% some of the factors of a monomial: the factor remains in the operator, as a zero in Z
mc1 = 3-s;      mc2 = 1+s+s^2;
kc1 = known_fun(mc1);       kc2 = known_fun(mc2);
fxx = build_f([1,0,1],{'x',m1,'x'},Id);
vx = vals;      vx.x = mc1;
vy = vals;      vy.y = mc2;
vxy = vals;     vxy.x = mc1;    vxy.y = mc2;
tests = {   subs(fm,{'x'},kc1),             [0,0,2],    prod_facs([1,0,2],{'x',m1,'y'},vx)
            subs(fm,1,kc1),                 [0,0,2],    prod_facs([1,0,2],{'x',m1,'y'},vx)
            subs(fm,{'y'},kc2),             [1,0,0],    prod_facs([1,0,2],{'x',m1,'y'},vy)
            subs(fm,3,kc2),                 [1,0,0],    prod_facs([1,0,2],{'x',m1,'y'},vy)
            subs(fm,{'x';'y'},[kc1;kc2]),   [0,0,0],    prod_facs([1,0,2],{'x',m1,'y'},vxy)
            subs(fxx,{'x'},kc1),            [0,0,0],    mc1*m1*mc1
            subs(fxx,3,kc1),                [1,0,0],    vals.x*m1*mc1
            subs(fm2,{'y'},kc1),            [1,0,0],    vals.x*mc1*m2   };
for k=1:size(tests,1)
    Zk = full(tests{k,1}.Z);    Ze = tests{k,2};
    n = max(numel(Zk),numel(Ze));   Zk(end+1:n) = 0;    Ze(end+1:n) = 0;
    if size(Zk,1)~=1 || ~isequal(Zk,Ze)
        error("Test_polyopvar_multiplier_factors failed (9.%d): Z = [%s], expected [%s].",k,num2str(Zk(:)'),num2str(Ze))
    end
    assert_close(eval_f(tests{k,1},vals),tests{k,3},tol,sprintf('(9.%d)',k));
    ncases = ncases+1;
end
% A known function and a linear function of w in the same substitute
valmix = 2*w+kc1;
g = subs(fm,{'x'},valmix);
assert_close(eval_f(g,vals),(2*vals.w+mc1)*m1*vals.y,tol,'(9) 2*w+m(s)');
% Operations on the result: products, sums
g = subs(fm,{'x'},kc1);
assert_close(eval_f(g.*fz,vals),mc1*m1*vals.y*vals.z,tol,'(9) times');
assert_close(eval_f(fz.*g,vals),vals.z*mc1*m1*vals.y,tol,'(9) times, reversed');
assert_close(eval_f(g+fm,vals),mc1*m1*vals.y+pm,tol,'(9) plus');
% A functional: the known function replaces one of the two factors in <x,y>
V = innerprod(fx,build_f(2,{'y'},Id));
tests = {   subs(V,{'y'},kc1),  int(vals.x*mc1,s,0,1)
            subs(V,2,kc1),      int(vals.x*mc1,s,0,1)
            subs(V,{'x'},kc2),  int(mc2*vals.y,s,0,1)   };
for k=1:size(tests,1)
    assert_close(eval_V(tests{k,1},vals),tests{k,2},tol,sprintf('(9) functional %d',k));
end
ncases = ncases+8+3;

fprintf('Test_polyopvar_multiplier_factors passed (%d cases).\n',ncases)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function f = build_f(Zrow,facs,Id,names)
% Polynomial f(x,y,z)(s) = Kop*Z(x,y,z), with Kop the pointwise product of the
% factors facs{j}: the identity operator on a state (name given), or the
% multiplier m(s) from R to L2[0,1] (polynomial given). Z = Zrow, with the
% states numbered as in names.
if nargin<4
    names = {'x','y','z'};
end
K = [];
for j=1:numel(facs)
    if ischar(facs{j})
        Fj = Id.(facs{j});
    else
        Fj = mat2tensopvar(facs{j},struct('in',{cell(1,0)},'out',{{'s'}}),struct('in',zeros(0,2),'out',[0,1]));
    end
    if isempty(K),  K = Fj;
    else,           K = otimes(K,Fj,[1,2]);     end
end
f = polyopvar();
f.varname = names(:);       f.Z = sparse(Zrow);     f.C = tensopmat(K);
f.dims.in = ones(numel(names),1);
f.vars.in = {'s'};          f.dom.in = [0,1];       f.depmat.in = true(numel(names),1);
f = initialize(f);
end

function g = known_fun(m)
% Degree-0 polyopvar representing the known function m(s) = polynomial in s on [0,1]:
% a multiplier from R to L2[0,1] on the dummy input 1
M = mat2tensopvar(m,struct('in',{cell(1,0)},'out',{{'s'}}),struct('in',zeros(0,2),'out',[0,1]));
g = polyopvar();
g.varname = cell(0,1);      g.Z = sparse(1,0);      g.C = tensopmat(M);
g = initialize(g);
end

function Y = eval_V(V,vals)
% Evaluate the functional V (intop coefficients) at the state values vals.(name): the
% dummy variables of s are numbered over the factors, s_j, as in apply_op of 'intop'
Y = 0;
Zf = full(V.Z);
for k=1:size(Zf,1)
    op = V.C.ops{1,k};
    if isempty(op),     continue,   end
    idcs = Zf(k,Zf(k,:)~=0);
    v = 1;
    for j=1:numel(idcs)
        xj = vals.(V.varname{idcs(j)});
        if numel(idcs)>1,   xj.varname(strcmp(xj.varname,'s')) = {['s_',num2str(j)]};   end
        v = kron(v,xj);
    end
    Y = Y + apply_op(op,v);
end
end

function Y = prod_facs(Zrow,facs,vals)
% Product of the values of the factors: state values for the states, m(s) for the known factors
names = {'x','y','z'};
Y = polynomial(1);
for j=1:numel(facs)
    if ischar(facs{j}),     Y = Y*vals.(names{Zrow(j)});
    else,                   Y = Y*facs{j};      end
end
end

function Y = eval_f(f,vals)
% Evaluate f at the state values vals.(name); unspecified states are set to 0
xv = cell(numel(f.varname),1);
for k=1:numel(xv)
    if isfield(vals,f.varname{k}),  xv{k} = vals.(f.varname{k});
    else,                           xv{k} = polynomial(0);          end
end
Y = eval_polyopvar(f,xv);
end

function assert_close(Ynew,Yref,tol,label)
Yc = polynomial(Yref).coeff;
if isempty(Yc) || max(abs(full(Yc(:))))<1e-6
    error("Test_polyopvar_multiplier_factors (%s): reference value is zero; test would be vacuous.",label)
end
err = polynomial(Ynew-Yref);
if ~isempty(err.coeff) && max(abs(full(err.coeff(:))))>tol
    error("Test_polyopvar_multiplier_factors failed (%s): residual %.3e.",label,max(abs(full(err.coeff(:)))))
end
end

function expect_error(fun,str)
try
    fun();
catch ME
    if ~contains(ME.message,str)
        error("Test_polyopvar_multiplier_factors failed: unexpected error message '%s' (expected to contain '%s').",ME.message,str)
    end
    return
end
error("Test_polyopvar_multiplier_factors failed: expected an error containing '%s', but none was thrown.",str)
end
