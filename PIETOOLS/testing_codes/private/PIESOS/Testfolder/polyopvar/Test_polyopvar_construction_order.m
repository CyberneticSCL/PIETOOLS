% Test_polyopvar_construction_order.m -- verify the convenience
% constructor polyopvar(name,vars,dom) sorts P.vars.in alphabetically and
% permutes P.dom.in to match, regardless of the order the spatial
% variables (and their domains) were supplied in; P.vars.out/P.dom.out
% always match the (independently, alphabetically sorted)
% P.C.vars.out/P.C.dom.out. Also confirms two states built over the same
% variables in opposite order agree exactly on vars.in/dom.in, and that
% the state's own action and '+' are unaffected. Oracle: the sorted
% names/domains written out by hand, and apply_poly values against the
% numeric polynomials.

% DJ (CC), 10/15/2026: Initial coding;

% NOTES:
% - P.vars.in (variables the state is defined in) differs in meaning from
%   P.C.vars.in (variables of the monomial basis): only vars.out must
%   match C.
% - Per-variable (non-uniform) domains are used so the dom.in permutation
%   is exercised.

tol = 1e-8;

for trial = 1:8
    rng(7000+trial);

    pvar sB sA
    x = polyopvar('x',[sB,sA],[2,3;0,1]);   % sB->[2,3], sA->[0,1]

    if ~isequal(x.vars.in,{'sA','sB'}) || ~isequal(x.dom.in,[0,1;2,3])
        error("Test_polyopvar_construction_order failed: vars.in/dom.in not sorted/permuted correctly (trial=%d).",trial)
    elseif ~isequal(x.vars.out,sort(x.vars.out)) || ~isequal(x.vars.out,x.C.vars.out)
        error("Test_polyopvar_construction_order failed: vars.out not sorted / does not match C.vars.out (trial=%d).",trial)
    elseif ~isequal(x.dom.out,x.C.dom.out)
        error("Test_polyopvar_construction_order failed: dom.out does not match C.dom.out (trial=%d).",trial)
    end

    xv = rand_poly([1,1],{'sA','sB'},2,5);
    Yx = apply_poly(x,{xv});
    err = cleanpoly(Yx-xv,1e-10);
    if max(max(abs(full(err.C))))>tol
        error("Test_polyopvar_construction_order failed: state's own action changed (trial=%d).",trial)
    end

    % same variables and domains in OPPOSITE order: agree on vars.in/dom.in
    y = polyopvar('y',[sA,sB],[0,1;2,3]);   % sA->[0,1], sB->[2,3] (same domains, forward order)
    if ~isequal(x.vars,y.vars) || ~isequal(x.dom,y.dom)
        error("Test_polyopvar_construction_order failed: states built with reversed variable order still disagree (trial=%d).",trial)
    end

    yv = rand_poly([1,1],{'sA','sB'},2,5);
    z = x+y;
    Yz = apply_poly(z,{xv,yv});
    err2 = cleanpoly(Yz-(xv+yv),1e-10);
    if max(max(abs(full(err2.C))))>tol
        error("Test_polyopvar_construction_order failed: x+y (reversed variable order) gave the wrong value (trial=%d).",trial)
    end
end

%%%% A genuine (non-pairwise) 3-variable permutation, with distinct domains.
pvar p1 p2 p3
u = polyopvar('u',[p3,p1,p2],[4,5;0,1;2,3]);   % p3->[4,5], p1->[0,1], p2->[2,3]
if ~isequal(u.vars.in,{'p1','p2','p3'}) || ~isequal(u.dom.in,[0,1;2,3;4,5])
    error("Test_polyopvar_construction_order failed: 3-variable vars.in/dom.in not sorted/permuted correctly.")
elseif ~isequal(u.vars.out,{'p1','p2','p3'}) || ~isequal(u.vars.out,u.C.vars.out)
    error("Test_polyopvar_construction_order failed: 3-variable vars.out not sorted/consistent with C.")
end
uv = rand_poly([1,1],{'p1','p2','p3'},2,5);
Yu = apply_poly(u,{uv});
erru = cleanpoly(Yu-uv,1e-10);
if max(max(abs(full(erru.C))))>tol
    error("Test_polyopvar_construction_order failed: 3-variable state's own action changed.")
end

disp('Test_polyopvar_construction_order passed.')
