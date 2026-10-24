% Test_tensopvar_apply_op.m -- tests apply_op of a 'tensopvar' with nested factors, in two
% parts, the descriptions of which follow (the first, formerly
% Test_tensopvar_apply_op_nested, and the second, formerly Test_tensopvar_apply_op_nested_tensor).
%
% ============ Part 1: check_nested ============
% nested factors (a 'tensopvar' as factor of another), applied to a cell
% with one input for each factor, counting those of the nested operators
% (num_factors(Cop)), as apply_poly/subs pass them for the operator acting
% on a monomial of a 'polyopvar', as well as to the 1 x d cell with the
% inputs of nested factors as cells.
%
% Such operators arise from the sum of monomials with the same Z and the
% same number of factors, with factors that depend on different variables
% (tensopvar/plus), and as the pointwise product of operators in a
% Kronecker product (nonlinear terms of a PIE).
% Verification is semantic and independent of the routine under test:
%   (1) 'polyopvar': the sum of monomials x*m1*y + x*m2*y with the m_k in
%       different variables, evaluated by apply_poly, vs. the product of
%       the polynomials x*(m1+m2)*y;
%   (2) 'tensopvar': the output of apply_op for the flat cell of leaf
%       inputs vs. (a) the value for the 1 x d form of the inputs, with
%       the inputs of nested factors as cells (built here by
%       nest_leaves), and (b) apply_nested, a recursion over the nested
%       operators that takes the inputs of the leaf operators one by one
%       from the flat list;
%   (3) operators of which the nesting differs between the terms
%       (apply_nested), which the 1 x d form cannot represent;
%   (4) wrong numbers or dimensions of inputs are rejected (the dimension
%       of the input of a factor after a nested one is checked against the
%       right factor), and non-nested operators are unaffected.
%
% Sums of three or more monomials are built from the left and from the
% right (cases 1.x).
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/22/2026: Also compare cases 2.1, 2.2 with apply_nested, as
%   apply_op now numbers the copies of the output variables of nested
%   tensor products (s_1,s_2,s_3); see Part 2 below;

% NOTES
% - Substitution of some, not all, states in a nested operator: see
%   Test_tensopvar_mtimes/case_nested.

% ============ Part 2: check_nested_tensor ============
% product (otype(1)=2) of which some factors are 'tensopvar' objects
% themselves (nested factors): a pointwise product A1.*A2, a tensor product
% N = A1 o B1 with several copies of the output variable, a single factor
% B1 (with otype [2,2] or [1,2]), and combinations of these, in up to four
% levels. The operators are built by hand (ops, otype, initialize), so
% that depmat has one row for each factor of each layer, with the
% dependency (deps) of the factor for a nested factor, as in otimes with
% operands of different types.
%
% Oracle (independent of apply_op): each of the 'sopvar' operators is
% applied to its own input by apply_sopvar; the output variable s is
% renamed to its k-th copy s_k, with k numbering the sopvar factors in the
% product in order (separately for each output variable; the sopvar
% factors of a pointwise product share the copy); the results are
% multiplied, summed over the terms.
%
% DJ (CC), 10/22/2026: Initial coding;

% NOTES
% - The inputs are passed as a flat cell, one for each sopvar factor in
%   order (num_factors), and, for some cases, as the 1 x d cell with the
%   inputs of nested factors as cells.
% - A pointwise product with both numbered (s_1,s_2) and unnumbered (s)
%   copies of a variable in a nested factor of a tensor product cannot be
%   renamed unambiguously: apply_op must raise an error (case 9).

% DJ (CC), 10/23/2026: Merge of Test_tensopvar_apply_op_nested and
%   Test_tensopvar_apply_op_nested_tensor;

check_nested();
check_nested_tensor();
disp('Test_tensopvar_apply_op passed.')


function check_nested()
tol = 1e-8;
ncases = 0;
pvar s s1 s2

%% (1) polyopvar: nested sums of monomials with known factors, evaluated by apply_poly
m1 = 1-s1^2;        m2 = 2+s2+s1*s2;        m3 = 3-s2;      m4 = 1+s1+s2^2;
polyopvar x(s1,s2) y(s1,s2)
Id = struct('x',x.C.ops{1,1},'y',y.C.ops{1,1});
xv = 1+s1;          yv = 2+s2-s1*s2;
% Each row: Z of the monomials, factors, coefficients, for which the sum is nested
specs = {   [1,0,2],    {{'x',m1,'y'},{'x',m3,'y'}},                  [1,1]
            [1,0,2],    {{'x',m1,'y'},{'x',m3,'y'}},                  [2,-1]
            [0,1,2],    {{m1,'x','y'},{m3,'x','y'}},                  [1,1]
            [1,2,0],    {{'x','y',m1},{'x','y',m3}},                  [1,1]
            [1,0],      {{'x',m1},{'x',m3}},                          [1,1]
            [1,0,0,2],  {{'x',m1,m3,'y'},{'x',m3,m1,'y'}},            [1,1]
            [1,0,2],    {{'x',m1,'y'},{'x',m3,'y'},{'x',m4,'y'}},     [1,1,1]
            [1,0,2],    {{'x',m1,'y'},{'x',m3,'y'},{'x',m4,'y'}},     [1,-2,3]   };
for k=1:size(specs,1)
    Zrow = specs{k,1};      facs = specs{k,2};      cf = specs{k,3};
    % Sum from the left, (f1+f2)+f3, and from the right, f1+(f2+f3)
    Yref = 0;       fl = cell(1,numel(facs));
    for l=1:numel(facs)
        fl{l} = cf(l)*build_f(Zrow,facs{l},Id);
        Yref = Yref + cf(l)*prod_facs(Zrow,facs{l},xv,yv);
    end
    f = fl{end};
    for l=numel(fl)-1:-1:1
        f = fl{l} + f;
    end
    fleft = fl{1};
    for l=2:numel(fl)
        fleft = fleft + fl{l};
    end
    for f = {combine_terms(f),combine_terms(fleft)}
        f = f{1};
        Top = f.C.ops{1,1};
        if ~isa(Top,'tensopvar') || ~any(cellfun('isclass',Top.ops(:),'tensopvar'))
            error("Test_tensopvar_apply_op failed (1.%d): setup: the operator should be nested.",k)
        end
        assert_close(apply_poly(f,{xv,yv}),Yref,tol,sprintf('(1.%d)',k));
        ncases = ncases+1;
    end
end

%% (2) tensopvar: flat inputs vs the 1 x d form of the inputs, and vs apply_nested
vars = struct('in',{{'s'}},'out',{{'s'}});
dom = struct('in',[0,1],'out',[0,1]);
dep1 = struct('in',true,'out',true);
pd = struct('in',1,'out',1);
mult = @() mat2tensopvar(randn+randn*s+randn*s^2,vars,dom).ops{1,1};      % multiplier on L2[0,1]
intl = @() rand_tensopvar([1,1],vars,dom,dep1,1,[2,2],pd,0.5).ops{1,1};    % random 3-PI operator
for variant=1:6
    rng(100+variant);
    switch variant
    case 1      % {{N1 o S1} o N2}, N1 = P1*P2 and N2 = P3*P4 pointwise products, S1 an integral: 5 leaves
        N1 = pointwise(mult(),mult());     S1 = intl();      N2 = pointwise(mult(),mult());
        R12 = kron_op({N1,S1});            Bop = kron_op({R12,N2});
    case 2      % {N1 o {S1 o N2}}, with an integral in N2: 5 leaves
        N1 = pointwise(mult(),mult());     S1 = mult();      N2 = pointwise(mult(),intl());
        R12 = kron_op({S1,N2});            Bop = kron_op({N1,R12});
    case 3      % one nested pointwise product per term
        Bop = kron_op({pointwise(mult(),mult());pointwise(mult(),mult())});
    case 4      % pointwise product of three operators as the first factor, then an integral, two terms
        Bop = kron_op({pointwise(pointwise(mult(),mult()),mult()),intl(); ...
                       pointwise(pointwise(mult(),mult()),mult()),intl()});
    case 5      % a nested operator with a single factor: as many inputs as factors
        Bop = kron_op({kron_op({mult()})});
    case 6      % nested operator with several terms, each a pointwise product, two factors in the outer operator
        N1 = pointwise_terms(mult(),mult(),mult(),mult());
        Bop = kron_op({N1,mult()});
    end
    nl = num_factors(Bop);
    pvars = polynomial(Bop.vars.in(:));
    leaves = cell(1,nl);
    for l=1:nl
        leaves{l} = rand_poly([1,1],pvars,2,10);
    end
    Yflat = apply_op(Bop,leaves);
    Ystruct = apply_op(Bop,nest_leaves(Bop,leaves));
    assert_close(Yflat,Ystruct,tol,sprintf('(2.%d) vs 1 x d form',variant));
%   if variant>=3                                                           % DJ (CC), 10/22/2026 (was)
        Yref = apply_nested(Bop,leaves,1,zeros(1,numel(Bop.vars.out)),Bop.vars.out);
        assert_close(Yflat,Yref,tol,sprintf('(2.%d) vs apply_nested',variant));
%   end                                                                     % DJ (CC), 10/22/2026 (was)
    ncases = ncases+1;
end

%% (3) The nesting differs between the terms: {N(2 leaves), S} and {S, N(2 leaves)}
% multipliers from R, so that the dependency of the factors on the variables is the same in each term
mR = @() mat2tensopvar(randn+randn*s+randn*s^2,struct('in',{cell(1,0)},'out',{{'s'}}),struct('in',zeros(0,2),'out',[0,1])).ops{1,1};
for variant=1:4
    rng(200+variant);
    switch variant
    case 1
        Bop = kron_op({pointwise(mR(),mR()),mR(); mR(),pointwise(mR(),mR())});
    case 2
        Bop = kron_op({pointwise(mR(),mR()),mR(); mR(),pointwise(mR(),mR()); pointwise(mR(),mR()),mR()});
    case 3
        Bop = kron_op({pointwise(mR(),mR()),mR(),mR(); mR(),mR(),pointwise(mR(),mR())});
    case 4  % nesting of depth two, in different factors
        Bop = kron_op({pointwise(pointwise(mR(),mR()),mR()),mR(); mR(),pointwise(pointwise(mR(),mR()),mR())});
    end
    nl = num_factors(Bop);
    leaves = num2cell(ones(1,nl));
    leaves = cellfun(@(c) c*(1+rand),leaves,'UniformOutput',false);      % column dimension 1: any scalar
    Yflat = apply_op(Bop,leaves);
    Yref = apply_nested(Bop,leaves,1,zeros(1,numel(Bop.vars.out)),Bop.vars.out);
    assert_close(Yflat,Yref,tol,sprintf('(3.%d)',variant));
    ncases = ncases+1;
end

%% (4) Wrong numbers or dimensions of inputs, and operators without nesting
% a factor after the nested one with a column dimension of 2: its input is the third of the flat list
mR2 = mat2tensopvar([1+s,2*s^2],struct('in',{cell(1,0)},'out',{{'s'}}),struct('in',zeros(0,2),'out',[0,1])).ops{1,1};
Bop = kron_op({pointwise(mR(),mR()),mR2});
leaves = {1.3,0.7,[1;-2]};
assert_close(apply_op(Bop,leaves),apply_nested(Bop,leaves,1,zeros(1,1),Bop.vars.out),tol,'(4) column dimension 2');
expect_error(@() apply_op(Bop,{1.3,0.7,1}),'Column dimension');
expect_error(@() apply_op(Bop,{1.3,[1;-2],0.7}),'Column dimension');
Bop = kron_op({pointwise(mR(),mR()),mR()});         % 2 factors, 3 leaves
expect_error(@() apply_op(Bop,{1,1}),'1xd cell');       % d inputs: the 1 x d form, and a number is not the cell of a nested factor
expect_error(@() apply_op(Bop,{1,1,1,1}),'input for each factor');
expect_error(@() apply_op(Bop,{1}),'input for each factor');
% the terms have different numbers of leaves: {N(2), S} and {N(3), S}
Bad = kron_op({pointwise(mR(),mR()),mR(); pointwise(pointwise(mR(),mR()),mR()),mR()});
expect_error(@() apply_op(Bad,{1,1,1}),'number of inputs must match');
% a 1 x d cell for an operator without nesting, and a bare input for a single factor, as before
T = kron_op({mult(),mult()});
v1 = rand_poly([1,1],polynomial({'s'}),2,10);       v2 = rand_poly([1,1],polynomial({'s'}),2,10);
Y = apply_op(T,{v1,v2});
Ya = apply_sopvar(T.ops{1,1},v1);       Yb = apply_sopvar(T.ops{1,2},v2);
Ya.varname(strcmp(Ya.varname,'s')) = {'s_1'};           Yb.varname(strcmp(Yb.varname,'s')) = {'s_2'};
assert_close(Y,Ya*Yb,tol,'(4) no nesting');
T1 = kron_op({mult()});
Y1 = apply_op(T1,v1);       Y1c = apply_op(T1,{v1});
assert_close(Y1,Y1c,tol,'(4) single factor');
expect_error(@() apply_op(T,v1),'1xd cell');
expect_error(@() apply_op(T,{v1,v2,v1}),'1xd cell');
ncases = ncases+11;
end

function check_nested_tensor()
tol = 1e-8;
ncases = 0;
rng(41);

% Output variable y (input x), and y2 (input x2)
dom01 = struct('in',[0,1],'out',[0,1]);
vx_y = struct('in',{{'x'}},'out',{{'y'}});
vx_y2 = struct('in',{{'x2'}},'out',{{'y2'}});
vx_yy2 = struct('in',{{'x'}},'out',{{'y','y2'}});
lf_y = @() rand_sopvar([1,1],vx_y,dom01,struct('in',2,'out',2),0.6);
lf_y2 = @() rand_sopvar([1,1],vx_y2,dom01,struct('in',2,'out',2),0.6);
lf_yy2 = @() rand_sopvar([1,1],vx_yy2,struct('in',[0,1],'out',[0,1;0,1]),struct('in',2,'out',[2,2]),0.6);
px = polynomial({'x'});     px2 = polynomial({'x2'});
rp = @() rand_poly([1,1],px,2,50);
rp2 = @() rand_poly([1,1],px2,2,50);
% Output of a sopvar operator R on its input v
ap = @(R,v) apply_sopvar(R,v);

%% (1) Pointwise product A1.*A2 and a single factor B1: Apw o B1 and B1 o Apw
for otB=[2,1]       % otype(1) of the single factor B1 has no effect
    a1 = lf_y();    a2 = lf_y();    b1 = lf_y();
    v1 = rp();      v2 = rp();      v3 = rp();
    Apw = mk({a1,a2},[1,2]);
    B1 = mk({b1},[otB,2]);
    pw = ap(a1,v1).*ap(a2,v2);
    C = mk({Apw,B1},[2,2]);
    check_rows(C,2,'(1) Apw o B1');
    check_equal(apply_op(C,{v1,v2,v3}),rn(pw,'y',1)*rn(ap(b1,v3),'y',2),'(1) Apw o B1');
    C = mk({B1,Apw},[2,2]);
    check_equal(apply_op(C,{v3,v1,v2}),rn(ap(b1,v3),'y',1)*rn(pw,'y',2),'(1) B1 o Apw');
    ncases = ncases+2;
end

%% (2) Tensor product N = A1 o B1 (2 copies of y) and a pointwise product: N o Apw, Apw o N
a1 = lf_y();    a2 = lf_y();    a3 = lf_y();    b1 = lf_y();
v1 = rp();      v2 = rp();      v3 = rp();      v4 = rp();
Apw = mk({a1,a2},[1,2]);
N = mk({a3,b1},[2,2]);
pw = ap(a1,v1).*ap(a2,v2);
C = mk({N,Apw},[2,2]);
check_equal(apply_op(C,{v3,v4,v1,v2}),rn(ap(a3,v3),'y',1)*rn(ap(b1,v4),'y',2)*rn(pw,'y',3),'(2) N o Apw');
C = mk({Apw,N},[2,2]);
check_equal(apply_op(C,{v1,v2,v3,v4}),rn(pw,'y',1)*rn(ap(a3,v3),'y',2)*rn(ap(b1,v4),'y',3),'(2) Apw o N');
ncases = ncases+2;

%% (3) N o N' (2+2 copies), and N o N' o B1
a4 = lf_y();    b2 = lf_y();    v5 = rp();      v6 = rp();
N2 = mk({a4,b2},[2,2]);
C = mk({N,N2},[2,2]);
check_equal(apply_op(C,{v3,v4,v5,v6}),rn(ap(a3,v3),'y',1)*rn(ap(b1,v4),'y',2)*rn(ap(a4,v5),'y',3)*rn(ap(b2,v6),'y',4),'(3) N o N');
B1 = mk({lf_y()},[2,2]);
B1op = B1.ops{1,1};     v7 = rp();
C = mk({N,B1,N2},[2,2]);
check_equal(apply_op(C,{v3,v4,v7,v5,v6}),rn(ap(a3,v3),'y',1)*rn(ap(b1,v4),'y',2)*rn(ap(B1op,v7),'y',3)*rn(ap(a4,v5),'y',4)*rn(ap(b2,v6),'y',5),'(3) N o B1 o N');
ncases = ncases+2;

%% (4) Three or more levels of nesting
% (a) { {Apw o N} o B1 }: the inner product has 3 copies of y
Inner = mk({Apw,N},[2,2]);
B1op = lf_y();      B1 = mk({B1op},[2,2]);      v7 = rp();
C = mk({Inner,B1},[2,2]);
check_equal(apply_op(C,{v1,v2,v3,v4,v7}),rn(pw,'y',1)*rn(ap(a3,v3),'y',2)*rn(ap(b1,v4),'y',3)*rn(ap(B1op,v7),'y',4),'(4a) {Apw o N} o B1');
% (b) { B1 o {N.*N'} o Apw }: pointwise product of 2 tensor products, with 2 copies of y
Pw2 = mk({N,N2},[1,2]);
pw2 = (rn(ap(a3,v3),'y',2)*rn(ap(b1,v4),'y',3)).*(rn(ap(a4,v5),'y',2)*rn(ap(b2,v6),'y',3));
C = mk({B1,Pw2,Apw},[2,2]);
check_rows(C,3,'(4b) B1 o {N.*N} o Apw');
check_equal(apply_op(C,{v7,v3,v4,v5,v6,v1,v2}),rn(ap(B1op,v7),'y',1)*pw2*rn(pw,'y',4),'(4b) B1 o {N.*N} o Apw');
% (c) { { {N o Apw} o B1 } o Apw' } : depth four
Inner2 = mk({N,Apw},[2,2]);
Inner3 = mk({Inner2,B1},[2,2]);
a5 = lf_y();    a6 = lf_y();    v8 = rp();      v9 = rp();
Apw2 = mk({a5,a6},[1,2]);
C = mk({Inner3,Apw2},[2,2]);
pw3 = ap(a5,v8).*ap(a6,v9);
check_equal(apply_op(C,{v3,v4,v1,v2,v7,v8,v9}),rn(ap(a3,v3),'y',1)*rn(ap(b1,v4),'y',2)*rn(pw,'y',3)*rn(ap(B1op,v7),'y',4)*rn(pw3,'y',5),'(4c) depth four');
% (d) pointwise product of pointwise products as nested factor: one copy
Pw4 = mk({Apw,mk({a5,a6},[1,2])},[1,2]);
C = mk({Pw4,B1},[2,2]);
check_equal(apply_op(C,{v1,v2,v8,v9,v7}),rn(pw.*pw3,'y',1)*rn(ap(B1op,v7),'y',2),'(4d) pointwise of pointwise');
ncases = ncases+4;

%% (5) Several terms
% (a) nested pointwise product with two terms: B1 o {A1.*A2 + A3.*A4}
a3b = lf_y();   a4b = lf_y();
Pt = mk({a1,a2;a3b,a4b},[1,2]);
C = mk({B1,Pt},[2,2]);
pt = pw + ap(a3b,v1).*ap(a4b,v2);
check_equal(apply_op(C,{v7,v1,v2}),rn(ap(B1op,v7),'y',1)*rn(pt,'y',2),'(5a) nested operator with two terms');
% (b) two terms in the outer product, each with a nested factor: Apw o B1 + Apw' o B1'
b1b = lf_y();   Apwb = mk({a3b,a4b},[1,2]);     B1b = mk({b1b},[2,2]);
C = mk({Apw,B1;Apwb,B1b},[2,2]);
check_equal(apply_op(C,{v1,v2,v3}),rn(pw,'y',1)*rn(ap(B1op,v3),'y',2) + rn(ap(a3b,v1).*ap(a4b,v2),'y',1)*rn(ap(b1b,v3),'y',2),'(5b) two outer terms');
% (c) nested tensor product with two terms, which both have 2 copies: {A3 o B1 + A4 o B2} o Apw
N2t = mk({a3,b1;a4,b2},[2,2]);
C = mk({N2t,Apw},[2,2]);
Nt = rn(ap(a3,v3),'y',1)*rn(ap(b1,v4),'y',2) + rn(ap(a4,v3),'y',1)*rn(ap(b2,v4),'y',2);
check_equal(apply_op(C,{v3,v4,v1,v2}),Nt*rn(pw,'y',3),'(5c) nested tensor product with two terms');
% (d) both: two terms, each a nested tensor product with two terms
C = mk({N2t,Apw;N2t,Apwb},[2,2]);
check_equal(apply_op(C,{v3,v4,v1,v2}),Nt*rn(pw,'y',3) + Nt*rn(ap(a3b,v1).*ap(a4b,v2),'y',3),'(5d) two terms of nested operators with two terms');
ncases = ncases+4;

%% (6) Different spatial variables: the copies are numbered for each output variable
% L1 (x -> y), Pw (x -> y,y2: pointwise product of two such operators), M2 (x2 -> y2), N2 (x2 -> y2, 2 copies)
L1 = lf_y();                l1 = lf_yy2();      l2 = lf_yy2();
M2 = lf_y2();               m1 = lf_y2();       m2 = lf_y2();
Pw = mk({l1,l2},[1,2]);     N2y2 = mk({m1,m2},[2,2]);
va = rp();      vb = rp();      vc = rp();      vz = rp2();     vz1 = rp2();    vz2 = rp2();
plw = ap(l1,vb).*ap(l2,vc);
% y copies: L1: 1, Pw: 2; y2 copies: Pw: 1, M2: 2, N2: 3,4
C = mk({mk({L1},[2,2]),Pw,mk({M2},[2,2]),N2y2},[2,2]);
check_rows(C,4,'(6) different variables');
Ey = rn(ap(L1,va),'y',1)*rn(rn(plw,'y',2),'y2',1)*rn(ap(M2,vz),'y2',2)*rn(ap(m1,vz1),'y2',3)*rn(ap(m2,vz2),'y2',4);
check_equal(apply_op(C,{va,vb,vc,vz,vz1,vz2}),Ey,'(6) different variables');
% N2 first, then Pw, M2 and L1: y2 copies: N2: 1,2, Pw: 3, M2: 4; y copies: Pw: 1, L1: 2
C = mk({N2y2,Pw,mk({M2},[2,2]),mk({L1},[2,2])},[2,2]);
Ey = rn(ap(m1,vz1),'y2',1)*rn(ap(m2,vz2),'y2',2)*rn(rn(plw,'y',1),'y2',3)*rn(ap(M2,vz),'y2',4)*rn(ap(L1,va),'y',2);
check_equal(apply_op(C,{vz1,vz2,vb,vc,vz,va}),Ey,'(6) different variables, other order');
% a nested tensor product depending on both variables, with 2 copies of each, after a factor depending on y only
Nyy2 = mk({l1,l2},[2,2]);
C = mk({mk({L1},[2,2]),Nyy2,mk({M2},[2,2])},[2,2]);
Ey = rn(ap(L1,va),'y',1)*rn(rn(ap(l1,vb),'y',2),'y2',1)*rn(rn(ap(l2,vc),'y',3),'y2',2)*rn(ap(M2,vz),'y2',3);
check_equal(apply_op(C,{va,vb,vc,vz}),Ey,'(6) tensor product in two variables');
ncases = ncases+3;

%% (7) Inputs as the 1 x d cell, with the inputs of nested factors as cells
C = mk({N,Apw},[2,2]);
Yflat = apply_op(C,{v3,v4,v1,v2});
check_equal(apply_op(C,{{v3,v4},{v1,v2}}),Yflat,'(7) 1 x d form, N o Apw');
C = mk({Apw,mk({a3},[2,2]),N},[2,2]);
check_equal(apply_op(C,{{v1,v2},v3,{v3,v4}}),apply_op(C,{v1,v2,v3,v3,v4}),'(7) 1 x d form, Apw o B1 o N');
check_equal(apply_op(C,{v1,v2,v3,v3,v4}),rn(pw,'y',1)*rn(ap(a3,v3),'y',2)*rn(ap(a3,v3),'y',3)*rn(ap(b1,v4),'y',4),'(7) Apw o B1 o N');
ncases = ncases+3;

%% (8) Operators without nested factors: unchanged
C = mk({a1,b1,a2},[2,2]);
check_equal(apply_op(C,{v1,v3,v2}),rn(ap(a1,v1),'y',1)*rn(ap(b1,v3),'y',2)*rn(ap(a2,v2),'y',3),'(8) flat tensor product');
C = mk({a1,a2},[1,2]);
check_equal(apply_op(C,{v1,v2}),ap(a1,v1).*ap(a2,v2),'(8) flat pointwise product');
C = mk({a1},[2,2]);
check_equal(apply_op(C,v1),rn(ap(a1,v1),'y',1),'(8) single factor, otype [2,2]');
C = mk({a1},[1,2]);
check_equal(apply_op(C,v1),ap(a1,v1),'(8) single factor, otype [1,2]');
% nested factor that has a constant output: no variables to rename
Ccst = mk({mk({a1,a2},[1,2]),mk({a3},[2,2])},[2,2]);
Ycst = polynomial(apply_op(Ccst,{0*v1,0*v2,0*v3}));
if any(abs(nonzeros(Ycst.coefficient))>tol)
    error("Test_tensopvar_apply_op failed (8): zero inputs give a nonzero output.")
end
ncases = ncases+5;

%% (9) A nested pointwise product with both numbered and unnumbered copies cannot be placed
Pmix = mk({a1,N},[1,2]);
% on its own, the pointwise product is evaluated as before
check_equal(apply_op(Pmix,{v1,v3,v4}),ap(a1,v1).*(rn(ap(a3,v3),'y',1)*rn(ap(b1,v4),'y',2)),'(9) pointwise product with N on its own');
C = mk({Pmix,B1},[2,2]);
expect_error(@() apply_op(C,{v1,v3,v4,v7}),'numbered and unnumbered');
ncases = ncases+2;
end

function N = pointwise(A,B)
% Pointwise product of the outputs of the operators A and B
N = tensopvar();        N.otype = [1,2];        N.ops = {A,B};      N = initialize(N);
end


function N = kron_op(ops)
% Kronecker product of the outputs of the operators in each row of the cell ops, summed over the rows
N = tensopvar();        N.otype = [2,2];        N.ops = ops;        N = initialize(N);
end


function N = pointwise_terms(A1,B1,A2,B2)
% Sum of two pointwise products, a nested operator with two terms
N = tensopvar();        N.otype = [1,2];        N.ops = {A1,B1; A2,B2};     N = initialize(N);
end


function X = nest_leaves(R,leaves)
% The 1 x d cell of inputs to the first term of R, with the inputs of nested factors as cells
pos = 1;        d = size(R.ops,2);      X = cell(1,d);
for j=1:d
    F = R.ops{1,j};
    if isa(F,'tensopvar')
        n = num_factors(F);
        X{j} = nest_leaves(F,leaves(pos:pos+n-1));
        pos = pos+n;
    else
        X{j} = leaves{pos};
        pos = pos+1;
    end
end
end


function [Y,pos,nslots] = apply_nested(R,leaves,pos,off,vars_out)
% Apply the (nested) tensor-PI operator R to the states in leaves(pos:end), one for each of the 'sopvar'
% operators in a term of R, in order, each term acting on the same states. The output variable vars_out{k}
% is replaced by vars_out{k}_(off(k)+1), vars_out{k}_(off(k)+2),... for the leaf operators in a Kronecker
% product of the outputs, whereas the outputs of the factors in a pointwise product share the variable.
M = numel(vars_out);
if isa(R,'sopvar')
    Y = apply_sopvar(R,leaves{pos});
    pos = pos+1;
    nslots = zeros(1,M);
    for k=1:M
        if ismember(vars_out{k},R.vars.out)
            is_v = strcmp(Y.varname,vars_out{k});
            Y.varname(is_v) = {[vars_out{k},'_',num2str(off(k)+1)]};
            nslots(k) = 1;
        end
    end
    return
end
Y = 0;
pos0 = pos;
for i=1:size(R.ops,1)
    Yi = 1;     pos = pos0;     off_j = off;    nslots = zeros(1,M);
    for j=1:size(R.ops,2)
        [Yij,pos,nslots_j] = apply_nested(R.ops{i,j},leaves,pos,off_j,vars_out);
        if R.otype(1)==1
            Yi = Yi.*Yij;
            nslots = max(nslots,nslots_j);
        else
            Yi = kron(Yi,Yij);
            off_j = off_j+nslots_j;
            nslots = nslots+nslots_j;
        end
    end
    Y = Y+Yi;
end
end


function f = build_f(Zrow,facs,Id)
% polyopvar with monomial Z=Zrow in the states x,y, acted on by the product of the factors facs: the
% identity operators of the states for 'x','y', multipliers from R for polynomials
K = [];
for j=1:numel(facs)
    if ischar(facs{j})
        Fj = Id.(facs{j});
    else
        vars = cell(1,0);
        if ismember('s1',facs{j}.varname),  vars = [vars,{'s1'}];   end
        if ismember('s2',facs{j}.varname),  vars = [vars,{'s2'}];   end
        Fj = mat2tensopvar(facs{j},struct('in',{cell(1,0)},'out',{vars}),struct('in',zeros(0,2),'out',repmat([0,1],numel(vars),1)));
    end
    if isempty(K),  K = Fj;
    else,           K = otimes(K,Fj,[1,2]);     end
end
f = polyopvar();
f.varname = {'x';'y'};      f.Z = sparse(Zrow);     f.C = tensopmat(K);
f.dims.in = [1;1];
f.vars.in = {'s1','s2'};    f.dom.in = [0,1;0,1];   f.depmat.in = true(2,2);
f = initialize(f);
end


function Y = prod_facs(Zrow,facs,xv,yv)
% Product of the values of the factors, as polynomials: x -> xv, y -> yv, polynomial -> itself
Y = 1;
for j=1:numel(facs)
    if ischar(facs{j})
        if strcmp(facs{j},'x'),     Y = Y*xv;
        else,                       Y = Y*yv;       end
    else
        Y = Y*facs{j};
    end
end
end


function assert_close(Y,Yref,tol,label)
% Output variables of Kronecker products are named s_1,s_2,..., and compared as they are named
Y = polynomial(Y);      Yref = polynomial(Yref);
D = Y-Yref;
if ~isequal(size(D),size(Yref))
    error("Test_tensopvar_apply_op failed %s: size %s, expected %s.",label,mat2str(size(Y)),mat2str(size(Yref)))
end
err = full(max([0;abs(nonzeros(D.coefficient))]));
if err>tol
    error("Test_tensopvar_apply_op failed %s: residual %.3e.",label,err)
end
if isempty(nonzeros(Yref.coefficient))
    error("Test_tensopvar_apply_op %s: reference value is zero; test would be vacuous.",label)
end
end


function expect_error(fun,str)
try
    fun();
catch ME
    if ~contains(ME.message,str)
        error("Test_tensopvar_apply_op failed: unexpected error message '%s' (expected to contain '%s').",ME.message,str)
    end
    return
end
error("Test_tensopvar_apply_op failed: expected an error containing '%s', but none was thrown.",str)
end

function T = mk(ops,otype)
% Tensor-PI operator with the factors in the cell ops, with otype(1) specifying a pointwise (1) or
% tensor (2) product of the outputs; a nested factor is a 'tensopvar' in ops
T = tensopvar();        T.otype = otype;        T.ops = ops;        T = initialize(T);
end


function p = rn(p,nm,k)
% The output variable nm of the polynomial p, renamed to its k-th copy nm_k
p.varname(strcmp(p.varname,nm)) = {[nm,'_',num2str(k)]};
end


function check_rows(C,d,label)
% depmat has one row for each of the d factors of C, also if a factor is nested
if size(C.depmat.out,1)~=d || size(C.depmat.in,1)~=d
    error("Test_tensopvar_apply_op failed (%s): depmat should have %d rows.",label,d)
end
end


function check_equal(Y,E,label)
err = cleanpoly(polynomial(Y)-polynomial(E),1e-10);
cf = nonzeros(err.coefficient);
if ~isempty(cf) && max(abs(cf))>1e-8
    error("Test_tensopvar_apply_op failed (%s): residual %.3e.",label,max(abs(cf)))
end
% the oracle must not be zero, or the comparison is vacuous
Ep = cleanpoly(polynomial(E),1e-10);
if isempty(nonzeros(Ep.coefficient))
    error("Test_tensopvar_apply_op (%s): reference value is zero; test would be vacuous.",label)
end
end
