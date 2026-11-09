% Test_tensopvar_combine_factors.m -- tests combine_factors of a 'tensopvar',
% which merges the factors that are operators from R^n to R^m (no spatial
% variables) into a neighboring factor by a Kronecker product, and the use of
% combine_factors in tensopvar/mtimes when the numbers of factors differ.
%
% Verification is semantic and independent of the routine under test:
%   (1) one factor R from R and one factor B on L2[0,1], in both orders:
%       the merged operator acting on kron(x,y) vs. kron of the outputs of
%       R and B applied to x and y by apply_sopvar (the oracle), for
%       matrices R and operators B of several sizes;
%   (2) R from R on both sides of B, and several operators from R: the
%       outputs of the operator before and after merging (apply_op on the
%       original factors, and on the Kronecker products of the inputs), the
%       number of factors, and the matrix kron(R1,R2,R3) if all factors are
%       from R, compared with the product of the matrices;
%   (3) operators with two terms, and with two factors on L2[0,1] around R;
%   (4) operators that are returned unchanged: no factor from R, a single
%       factor, a pointwise product of the inputs (otype [2,1], [1,1]);
%   (5) mtimes: a single-factor operator times a multi-factor operator with
%       factors from R composes, as the composition of the merged operators
%       (apply_sopvar of the single operator on the output of the multi-factor
%       operator by the oracle in (1)), and an error is thrown if the numbers
%       of factors still differ;
%   (6) the same for a pointwise product of the outputs (otype [1,2]): the
%       merged operator vs. the elementwise product of the outputs of R and
%       B (R: R^n to R^m, B with m rows), in both orders, around B, with two
%       terms, for all factors from R (row-wise Kronecker product of the
%       matrices), and in mtimes.
%
% DJ (CC), 11/08/2026: Initial coding;
% DJ (CC), 11/08/2026: Add case (6), the pointwise product of the outputs;

rng(11);
tol = 1e-9;
pvar s
vars = struct('in',{{'s'}},'out',{{'s'}});
dom = struct('in',[0,1],'out',[0,1]);
degs = struct('in',2,'out',2);
ncases = 0;

%% (1) R and B, both orders, several sizes
for sz = [1 1 1 1; 2 3 1 1; 2 1 2 3; 3 2 2 2; 1 2 2 1]'      % [mR,nR,pB,qB]
    mR = sz(1);     nR = sz(2);     pB = sz(3);     qB = sz(4);
    Rop = mat2tensopvar(randn(mR,nR)).ops{1,1};
    Bop = rand_sopvar([pB,qB],vars,dom,degs,0.6);
    xr = randn(nR,1);
    xb = rand_poly_vec(qB,s);
    Yr = apply_sopvar(Rop,polynomial(xr));
    Yb = apply_sopvar(Bop,xb);
    % R o B
    C = make_op({Rop,Bop});
    Cc = combine_factors(C);
    assert(size(Cc.ops,2)==1,'(1) R o B: expected one factor.')
    assert_close(apply_op(Cc,kron_vec(xr,xb)),rn(kron_vec(Yr,Yb)),tol,sprintf('(1) R o B, size [%d %d %d %d]',sz));
    % B o R
    C = make_op({Bop,Rop});
    Cc = combine_factors(C);
    assert(size(Cc.ops,2)==1,'(1) B o R: expected one factor.')
    assert_close(apply_op(Cc,kron_vec(xb,xr)),rn(kron_vec(Yb,Yr)),tol,sprintf('(1) B o R, size [%d %d %d %d]',sz));
    ncases = ncases+2;
end

%% (2) Several operators from R, around B
R1 = mat2tensopvar(randn(2,3)).ops{1,1};
R2 = mat2tensopvar(randn(1,2)).ops{1,1};
R3 = mat2tensopvar(randn(2,2)).ops{1,1};
Bop = rand_sopvar([2,1],vars,dom,degs,0.6);
x1 = randn(3,1);    x2 = randn(2,1);    x3 = randn(2,1);    xb = rand_poly_vec(1,s);
% R1 o B o R3
C = make_op({R1,Bop,R3});
Cc = combine_factors(C);
assert(size(Cc.ops,2)==1,'(2) R1 o B o R3: expected one factor.')
Yref = kron_vec(kron_vec(apply_sopvar(R1,polynomial(x1)),apply_sopvar(Bop,xb)),apply_sopvar(R3,polynomial(x3)));
assert_close(apply_op(Cc,kron_vec(kron_vec(x1,xb),x3)),rn(Yref),tol,'(2) R1 o B o R3');
% R1 o R2 o B o R3 o R2
C = make_op({R1,R2,Bop,R3,R2});
Cc = combine_factors(C);
assert(size(Cc.ops,2)==1,'(2) five factors: expected one factor.')
Yref = kron_vec(kron_vec(kron_vec(kron_vec(apply_sopvar(R1,polynomial(x1)),apply_sopvar(R2,polynomial(x2))), ...
        apply_sopvar(Bop,xb)),apply_sopvar(R3,polynomial(x3))),apply_sopvar(R2,polynomial(x2)));
assert_close(apply_op(Cc,kron_vec(kron_vec(kron_vec(kron_vec(x1,x2),xb),x3),x2)),rn(Yref),tol,'(2) R1 o R2 o B o R3 o R2');
% All factors from R: a single matrix operator, the Kronecker product of the matrices
M1 = randn(2,3);    M2 = randn(1,2);    M3 = randn(2,2);
C = make_op({mat2tensopvar(M1).ops{1,1},mat2tensopvar(M2).ops{1,1},mat2tensopvar(M3).ops{1,1}});
Cc = combine_factors(C);
assert(size(Cc.ops,2)==1,'(2) all from R: expected one factor.')
assert(isequal(size(Cc.ops{1,1}.vars.in),size(cell(1,0))) || isempty(Cc.ops{1,1}.vars.in),'(2) all from R: expected no variables.')
assert(norm(full(Cc.ops{1,1}.params{1})-kron(kron(M1,M2),M3),'fro')<tol,'(2) all from R: wrong matrix.')
assert(isequal(Cc.ops{1,1}.dims,[4,12]),'(2) all from R: wrong dimensions.')
ncases = ncases+4;

%% (3) Two terms, and two factors on L2[0,1] around R
Ra = mat2tensopvar(randn(2,2)).ops{1,1};    Rb = mat2tensopvar(randn(2,2)).ops{1,1};
Ba = rand_sopvar([1,1],vars,dom,degs,0.6);  Bb = rand_sopvar([1,1],vars,dom,degs,0.6);
xr = randn(2,1);    xb = rand_poly_vec(1,s);
C = make_op({Ra,Ba;Rb,Bb});
Cc = combine_factors(C);
assert(isequal(size(Cc.ops),[2,1]),'(3) two terms: expected 2 x 1 operators.')
Yref = kron_vec(apply_sopvar(Ra,polynomial(xr)),apply_sopvar(Ba,xb)) + kron_vec(apply_sopvar(Rb,polynomial(xr)),apply_sopvar(Bb,xb));
assert_close(apply_op(Cc,kron_vec(xr,xb)),rn(Yref),tol,'(3) two terms, R o B');
% Two factors on L2[0,1] around R: R merges into the second, compare apply_op before and after
B1 = rand_sopvar([1,1],vars,dom,degs,0.6);  B2 = rand_sopvar([1,1],vars,dom,degs,0.6);
xa = rand_poly_vec(1,s);    xc = rand_poly_vec(1,s);
C = make_op({B1,Ra,B2});
Cc = combine_factors(C);
assert(size(Cc.ops,2)==2,'(3) B1 o R o B2: expected two factors.')
assert_close(apply_op(Cc,{xa,kron_vec(xr,xc)}),apply_op(C,{xa,xr,xc}),tol,'(3) B1 o R o B2');
ncases = ncases+2;

%% (4) Operators that are returned unchanged
C = make_op({B1,B2});
assert(isequal(combine_factors(C),C),'(4) no factor from R: expected unchanged.')
C = make_op({Ra});
assert(isequal(combine_factors(C),C),'(4) single factor: expected unchanged.')
% pointwise product of the inputs: (A z).*(B z) is not linear in z, so that the factors cannot be merged
for ot = {[2,1],[1,1]}
    C = make_op({mat2tensopvar(randn(1,1)).ops{1,1},B1},ot{1});
    Cc = combine_factors(C);
    assert(isequal(Cc.otype,ot{1}) && size(Cc.ops,2)==2,'(4) otype %s: expected unchanged.',mat2str(ot{1}))
end
ncases = ncases+4;

%% (6) Pointwise product of the outputs, otype [1,2]: row r of the merged operator is the
% Kronecker product of the rows r of the operators (same numbers of rows)
for sz = [1 3 1; 2 3 1; 2 1 2; 3 2 2; 2 2 3]'       % [m,nR,qB]: R: R^nR to R^m, B: L2^qB to L2^m
    m = sz(1);     nR = sz(2);     qB = sz(3);
    Rop = mat2tensopvar(randn(m,nR)).ops{1,1};
    Bop = rand_sopvar([m,qB],vars,dom,degs,0.6);
    xr = randn(nR,1);       xb = rand_poly_vec(qB,s);
    Yr = apply_sopvar(Rop,polynomial(xr));      Yb = apply_sopvar(Bop,xb);
    Cc = combine_factors(make_op({Rop,Bop},[1,2]));
    assert(size(Cc.ops,2)==1 && isequal(Cc.otype,[1,2]),'(6) R o B: expected one factor of type [1,2].')
    assert_close(apply_op(Cc,kron_vec(xr,xb)),Yr.*Yb,tol,sprintf('(6) R o B, pointwise, size [%d %d %d]',sz));
    Cc = combine_factors(make_op({Bop,Rop},[1,2]));
    assert(size(Cc.ops,2)==1,'(6) B o R: expected one factor.')
    assert_close(apply_op(Cc,kron_vec(xb,xr)),Yb.*Yr,tol,sprintf('(6) B o R, pointwise, size [%d %d %d]',sz));
    ncases = ncases+2;
end
% R1, B, R2 with two rows, and two terms
R1 = mat2tensopvar(randn(2,3)).ops{1,1};    R2 = mat2tensopvar(randn(2,2)).ops{1,1};
Bop = rand_sopvar([2,2],vars,dom,degs,0.6);
x1 = randn(3,1);    x2 = randn(2,1);    xb = rand_poly_vec(2,s);
Cc = combine_factors(make_op({R1,Bop,R2},[1,2]));
assert(size(Cc.ops,2)==1,'(6) R1 o B o R2: expected one factor.')
Yref = apply_sopvar(R1,polynomial(x1)).*apply_sopvar(Bop,xb).*apply_sopvar(R2,polynomial(x2));
assert_close(apply_op(Cc,kron_vec(kron_vec(x1,xb),x2)),Yref,tol,'(6) R1 o B o R2, pointwise');
Ra2 = mat2tensopvar(randn(2,3)).ops{1,1};   Bb2 = rand_sopvar([2,2],vars,dom,degs,0.6);
Cc = combine_factors(make_op({R1,Bop;Ra2,Bb2},[1,2]));
Yref = apply_sopvar(R1,polynomial(x1)).*apply_sopvar(Bop,xb) + apply_sopvar(Ra2,polynomial(x1)).*apply_sopvar(Bb2,xb);
assert(isequal(size(Cc.ops),[2,1]),'(6) two terms: expected 2 x 1 operators.')
assert_close(apply_op(Cc,kron_vec(x1,xb)),Yref,tol,'(6) two terms, pointwise');
% All factors from R: the row-wise Kronecker product of the matrices
M1 = randn(2,3);    M2 = randn(2,2);
Cc = combine_factors(make_op({mat2tensopvar(M1).ops{1,1},mat2tensopvar(M2).ops{1,1}},[1,2]));
Mref = [kron(M1(1,:),M2(1,:)); kron(M1(2,:),M2(2,:))];
assert(size(Cc.ops,2)==1 && norm(full(Cc.ops{1,1}.params{1})-Mref,'fro')<tol,'(6) all from R: wrong matrix.')
% mtimes: a single-factor operator times R o B with pointwise outputs
Zop = make_op({rand_sopvar([2,2],vars,dom,degs,0.6)},[1,2]);
Rop = mat2tensopvar(randn(2,3)).ops{1,1};   Bop = rand_sopvar([2,2],vars,dom,degs,0.6);
xr = randn(3,1);    xb = rand_poly_vec(2,s);
Cm = Zop*make_op({Rop,Bop},[1,2]);
assert(size(Cm.ops,2)==1,'(6) mtimes: expected one factor after composition.')
Yref = apply_sopvar(Zop.ops{1,1},apply_sopvar(Rop,polynomial(xr)).*apply_sopvar(Bop,xb));
assert_close(apply_op(Cm,kron_vec(xr,xb)),Yref,tol,'(6) mtimes, pointwise');
ncases = ncases+5;

%% (5) mtimes with a different number of factors
% Single-factor Z on L2[0,1], times B1 o R o ... : combine, then compose
Zop = make_op({rand_sopvar([2,2],vars,dom,degs,0.6)});     % acts on R^2 o L2 = L2^2: kron input of dimension 2
Bop = rand_sopvar([1,1],vars,dom,degs,0.6);
Rc = mat2tensopvar(eye(2)).ops{1,1};
C = make_op({Rc,Bop});            % identity on R^2 o B: maps L2[0,1] o R^2 to L2^2
xr = randn(2,1);    xb = rand_poly_vec(1,s);
Cm = Zop*C;
assert(size(Cm.ops,2)==1,'(5) expected one factor after composition.')
Yref = apply_sopvar(Zop.ops{1,1},kron_vec(apply_sopvar(Rc,polynomial(xr)),apply_sopvar(Bop,xb)));
assert_close(apply_op(Cm,kron_vec(xr,xb)),rn(Yref),tol,'(5) single factor times R o B');
% Still different numbers of factors: error
B3 = rand_sopvar([1,1],vars,dom,degs,0.6);
Cthree = make_op({B1,B2,B3});
Zone = make_op({rand_sopvar([1,1],vars,dom,degs,0.6)});
expect_error(@() Zone*Cthree,'number of factors');
ncases = ncases+2;

fprintf('Test_tensopvar_combine_factors passed (%d cases).\n',ncases)


function C = make_op(ops,otype)
% Tensor-PI operator with the factors ops (K x d cell), with proper tensor products,
% or of the type otype
if nargin<2,    otype = [2,2];  end
C = tensopvar();
C.ops = ops;
C.otype = otype;
C = initialize(C);
end

function x = rand_poly_vec(n,s)
% n x 1 vector of random polynomials in s
x = polynomial(zeros(n,1));
for i=1:n
    x(i) = randn + randn*s + randn*s^2;
end
end

function z = kron_vec(x,y)
% Kronecker product kron(x,y) of column vectors x and y, either of which may be
% a 'polynomial' array: x(1)*y; x(2)*y; ...
z = [];
for i=1:numel(x)
    xi = x(i);
    z = [z; xi*y];
end
end

function p = rn(p)
% The output variable s of the polynomial p, renamed to its first copy s_1, as apply_op names the
% output variables of Kronecker products
p = polynomial(p);
p.varname(strcmp(p.varname,'s')) = {'s_1'};
end

function assert_close(Y1,Y2,tol,label)
% Assert that the polynomial arrays Y1, Y2 are equal up to tol
D = Y1-Y2;
err = full(max([0;abs(nonzeros(D.coefficient))]));
assert(~isempty(nonzeros(polynomial(Y2).coefficient)),'%s: reference value is zero, test is vacuous.',label);
assert(size(Y1,1)==size(Y2,1) && err<=tol*(1+max(abs(full(Y2.coefficient(:))))),'%s: mismatch (size %s vs %s, error %g).',label,mat2str(size(Y1)),mat2str(size(Y2)),err);
end

function expect_error(f,needle)
% Assert that f throws an error with the message containing needle
try
    f();
catch err
    assert(contains(err.message,needle),'Wrong error message: "%s" (expected "%s").',err.message,needle);
    return
end
error('Expected an error containing "%s", none thrown.',needle);
end



