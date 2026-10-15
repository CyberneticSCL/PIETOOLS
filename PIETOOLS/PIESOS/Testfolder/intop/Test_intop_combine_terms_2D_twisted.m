% Test_intop_combine_terms_2D_twisted.m -- Tests @intop/combine_terms.m
% (and @intop/reorder_states.m, @intop/reorder_vars.m, which it calls) with
% a "twisted" nesting order: the two spatial directions s1, s2 integrate
% the SAME pair of repeated factor-slots in OPPOSITE order (slot1 inner/
% slot2 outer along s1, but slot1 outer/slot2 inner along s2). Parts 4/5
% are multi-group cases (two distinct nesting orders along s1, so
% combine_terms splits into two groups, canonicalizes each with a
% different permutation, then sums). Oracle: apply_op of the functional
% before vs. after combine_terms on x(s1,s2) o x(s1,s2) (o x(s1,s2) for
% degree 3), exact polynomial arithmetic.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - Every copy of the state must be the same base formula, renamed by
%   subs() to its own dummy variables (s1_1/s2_1, s1_2/s2_2, ...):
%   apply_op matches dummy variables to slots by NAME, not by kron()
%   position, and combine_terms assumes repeated slots hold the same
%   function (one state to a tensor power). Independent functions per slot
%   violate that precondition.
% - Value check only; the canonical form is tested in
%   Test_intop_combine_terms_canonical.m.

tol = 1e-8;
vars = struct('in',{{'s1','s2'}},'out',{cell(1,0)});
dom = struct('in',[0,1;0,1],'out',zeros(0,2));
depmat = true(1,2);   % the (single) state depends on both s1 and s2

%% Part 1: single-group, SCALAR state, degree 2, twisted order
I1 = [0,2, 0,3];   % s1: slot1 inner [0,dum2], slot2 outer [0,b]
I2 = [0,3, 0,1];   % s2: slot1 outer [0,b], slot2 inner [0,dum1]
intmat = {I1,I2};

pvar s1_dum1 s1_dum2 s2_dum1 s2_dum2
K = s1_dum1*s2_dum1^2 + 5*s1_dum2*s2_dum2^3 - 2*s1_dum2^2*s2_dum1;
Kop = intop({K},intmat,vars,dom,false);
Kop2 = combine_terms(Kop,2,1,depmat);

pvar s t s1_1 s2_1 s1_2 s2_2
vx = s^2*t + 3*t^3 - s*t^2 + 1;
xval = subs(vx,[s;t],[s1_1;s2_1]) * subs(vx,[s;t],[s1_2;s2_2]);

Y0 = apply_op(Kop,xval);
Y1 = apply_op(Kop2,xval);
err = cleanpoly(Y0-Y1,1e-10);
if ~isdouble(err) || abs(double(err))>tol
    error("Test_intop_combine_terms_2D_twisted failed: Part 1 (scalar, single-group).")
end

%% Part 2: single-group, VECTOR state (dims=2), degree 2, twisted order
pvar s1_dum1 s1_dum2 s2_dum1 s2_dum2
Kvec = [s1_dum1*s2_dum1^2, 5*s1_dum2*s2_dum2^3, -2*s1_dum2^2*s2_dum1, s1_dum1+s2_dum2];
KopB = intop({Kvec},intmat,vars,dom,false);
KopB2 = combine_terms(KopB,2,2,depmat);

vxvec = [s^2*t+1; 3*t^3-s*t^2+2];
v1vec = subs(vxvec,[s;t],[s1_1;s2_1]);
v2vec = subs(vxvec,[s;t],[s1_2;s2_2]);
xvalB = kron(v1vec,v2vec);

Y0B = apply_op(KopB,xvalB);
Y1B = apply_op(KopB2,xvalB);
errB = cleanpoly(Y0B-Y1B,1e-10);
if ~isdouble(errB) || abs(double(errB))>tol
    error("Test_intop_combine_terms_2D_twisted failed: Part 2 (vector, single-group).")
end

%% Part 3: single-group, SCALAR state, degree 3, twisted order
pvar s1_dum1 s1_dum2 s1_dum3 s2_dum1 s2_dum2 s2_dum3
I1c = [0,3, 0,1, 0,4];   % s1: slot2 inner[0,slot1], slot1 mid[0,slot3], slot3 outer
I2c = [0,4, 0,1, 0,2];   % s2: slot1 outer, slot2 mid[0,slot1], slot3 inner[0,slot2]
intmatC = {I1c,I2c};
Kc = s1_dum1*s2_dum1 + 2*s1_dum2*s2_dum2^2 - 3*s1_dum3^2*s2_dum3 + s1_dum1*s1_dum2*s2_dum3;
KopC = intop({Kc},intmatC,vars,dom,false);
KopC2 = combine_terms(KopC,3,1,depmat);

pvar s1_3 s2_3
vxC1 = subs(vx,[s;t],[s1_1;s2_1]);
vxC2 = subs(vx,[s;t],[s1_2;s2_2]);
vxC3 = subs(vx,[s;t],[s1_3;s2_3]);
xvalC = vxC1*vxC2*vxC3;

Y0C = apply_op(KopC,xvalC);
Y1C = apply_op(KopC2,xvalC);
errC = cleanpoly(Y0C-Y1C,1e-10);
if ~isdouble(errC) || abs(double(errC))>tol
    error("Test_intop_combine_terms_2D_twisted failed: Part 3 (scalar, degree 3, single-group).")
end

%% Part 4: MULTI-group, SCALAR state -- the key stress test. Two terms
% along s1 with OPPOSITE nesting order of the same factor-slot pair
% (forcing combine_terms to split into 2 groups, apply a DIFFERENT
% canonicalizing permutation to each, then sum), plus a twisted order
% along s2.
I1m = [0,2,0,3;  0,3,0,1];   % s1: term1 slot1-inner/slot2-outer, term2 REVERSED
I2m = [0,3,0,1];             % s2: slot1 outer, slot2 inner (one term)
intmatM = {I1m,I2m};
pvar s1_dum1 s1_dum2 s2_dum1 s2_dum2
K1m = s1_dum1*s2_dum1^2 + 5*s1_dum2*s2_dum2^3 - 2*s1_dum2^2*s2_dum1;
K2m = 3*s1_dum1^2*s2_dum2 - s1_dum2*s2_dum1^2 + 7;
KopM = intop({K1m;K2m},intmatM,vars,dom,false);
KopM2 = combine_terms(KopM,2,1,depmat);
% Confirm this genuinely exercised the multi-group path (not a no-op):
if size(KopM.intmat{1},1)<2
    error("Test_intop_combine_terms_2D_twisted failed: Part 4 setup assumption violated (need 2+ terms along s1).")
end

xvalM = subs(vx,[s;t],[s1_1;s2_1]) * subs(vx,[s;t],[s1_2;s2_2]);
Y0M = apply_op(KopM,xvalM);
Y1M = apply_op(KopM2,xvalM);
errM = cleanpoly(Y0M-Y1M,1e-10);
if ~isdouble(errM) || abs(double(errM))>tol
    error("Test_intop_combine_terms_2D_twisted failed: Part 4 (scalar, MULTI-group, twisted).")
end

%% Part 5: MULTI-group, VECTOR state (dims=2), same twisted multi-group setup.
pvar s1_dum1 s1_dum2 s2_dum1 s2_dum2
K1mv = [s1_dum1*s2_dum1^2, 5*s1_dum2*s2_dum2^3, -2*s1_dum2^2*s2_dum1, s1_dum1+s2_dum2];
K2mv = [3*s1_dum1^2*s2_dum2, -s1_dum2*s2_dum1^2+7, s1_dum1-s1_dum2, s2_dum1*s2_dum2];
KopMv = intop({K1mv;K2mv},intmatM,vars,dom,false);
KopMv2 = combine_terms(KopMv,2,2,depmat);

xvalMv = kron(subs(vxvec,[s;t],[s1_1;s2_1]), subs(vxvec,[s;t],[s1_2;s2_2]));
Y0Mv = apply_op(KopMv,xvalMv);
Y1Mv = apply_op(KopMv2,xvalMv);
errMv = cleanpoly(Y0Mv-Y1Mv,1e-10);
if ~isdouble(errMv) || abs(double(errMv))>tol
    error("Test_intop_combine_terms_2D_twisted failed: Part 5 (vector, MULTI-group, twisted).")
end

disp('Test_intop_combine_terms_2D_twisted passed.')
