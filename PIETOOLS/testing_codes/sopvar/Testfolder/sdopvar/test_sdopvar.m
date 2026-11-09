clc; clear;
% TEST_SDOPVAR exercises sdopvar operations using random sdopvar objects.

rng(1);
tol = 1e-10;

%% Random sdopvar fixtures
vars = struct('in',{{'s1','s2','s_3'}},'out',{{'s1','s2','s_3'}});
dom = struct('in',[0 1; -1 1; 0 1],'out',[0 1; -1 1; 0 1]);
deg = 2;
ndec = 6;
density = 0.35;

P = rand_sdopvar([2,2],vars,dom,deg,ndec,density);
Q = rand_sdopvar([2,2],vars,dom,deg,ndec,density);

nL = P.dims(1)*prod(cellfun(@numel,P.ZL));
nR = P.dims(2)*prod(cellfun(@numel,P.ZR));
nC = nL*nR;
for i = 1:numel(P.params.A)
    if ~isequal(size(P.params.A{i}),[nC,1])
        error('rand_sdopvar test failed: A{%d} has the wrong size.',i);
    end
    if ~isequal(size(P.params.B{i}),[ndec,nC])
        error('rand_sdopvar test failed: B{%d} has the wrong size.',i);
    end
end
disp('rand_sdopvar storage test passed.');

%% Addition and equality
if ~eq(P+Q,Q+P,tol)
    error('sdopvar addition test failed: P+Q differs from Q+P.');
end
if ~eq(P+P,2*P,tol)
    error('sdopvar scalar/addition test failed: P+P differs from 2*P.');
end
Z = 0*P;
if ~eq(P+Z,P,tol)
    error('sdopvar addition identity test failed: P+0 differs from P.');
end
disp('sdopvar addition tests passed.');

%% Equality with different decision bases
P_dec = P;
Zd_perm = P.Zd.';
P_dec.Zd = [Zd_perm(end); Zd_perm(1:end-1)];
for i = 1:numel(P_dec.params.B)
    P_dec.params.B{i} = P_dec.params.B{i}([end,1:end-1],:);
end
if ~eq(P,P_dec,tol)
    error('sdopvar eq test failed: reordered decision basis compares unequal.');
end

P_dec_extra = P;
P_dec_extra.Zd = [string(P_dec_extra.Zd(:)); "unused_decision"];
for i = 1:numel(P_dec_extra.params.B)
    P_dec_extra.params.B{i} = [P_dec_extra.params.B{i}; sparse(1,size(P_dec_extra.params.B{i},2))];
end
if ~eq(P,P_dec_extra,tol)
    error('sdopvar eq test failed: extra zero decision row changes equality.');
end
disp('sdopvar decision-basis eq tests passed.');

%% Equality with different monomial bases
vars_eq = struct('in',{{'s'}},'out',{{'s'}});
dom_eq = struct('in',[0 1],'out',[0 1]);
params_eq_1 = struct('A',{{sparse([1;2])}},'B',{{sparse([3 4])}});
params_eq_2 = struct('A',{{sparse([1;2;0])}},'B',{{sparse([3 4 0])}});
P_eq = sdopvar(params_eq_1,vars_eq,{'d1'},{[0;1]},{0},dom_eq,[1 1]);
Q_eq = sdopvar(params_eq_2,vars_eq,{'d1'},{[0;1;2]},{0},dom_eq,[1 1]);
if ~eq(P_eq,Q_eq,tol)
    error('sdopvar eq test failed: equal sdopvars with different monomial bases.');
end
if eq(P_eq,2*Q_eq,tol)
    error('sdopvar eq test failed: different sdopvars.');
end
disp('sdopvar monomial-basis eq tests passed.');

%% matrix multiplication
if ~eq(eye(P.dims(1))*P,P,tol)
    error('sdopvar mtimes test failed: left identity.');
end
if ~eq(P*eye(P.dims(2)),P,tol)
    error('sdopvar mtimes test failed: right identity.');
end
L = [1 2; -1 0.5];
R = [2 -0.25; 0 1];
if ~eq((L*P)*R,L*(P*R),tol)
    error('sdopvar mtimes test failed: left/right multiplication is not associative.');
end
disp('sdopvar numeric multiplication tests passed.');

%% Transpose
Pt = P';
Ptt = Pt';
if ~eq(Ptt,P,tol)
    error('sdopvar transpose test failed: double transpose differs from original sdopvar.');
end
disp('sdopvar transpose tests passed.');

disp('All sdopvar tests passed.');

