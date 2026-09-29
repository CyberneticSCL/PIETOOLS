function test_tier1_fixes_sop()
% TEST_TIER1_FIXES_SOP checks the three defects the Tier 1 review
% reproduced (09/29/2026), each against semantics, never a round trip:
%
%   (1) gam*P for a block P holding the scalar-0 zero-cell shorthand built
%       an sdopvar with 1 x 1 A and nd x 1 B cells ('dpvar_op_copvar',
%       scale_block), which getsol and composition then refused;
%   (2) getsol refused an all-zero B cell that is too narrow, as
%       sopvar2sdopvar produces from a scalar-0 cell ('getsol_lpivar_sop');
%   (3) an 'opvar' holding dpvar-valued fields (legacy [dpvar, opvar])
%       came back from lpigetsol_sop / subs_dvar_sop UNSUBSTITUTED, and a
%       nested cell {{cdopvar}} errored ('lpigetsol_sop').
%
% (1)-(2): the substituted operator is applied to polynomial test functions
% by quadrature from the class definition ('heatNd_apply', no class method
% in the loop) and compared with the scalar times the fixed block, applied
% the same way. (3): the substituted opvar is compared with the opvar built
% by legacy concatenation of the substituted VALUES (a different path), and
% a legacy opvar without dpvar fields must come back bit-identical.
%
% Initial coding MMP, 09/29/2026

warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
nchk = 0;

% ------------------------------------------ a block with a scalar-0 cell
rng(1);
b0 = rand_copvar(struct('out',{{{'s'}}},'in',{{{'s'}}}),struct('out',2,'in',2),[0 1],1,0.9);
bb = b0.C{1,1};     prm = bb.params;    prm{2} = 0;
bz = sopvar(prm,bb.vars,bb.ZL,bb.ZR,bb.dom,bb.dims);
assert(isscalar(bz.params{2}) && nnz(bz.params{1}) && nnz(bz.params{3}),...
    'setup: the constructor must keep the scalar-0 shorthand');
gam = dpvar('gam');     g0 = 0.5;
f = @(X) [1+X(:,1), X(:,1).^2];     S = linspace(0.05,0.95,7)';
yb = heatNd_apply(bz,f,S);

% (1) scalar forms: subs(form, gam = g0) = c*bz by action
forms = {@() gam*bz, g0; @() bz*gam, g0; @() -gam*bz, -g0; @() (gam*bz)', NaN};
for k = 1:size(forms,1)
    X = forms{k,1}();
    Y = subs_dvar_sop(X,{'gam'},g0);
    assert(isa(Y,'sopvar'),'(1) form %d: %s after substitution',k,class(Y));
    if isnan(forms{k,2})    % adjoint: compare with g0*bf' applied, bf = bz
        % with the zero cell stored at full size (@sopvar/ctranspose errors
        % on the scalar-0 shorthand itself, a pre-existing defect)
        prf = bb.params;    prf{2} = sparse(size(prf{2},1),size(prf{2},2));
        bf0 = sopvar(prf,bb.vars,bb.ZL,bb.ZR,bb.dom,bb.dims);
        yr = g0*heatNd_apply(bf0',f,S);
    else
        yr = forms{k,2}*yb;
    end
    e = norm(heatNd_apply(Y,f,S)-yr)/norm(yr);
    assert(e<1e-12,'(1) form %d: relative error %.2e',k,e);
    nchk = nchk+1;
end
% composition with the scaled block: subs((gam*bz)*bf) = g0*(bf*bf), bf
% the full-size copy of bz (the class compositions error on a shorthand
% OPERAND, pre-existing; here only the scaled factor is under test)
prf = bb.params;    prf{2} = sparse(size(prf{2},1),size(prf{2},2));
bf = sopvar(prf,bb.vars,bb.ZL,bb.ZR,bb.dom,bb.dims);
assert(norm(heatNd_apply(bf,f,S)-yb)<1e-14*norm(yb),'setup: bf is not bz');
Y = subs_dvar_sop((gam*bz)*bf,{'gam'},g0);
yr = g0*heatNd_apply(bf*bf,f,S);
e = norm(heatNd_apply(Y,f,S)-yr)/norm(yr);
assert(e<1e-12,'(1) (gam*bz)*bz: relative error %.2e',e);
nchk = nchk+1;

% (2) sopvar2sdopvar of the shorthand block, substituted, is the block
for a0 = [0, 2.7]
    Y = subs_dvar_sop(sopvar2sdopvar(bz,{'a'}),{'a'},a0);
    e = norm(heatNd_apply(Y,f,S)-yb)/norm(yb);
    assert(isa(Y,'sopvar') && e<1e-14,'(2) a = %g: relative error %.2e',a0,e);
    nchk = nchk+1;
end
% a nonzero mis-sized cell is still refused (not silently zeroed)
Xb = gam*bb;        prmb = Xb.params;   % B = b*vec(C)', nonzero
assert(nnz(prmb.B{1})>0,'setup: the control needs a nonzero B');
prmb.B{1} = prmb.B{1}(:,1:end-1);
Xb.params = prmb;
e = get_err(@() subs_dvar_sop(Xb,{'gam'},1));
assert(strcmp(e.identifier,'getsol_lpivar_sop:badB'),'(2) mis-sized nonzero B: %s',e.identifier);
nchk = nchk+1;

% --------------------------------------- (3) legacy opvar with dpvar fields
pvar s s_dum
T5 = rand_opvar([0 1; 2 2],1,s,s_dum,[0 1]);
h = dpvar('h');
L = [[gam; 1-h], T5];
assert(isa(L,'opvar') && isa(L.Q2,'dpvar'),'setup: legacy [dpvar, opvar] should hold a dpvar Q2');
Sx = subs_dvar_sop(L,{'gam';'h'},[0.7;-1.3]);
Ref = [[0.7; 2.3], T5];
fl = {'P','Q1','Q2'};   fr = {'R0','R1','R2'};
for k = 1:3
    assert(~isa(Sx.(fl{k}),'dpvar') && ~isa(Sx.R.(fr{k}),'dpvar'),'(3) a dpvar field is left');
end
D = Sx - Ref;
m = max([maxc(D.P),maxc(D.Q1),maxc(D.Q2),maxc(D.R.R0),maxc(D.R.R1),maxc(D.R.R2)]);
assert(m<1e-14,'(3) substituted opvar differs from the legacy value by %.2e',m);
prog = struct('decvartable',{{'gam';'h'}},'solinfo',struct('RRx',[0.7;-1.3],'x',[],'info',1));
% (isequal(T5,T5) is false for a legacy opvar - its R struct fails isequal
% against itself - so the legacy route is checked by an exact difference)
Sa = lpigetsol_sop(prog,T5);    Sb = lpigetsol(prog,T5);
for X = {Sa - T5, Sb - T5}
    Dx = X{1};
    m = max([maxc(Dx.P),maxc(Dx.Q1),maxc(Dx.Q2),maxc(Dx.R.R0),maxc(Dx.R.R1),maxc(Dx.R.R2)]);
    assert(isa(Sa,'opvar') && m==0,'(3) opvar without dpvar fields is not returned as is (%.2e)',m);
end
nchk = nchk+3;

% nested cells of container decision operators
prog = lpiprogram_sop({'s'},[0 1]);
[prog,Q] = lpivar_cdopvar(prog,1,{{'s'}},[0 1],1);
prog.solinfo = struct('RRx',randn(numel(prog.decvartable),1),'x',[],'info',1);
Qs = lpigetsol_sop(prog,Q);
r = lpigetsol_sop(prog,{{Q},Q});
assert(isa(r{1}{1},'copvar') && isequal(r{1}{1},Qs) && isequal(r{2},Qs),'(3) nested cell {{Q},Q}');
nchk = nchk+1;

fprintf('test_tier1_fixes_sop passed (%d checks).\n',nchk);
end


function e = get_err(f)
e = struct('identifier','','message','');
try,        f();
catch ME,   e = ME;
end
end

function m = maxc(X)
% Largest absolute coefficient of a double or polynomial.
if isnumeric(X)
    m = full(max(abs(X(:))));
else
    cf = polynomial(X).coefficient;     m = full(max(abs(cf(:))));
end
if isempty(m),  m = 0;  end
end
