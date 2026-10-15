% Test_mtimes_id_bypass.m -- verify the identity bypass in mtimes of
% 'intop', 'tensopvar' and 'tensopmat' (and in mtimes_tensopvar).
%
% Whenever one factor of a product is an identity operator, the product
% must return the other factor unchanged. Checked three ways, none of
% which uses the bypass's own conventions:
%  (i)   semantics: apply(A*B,X) against apply(A,apply(B,X)) with the
%        apply_* helpers (independent of mtimes), for identity and
%        non-identity factors, so both the bypass and the full path are
%        hit;
%  (ii)  structure: the bypass fires (C equals the other factor) exactly
%        for identity factors, incl. non-canonical identities (padded
%        bases);
%  (iii) guards: mismatched function spaces/dimensions still error, and
%        perturbed identities (beyond ztol) take the full path.
% Requires PIESOS_2D (not polyopvar) and PIETOOLS on the MATLAB path.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES:
% - A non-identity multiplier that only scales zero rows of the 'intop'
%   kernels leaves the coefficients equal to those of B (e.g.
%   diag(1.5,1,1) on a kernel with a zero first row), so C=A*B must
%   differ from B only if A*B(X) differs from B(X) (is_changed).
% - The full path of mtimes_tensopvar errors (reshape in get_params) when
%   a variable only passes through B; non-identity products there may
%   error, or must differ from the 'intop' operand.
% - eqop ignores 'chkval': isequal(B,B) is false for objects holding
%   'polynomial' parameters.
% - Multi-factor products with otype [1,2] and pointwise A*A are
%   unsupported and must error.

warning('off','sopvar:noncanonicalMultiplier');
rng(4242);
tol = 1e-8;
nchk = 0;
pvar s1 s2
dpvar dv1

% 1D building blocks on s1, s2 in [0,1]
va = struct('in',{{'s1'}},'out',{{'s1'}});  vb = struct('in',{{'s2'}},'out',{{'s2'}});
dom = struct('in',[0 1],'out',[0 1]);
Ia = make_id(1,va,dom);   Ib = make_id(1,vb,dom);   I2a = make_id(2,va,dom);
Tv = @(P) tensopvar(P);

%% ============================================================ intop mtimes
% A*B, A a multiplier, B a random functional on N variables
for trial=1:6
    m = randi(3);   n = randi(3);   N = randi(2);
    vars = arrayfun(@(k) ['x',num2str(k)],1:N,'UniformOutput',false);
    domN = repmat([0 1],N,1);
    Bop = rand_intop([m,n],vars,domN,randi(2,[1,N]),2);
    X = rand_poly([n,1],polynomial(vars(:)),2,50);
    BX = apply_op(Bop,X);
    ids = {eye(m),polynomial(eye(m)),dpvar(eye(m))};
    if m==1,  ids = [ids,{1}];  end
    for k=1:numel(ids)
        C = ids{k}*Bop;
        nchk = nchk+1;  chk(eqop(C,Bop),true,"intop id*B returns B, type "+class(ids{k})+" m="+m);
        nchk = nchk+1;  chk(small(apply_op(C,X)-BX,tol),true,"intop id*B semantics m="+m);
    end
    % non-identity multipliers take the full path and stay correct
    Amat = eye(m);  Amat(1,1) = 1.5;
    nonids = {Amat,2*eye(m),polynomial(Amat),eye(m)+1e-9*ones(m)};
    for k=1:numel(nonids)
        C = nonids{k}*Bop;
        dBX = nonids{k}*BX-BX;
        is_changed = nnz(dBX.C)>0;          % A*B(X) differs from B(X)
        nchk = nchk+1;  chk(eqop(C,Bop),~is_changed,"intop nonid*B not bypassed, k="+k+" m="+m);
        nchk = nchk+1;  chk(small(apply_op(C,X)-nonids{k}*BX,tol),true,"intop nonid*B semantics, k="+k+" m="+m);
    end
    % scalar identity against matrix-valued B; scalar on the right of an intop
    nchk = nchk+1;  chk(eqop(1*Bop,Bop),true,"intop 1*B");
    nchk = nchk+1;  chk(eqop(Bop*1,Bop),true,"intop B*1");
    nchk = nchk+1;  chk(throws(@() eye(m+1)*Bop),true,"intop eye(m+1)*B must error");
end
% functional with decision variables in the kernel
Bq = rand_intop([2,2],{'x1'},[0 1],1,2,3);
nchk = nchk+1;  chk(eqop(eye(2)*Bq,Bq),true,"intop eye*B with decision variables");
nchk = nchk+1;  chk(eqop(dpvar(eye(2))*Bq,Bq),true,"intop dpvar eye*B with decision variables");
nchk = nchk+1;  chk(eqop((eye(2)+0*dv1*eye(2))*Bq,Bq),true,"intop (eye+0*dvar)*B");
% a multiplier with decision variables is not an identity: dpvar*dpvar is not defined, so it errors as before
nchk = nchk+1;  chk(throws(@() (eye(2)+dv1*eye(2))*Bq),true,"intop (eye+dvar*eye)*B must not bypass");

%% ========================================= tensopvar mtimes (and sopvar)
Pa = rand_sop(1,va,dom);      Pb = rand_sop(1,vb,dom);     P2a = rand_sop(2,va,dom);
Ia_pad = pad_id(1,va,dom);
Xa = rand_poly([1,1],polynomial({'s1'}),3,50);
Xb = rand_poly([1,1],polynomial({'s2'}),3,50);
X2a = rand_poly([2,1],polynomial({'s1'}),3,50);

% single factor: identity on either side, canonical and padded
for kind={'id','pad'}
    Id = Ia;  if strcmp(kind{1},'pad'),  Id = Ia_pad;  end
    A = Tv(Pa);
    C1 = A*Tv(Id);      C2 = Tv(Id)*A;
    nchk = nchk+1;  chk(eqop(C1,A),true,"tensopvar A*id returns A ("+kind{1}+")");
    nchk = nchk+1;  chk(eqop(C2,A),true,"tensopvar id*A returns A ("+kind{1}+")");
    ref = apply_op(A,{apply_op(Tv(Id),{Xa})});
    nchk = nchk+1;  chk(small(apply_op(C1,{Xa})-ref,tol),true,"tensopvar A*id semantics ("+kind{1}+")");
    nchk = nchk+1;  chk(small(apply_op(C2,{Xa})-ref,tol),true,"tensopvar id*A semantics ("+kind{1}+")");
    nchk = nchk+1;  chk(eqop(A*Id,A),true,"tensopvar A*sopvar-id ("+kind{1}+")");
    nchk = nchk+1;  chk(is_id_op(Tv(Id)*Tv(Id)),true,"tensopvar id*id is identity ("+kind{1}+")");
end
for sc = [2,1+1e-6]
    A = Tv(Pa);   B = Tv(sc*Ia);
    C = A*B;
    nchk = nchk+1;  chk(eqop(C,A),false,"tensopvar A*(sc*id) not bypassed sc="+sc);
    ref = apply_op(A,{apply_op(B,{Xa})});
    nchk = nchk+1;  chk(small(apply_op(C,{Xa})-ref,tol),true,"tensopvar A*(sc*id) semantics sc="+sc);
    C = B*A;
    nchk = nchk+1;  chk(eqop(C,A),false,"tensopvar (sc*id)*A not bypassed sc="+sc);
    ref = apply_op(B,{apply_op(A,{Xa})});
    nchk = nchk+1;  chk(small(apply_op(C,{Xa})-ref,tol),true,"tensopvar (sc*id)*A semantics sc="+sc);
end
A = Tv(Pa);   B = Tv(rand_sop(1,va,dom));
C = A*B;     ref = apply_op(A,{apply_op(B,{Xa})});
nchk = nchk+1;  chk(small(apply_op(C,{Xa})-ref,tol),true,"tensopvar random*random semantics");
% multi-term operand, identity as one of two terms is NOT an identity
A2 = Tv(Pa) + Tv(Ia);
nchk = nchk+1;  chk(eqop(A*Tv(Ia+Ia),A),false,"tensopvar A*(id+id) not bypassed");
nchk = nchk+1;  chk(eqop(Tv(P2a)*Tv(I2a),Tv(P2a)),true,"tensopvar p=2 A*id");
nchk = nchk+1;  chk(eqop(Tv(I2a)*Tv(P2a),Tv(P2a)),true,"tensopvar p=2 id*A");
% function-space mismatches must still error, also for identity factors
nchk = nchk+1;  chk(throws(@() Tv(Pb)*Tv(Ia)),true,"tensopvar id on s1 against op on s2 must error");
nchk = nchk+1;  chk(throws(@() Tv(Ib)*Tv(Pa)),true,"tensopvar op on s1 into id on s2 must error");
nchk = nchk+1;  chk(throws(@() Tv(P2a)*Tv(Ia)),true,"tensopvar dimension mismatch must error");
nchk = nchk+1;  chk(throws(@() Tv(I2a)*Tv(Pa)),true,"tensopvar id(2) times 1-dim op must error");
% identity on a different interval is not composable: not bypassed
vdom2 = struct('in',[0 2],'out',[0 2]);
Ia2 = make_id(1,va,vdom2);
try
    Cd = Tv(Ia2)*Tv(Pa);
    nchk = nchk+1;  chk(eqop(Cd,Tv(Pa)),false,"tensopvar id on [0,2] times op on [0,1] not bypassed");
catch
    nchk = nchk+1;  % erroring is also acceptable: the bypass must not hide the mismatch
end

A = Tv(Pa);   A2d = Tv(P2a);
nchk = nchk+1;  chk(eqop(eye(1)*A,A),true,"tensopvar eye(1)*A");
nchk = nchk+1;  chk(eqop(1*A,A),true,"tensopvar 1*A");
nchk = nchk+1;  chk(eqop(A*1,A),true,"tensopvar A*1");
nchk = nchk+1;  chk(eqop(eye(2)*A2d,A2d),true,"tensopvar eye(2)*A (p=2)");
nchk = nchk+1;  chk(eqop(polynomial(eye(2))*A2d,A2d),true,"tensopvar polynomial eye(2)*A");
nchk = nchk+1;  chk(eqop(dpvar(eye(2))*A2d,A2d),true,"tensopvar dpvar eye(2)*A");
nchk = nchk+1;  chk(eqop(2*eye(2)*A2d,A2d),false,"tensopvar 2*eye(2)*A not bypassed");
nchk = nchk+1;  chk(small(apply_op(2*eye(2)*A2d,{X2a})-2*apply_op(A2d,{X2a}),tol),true,"tensopvar 2*eye(2)*A semantics");
nchk = nchk+1;  chk(throws(@() eye(3)*A2d),true,"tensopvar eye(3)*A (2 rows) must error");
nchk = nchk+1;  chk(throws(@() eye(2)*A),true,"tensopvar eye(2)*A (1 row) must error");

% multi-factor, Kronecker rows/cols ([2,2]), and pointwise products
mkT = @(c,type) sopvar2tensopvar(c,type);
for type = {[2,2],[1,2]}
    Tid = mkT({Ia,Ib},type{1});    Tab = mkT({Pa,Pb},type{1});
    nchk = nchk+1;  chk(is_id_op(Tid),type{1}(1)==type{1}(2),"tensopvar two-factor identity is identity iff otype(1)==otype(2), type="+mat2str(type{1}));
    if ~isequal(type{1},[2,2])
        % not an identity (otype differs), pointwise rows: composition still unsupported
        nchk = nchk+1;  chk(throws(@() Tab*Tid),true,"tensopvar [1,2] multi-factor product still errors");
    else
        C = Tab*Tid;
        nchk = nchk+1;  chk(eqop(C,Tab),true,"tensopvar multi-factor A*id, type="+mat2str(type{1}));
        C = Tid*Tab;
        nchk = nchk+1;  chk(eqop(C,Tab),true,"tensopvar multi-factor id*A, type="+mat2str(type{1}));
        Cs = mkT({2*Ia,Ib},type{1});
        C = Tab*Cs;
        nchk = nchk+1;  chk(eqop(C,Tab),false,"tensopvar multi-factor A*(2id (x) id) not bypassed");
        % semantics of the full path next to the bypass: (A*(2id (x) id))(Xa,Xb) = A(2*Xa,Xb) (Kronecker factors)
        nchk = nchk+1;  chk(small(apply_op(C,{Xa,Xb})-apply_op(Tab,{apply_sopvar(2*Ia,Xa),Xb}),tol),true,"tensopvar multi-factor A*(2id (x) id) semantics");
    end
    nchk = nchk+1;  chk(eqop(1*Tab,Tab),true,"tensopvar 1*multi-factor, type="+mat2str(type{1}));
    nchk = nchk+1;  chk(eqop(eye(size(Tab,1))*Tab,Tab),true,"tensopvar eye(size)*multi-factor, type="+mat2str(type{1}));
    nchk = nchk+1;  chk(small(apply_op(3*Tab,{Xa,Xb})-3*apply_op(Tab,{Xa,Xb}),tol),true,"tensopvar 3*multi-factor semantics");
end
nchk = nchk+1;  chk(throws(@() mkT({Pa,Pb},[2,2])*Tv(Ia)),true,"tensopvar 2-factor times 1-factor identity must error");
nchk = nchk+1;  chk(throws(@() mkT({Ia,Ib},[2,2])*Tv(Pa)),true,"tensopvar 2-factor identity times 1-factor must error");
% pointwise product, identity: no mixed product property needed
Tp_id = mkT({Ia,Ia},[1,1]);    Tp = mkT({Pa,Pa},[1,1]);
nchk = nchk+1;  chk(eqop(Tp*Tp_id,Tp),true,"tensopvar pointwise A*id");
nchk = nchk+1;  chk(eqop(Tp_id*Tp,Tp),true,"tensopvar pointwise id*A");
nchk = nchk+1;  chk(throws(@() Tp*Tp),true,"tensopvar pointwise A*A still errors");
% identity whose otype disagrees with that of the other operand: not composable
nchk = nchk+1;  chk(throws(@() mkT({Pa,Pa},[2,2])*mkT({Ia,Ia},[1,1])),true,"tensopvar otype mismatch with identity must error");

%% ================================================ mtimes_tensopvar (intop)
% functional on L2[s1] -> R^m times tensopvar
for trial=1:6
    m = randi(3);   n = 1;
    Kop = rand_intop([m,n],{'s1'},[0 1],1,2);
    Xk = rand_poly([n,1],polynomial({'s1'}),2,50);
    for kind={'id','pad'}
        Id = Ia;  if strcmp(kind{1},'pad'),  Id = Ia_pad;  end
        C = Kop*Tv(Id);
        nchk = nchk+1;  chk(eqop(C,Kop),true,"intop*tensopvar id returns K ("+kind{1}+")");
        ref = apply_op(Kop,apply_op(Tv(Id),{Xk}));
        nchk = nchk+1;  chk(small(apply_op(C,Xk)-ref,tol),true,"intop*tensopvar id semantics ("+kind{1}+")");
        Cm = mtimes_tensopvar(Kop,Tv(Id));
        nchk = nchk+1;  chk(eqop(Cm,Kop),true,"mtimes_tensopvar K*id returns K ("+kind{1}+")");
    end
    % Non-identity operands must not be bypassed; the full path may error (see NOTES).
    for Bs = {Tv(2*Ia),Tv(rand_mult(1,va,dom))}
        try
            C = Kop*Bs{1};
            nchk = nchk+1;  chk(eqop(C,Kop),false,"intop*tensopvar non-identity not bypassed");
        catch
            nchk = nchk+1;
        end
    end
    nchk = nchk+1;  chk(throws(@() Kop*Tv(Ib)),true,"intop*tensopvar id on other variable must error");
end
% tensor-product identity on two variables: a functional in both, and a mismatched one
vn2 = {'s1','s2'};   dom2 = [0 1;0 1];
Tid2 = sopvar2tensopvar({Ia,Ib},[2,2]);
for trial=1:4
    K2 = rand_intop([2,1],vn2,dom2,[1 1],2);
    X2 = rand_poly([1,1],polynomial(vn2(:)),2,50);
    C = K2*Tid2;
    nchk = nchk+1;  chk(eqop(C,K2),true,"intop*tensopvar id (x) id returns K, two variables");
    nchk = nchk+1;  chk(small(apply_op(C,X2)-apply_op(K2,X2),tol),true,"intop*tensopvar id (x) id semantics, two variables");
end
nchk = nchk+1;  chk(throws(@() rand_intop([1,1],{'s1'},[0 1],1,2)*Tid2),true,"intop*tensopvar id (x) id with mismatched input space must error");

%% ================================================== tensopmat mtimes
% Blocks are 1x1 operators on s1, so every block dimension is 1
mk = @() tensopvar(rand_sop(1,va,dom));
Ia_ = Tv(Ia);
for nb = 1:3
    % B: nb x nb blocks; A: nb x nb identity (empty off-diagonal, or explicit zeros)
    Bm = cell(nb,nb);   Am_empty = cell(nb,nb);   Am_zero = cell(nb,nb);   Am_pad = cell(nb,nb);
    for i=1:nb
        for j=1:nb
            Bm{i,j} = mk();
            Am_zero{i,j} = Tv(0*Ia);
        end
        Am_empty{i,i} = Ia_;   Am_zero{i,i} = Ia_;   Am_pad{i,i} = Tv(Ia_pad);
    end
    Bmat = tensopmat(Bm);
    X = cell(1,nb);
    for j=1:nb,  X{j} = rand_poly([1,1],polynomial({'s1'}),3,50);  end
    Yb = apply_op(Bmat,X);
    for Aname={'empty','zero','pad'}
        switch Aname{1}
            case 'empty',  Aid = tensopmat(Am_empty);
            case 'zero',   Aid = tensopmat(Am_zero);
            case 'pad',    Aid = tensopmat(Am_pad);
        end
        nchk = nchk+1;  chk(is_id_op(Aid),true,"tensopmat test identity is identity, nb="+nb+" "+Aname{1});
        C1 = Aid*Bmat;       C2 = Bmat*Aid;
        % the bypass returns the cleaned operand
        Bc = clean_op(Bmat);
        nchk = nchk+1;  chk(eqop(C1,Bc),true,"tensopmat id*B returns B, nb="+nb+" "+Aname{1});
        nchk = nchk+1;  chk(eqop(C2,Bc),true,"tensopmat B*id returns B, nb="+nb+" "+Aname{1});
        % semantics: apply(id*B,X) against apply(id,apply(B,X)), and apply(B*id,X) against apply(B,apply(id,X))
        Yid = apply_op(Aid,num2cell(Yb(:).'));
        nchk = nchk+1;  chk(small(apply_op(C1,X)-Yid,tol),true,"tensopmat id*B semantics, nb="+nb+" "+Aname{1});
        Xid = num2cell(apply_op(Aid,X).');
        nchk = nchk+1;  chk(small(apply_op(C2,X)-apply_op(Bmat,Xid),tol),true,"tensopmat B*id semantics, nb="+nb+" "+Aname{1});
    end
    % perturbed identities beyond ztol: full path
    Apert = tensopmat(Am_empty);   Apert.ops{1,1} = Tv(2*Ia);
    C1 = Apert*Bmat;
    nchk = nchk+1;  chk(eqop(C1,clean_op(Bmat)),false,"tensopmat (2,1,..)*B not bypassed, nb="+nb);
    Yp = apply_op(Apert,num2cell(Yb(:).'));
    nchk = nchk+1;  chk(small(apply_op(C1,X)-Yp,tol),true,"tensopmat (2,1,..)*B semantics, nb="+nb);
    Apert = tensopmat(Am_zero);    Apert.ops{1,nb} = Tv(1e-3*Ia);
    if nb>1,  Apert.ops{1,nb} = Tv(1e-3*Ia);  end
    C1 = Apert*Bmat;
    if nb>1
        nchk = nchk+1;  chk(eqop(C1,clean_op(Bmat)),false,"tensopmat nonzero off-diagonal not bypassed, nb="+nb);
        Yp = apply_op(Apert,num2cell(Yb(:).'));
        nchk = nchk+1;  chk(small(apply_op(C1,X)-Yp,tol),true,"tensopmat nonzero off-diagonal semantics, nb="+nb);
    end
    if nb==1
        nchk = nchk+1;  chk(eqop(eye(1)*Bmat,clean_op(Bmat)) || eqop(eye(1)*Bmat,Bmat),true,"tensopmat eye(1)*B");
    end
    nchk = nchk+1;  chk(eqop(1*Bmat,Bmat),true,"tensopmat 1*B, nb="+nb);
    nchk = nchk+1;  chk(eqop(Bmat*1,Bmat),true,"tensopmat B*1, nb="+nb);
    nchk = nchk+1;  chk(eqop(polynomial(1)*Bmat,Bmat),true,"tensopmat polynomial(1)*B, nb="+nb);
    nchk = nchk+1;  chk(eqop(dpvar(1)*Bmat,Bmat),true,"tensopmat dpvar(1)*B, nb="+nb);
    nchk = nchk+1;  chk(eqop(2*Bmat,Bmat),false,"tensopmat 2*B not bypassed, nb="+nb);
    nchk = nchk+1;  chk(small(apply_op(2*Bmat,X)-2*Yb,tol),true,"tensopmat 2*B semantics, nb="+nb);
    nchk = nchk+1;  chk(throws(@() eye(nb+1)*Bmat),true,"tensopmat eye(nb+1)*B must error, nb="+nb);
end

% single row of blocks: eye(m)*B returns B (block structure preserved)
B1 = tensopmat({mk(),mk(),mk()});
X = {rand_poly([1,1],polynomial({'s1'}),3,50),rand_poly([1,1],polynomial({'s1'}),3,50),rand_poly([1,1],polynomial({'s1'}),3,50)};
nchk = nchk+1;  chk(eqop(eye(1)*B1,B1),true,"tensopmat eye(1)*(1x3 blocks) returns B");
nchk = nchk+1;  chk(small(apply_op(2*eye(1)*B1,X)-2*apply_op(B1,X),tol),true,"tensopmat 2*eye(1)*(1x3) semantics");
% several rows of blocks: eye(m)*B keeps the (existing) single-row-block output
B2 = tensopmat({mk(),mk();mk(),mk()});
X2 = {rand_poly([1,1],polynomial({'s1'}),3,50),rand_poly([1,1],polynomial({'s1'}),3,50)};
C = eye(2)*B2;
nchk = nchk+1;  chk(size(C.ops,1),1,"tensopmat eye(2)*(2x2 blocks) is a single row-block");
nchk = nchk+1;  chk(small(apply_op(C,X2)-apply_op(B2,X2),tol),true,"tensopmat eye(2)*(2x2 blocks) semantics");
va_s2 = Tv(Ib);
nchk = nchk+1;  chk(throws(@() tensopmat({Ib})*Bmat),true,"tensopmat id on s2 times op on s1 must error");
nchk = nchk+1;  chk(throws(@() tensopmat({Tv(I2a)})*B1),true,"tensopmat id with mismatched dimension must error");
% blocks of different types: intop blocks next to tensopvar blocks
Kid = intop(eye(2));
Bmix = tensopmat({Tv(Pa),[];[],Kid});     % block-diagonal, L2[s1] x R^2
Amix = tensopmat({Ia_,[];[],intop(eye(2))});
C = Amix*Bmix;
nchk = nchk+1;  chk(eqop(C,clean_op(Bmix)),true,"tensopmat mixed-type id*B returns B");
C = Bmix*Amix;
nchk = nchk+1;  chk(eqop(C,clean_op(Bmix)),true,"tensopmat mixed-type B*id returns B");

fprintf('Test_mtimes_id_bypass passed (%d checks).\n',nchk);


%% -------------------------------------------------------------- helpers
function chk(val,expected,name)
if ~isequal(val,expected) && ~isequal(logical(val),logical(expected))
    error("Test_mtimes_id_bypass failed: %s (got %s, expected %s).",name,num2str(val),num2str(expected));
end
end

function tf = small(D,tol)
% true if the (constant) polynomial D is zero up to a relative tolerance
D = cleanpoly(D,1e-10);
if ~isdouble(D),  tf = false;  return,  end
tf = max(abs(full(double(D(:)))))<=tol;
end

function tf = eqop(a,b)
% Structural operator equality, ignoring the hidden 'chkval' flag. Needed since
% isequal(B,B) is false for objects holding 'polynomial' parameters.
if isobject(a)
    tf = strcmp(class(a),class(b));
    if tf
        w = warning('off','MATLAB:structOnObject');
        sa = struct(a);   sb = struct(b);
        warning(w);
        if isfield(sa,'chkval'),  sa = rmfield(sa,'chkval');  end
        if isfield(sb,'chkval'),  sb = rmfield(sb,'chkval');  end
        tf = eqop(sa,sb);
    end
elseif iscell(a)
    tf = iscell(b) && isequal(size(a),size(b));
    for k=1:numel(a)
        if ~tf,  break,  end
        tf = eqop(a{k},b{k});
    end
elseif isstruct(a)
    tf = isstruct(b) && isequal(sort(fieldnames(a)),sort(fieldnames(b)));
    fn = fieldnames(a);
    for k=1:numel(fn)
        if ~tf,  break,  end
        tf = eqop(a.(fn{k}),b.(fn{k}));
    end
else
    tf = isequal(a,b);
end
end
function tf = throws(f)
tf = false;
try
    f();
catch
    tf = true;
end
end

function Iop = make_id(p,vars,dom)
% Canonical identity I_p on L2^p[dom]: multiplier in every variable.
N = numel(vars.in);
params = [{speye(p)},repmat({sparse(p,p)},1,3^N-1)];
Iop = sopvar(params,vars,repmat({0},1,N),repmat({0},1,N),dom,[p p]);
end

function Iop = pad_id(p,vars,dom)
% Identity written on an enlarged monomial basis {1,s} (non-canonical form)
Iop = make_id(p,vars,dom);
Q = sopvar({sparse(p*2,p*2),sparse(p*2,p*2),sparse(p*2,p*2)},vars,{[0;1]},{[0;1]},dom,[p p]);
Iop = Iop+Q;
end

function P = rand_sop(p,vars,dom)
% Random p x p operator on one variable, monomials {1,s} in all kernels
params = cell(1,3);
for c=1:3,  params{c} = sprandn(p*2,p*2,0.7);  end
P = sopvar(params,vars,{[0;1]},{[0;1]},dom,[p p]);
end

function P = rand_mult(p,vars,dom)
% Random p x p multiplier operator on one variable (monomials {1,s}), no integral kernels
P = sopvar({sprandn(p*2,p*2,0.9),sparse(p*2,p*2),sparse(p*2,p*2)},vars,{[0;1]},{[0;1]},dom,[p p]);
end
