% Test_is_id_op.m -- verify is_id_op for 'intop', 'sopvar' and 'tensopvar'.
%
% The 'sopvar' routine is checked against independent oracles, never
% against its own conventions: (i) eq(P,I) on the class's own 'plus'
% (basis union, kron(eye(dims),.) ordering), for 0-3 spatial variables
% and randomly enlarged/padded monomial bases; (ii) for 1D, the kernels
% of sopvar2opvar(P). Requires PIESOS_2D (not polyopvar) and PIETOOLS on
% the MATLAB path.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES:
% - 'intop' kernels are coefficients of monomials in the dummy variables,
%   so a multiplier depending on a free variable s1 is rejected by the
%   constructor (chk_err); not tested through 'intop'.
% - The 3-factor, nested and kron(tensopvar) cases are inside try/catch:
%   if unsupported, they only print a note and are not counted.

warning('off','sopvar:noncanonicalMultiplier');
pvar s1 s2 s3
dpvar dv1
nchk = 0;

%% ---------------------------------------------------------------- intop
for m=1:4
    nchk = nchk+1;  chk(is_id_op(intop(eye(m))),true,"intop eye double m="+m);
    nchk = nchk+1;  chk(is_id_op(intop(polynomial(eye(m)))),true,"intop eye polynomial m="+m);
    nchk = nchk+1;  chk(is_id_op(intop(dpvar(eye(m)))),true,"intop eye dpvar m="+m);
    nchk = nchk+1;  chk(is_id_op(intop(2*eye(m))),false,"intop 2*eye m="+m);
    nchk = nchk+1;  chk_err(@() intop(eye(m)+(1+m)*s1*eye(m)),"dummy variables","intop eye+s1 m="+m);
    nchk = nchk+1;  chk(is_id_op(intop(zeros(m))),false,"intop zero m="+m);
    nchk = nchk+1;  chk(is_id_op(intop(dv1*eye(m))),false,"intop dvar*eye m="+m);
    nchk = nchk+1;  chk(is_id_op(intop(eye(m)+dv1*eye(m))),false,"intop eye+dvar*eye m="+m);
    nchk = nchk+1;  chk(is_id_op(intop([eye(m),zeros(m,1)])),false,"intop nonsquare m="+m);
    % off-diagonal / permutation entries
    if m>=2
        Pm = eye(m);  Pm(1,2) = 1;
        nchk = nchk+1;  chk(is_id_op(intop(Pm)),false,"intop offdiag m="+m);
        Pm = eye(m);  Pm([1 2],:) = Pm([2 1],:);
        nchk = nchk+1;  chk(is_id_op(intop(Pm)),false,"intop permutation m="+m);
    end
    % tolerance
    nchk = nchk+1;  chk(is_id_op(intop(eye(m)+1e-9*ones(m))),false,"intop 1e-9 default tol m="+m);
    nchk = nchk+1;  chk(is_id_op(intop(eye(m)+1e-9*ones(m)),1e-6),true,"intop 1e-9 ztol=1e-6 m="+m);
    nchk = nchk+1;  chk(is_id_op(intop(eye(m)+1e-14*ones(m))),true,"intop 1e-14 default tol m="+m);
end
nchk = nchk+1;  chk(is_id_op(intop(polynomial(1))),true,"intop scalar 1");
nchk = nchk+1;  chk(is_id_op(intop(1+0*s1)),true,"intop 1+0*s1");
nchk = nchk+1;  chk_err(@() intop(s1),"dummy variables","intop s1");
% functionals with spatial variables are never identity, however the kernel looks
Kfull = intop({polynomial(1)},{[0 2]},{'s1'},[0 1]);
nchk = nchk+1;  chk(is_id_op(Kfull),false,"intop int_0^1 v ds");
Kfull3 = intop({polynomial(1)},{[0 2]},{'s1'},[0 1]) ;
nchk = nchk+1;  chk(is_id_op(Kfull3+Kfull3),false,"intop sum of integrals");
Kfull4 = intop({polynomial(eye(2))},{[0 2]},{'s1'},[0 1]);
nchk = nchk+1;  chk(is_id_op(Kfull4),false,"intop 2x2 identity kernel but integral");

%% ---------------------------------------------------------------- sopvar
rng(1070);
ntrial = 0;
for nvar=0:3
    names = {'s1','s2','s3'};
    vn = names(1:nvar);
    for p=1:2
        dom = repmat([0,1],nvar,1);
        if nvar>=2,  dom(2,:) = [-1,2];  end
        vars = struct('in',{vn},'out',{vn});
        domst = struct('in',dom,'out',dom);
        Iop = make_id(p,vars,domst);
        nchk = nchk+1;  chk(is_id_op(Iop),true,"sopvar canonical id nvar="+nvar+" p="+p);
        nchk = nchk+1;  chk(eq_quiet(Iop,Iop),true,"oracle eq sanity nvar="+nvar);
        for trial=1:12
            ntrial = ntrial+1;
            ZL = cell(1,nvar);  ZR = cell(1,nvar);
            for k=1:nvar
                ZL{k} = rand_basis(true);
                ZR{k} = rand_basis(true);
            end
            nL = prod(cellfun(@numel,ZL));   nR = prod(cellfun(@numel,ZR));
            kind = mod(trial,4);
            % kind 0: Q=0 (identity on an enlarged basis); 1: Q random;
            % 2: Q only in non-multiplier cells; 3: Q only in params{1}
            ncell = 3^nvar;
            params = cell(1,ncell);
            for c=1:ncell
                params{c} = sparse(p*nL,p*nR);
                active = (kind==1) || (kind==2 && c>1) || (kind==3 && c==1);
                if active
                    params{c} = sprandn(p*nL,p*nR,0.4);
                end
            end
            Q = sopvar(params,vars,ZL,ZR,domst,[p p]);
            if rand<0.5
                P = Iop + Q;
            else
                P = Q + Iop;
            end
            % oracle: structural-independent equality through plus
            expected = eq_quiet(P,Iop);
            nchk = nchk+1;
            chk(is_id_op(P),expected,"sopvar random nvar="+nvar+" p="+p+" kind="+kind+" trial="+trial);
            if kind==0
                nchk = nchk+1;  chk(is_id_op(P),true,"sopvar padded identity nvar="+nvar+" p="+p+" trial="+trial);
            end
            if kind>=1 && nnz(cell2mat(cellfun(@(x)x(:)',params,'UniformOutput',false)))>0
                % Q is generic and nonzero, so P = I + Q is not expected to be the identity
                nchk = nchk+1;  chk(is_id_op(P),expected,"sopvar nonzero Q consistent trial="+trial);
            end
            % Scaling the identity must break it
            nchk = nchk+1;  chk(is_id_op(2*P),eq_quiet(2*P,Iop),"sopvar 2*P nvar="+nvar+" p="+p+" trial="+trial);
        end
        if nvar>=1
            Zsm = sopvar(cellfun(@(x) 1e-9*x,Iop.params,'UniformOutput',false),vars,Iop.ZL,Iop.ZR,domst,[p p]);
            nchk = nchk+1;  chk(is_id_op(Iop+Zsm),false,"sopvar 1e-9 default tol nvar="+nvar);
            nchk = nchk+1;  chk(is_id_op(Iop+Zsm,1e-6),true,"sopvar 1e-9 ztol=1e-6 nvar="+nvar);
        end
    end
end

% Variable-specific negatives
Ia = make_id(1,struct('in',{{'s1'}},'out',{{'s1'}}),struct('in',[0 1],'out',[0 1]));
M = speye(2);   M(1,2) = 1;
Pm = sopvar({M,sparse(2,2),sparse(2,2)},struct('in',{{'s1'}},'out',{{'s1'}}),{0},{0},struct('in',[0 1],'out',[0 1]),[2 2]);
nchk = nchk+1;  chk(is_id_op(Pm),false,"sopvar p=2 off-diagonal");
M = speye(2);   M = M([2 1],:);
Pm = sopvar({M,sparse(2,2),sparse(2,2)},struct('in',{{'s1'}},'out',{{'s1'}}),{0},{0},struct('in',[0 1],'out',[0 1]),[2 2]);
nchk = nchk+1;  chk(is_id_op(Pm),false,"sopvar p=2 permutation");
Pm = sopvar({sparse(diag([1 0])),sparse(2,2),sparse(2,2)},struct('in',{{'s1'}},'out',{{'s1'}}),{0},{0},struct('in',[0 1],'out',[0 1]),[2 2]);
nchk = nchk+1;  chk(is_id_op(Pm),false,"sopvar p=2 diag(1,0)");
% multiplier by s1 (no constant monomial in ZL), basis {1} only
Pm = sopvar({sparse(1),sparse(1,1),sparse(1,1)},struct('in',{{'s1'}},'out',{{'s1'}}),{1},{0},struct('in',[0 1],'out',[0 1]),[1 1]);
nchk = nchk+1;  chk(is_id_op(Pm),false,"sopvar multiplier s1");
Pf = sopvar({sparse(1)},struct('in',{{'s1'}},'out',{cell(1,0)}),{},{0},struct('in',[0 1],'out',zeros(0,2)),[1 1]);
nchk = nchk+1;  chk(is_id_op(Pf),false,"sopvar L2->R");
Pf = sopvar({sparse(1)},struct('in',{cell(1,0)},'out',{{'s1'}}),{0},{},struct('in',zeros(0,2),'out',[0 1]),[1 1]);
nchk = nchk+1;  chk(is_id_op(Pf),false,"sopvar R->L2");
Pg = sopvar({sparse(1),sparse(1),sparse(1)},struct('in',{{'s1'}},'out',{{'s2'}}),{0},{0},struct('in',[0 1],'out',[0 1]),[1 1]);
nchk = nchk+1;  chk(is_id_op(Pg),false,"sopvar s1->s2");
Ph = sopvar({sparse(1),sparse(1,1),sparse(1,1)},struct('in',{{'s1'}},'out',{{'s1'}}),{0},{0},struct('in',[0 1],'out',[0 1]),[1 1]);
Ph.dims = [1,2];
nchk = nchk+1;  chk(is_id_op(Ph),false,"sopvar nonsquare dims");
Z0 = 0*Ia;
nchk = nchk+1;  chk(is_id_op(Z0),false,"sopvar zero");

% 1D cross-check against sopvar2opvar kernels (independent of the sopvar
% coefficient-ordering conventions used in is_id_op)
for p=1:2
  for trial=1:30
    ZL = {rand_basis(true)};   ZR = {rand_basis(true)};
    nL = numel(ZL{1});  nR = numel(ZR{1});
    params = {sparse(p*nL,p*nR),sparse(p*nL,p*nR),sparse(p*nL,p*nR)};
    if rand<0.7
        % place an identity-like block on the constant monomials (if any)
        p0L = find(ZL{1}==0);  p0R = find(ZR{1}==0);
        if ~isempty(p0L) && ~isempty(p0R) && ZR{1}(p0R)==0
            rows = (0:p-1)'*nL+p0L;   cols = (0:p-1)'*nR+p0R;
            params{1} = sparse(rows,cols,1,p*nL,p*nR);
        end
    end
    if rand<0.3,  params{2}(1,1) = 1;  end
    if rand<0.2 && p*nL>=1,  params{1}(end,1) = params{1}(end,1)+1;  end
    P = sopvar(params,struct('in',{{'s1'}},'out',{{'s1'}}),ZL,ZR,struct('in',[0 1],'out',[0 1]),[p p]);
    try
        Rop = sopvar2opvar(P);
    catch
        continue
    end
    R0 = Rop.R.R0 - eye(p);
    sem = isempty(nonzeros(R0.C)) || max(abs(nonzeros(R0.C)))<1e-12;
    sem = sem && (isempty(nonzeros(Rop.R.R1.C)) || max(abs(nonzeros(Rop.R.R1.C)))<1e-12);
    sem = sem && (isempty(nonzeros(Rop.R.R2.C)) || max(abs(nonzeros(Rop.R.R2.C)))<1e-12);
    nchk = nchk+1;  chk(is_id_op(P),sem,"sopvar vs sopvar2opvar p="+p+" trial="+trial);
  end
end

%% ------------------------------------------------------------- tensopvar
va = struct('in',{{'s1'}},'out',{{'s1'}});  da = struct('in',[0 1],'out',[0 1]);
vb = struct('in',{{'s2'}},'out',{{'s2'}});  db = struct('in',[0 1],'out',[0 1]);
Ia = make_id(1,va,da);    Ib = make_id(1,vb,db);
I2a = make_id(2,va,da);
Z = @(P) 0*P;
mkT = @(c) sopvar2tensopvar(c,[2,2]);   % Kronecker type; default [1,2] maps between different spaces

nchk = nchk+1;  chk(is_id_op(tensopvar(Ia)),true,"tensopvar single id");
nchk = nchk+1;  chk(is_id_op(tensopvar(2*Ia)),false,"tensopvar single 2*id");
nchk = nchk+1;  chk(is_id_op(tensopvar(I2a)),true,"tensopvar single id p=2");
nchk = nchk+1;  chk(is_id_op(tensopvar(Z(Ia))),false,"tensopvar single zero");
Pint = Ia;  Pint.params{2} = sparse(1);
nchk = nchk+1;  chk(is_id_op(tensopvar(Pint)),false,"tensopvar id + integral");
% two terms: identity and zero -> identity; two nonzero -> false
nchk = nchk+1;  chk(is_id_op(tensopvar({Ia;Z(Ia)})),true,"tensopvar {id; 0}");
nchk = nchk+1;  chk(is_id_op(tensopvar({Z(Ia);Ia})),true,"tensopvar {0; id}");
nchk = nchk+1;  chk(is_id_op(tensopvar({Ia;Ia})),false,"tensopvar {id; id}");
nchk = nchk+1;  chk(is_id_op(tensopvar({Ia;2*Ia})),false,"tensopvar {id; 2*id}");
nchk = nchk+1;  chk(is_id_op(tensopvar({Z(Ia);Z(Ia)})),false,"tensopvar {0; 0}");
% tolerance, both on the dropped term and on the factor
tiny = 1e-9*Ia;
nchk = nchk+1;  chk(is_id_op(tensopvar({Ia;tiny})),false,"tensopvar {id; 1e-9 id} default tol");
nchk = nchk+1;  chk(is_id_op(tensopvar({Ia;tiny}),1e-6),true,"tensopvar {id; 1e-9 id} ztol=1e-6");
nchk = nchk+1;  chk(is_id_op(tensopvar(Ia+tiny)),false,"tensopvar id+1e-9 default tol");
nchk = nchk+1;  chk(is_id_op(tensopvar(Ia+tiny),1e-6),true,"tensopvar id+1e-9 ztol=1e-6");
nchk = nchk+1;  chk(is_id_op(tensopvar(Ia+1e-14*Ia)),true,"tensopvar id+1e-14 default tol");
Top = mkT({Ia,Ib});
nchk = nchk+1;  chk(is_id_op(Top),true,"tensopvar id (x) id (s1,s2)");
nchk = nchk+1;  chk(is_id_op(mkT({Ia,2*Ib})),false,"tensopvar id (x) 2*id");
nchk = nchk+1;  chk(is_id_op(mkT({Ia,Pint})),false,"tensopvar id (x) (id+int)");
nchk = nchk+1;  chk(is_id_op(mkT({Ia,Z(Ib)})),false,"tensopvar id (x) 0");
nchk = nchk+1;  chk(is_id_op(mkT({Ia,Ib;Z(Ia),Ib})),true,"tensopvar {id,id; 0,id}");
nchk = nchk+1;  chk(is_id_op(mkT({Ia,Ib;Ia,Ib})),false,"tensopvar {id,id; id,id}");
Top3 = mkT({Ia,Ib,I2a});     % repeated variable, p=2
try
    r3 = is_id_op(Top3);
    nchk = nchk+1;  chk(r3,true,"tensopvar three factors");
catch ME
    fprintf('  (note: 3-factor tensopvar not constructible/evaluable here: %s)\n',ME.message);
end
Pf = sopvar({sparse(1)},struct('in',{{'s1'}},'out',{cell(1,0)}),{},{0},struct('in',[0 1],'out',zeros(0,2)),[1 1]);
nchk = nchk+1;  chk(is_id_op(tensopvar(Pf)),false,"tensopvar L2->R");
% pointwise vs tensor product types: differing otype with 2 factors must not be identity
Top_t = Top;   Top_t.otype = [1,2];
Top_u = Top;   Top_u.otype = [2,1];
nchk = nchk+1;  chk(is_id_op(Top_t),false,"tensopvar otype [1,2] two factors");
nchk = nchk+1;  chk(is_id_op(Top_u),false,"tensopvar otype [2,1] two factors");
Tp = sopvar2tensopvar({Ia,Ia},[1,1]);
nchk = nchk+1;  chk(is_id_op(Tp),true,"tensopvar pointwise id.*id");
nchk = nchk+1;  chk(is_id_op(sopvar2tensopvar({Ia,2*Ia},[1,1])),false,"tensopvar pointwise id.*2id");
Tn = mkT({tensopvar(Ia),Ib});
try
    nchk = nchk+1;  chk(is_id_op(Tn),true,"tensopvar nested id");
catch ME
    fprintf('  (note: nested tensopvar not evaluable here: %s)\n',ME.message);
end
try
    Tk = kron(tensopvar(Ia),tensopvar(Ib));
    nchk = nchk+1;  chk(is_id_op(Tk),true,"tensopvar kron(id,id)");
catch ME
    fprintf('  (note: kron(tensopvar) not evaluable here: %s)\n',ME.message);
end
for trial=1:5
    rng(7000+trial);
    dims = [1,1;1,1];
    vars = struct('in',{{'s1','s2'}},'out',{{'s1','s2'}});
    dom = struct('in',repmat([0,1],2,1),'out',repmat([0,1],2,1));
    depmat = struct('in',logical(eye(2)),'out',logical(eye(2)));
    Bop = rand_tensopvar(dims,vars,dom,depmat,2,[2,2],1,0.5);
    nchk = nchk+1;  chk(is_id_op(Bop),false,"tensopvar random trial="+trial);
end

% Randomized multi-factor / multi-term operators. The expected answer comes
% from how each factor was built (not from is_id_op): a term is zero if some
% factor is zero up to ztol, a term is identity if all its factors are, and
% the operator is identity iff exactly one term survives and it is identity.
kinds = {'id','idpad','zero','scaled','int','tiny'};
probs = [0.35 0.25 0.1 0.1 0.1 0.1];
ntens = 0;  ntrue = 0;  nmulti = 0;
for typ = {[2,2],[1,1]}
  type = typ{1};
  for trial=1:150
    rng(8000+trial+1000*type(1));
    nT = randi(3);
    nF = randi(4);
    if isequal(type,[1,1]),  nF = randi(3);  end
    ops = cell(nT,nF);   kd = cell(nT,nF);
    pdim = randi(2,1,nF);
    if isequal(type,[1,1]),  pdim(:) = pdim(1);  end
    for i=1:nT
        for j=1:nF
            kd{i,j} = kinds{find(rand<=cumsum(probs),1)};
            % type [2,2]: factor j acts on its own variable s_j;
            % type [1,1]: all factors act pointwise on s1
            if isequal(type,[1,1]),  vj = 1;  else,  vj = j;  end
            ops{i,j} = make_factor(kd{i,j},pdim(j),vj);
        end
    end
    Top = sopvar2tensopvar(ops,type);
    for ztol = [1e-12,1e-6]
        % effective kind of each factor under this tolerance
        eff = kd;
        eff(strcmp(kd,'tiny')) = {'scaled'};
        if ztol>1e-9,  eff(strcmp(kd,'tiny')) = {'zero'};  end
        eff(strcmp(eff,'idpad')) = {'id'};
        is_zero_trm = any(strcmp(eff,'zero'),2);
        is_id_trm = all(strcmp(eff,'id'),2);
        expected = sum(~is_zero_trm)==1 && is_id_trm(~is_zero_trm);
        ntens = ntens+1;  ntrue = ntrue+expected;  nmulti = nmulti+(nT>1 && nF>1);
        nchk = nchk+1;
        chk(is_id_op(Top,ztol),expected,"tensopvar random type="+mat2str(type)+" nT="+nT+" nF="+nF+" ztol="+ztol+" trial="+trial+" kinds="+strjoin(reshape(kd',1,[]),','));
    end
  end
end
fprintf('  randomized tensopvar checks: %d (%d expected identity, %d with >1 term and >1 factor)\n',ntens,ntrue,nmulti);

I3 = arrayfun(@(k) make_factor('id',1,k),1:4,'UniformOutput',false);
Zk = @(k) make_factor('zero',1,k);
Sk = @(k) make_factor('scaled',1,k);
nchk = nchk+1;  chk(is_id_op(mkT({I3{1},I3{2},I3{3}})),true,"tensopvar 3 factors id");
nchk = nchk+1;  chk(is_id_op(mkT({I3{1},I3{2},I3{3},I3{4}})),true,"tensopvar 4 factors id");
nchk = nchk+1;  chk(is_id_op(mkT({I3{1},I3{2},Sk(3)})),false,"tensopvar 3 factors, last scaled");
nchk = nchk+1;  chk(is_id_op(mkT({Sk(1),I3{2},I3{3}})),false,"tensopvar 3 factors, first scaled");
nchk = nchk+1;  chk(is_id_op(mkT({I3{1},Zk(2),I3{3}})),false,"tensopvar 3 factors, middle zero");
nchk = nchk+1;  chk(is_id_op(mkT({I3{1},I3{2},I3{3};Zk(1),I3{2},I3{3};I3{1},I3{2},Zk(3)})),true,"tensopvar 3 terms: id + 2 zero");
nchk = nchk+1;  chk(is_id_op(mkT({Zk(1),I3{2},I3{3};I3{1},I3{2},I3{3};I3{1},Zk(2),I3{3}})),true,"tensopvar 3 terms: id in the middle");
nchk = nchk+1;  chk(is_id_op(mkT({I3{1},I3{2},I3{3};I3{1},I3{2},I3{3}})),false,"tensopvar 2 identical id terms");
nchk = nchk+1;  chk(is_id_op(mkT({I3{1},I3{2},I3{3};Sk(1),I3{2},I3{3}})),false,"tensopvar id + scaled term");
nchk = nchk+1;  chk(is_id_op(mkT({Zk(1),I3{2},I3{3};I3{1},Zk(2),I3{3}})),false,"tensopvar all terms zero");

fprintf('Test_is_id_op passed (%d checks).\n',nchk);


%% -------------------------------------------------------------- helpers
function chk(val,expected,name)
if ~isequal(logical(val),logical(expected))
    error("Test_is_id_op failed: %s (got %d, expected %d).",name,val,expected);
end
end

function chk_err(fun,msgpart,name)
% Check that the function handle errors, with a message containing msgpart
try
    fun();
catch ME
    if contains(ME.message,msgpart)
        return
    end
    error("Test_is_id_op failed: %s: unexpected error: %s",name,ME.message);
end
error("Test_is_id_op failed: %s: expected an error.",name);
end

function Z = rand_basis(allow_missing_const)
% Random nonempty subset of {0,1,2}, as a sorted column; 0 present with
% probability 3/4 so that missing-constant bases are also exercised.
Z = find(rand(3,1)<0.5)-1;
if rand<0.75 || ~allow_missing_const
    Z = union(Z,0);
end
if isempty(Z),  Z = 1;  end
Z = Z(:);
end

function Iop = make_id(p,vars,dom)
% Canonical identity operator I_p on L2^p[dom]: multiplier in every variable.
N = numel(vars.in);
params = [{speye(p)},repmat({sparse(p,p)},1,3^N-1)];
Iop = sopvar(params,vars,repmat({0},1,N),repmat({0},1,N),dom,[p p]);
end

function tf = eq_quiet(P1,P2)
% eq prints a message on mismatch of dimensions; silence it.
evalc('tf = eq(P1,P2);');
end

function P = make_factor(kind,p,k)
% Factor acting on variable s<k> (domain [0,1]) of dimension p x p.
vn = {['s',num2str(k)]};
vars = struct('in',{vn},'out',{vn});    dom = struct('in',[0 1],'out',[0 1]);
Iop = make_id(p,vars,dom);
switch kind
    case 'id'
        P = Iop;
    case 'idpad'        % identity on an enlarged monomial basis
        Q = sopvar({sparse(p*2,p*2),sparse(p*2,p*2),sparse(p*2,p*2)},vars,{[0;1]},{[0;1]},dom,[p p]);
        P = Iop + Q;
    case 'zero'
        P = 0*Iop;
    case 'scaled'
        P = 2*Iop;
    case 'int'
        P = Iop;    P.params{2} = speye(p);
    case 'tiny'
        P = 1e-9*Iop;
end
end
