function test_sdopvar_tier2()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_SDOPVAR_TIER2 checks four changes of 09/29/2026 to 'sdopvar' against
% the class definition, never one routine against another (CLAUDE.md s4):
%
%   SD-2  private 'ChangeDecVar', reached through 'setdvars': new
%         coefficient Pi*B with Pi = T1' (spec sec. 8.1.4);
%   SD-3  'sync_basis' reads operand positions in the merged monomial
%         bases directly ('basis_positions');
%   SD-5  'plus' delegates to 'plus_batch';
%   SD-6  'plus_batch' skips a summand without entries before building
%         its zeros.
%
% ChangeDecVar: for random operators and random supersets / permutations
% of their decision list, every decision value is drawn once per NAME and
% read through each object's own Zd. Then C_gam(d) = A_gam + B_gam'*d must
% be unchanged for every parameter, and so must the kernels
% ('pi_sdopvar_kernels'); A is untouched, B is sparse (also from a full B),
% has one row per new name, and the LOC path gives the same C(d). A name
% missing from the new list, and a LOC of the wrong length, must error.
%
% plus / plus_batch: A+B, B+A and plus_batch(P1,...,PN) must equal the sum
% of the operands, decided two ways from the kernel definition, with the
% decision values drawn per NAME:
%   (i)  kernels of the result against the sum of the operands' kernels
%        ('pi_blk_kernels'), compared as polynomials. Valid because the
%        canonical multiplier form (sdopvar.m) makes the stored kernel
%        unique for a given operator;
%   (ii) the operators applied to a random polynomial test function at
%        random points by Gauss quadrature straight from the kernel
%        definition ('heatNd_apply'), which does not rely on (i)'s
%        uniqueness.
% Operand pairs: overlapping, disjoint and identical decision lists;
% sopvar+sdopvar in both orders; a summand without entries on a strict
% subset of the other operand's bases (so its basis map is not the
% identity), once with a sparse and once with a full running sum; N = 3..5
% mixed operands; and N operands on one list with 'shared_Zd'.
%
% Spatial variables: 1 (shared), 2 (shared), 3 (shared), and out = {s1,s2}
% against in = {s2,s3} (one shared, one output-only, one input-only).
%
% Requires 'heatNd_apply' (PIETOOLS_demos/sopvar_demos) and the helpers of
% this folder and of Testfolder/sdopvar ('rand_sdopvar', 'rand_sopvar').
%
% Initial coding MMP, 09/29/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(20260929);
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
assert(exist('heatNd_apply','file')==2,'heatNd_apply (PIETOOLS_demos/sopvar_demos) is not on the path');
TOL = 1e-10;    % relative; every check below is exact up to rounding
NQ  = 6;        % Gauss nodes per direction: exact to degree 11, and the
                % integrands have degree <= 6 per variable (bases <= 2,
                % test function <= 2, one integration)
nchk = 0;   nfail = 0;

cfgs = { struct('lbl','1 var',  'out',{{'s1'}},          'in',{{'s1'}}), ...
         struct('lbl','2 vars', 'out',{{'s1','s2'}},     'in',{{'s1','s2'}}), ...
         struct('lbl','3 vars', 'out',{{'s1','s2','s3'}},'in',{{'s1','s2','s3'}}), ...
         struct('lbl','mixed',  'out',{{'s1','s2'}},     'in',{{'s2','s3'}}) };

% ================================================== SD-2: ChangeDecVar
nfullB = 0;
for c = 1:numel(cfgs)
    cf = cfgs{c};   vars = struct('in',{cf.in},'out',{cf.out});
    for rep = 1:6
        dims = randi(2,1,2);    nd = randi([3 12]);
        P0 = rand_sdopvar(dims,vars,[0 1],randi(2),nd,0.5);
        Zd = compose('a%d',randperm(3*nd,nd).');         % cellstr column
        if rep==2,  Zd = string(Zd).';  end               % string row
        P = with_dvars(P0,Zd,rep==3);                     % rep 3: full B
        if rep==3
            nfullB = nfullB + ~issparse(P.params.B{1});
        end
        Zo = names(P.Zd);
        switch rep
            case 4,     Zn = Zo(randperm(nd));            % permutation only
            case 5,     Zn = Zo;                          % the same list
            otherwise
                Zn = [Zo; compose('e%d',(1:randi([1 6])).')];
                Zn = Zn(randperm(numel(Zn)));
        end
        dmap = draw(containers.Map('KeyType','char','ValueType','double'),Zn);
        lbl = sprintf('ChangeDecVar %s rep %d',cf.lbl,rep);
        Q  = setdvars(P,Zn);
        [~,loc] = ismember(Zo,Zn);
        Q2 = setdvars(P,Zn,loc);
        Q3 = setdvars(P,Zn,int32(loc));     % an integer LOC, as HEAD accepted
        for Qx = {Q,Q2,Q3}
            R = Qx{1};
            ck(isequal(names(R.Zd),Zn(:)),[lbl ': Zd is the new list']);
            ck(isequal(R.params.A,P.params.A),[lbl ': A untouched']);
            ck(all(cellfun(@issparse,R.params.B(:))),[lbl ': B sparse']);
            ck(all(cellfun(@(b) size(b,1),R.params.B(:))==numel(Zn)),[lbl ': one B row per name']);
            ck(cdiff(R,P,dmap)<=TOL,sprintf('%s: C(d) by name, rel err %.2e',lbl,cdiff(R,P,dmap)));
        end
        if rep<=3
            e = kdiff({pkern(Q,dmap)},{pkern(P,dmap)});
            ck(e<=TOL,sprintf('%s: kernels, rel err %.2e',lbl,e));
        end
        % a name missing from the new list, and a LOC for another list
        Zbad = Zn(~strcmp(Zn,Zo{1}));
        ck(threw(@() setdvars(P,Zbad)),[lbl ': missing name errors']);
        ck(strcmp(errid(@() setdvars(P,[Zn;{'zz'}],loc(1:end-1))),'ChangeDecVar:badLoc'), ...
           [lbl ': short LOC errors']);
    end
end
ck(nfullB>0,'a full B reached ChangeDecVar (constructor kept it full)');

% no old decision variables: all-zero B rows on the new list
P0 = rand_sdopvar([1 2],struct('in',{{'s1','s2'}},'out',{{'s1','s2'}}),[0 1],1,2,0.5);
prm = P0.params;
for k = 1:numel(prm.B),     prm.B{k} = sparse(0,numel(prm.A{k}));    end
P = sdopvar(prm,P0.vars,cell(0,1),P0.ZL,P0.ZR,P0.dom,P0.dims);
Q = setdvars(P,{'z1';'z2';'z3'});
ck(all(cellfun(@(b) isequal(size(b),[3 numel(b)/3]) && issparse(b) && nnz(b)==0, ...
       Q.params.B(:))),'ChangeDecVar from an empty list: zero sparse B, 3 rows');
ck(isequal(Q.params.A,P.params.A),'ChangeDecVar from an empty list: A untouched');

% ============================================ SD-3/5/6: plus, plus_batch
nfullsum = 0;
for c = 1:numel(cfgs)
    cf = cfgs{c};   vars = struct('in',{cf.in},'out',{cf.out});
    dims = randi(2,1,2);
    for rep = 1:2
        lblc = sprintf('%s rep %d',cf.lbl,rep);
        na = randi([3 8]);  nb = randi([3 8]);
        ZA = compose('p%d',randperm(20,na).');
        ZB = compose('p%d',randperm(20,nb).');                    % overlap
        ZC = compose('q%d',(1:nb).');                              % disjoint
        A  = with_dvars(rand_sdopvar(dims,vars,[0 1],2,na,0.5),ZA,false);
        B  = with_dvars(rand_sdopvar(dims,vars,[0 1],2,nb,0.5),ZB,false);
        Bc = with_dvars(rand_sdopvar(dims,vars,[0 1],2,nb,0.5),ZC,false);
        Bs = with_dvars(rand_sdopvar(dims,vars,[0 1],2,na,0.5),A.Zd,false);   % same list
        S  = rand_sopvar(dims,vars,[0 1],2,0.5);
        E  = zero_sub(A,compose('r%d',(1:3).'),false);            % no entries
        Af = full_A(A);                                            % full A
        Ef = zero_sub(Af,compose('r%d',(1:3).'),false);
        ck(~issparse(Af.params.A{1}),[lblc ': full A kept by the constructor']);
        dmap = containers.Map('KeyType','char','ValueType','double');
        dmap = draw(dmap,[ZA;ZB;ZC;{'r1';'r2';'r3'}]);
        pts  = rand(3,numel(A.vars.out));
        f    = rand_testfun(A.dims(2),numel(A.vars.in));
        pairs = { {A,B,'overlapping lists'}, {A,Bc,'disjoint lists'}, ...
                  {A,Bs,'identical lists'}, {A,S,'sdopvar+sopvar'}, ...
                  {S,A,'sopvar+sdopvar'}, {A,E,'empty summand, sparse sum'}, ...
                  {E,A,'empty summand first'}, {Af,Ef,'empty summand, full sum'} };
        for p = 1:numel(pairs)
            X = pairs{p}{1};    Y = pairs{p}{2};
            lbl = sprintf('%s, %s',lblc,pairs{p}{3});
            chk_sum(X+Y,{X,Y},dmap,f,pts,[lbl ': X+Y']);
            chk_sum(Y+X,{X,Y},dmap,f,pts,[lbl ': Y+X']);
            chk_sum(plus_batch(X,Y),{X,Y},dmap,f,pts,[lbl ': plus_batch(X,Y)']);
        end
        nfullsum = nfullsum + 1;
        % N = 3..5 mixed operands, in a random order
        ops = {A,B,Bc,S,E};
        ops = ops(randperm(5,randi([3 5])));
        if ~any(cellfun(@(x) isa(x,'sdopvar'),ops)),   ops{end+1} = A;     end %#ok<AGROW>
        chk_sum(plus_batch(ops{:}),ops,dmap,f,pts,sprintf('%s: plus_batch of %d mixed',lblc,numel(ops)));
        % N operands built on one list, with the caller's 'shared_Zd'
        ops = {A,Bs,with_dvars(rand_sdopvar(dims,vars,[0 1],1,na,0.5),A.Zd,false),E};
        ops{4} = zero_sub(A,A.Zd,false);
        chk_sum(plus_batch(ops{:},'shared_Zd'),ops,dmap,f,pts,[lblc ': plus_batch shared_Zd']);
        % incompatible summands must not return
        D = with_dvars(rand_sdopvar(dims+[1 0],vars,[0 1],1,na,0.5),ZA,false);
        ck(threw(@() A+D),[lblc ': dims mismatch errors']);
        ck(threw(@() A+3),[lblc ': numeric summand errors']);
        ck(threw(@() plus(A,'shared_Zd')) && threw(@() A+"x"),[lblc ': text summand errors']);
    end
end
ck(nfullsum>0,'full running sum case exercised');

fprintf('test_sdopvar_tier2: %d checks, %d failed\n',nchk,nfail);
assert(nfail==0,'test_sdopvar_tier2: %d of %d checks failed',nfail,nchk);


    % ------------------------------------------------------------ nested
    function ck(ok,msg)
        nchk = nchk + 1;
        if ~ok
            nfail = nfail + 1;
            fprintf('FAIL: %s\n',msg);
        end
    end

    function chk_sum(R,ops,dmap,f,pts,lbl)
    % R must be an sdopvar equal to ops{1}+...+ops{N}: (i) kernels, (ii)
    % quadrature, both at the decision values of dmap (by NAME).
        ck(isa(R,'sdopvar'),[lbl ': class sdopvar']);
        if ~isa(R,'sdopvar'),   return,     end
        zs = {};
        for kk = 1:numel(ops)
            if isa(ops{kk},'sdopvar'),  zs = [zs; names(ops{kk}.Zd)]; end %#ok<AGROW>
        end
        ck(isempty(setxor(names(R.Zd),zs)),[lbl ': Zd is the union of the lists']);
        ck(numel(unique(names(R.Zd)))==numel(R.Zd),[lbl ': Zd has no repeated name']);
        Ks = cellfun(@(x) pkern(x,dmap),ops,'UniformOutput',false);
        e1 = kdiff({pkern(R,dmap)},{ksum(Ks)});
        ck(e1<=TOL,sprintf('%s: kernels, rel err %.2e',lbl,e1));
        y = heatNd_apply(fixblock(R,dmap),f,pts,NQ);
        ys = zeros(size(y));
        for kk = 1:numel(ops)
            ys = ys + heatNd_apply(fixblock(ops{kk},dmap),f,pts,NQ);
        end
        e2 = max(abs(y(:)-ys(:)))/max(1,max(abs(ys(:))));
        ck(e2<=TOL,sprintf('%s: quadrature, rel err %.2e',lbl,e2));
    end
end


% ================================================================ helpers
function z = names(Zd)
% Decision names as a cellstr column, whatever the stored orientation/type.
z = cellstr(string(Zd(:)));
end

function dmap = draw(dmap,nm)
% One random value per decision variable NAME, kept across calls.
for i = 1:numel(nm)
    if ~isKey(dmap,nm{i}),  dmap(nm{i}) = randn;    end
end
end

function d = dvals(P,dmap)
% Values of P's decision variables, read by name from P's own list.
z = names(P.Zd);    d = zeros(numel(z),1);
for i = 1:numel(z),     d(i) = dmap(z{i});      end
end

function e = cdiff(R,P,dmap)
% max over parameters of |C_R(d) - C_P(d)| / max(1,|C_P(d)|), where
% C(d) = A + B'*d (spec sec. 8.1.1) with each object's own list.
dR = dvals(R,dmap);     dP = dvals(P,dmap);    e = 0;
for k = 1:numel(P.params.A)
    cR = full(R.params.A{k}(:) + R.params.B{k}.'*dR);
    cP = full(P.params.A{k}(:) + P.params.B{k}.'*dP);
    e = max(e, max(abs(cR-cP),[],'all')/max(1,max(abs(cP),[],'all')));
end
end

function K = pkern(B,dmap)
% Kernels of an sopvar/sdopvar, a decision operator at the drawn values.
if isa(B,'sdopvar'),    K = pi_blk_kernels(B,dvals(B,dmap));
else,                   K = pi_blk_kernels(B,[]);
end
end

function K = ksum(Ks)
K = Ks{1};
for i = 2:numel(Ks)
    for k = 1:numel(K),     K{k} = K{k} + Ks{i}{k};     end
end
end

function e = kdiff(K1,K2)
% Largest coefficient of K1 - K2 over all cells, relative to K2.
K1 = K1{1};     K2 = K2{1};     e = 0;
for k = 1:numel(K2)
    e = max(e, maxcoef(K1{k}-K2{k})/max(1,maxcoef(K2{k})));
end
end

function m = maxcoef(X)
if isnumeric(X)
    m = full(max(abs(X(:))));
else
    X = polynomial(X);  cf = X.coefficient;
    m = full(max(abs(cf(:))));
end
if isempty(m),  m = 0;  end
end

function P = with_dvars(P0,Zd,mkfull)
% P0's operator coefficients on the decision list Zd (numel(P0.Zd) names),
% optionally with every B stored full.
prm = P0.params;
if mkfull
    for k = 1:numel(prm.B),     prm.B{k} = full(prm.B{k});     end
end
P = sdopvar(prm,P0.vars,Zd,P0.ZL,P0.ZR,P0.dom,P0.dims);
end

function P = full_A(P0)
% P0 with every A stored full.
prm = P0.params;
for k = 1:numel(prm.A),     prm.A{k} = full(prm.A{k});     end
P = sdopvar(prm,P0.vars,P0.Zd,P0.ZL,P0.ZR,P0.dom,P0.dims);
end

function E = zero_sub(P0,Zd,~)
% An operator without entries on the list Zd, whose bases are a strict
% subset of P0's in every direction that has more than one monomial, so
% that its map onto the merged bases is not the identity.
ZL = P0.ZL;     ZR = P0.ZR;
for i = 1:numel(ZL),    if numel(ZL{i})>1, ZL{i} = ZL{i}(1:end-1); end,  end
for i = 1:numel(ZR),    if numel(ZR{i})>1, ZR{i} = ZR{i}(2:end);   end,  end
nC = P0.dims(1)*prod([cellfun(@numel,ZL),1])*P0.dims(2)*prod([cellfun(@numel,ZR),1]);
prm.A = cell(size(P0.params.A));    prm.B = cell(size(P0.params.B));
for k = 1:numel(prm.A)
    prm.A{k} = sparse(nC,1);    prm.B{k} = sparse(numel(Zd),nC);
end
E = sdopvar(prm,P0.vars,Zd,ZL,ZR,P0.dom,P0.dims);
end

function f = rand_testfun(p,nin)
% Random polynomial x: R^nin -> R^p of degree <= 2 in each variable.
c0 = randn(1,p);    c1 = randn(nin,p);    c2 = randn(nin,p);    c3 = randn(1,p);
f = @(X) ones(size(X,1),1)*c0 + X*c1 + (X.^2)*c2 + prod(X,2)*c3;
end

function B = fixblock(B,dmap)
% An sdopvar at the decision values of dmap, as an sopvar with
% C_gam = unvec(A_gam + B_gam'*d); an sopvar is returned as it is.
if ~isa(B,'sdopvar'),   return,     end
d = dvals(B,dmap);
nr = B.dims(1)*prod([cellfun(@numel,B.ZL),1]);
nc = B.dims(2)*prod([cellfun(@numel,B.ZR),1]);
prm = cell(size(B.params.A));
for k = 1:numel(prm)
    v = zeros(nr*nc,1);
    Ak = B.params.A{k};     Bk = B.params.B{k};
    if numel(Ak)==nr*nc,    v = v + full(Ak(:));    end
    if size(Bk,2)==nr*nc,   v = v + full(Bk.'*d);   end
    prm{k} = reshape(v,nr,nc);
end
B = sopvar(prm,B.vars,B.ZL,B.ZR,B.dom,B.dims);
end

function tf = threw(f)
try,    f();    tf = false;
catch,  tf = true;
end
end

function id = errid(f)
try,    f();    id = '';
catch ME,       id = ME.identifier;
end
end
