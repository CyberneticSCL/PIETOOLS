function R = test_heatNd_poincare()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = TEST_HEATND_POINCARE() asserting tests of the Poincare diagnostic
% (HEATND_POINCARE_OPS, HEATND_POINCARE), each object against its own
% definition (CLAUDE.md sec. 4); N = 1, 2 for anything with a solve, N = 3
% for the SDP-free parts only:
%
%  (1) I0 passes for N = 1, 2, 3 (DD, DN, ND; unit and non-unit boxes);
%  (2) negative controls - I0 must FAIL, and in the layer expected:
%      (a) BC labels swapped (the T kernels are DD x DN, the labels say DN
%          x DD): kernel check (k) and identity (c) fail;
%      (b) lower/upper cells swapped in D_1T's direction 1: (k), (c), (q);
%      (c) sign of D_2T flipped (DN): (k) and (q) fail, (c) PASSES - the
%          identity T*A_i = -(D_iT)*(D_iT) is blind to the sign, which is
%          why I0 has three layers;
%  (3) SDP shapes against the reference pricing formulas (c = (d+1)(5d+3),
%      c_p = (d+1)(5d+7)): L2 none m = 2^(N-1) c^N, product 2^(N-1) c_p^N,
%      faces 2^(N-1) [c^N + 2N(d+1) c^(N-1)], Gram (2(d+1)^2)^N; energy
%      faces 1-D m 24 [10], 2-D full 920 [96], per direction 844 [80]
%      (stock, measured); At identical at 3 lam and b affine (<= 1e-12);
%      the multi-indices equal 'eq_opts_sopvar''s support rule;
%  (4) 'op' and 'grad' targets give the same SDP (At equal, b to 1e-12);
%  (5) values, 'objective' mode (one SDP + a fixed-lam re-certification):
%      L2 faces d = 1 per direction 1.0000 (+-2e-4) in 1-D and 2-D, full 2-D
%      2.0000 (the sum rule of the tensor bound B2); energy faces d = 1 1-D
%      DD equals 2-D direction 1 (bounds B1 + B2: equality) to 1e-5;
%  (6) soundness: at 1.001 x exact no verdict is +1 (L2 faces d = 2, energy
%      faces d = 1; DD and DN);
%  (7) coverage: tensor d = 0 misses the DD monomial s th (reference):
%      uncovered rows > 0 and a lam > 0 is -1 'deficient';
%  (8) bisection reuse (HEATND_BISECT through the family's solve hook): the
%      certified bracket of 1-D DN L2 faces d = 1 contains 1.0 within 2e-4;
%  (9) a certified point is a certificate SEMANTICALLY: at the returned
%      decision values, G = Z*MZ + sum_g Z*gM_gZ (HEATND_POINCARE aux,
%      substituted from the class definition, sdopvar.m) matches Q(lam) in
%      every coefficient; on the tensor Legendre space v (degree <= P per
%      direction, P = 4 in 1-D, 3 in 2-D) its Galerkin matrix <v_a,Gv_b>
%      (quadrature, HEATND_APPLY) equals that of the class target and of
%      the ANALYTIC form (||d_sel u||^2 - lam||u||^2, or ||Lap_sel u||^2 -
%      lam||d_sel u||^2, u = Tv) - no SDP data in that comparison; and the
%      Rayleigh-Ritz minimum of <v,Gv>/nrm(v) over the whole space (nrm =
%      ||u||^2 L2, ||d_sel u||^2 energy; lam units) is > 0 and equals the
%      analytic one to 1e-4 mu, while the analytic form at 1.001 mu is
%      NEGATIVE there (sensitivity control: the space would expose a false
%      certificate 1e-3 above mu).
%
% Initial coding MMP, 09/27/2026
% MMP, 09/27/2026 (final reviews): (9) sign test by the Rayleigh-Ritz
%   minimum over the tensor polynomial space, with the 1.001 mu control (it
%   was 6 functions: the lowest mode and 5 random polynomials).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info = heatNd_path();
fprintf('test_heatNd_poincare: maxNumCompThreads = %d (MOSEK uses its own count)\n',info.threads);
rng(20260927);
R = struct('check',{},'val',{},'ok',{});
% % % (1) I0.
cases = { {1,{'DD'},[0 1]}, {1,{'DN'},[0 1]}, {1,{'ND'},[-1 2]}, {2,{'DD','DN'},[0 1;0 1]}, ...
          {2,{'ND','DD'},[0 2;-1 1]}, {3,{'DD','DN','DN'},[0 1;0 1;0 1]}, {3,{'DN','ND','DD'},[0 1;1 2;0 3]} };
for c = 1:numel(cases)
    [N,bc,dom] = cases{c}{:};
    O = heatNd_poincare_ops(heatNd_pie(N,3.7,bc,dom),true);
    e = max([O.I0.err]./[O.I0.tol]);
    R(end+1) = chk(sprintf('(1) I0 N=%d %s %s: %d checks, max err/tol %.1e',N,strjoin(bc,'x'),mat2str(dom), ...
                   numel(O.I0),e),e,O.I0ok);                                    %#ok<AGROW>
end
% % % (2) negative controls.
pie = heatNd_pie(2,0);
pw = pie;   pw.bc = {'DN','DD'};                        % (a) wrong labels
Oa = heatNd_poincare_ops(pw,true);
R(end+1) = chk('(2a) swapped BC labels: I0 fails in (k) and (c)',0, ...
               ~Oa.I0ok && failed(Oa,'(k) T kernels') && failed(Oa,'(c) T*A_1 + '));  %#ok<AGROW>
O = heatNd_poincare_ops(pie);
Ob = O;     B = Ob.DT{1}.C{1,1};    prm = B.params;     prm = prm([1 3 2],:);  % (b) swap in s1
Ob.DT{1} = copvar({sopvar(prm,B.vars,B.ZL,B.ZR,B.dom,B.dims)});
Ob = heatNd_poincare_ops(pie,true,Ob);
R(end+1) = chk('(2b) D_1T lower/upper swapped: I0 fails in (k), (c) and (q)',0, ...
               ~Ob.I0ok && failed(Ob,'(k) D_1T') && failed(Ob,'(c) T*A_1 + ') && failed(Ob,'(q) D_1T v'));  %#ok<AGROW>
Oc = O;     Oc.DT{2} = -1*Oc.DT{2};                     % (c) sign of the DN factor
Oc = heatNd_poincare_ops(pie,true,Oc);
R(end+1) = chk('(2c) D_2T sign flipped: (k), (q) fail, (c) passes (sign-blind)',0, ...
               ~Oc.I0ok && failed(Oc,'(k) D_2T') && failed(Oc,'(q) D_2T v') && ~failed(Oc,'(c) T*A_2 + '));  %#ok<AGROW>
% % % (3) shapes, affinity, multi-indices.
nof = struct('verbose',false);
for N = 1:2
    O = heatNd_poincare_ops(heatNd_pie(N,0));
    for d = 1:2
        c = (d+1)*(5*d+3);  cp = (d+1)*(5*d+7);     nG = (2*(d+1)^2)^N;
        mx = {'none',2^(N-1)*c^N,nG; 'product',2^(N-1)*cp^N,[nG nG]; ...
              'faces',2^(N-1)*(c^N + 2*N*(d+1)*c^(N-1)),nG*ones(1,1+2*N)};
        for g = 1:3
            Rg = heatNd_poincare(O,struct('form','L2','d',d,'gen',mx{g,1}),nof);
            ok = Rg.shape.m==mx{g,2} && isequal(Rg.shape.Ks,mx{g,3}) && Rg.affine.At_same && ...
                 Rg.affine.b_err<=1e-12 && Rg.include_eqopts && Rg.shape.uncovered==0;
            R(end+1) = chk(sprintf('(3) N=%d L2 d=%d %s: m %d (formula %d), Ks %s, b affine %.0e, eq_opts %d', ...
                           N,d,mx{g,1},Rg.shape.m,mx{g,2},mat2str(Rg.shape.Ks),Rg.affine.b_err,Rg.include_eqopts),Rg.shape.m,ok); %#ok<AGROW>
        end
    end
    sels = [{'full'}, num2cell(1:N)];
    for k = 1:numel(sels)
        Rg = heatNd_poincare(O,struct('form','energy','sel',sels{k},'d',1,'gen','faces'),nof);
        if N==1,                        mk = 24;    gk = 10;
        elseif strcmp(num2str(sels{k}),'full'),  mk = 920;   gk = 96;
        else,                           mk = 844;   gk = 80;
        end
        ok = Rg.shape.m==mk && isequal(Rg.shape.Ks,gk*ones(1,1+2*N)) && Rg.affine.At_same && ...
             Rg.affine.b_err<=1e-12 && Rg.include_eqopts && Rg.shape.uncovered==0;
        R(end+1) = chk(sprintf('(3) N=%d energy sel=%s d=1 faces: m %d (stock %d), Ks %s',N,num2str(sels{k}), ...
                       Rg.shape.m,mk,mat2str(Rg.shape.Ks)),Rg.shape.m,ok);        %#ok<AGROW>
    end
end
% % % (4) 'op' = 'grad'.
O2 = heatNd_poincare_ops(heatNd_pie(2,0));
for f = {'L2','energy'}
    sp = struct('form',f{1},'sel',1,'d',1,'gen','faces');
    Rg = heatNd_poincare(O2,sp,struct('verbose',false,'keepF',true));
    sp.target = 'op';   Ro = heatNd_poincare(O2,sp,struct('verbose',false,'keepF',true));
    e = max(norm(Rg.F.b1-Ro.F.b1)/norm(Rg.F.b1),norm(Rg.F.db-Ro.F.db)/norm(Rg.F.db));
    R(end+1) = chk(sprintf('(4) 2-D %s dir 1: op vs grad, At equal %d, b rel diff %.1e',f{1}, ...
                   isequal(Rg.F.At,Ro.F.At),e),e,isequal(Rg.F.At,Ro.F.At) && e<=1e-12);  %#ok<AGROW>
end
% % % (5) values by the objective mode.
ob = struct('mode','objective','verbose',false);
O1 = {heatNd_poincare_ops(heatNd_pie(1,0,{'DD'})), heatNd_poincare_ops(heatNd_pie(1,0,{'DN'}))};
for b = 1:2
    Ro = heatNd_poincare(O1{b},struct('form','L2','d',1,'gen','faces'),ob);
    R(end+1) = chk(sprintf('(5) 1-D %s L2 faces d=1: lam-hat %.6f (1.0000), recert %+d',O1{b}.bc{1}, ...
                   Ro.obj.lamhat,Ro.obj.recert.st),Ro.obj.lamhat,abs(Ro.obj.lamhat-1)<=2e-4 && Ro.obj.recert.st==1); %#ok<AGROW>
end
sels = {'full',1,2};    pred = [2 1 1];
for k = 1:3
    Ro = heatNd_poincare(O2,struct('form','L2','sel',sels{k},'d',1,'gen','faces'),ob);
    R(end+1) = chk(sprintf('(5) 2-D L2 faces d=1 sel=%s: lam-hat %.6f (%g), recert %+d',num2str(sels{k}), ...
                   Ro.obj.lamhat,pred(k),Ro.obj.recert.st),Ro.obj.lamhat, ...
                   abs(Ro.obj.lamhat-pred(k))<=2e-4*pred(k) && Ro.obj.recert.st==1);  %#ok<AGROW>
end
E1 = heatNd_poincare(O1{1},struct('form','energy','d',1,'gen','faces'),ob);
E2 = heatNd_poincare(O2,struct('form','energy','sel',1,'d',1,'gen','faces'),ob);
e = abs(E1.obj.lamhat-E2.obj.lamhat)/E1.obj.lamhat;
R(end+1) = chk(sprintf('(5) energy faces d=1 DD: 1-D %.7f = 2-D dir 1 %.7f (rel %.1e), recert %+d %+d', ...
               E1.obj.lamhat,E2.obj.lamhat,e,E1.obj.recert.st,E2.obj.recert.st),e, ...
               e<=1e-5 && E1.obj.recert.st==1 && E2.obj.recert.st==1);           %#ok<AGROW>
% % % (6) soundness at 1.001 x exact.
to = struct('mode','test','verbose',false,'bisect',struct('retry',{{}},'verbose',false));
for b = 1:2
    for f = {{'L2',2},{'energy',1}}
        mu = O1{b}.mu;  to.lams = 1.001*mu;
        Rt = heatNd_poincare(O1{b},struct('form',f{1}{1},'d',f{1}{2},'gen','faces'),to);
        st = Rt.test.trace(end,2);
        R(end+1) = chk(sprintf('(6) 1-D %s %s faces d=%d at 1.001 mu: st %+d (never +1)',O1{b}.bc{1}, ...
                       f{1}{1},f{1}{2},st),st,st~=1);                            %#ok<AGROW>
    end
end
% % % (7) coverage at tensor d = 0 (DD).
Rc = heatNd_poincare(O1{1},struct('form','L2','d',0,'gen','faces'),struct('mode','test','lams',1, ...
                     'verbose',false,'bisect',struct('retry',{{}},'verbose',false)));
R(end+1) = chk(sprintf('(7) 1-D DD L2 d=0: uncovered rows %d, st at lam 1 %+d (%s)',Rc.shape.uncovered, ...
               Rc.test.trace(end,2),Rc.test.why{end}),Rc.shape.uncovered, ...
               Rc.shape.uncovered>0 && Rc.test.trace(end,2)==-1 && contains(Rc.test.why{end},'deficient'));  %#ok<AGROW>
% % % (8) bisection through HEATND_BISECT.
Rb = heatNd_poincare(O1{2},struct('form','L2','d',1,'gen','faces'),struct('mode','bisect','verbose',false, ...
                     'bisect',struct('rtol',1e-4,'verbose',false,'kappa0',[0.5 1.5])));
R(end+1) = chk(sprintf('(8) 1-D DN L2 faces d=1 bisect: [%.6f, %.6f) (%s)',Rb.B.lo,Rb.B.hi,Rb.B.stop), ...
               [Rb.B.lo Rb.B.hi],Rb.B.lo>=1-2e-4 && Rb.B.lo<=1+2e-4 && Rb.B.hi>Rb.B.lo && Rb.B.hi<=1.01); %#ok<AGROW>
% % % (9) semantic certificates.
sc = { {O1{1},struct('form','L2','sel','full','d',1,'gen','faces'),0.9}, ...
       {O2,struct('form','energy','sel',1,'d',1,'gen','faces'),9.0}, ...
       {O2,struct('form','L2','sel','full','d',1,'gen','faces'),1.8} };
for k = 1:numel(sc)
    [Ok,sp,lam] = sc{k}{:};
    Rs = heatNd_poincare(Ok,sp,struct('mode','sdp','lams',lam,'verbose',false,'keepF',true));
    v = heatNd_solve(Rs.D,struct('tight',true,'keepx',true));
    [res,qerr,aerr,rG,rA,rC] = semantic(Ok,Rs,v,lam,sp);
    R(end+1) = chk(sprintf(['(9) N=%d %s sel=%s lam=%g: st %+d, coef res %.1e, Gal G vs Q %.1e, vs analytic %.1e; ' ...
                   'Ritz min G %.6g, analytic %.6g, at 1.001 mu %.2e'],Ok.N,sp.form,num2str(sp.sel),lam,v.st, ...
                   res,qerr,aerr,rG,rA,rC),[res qerr aerr rG rA rC], ...
                   v.st==1 && res<=1e-6 && qerr<=1e-6 && aerr<=1e-6 && rG>0 && abs(rG-rA)<=1e-4*Rs.mu && rC<0);  %#ok<AGROW>
end
ok = true;
fprintf('\n');
for i = 1:numel(R)
    ok = ok && R(i).ok;
    fprintf('%-120s %s\n',R(i).check,char(string(R(i).ok)));
end
assert(ok,'test_heatNd_poincare: a check failed.');
fprintf('test_heatNd_poincare: all %d checks pass.\n',numel(R));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = chk(c,v,ok),   s = struct('check',c,'val',v,'ok',logical(ok));   end

function tf = failed(O,prefix)
% True if some I0 row whose name starts with PREFIX failed.
k = startsWith({O.I0.check},prefix);
tf = any(k) && any(~[O.I0(k).ok]);
end

function [res,qerr,aerr,rG,rA,rC] = semantic(O,Rs,v,lam,sp)
% Certificate G at the solve's decision values, from the class definition;
% Galerkin matrices on the tensor Legendre space (header (9)) by quadrature;
% Rayleigh-Ritz minima of G, of the analytic form at lam, and at 1.001 mu.
A = Rs.aux;     dv = full(A.RR*v.x)*Rs.D.bscl;          % decvartable order
nm = cellstr(string(A.prog.decvartable));   nm = nm(1:numel(dv));
Qt = A.Q0 + lam*A.Q1;
res = coefres(A.G - Qt,nm,dv);
Gv = subsd(A.G,nm,dv);
% Degree <= P per direction. At d = 1 (Z: a, c <= 1, face weight degree 1)
% G's kernels have joint degree <= 6 per direction (one integration), so
% every integrand has degree <= 2P+7 per direction <= 2nq-1 with nq = 2P+2:
% Gauss-exact (DERIVED; the (9) cases are all d = 1).
N = O.N;    dom = O.dom;    P = 5-N;    nq = 2*P+2;     [Xq,wq] = qgrid(dom,nq);
if ischar(sp.sel),  I = 1:N;    else,   I = sp.sel;     end
E = dec2base(0:(P+1)^N-1,P+1,N)-'0';    nb = size(E,1);  % multi-index per row
nx = numel(wq);     Vq = zeros(nx,nb);  GV = Vq;    QV = Vq;    UV = Vq;    AV = Vq;
DV = zeros(nx,nb,numel(I));
for b = 1:nb
    f = @(X) legt(X,E(b,:),dom);
    Vq(:,b) = f(Xq);
    GV(:,b) = heatNd_apply(Gv,f,Xq,nq);     QV(:,b) = heatNd_apply(Qt,f,Xq,nq);
    % Analytic pieces: d_i u = D_iT v, d_i^2 u = A_i v, u = T v.
    for k = 1:numel(I)
        DV(:,b,k) = heatNd_apply(O.DT{I(k)},f,Xq,nq);
        if strcmp(sp.form,'energy'),    AV(:,b) = AV(:,b) + heatNd_apply(O.A{I(k)},f,Xq,nq);   end
    end
    if strcmp(sp.form,'L2'),    UV(:,b) = heatNd_apply(O.T,f,Xq,nq);   end
end
ip = @(X,Y) X'*(wq.*Y);                                 % <x_a, y_b>
MG = ip(Vq,GV);     MQ = ip(Vq,QV);
GG = 0;     for k = 1:numel(I),     GG = GG + ip(DV(:,:,k),DV(:,:,k));  end   % ||d_sel u||^2
if strcmp(sp.form,'L2'),    Nr = ip(UV,UV);     Qa = @(l) GG - l*Nr;           % nrm ||u||^2
else,                       Nr = GG;            AA = ip(AV,AV);     Qa = @(l) AA - l*GG;
end
MA = Qa(lam);
qerr = max(abs(MG(:)-MQ(:)))/max(abs(MQ(:)));   aerr = max(abs(MG(:)-MA(:)))/max(abs(MA(:)));
rG = ritz(MG,Nr);   rA = ritz(MA,Nr);   rC = ritz(Qa(1.001*Rs.mu),Nr);
end

function y = legt(X,e,dom)
% Tensor Legendre polynomial prod_i P_{e_i}, direction i mapped to [-1,1].
y = ones(size(X,1),1);
for i = 1:numel(e)
    x = 2*(X(:,i)-dom(i,1))/(dom(i,2)-dom(i,1)) - 1;
    p0 = ones(size(x));     p1 = x;
    if e(i)==0,     p1 = p0;    end
    for n = 1:e(i)-1                                    % Bonnet recurrence
        p2 = ((2*n+1)*x.*p1 - n*p0)/(n+1);  p0 = p1;    p1 = p2;
    end
    y = y.*p1;
end
end

function r = ritz(A,B)
% min <v,Av>/<v,Bv> over the space (B > 0 on it): least eigenvalue of the
% symmetric parts after whitening B.
A = (A+A')/2;   B = (B+B')/2;
[V,d] = eig(B,'vector');    k = d>1e-13*max(d);
W = V(:,k)./sqrt(d(k))';    S = W'*A*W;
r = min(eig((S+S')/2));
end

function amax = coefres(X,names,vals)
% Largest kernel coefficient of the 1x1 'cdopvar' X at the decision values,
% relative to the no-cancellation scale max(|A| + |B||d|) (as TEST_HEATND_LPI).
B = X.C{1};     zd = cellstr(string(B.Zd));
[tf,loc] = ismember(zd(:),names(:));    d = zeros(numel(zd),1);    d(tf) = vals(loc(tf));
amax = 0;   scale = 0;
for g = 1:numel(B.params.A)
    a = full(B.params.A{g}(:));     Bg = B.params.B{g};
    if isempty(Bg),     bd = zeros(size(a));    ad = bd;
    else,               bd = full(Bg'*d);       ad = full(abs(Bg)'*abs(d));
    end
    if isempty(a),      a = zeros(size(bd));    end
    if isempty(a),      continue,   end
    amax = max(amax,max(abs(a+bd)));    scale = max([scale;abs(a)+ad]);
end
amax = amax/max(scale,eps);
end

function Y = subsd(X,names,vals)
% Fixed operator X(d) of a 1x1 'cdopvar' (params = unvec(A + B'd), sdopvar.m).
B = X.C{1};     zd = cellstr(string(B.Zd));
[tf,loc] = ismember(zd(:),names(:));    d = zeros(numel(zd),1);    d(tf) = vals(loc(tf));
nL = prod(cellfun(@numel,B.ZL));    nR = prod(cellfun(@numel,B.ZR));
m = B.dims(1)*nL;   n = B.dims(2)*nR;
prm = cell(size(B.params.A));
for g = 1:numel(prm)
    a = B.params.A{g};  Bg = B.params.B{g};
    if isempty(a),  a = sparse(m*n,1);  end
    if ~isempty(Bg),    a = a + Bg'*d;  end
    prm{g} = reshape(full(a),m,n);
end
Y = copvar({sopvar(prm,B.vars,B.ZL,B.ZR,B.dom,B.dims)});
end

function [X,w] = qgrid(dom,n)
b = (1:n-1)./sqrt(4*(1:n-1).^2-1);
[V,L] = eig(diag(b,1)+diag(b,-1));
[x,ix] = sort(diag(L));     wg = 2*V(1,ix)'.^2;     x = (x+1)/2;    wg = wg/2;
N = size(dom,1);    X = zeros(1,0);     w = 1;
for d = 1:N
    L = dom(d,2)-dom(d,1);  p = dom(d,1)+L*x;   q = L*wg;   m = size(X,1);
    X = [repmat(X,n,1), repelem(p,m,1)];    w = repmat(w,n,1).*repelem(q,m,1);
end
end
