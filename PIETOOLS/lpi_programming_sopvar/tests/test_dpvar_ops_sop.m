function test_dpvar_ops_sop()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_DPVAR_OPS_SOP checks Tier 1b/1c of the container parity map: the
% constructors eye_copvar_sop, zeros_copvar_sop, mat2copvar_sop, and a dpvar
% (scalar, matrix, or constant) acting as an operator on copvar, cdopvar,
% sopvar and sdopvar objects through mtimes, plus, minus, horzcat and
% vertcat. Numeric operands keep their old behavior (checked).
%
% Every result is checked SEMANTICALLY, never against another class routine:
%   1. its decision variables are replaced by random values, A + B'*d per
%      gamma cell (written here, 'fix_op'), and the dpvar by the same values
%      ('dp_val');
%   2. the fixed operator is applied to random polynomial test functions at
%      random points by Gauss-Legendre quadrature of its kernel definition
%      (Sec. 4-5 of the sopvar document; 'act'), and so are the operands;
%   3. the result must equal what the legacy semantics say it is: the
%      scalar or matrix times the operands' outputs, the operands' outputs
%      plus scalar*x, a constant block's action by its definition
%      ('const_act': multiplier / constant extension / integral).
% In 1-D each form is also compared with the LEGACY opvar/dopvar result on
% the same data, both evaluated at the same decision values and compared
% through opvar2copvar with the container eq (which works across block
% partitions). Spatial dimension 1, 2 and 3; verify and the canonical
% multiplier form on every result; errors for bilinear forms, shape
% mismatches and the refused cases.
%
% Uses: rand_copvar, rand_cdopvar, opvar2copvar (sopvar/Testfolder),
% pi_dpvar_eval (claude_tests), rand_opvar, is_canonical_multiplier.
%
% Initial coding MMP, 09/29/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(20260929,'twister');
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
npass = 0;  nfail = 0;
TOL = 1e-9;

    function ck(c,msg)
        if c,   npass = npass+1;
        else,   nfail = nfail+1;    fprintf('  FAIL: %s\n',msg);
        end
    end
    function ckerr(f,id,msg)
        try
            f();
            ck(false,[msg ': no error']);
        catch ME
            ck(strcmp(ME.identifier,id),sprintf('%s: got %s (%s)',msg,ME.identifier,ME.message));
        end
    end
    function ckval(Y,Yx,msg)
        e = 0;  sc = 1;
        for ii = 1:numel(Y)
            e = max(e,max(abs(Y{ii}(:)-Yx{ii}(:)),[],'all'));
            sc = max(sc,max(abs(Yx{ii}(:)),[],'all'));
        end
        ck(e/sc<=TOL,sprintf('%s: rel err %.2e',msg,e/sc));
    end
    function ckok(P,msg)
        % 'verify' (containers), and the canonical multiplier form of every
        % block AS RETURNED, A and B for a decision block (checking a fixed
        % copy would be vacuous: the sopvar constructor canonicalizes).
        if isa(P,'copvar') || isa(P,'cdopvar')
            info = verify(P);
            ck(info.true==1,[msg ': verify ' strjoin(info.flags',' | ')]);
            B = P.C;
        else
            B = {P};
        end
        for kk = 1:numel(B)
            b = B{kk};
            if isempty(b),  continue,   end
            tf = is_canonical_multiplier(b.params,b.vars,b.ZL,b.ZR,b.dims);
            ck(tf,sprintf('%s: block %d not canonical',msg,kk));
        end
    end

gam = dpvar('gam');     h = dpvar({'h1','h2'});
NAMES = {'gam','h1','h2'};      VALS = randn(1,3);
gv = VALS(1);
dom1 = [0 1];

% % % ==================================================== the evaluator
% 'act' is checked first against an independent evaluator, the symbolic
% 'apply_sopvar' (DJ, exact polynomial integration), on random blocks of
% every variable role (shared, output-only, input-only), so that a check
% below cannot pass because the evaluator is wrong or returns zeros.
fprintf('--- evaluator against apply_sopvar\n');
for nd = 1:3
    v = arrayfun(@(k) sprintf('s%d',k),1:nd,'UniformOutput',false);
    roles = {{v,v},{{'s1'},v},{v,{'s1'}},{v,{}},{{},v}};
    for r = 1:numel(roles)
        so = roles{r}{1};   si = roles{r}{2};
        P = rand_copvar(struct('out',{{so}},'in',{{si}}),struct('out',2,'in',2),dom1,2,0.8);
        mb = metadata(P);   xf = rand_funs(mb,2);   T = rand_pts(mb);
        Y = act(P,xf,T);
        if isempty(xf(1).vars),     xp = polynomial(xf(1).coef.');
        else,                       xp = polynomial(xf(1).coef,xf(1).exps,xf(1).vars,[2 1]);
        end
        Ps = apply_sopvar(P.C{1,1},xp);
        vo = mb.vars(mb.space_out(1,:));
        if isempty(vo),     Ys = double(Ps).';
        else,               Ys = pi_poly_grid(Ps,vo,T{1});
        end
        e = max(abs(Y{1}(:)-Ys(:)))/max(1,max(abs(Ys(:))));
        ck(e<=TOL && max(abs(Ys(:)))>1e-3,sprintf('nd=%d role %d: act vs apply_sopvar %.2e',nd,r,e));
    end
end

% % % ================================================================ 1c
fprintf('--- constructors\n');
for nd = 1:3
    v = arrayfun(@(k) sprintf('s%d',k),1:nd,'UniformOutput',false);
    dstr = struct('vars',{v},'dom',[0 1; -1 2; 0 0.5]);
    dstr.dom = dstr.dom(1:nd,:);
    % Square: R^1, L2^2[all], L2^1[s1] (a proper subspace for nd >= 2).
    if nd==1,   sp = {{},v};        dm = [1;2];
    else,       sp = {{},v,{'s1'}}; dm = [1;2;1];
    end
    I = eye_copvar_sop(dm,sp,dstr);
    ck(isa(I,'copvar'),sprintf('nd=%d eye class',nd));
    ckok(I,sprintf('nd=%d eye',nd));
    meta = metadata(I);
    xf = rand_funs(meta,2);     T = rand_pts(meta);
    ckval(act(I,xf,T),ident_act(meta,xf,T,1),sprintf('nd=%d eye action',nd));
    ck(nnz(~cellfun(@isempty,I.C))==numel(dm),sprintf('nd=%d eye diagonal only',nd));
    % Rectangular: out R^2, L2^1[all]; in L2^2[all], R^1, L2^1[s1].
    so = {{},v};    si = {v,{},{'s1'}};
    ds = struct('out',[2;1],'in',[2;1;1]);
    ss = struct('out',{so},'in',{si});
    Z = zeros_copvar_sop(ds,ss,dstr);
    ckok(Z,sprintf('nd=%d zeros',nd));
    ck(Z==0,sprintf('nd=%d zeros == 0',nd));
    ck(nnz(~cellfun(@isempty,Z.C))<=numel(so)+numel(si)-1,sprintf('nd=%d zeros fill count',nd));
    mz = metadata(Z);   xz = rand_funs(mz,2);   Tz = rand_pts(mz);
    Yz = act(Z,xz,Tz);
    ck(all(cellfun(@(y) all(y(:)==0),Yz)),sprintf('nd=%d zeros action',nd));
    % Constant matrix, integral blocks included (L2 -> R^2).
    Mn = randn(3,4);
    Mc = mat2copvar_sop(Mn,ds,ss,dstr);
    ck(isa(Mc,'copvar'),sprintf('nd=%d mat class',nd));
    ckok(Mc,sprintf('nd=%d mat',nd));
    ckval(act(Mc,xz,Tz),const_act(Mn,mz,xz,Tz),sprintf('nd=%d mat action',nd));
    % dpvar matrix: affine in gam, h1, h2.
    Dm = rand_dp(3,4,gam,h);
    Md = mat2copvar_sop(Dm,ds,ss,dstr);
    ck(isa(Md,'cdopvar'),sprintf('nd=%d dpvar mat class',nd));
    ck(isempty(setxor(Md.Zd,Dm.dvarname)),sprintf('nd=%d dpvar mat Zd',nd));
    ckok(Md,sprintf('nd=%d dpvar mat',nd));
    ckval(act(fix_op(Md,NAMES,VALS),xz,Tz),const_act(dp_val(Dm,NAMES,VALS),mz,xz,Tz),...
        sprintf('nd=%d dpvar mat action',nd));
    % Refusals.
    ckerr(@() mat2copvar_sop(Mn,ds,ss,dstr,struct('mult_only',true)),...
        'mat2copvar_grid:integral',sprintf('nd=%d mult_only refuses Q1',nd));
    ckerr(@() mat2copvar_sop(randn(3,3),ds,ss,dstr),'mat2copvar_grid:size',...
        sprintf('nd=%d mat size',nd));
    pv = polynomial(v(1));
    ckerr(@() mat2copvar_sop(gam*pv*ones(3,4),ds,ss,dstr),'mat2copvar_grid:spatial',...
        sprintf('nd=%d spatial dpvar',nd));
    ckerr(@() eye_copvar_sop(ds,ss,dstr),'eye_copvar_sop:notSquare',sprintf('nd=%d eye nonsquare',nd));
end
% Against the legacy detour in 1-D: opvar2copvar(mat2opvar(...)).
s1 = polynomial({'s1'});  th1 = polynomial({'s1_dum'});   % pvar would assignin into a static workspace
Mo = randn(3,3);    Mo(1,2:3) = 0;      % mat2opvar refuses a nonzero Q1
Lg = opvar2copvar(mat2opvar(Mo,[1 1;2 2],[s1,th1],[0 1]));
Mc = mat2copvar_sop(Mo,[1;2],{{},{'s1'}},[0 1],struct('mult_only',true));
ck(eq(Lg,Mc,1e-12),'1-D mat2copvar_sop == opvar2copvar(mat2opvar)');
Ic = eye_copvar_sop([1;2],{{},{'s1'}},[0 1]);
ck(eq(opvar2copvar(mat2opvar(eye(3),[1 1;2 2],[s1,th1],[0 1])),Ic,1e-12),...
    '1-D eye_copvar_sop == opvar2copvar(mat2opvar(eye))');

% % % ================================================================ 1b
for nd = 1:3
    fprintf('--- dpvar operators, %d spatial variable(s)\n',nd);
    v = arrayfun(@(k) sprintf('s%d',k),1:nd,'UniformOutput',false);
    sp = struct('out',{{{},v}},'in',{{{},v}});
    dm = struct('out',[1;2],'in',[1;2]);
    X  = rand_copvar(sp,dm,dom1,1,0.7);
    Xd = rand_cdopvar(sp,dm,dom1,1,3,0.7);
    NAMES = [{'gam','h1','h2'}, reshape(Xd.Zd,1,[])];
    VALS  = [VALS(1:3), randn(1,numel(Xd.Zd))];
    meta = metadata(X);     xf = rand_funs(meta,2);     T = rand_pts(meta);
    YX  = act(X,xf,T);
    YXd = act(fix_op(Xd,NAMES,VALS),xf,T);
    Ix  = ident_act(meta,xf,T,1);
    a0 = -0.5;  a1 = 2;     av = a1*gv + a0;
    lin = @(c,Ya,Yb) cellfun(@(p,q) c(1)*p + c(2)*q,Ya,Yb,'UniformOutput',false);
    tag = @(s) sprintf('nd=%d %s',nd,s);

    % Scalar products.
    R = gam*X;      ck(isa(R,'cdopvar'),tag('gam*X class'));   ckok(R,tag('gam*X'));
    ckval(act(fix_op(R,NAMES,VALS),xf,T),lin([gv 0],YX,YX),tag('gam*X'));
    R = X*gam;      ckok(R,tag('X*gam'));
    ckval(act(fix_op(R,NAMES,VALS),xf,T),lin([gv 0],YX,YX),tag('X*gam'));
    R = -gam*X;     ckok(R,tag('-gam*X'));
    ckval(act(fix_op(R,NAMES,VALS),xf,T),lin([-gv 0],YX,YX),tag('-gam*X'));
    R = (a1*gam + a0)*X;    ckok(R,tag('(2gam-0.5)*X'));
    ckval(act(fix_op(R,NAMES,VALS),xf,T),lin([av 0],YX,YX),tag('(2gam-0.5)*X'));
    R = dpvar(2.5)*X;       ck(isa(R,'copvar'),tag('const dpvar*X is copvar'));
    ckval(act(R,xf,T),lin([2.5 0],YX,YX),tag('const dpvar*X'));
    R = dpvar(2.5)*Xd;      ck(isa(R,'cdopvar'),tag('const dpvar*Xd class'));
    ckval(act(fix_op(R,NAMES,VALS),xf,T),lin([2.5 0],YXd,YXd),tag('const dpvar*Xd'));
    ckerr(@() gam*Xd,'cdopvar:decisionTimesDecision',tag('gam*Xd refused'));
    ckerr(@() Xd*gam,'cdopvar:decisionTimesDecision',tag('Xd*gam refused'));

    % Sums: scalar is scalar*I on the square operator.
    R = gam + X;    ck(isa(R,'cdopvar'),tag('gam+X class'));   ckok(R,tag('gam+X'));
    ckval(act(fix_op(R,NAMES,VALS),xf,T),lin([1 gv],YX,Ix),tag('gam+X'));
    R = X - gam;    ckok(R,tag('X-gam'));
    ckval(act(fix_op(R,NAMES,VALS),xf,T),lin([1 -gv],YX,Ix),tag('X-gam'));
    R = gam - X;    ckok(R,tag('gam-X'));
    ckval(act(fix_op(R,NAMES,VALS),xf,T),lin([-1 gv],YX,Ix),tag('gam-X'));
    R = Xd + gam;   ckok(R,tag('Xd+gam'));
    ck(isempty(setxor(R.Zd,[Xd.Zd(:);{'gam'}])),tag('Xd+gam Zd union'));
    ck(isequal(R.Zd(1:numel(Xd.Zd)),Xd.Zd(:)),tag('Xd+gam keeps Xd list first'));
    ckval(act(fix_op(R,NAMES,VALS),xf,T),lin([1 gv],YXd,Ix),tag('Xd+gam'));
    R = gam - Xd;   ckok(R,tag('gam-Xd'));
    ckval(act(fix_op(R,NAMES,VALS),xf,T),lin([-1 gv],YXd,Ix),tag('gam-Xd'));
    R = dpvar(1e-3) + X;    ck(isa(R,'copvar'),tag('dpvar(1e-3)+X is copvar'));
    ckval(act(R,xf,T),lin([1 1e-3],YX,Ix),tag('dpvar(1e-3)+X'));
    ckerr(@() 1e-3 + X,'plus:badInput',tag('numeric summand refused, as before'));
    % Full-size matrix: a multiplier split by the dimensions; the
    % R^1 <- L2 block (row 1, cols 2:3) would be an integral and is zero.
    Dm = rand_dp(3,3,gam,h);    Dm(1,2:3) = 0;
    R = X + Dm;     ckok(R,tag('X+Dm'));
    Yc = const_act(dp_val(Dm,NAMES,VALS),meta,xf,T);
    ckval(act(fix_op(R,NAMES,VALS),xf,T),lin([1 1],YX,Yc),tag('X+Dm'));
    R = Dm - X;     ckok(R,tag('Dm-X'));
    ckval(act(fix_op(R,NAMES,VALS),xf,T),lin([-1 1],YX,Yc),tag('Dm-X'));
    Dq = rand_dp(3,3,gam,h);    Dq(1,2) = gam;
    ckerr(@() X + Dq,'mat2copvar_grid:integral',tag('X+Dm with an integral block refused'));
    ckerr(@() X + rand_dp(2,2,gam,h),'plus:dimMismatch',tag('X+Dm size mismatch'));
    Xr = rand_copvar(struct('out',{{{},v}},'in',{{v}}),struct('out',[1;2],'in',2),dom1,1,0.7);
    ckerr(@() gam + Xr,'plus:dimMismatch',tag('gam+nonsquare refused'));

    % Matrix products: D*X needs one output space, X*D one input space.
    X1 = rand_copvar(struct('out',{{v}},'in',{{{},v}}),struct('out',2,'in',[1;1]),dom1,1,0.7);
    m1 = metadata(X1);  x1 = rand_funs(m1,2);   T1 = rand_pts(m1);  Y1 = act(X1,x1,T1);
    Dm = rand_dp(3,2,gam,h);    Dv = dp_val(Dm,NAMES,VALS);
    R = Dm*X1;      ck(isa(R,'cdopvar'),tag('Dm*X1 class'));  ckok(R,tag('Dm*X1'));
    T1o = T1;       % the output space is X1's, with 3 components
    ckval(act(fix_op(R,NAMES,VALS),x1,T1o),{Y1{1}*Dv.'},tag('Dm*X1'));
    Mn = randn(3,2);
    R = dpvar(Mn)*X1;   ck(isa(R,'copvar'),tag('dpvar(M)*X1 is copvar'));
    ckval(act(R,x1,T1),{Y1{1}*Mn.'},tag('dpvar(M)*X1'));
    ckerr(@() Mn*X1,'mtimes:nonScalarNumeric',tag('numeric matrix factor refused, as before'));
    X2 = rand_copvar(struct('out',{{{},v}},'in',{{v}}),struct('out',[1;1],'in',2),dom1,1,0.7);
    m2 = metadata(X2);  T2 = rand_pts(m2);
    D2 = rand_dp(2,3,gam,h);    D2v = dp_val(D2,NAMES,VALS);
    x3 = rand_funs(struct('vars',{m2.vars},'space_in',m2.space_in,'dim_in',3),2);
    x3D = x3;       x3D(1).mix = D2v;       % the test function D2v*x3
    R = X2*D2;      ckok(R,tag('X2*D2'));
    ckval(act(fix_op(R,NAMES,VALS),x3,T2),act(X2,x3D,T2),tag('X2*D2'));
    ckerr(@() Dm*X,'mtimes:ambiguousMatrix',tag('Dm*X with two row spaces'));
    ckerr(@() rand_dp(3,3,gam,h)*X1,'mtimes:dimMismatch',tag('Dm*X1 inner dims'));
    Xd1 = rand_cdopvar(struct('out',{{v}},'in',{{{},v}}),struct('out',2,'in',[1;1]),dom1,1,2,0.7);
    ckerr(@() Dm*Xd1,'cdopvar:decisionTimesDecision',tag('Dm*Xd1 refused'));

    % Concatenation, the DEMO5 shape [-gam, D', B'; D, -gam, C; B, C', A].
    sR = {{}};  sL = {v};
    Dp = rand_copvar(struct('out',{sR},'in',{sR}),struct('out',1,'in',1),dom1,1,1);
    Bp = rand_copvar(struct('out',{sL},'in',{sR}),struct('out',2,'in',1),dom1,1,1);
    Cp = rand_copvar(struct('out',{sR},'in',{sL}),struct('out',1,'in',2),dom1,1,1);
    Ap = rand_copvar(struct('out',{sL},'in',{sL}),struct('out',2,'in',2),dom1,1,0.7);
    Dpt = Dp';  Bpt = Bp';  Cpt = Cp';
    Q = [-gam, Dpt, Bpt; Dp, -gam, Cp; Bp, Cpt, Ap];
    ck(isa(Q,'cdopvar') && isequal(size(Q),[3 3]),tag('DEMO5 shape class and grid'));
    ckok(Q,tag('DEMO5 shape'));
    O = {-gam, Dpt, Bpt; Dp, -gam, Cp; Bp, Cpt, Ap};
    mq = metadata(Q);   xq = rand_funs(mq,2);   Tq = rand_pts(mq);
    ckval(act(fix_op(Q,NAMES,VALS),xq,Tq),grid_act(O,xq,Tq,NAMES,VALS),tag('DEMO5 shape'));
    % A decision block beside a dpvar entry: Zd merge.
    Pd = rand_cdopvar(struct('out',{sL},'in',{sL}),struct('out',2,'in',2),dom1,1,3,0.7);
    NAMES = [NAMES, reshape(Pd.Zd,1,[])];   VALS = [VALS, randn(1,numel(Pd.Zd))];
    Q2 = [-gam, Bpt; Bp, Pd];
    ckok(Q2,tag('[-gam,B''; B,Pd]'));
    ck(isempty(setxor(Q2.Zd,[{'gam'};Pd.Zd(:)])),tag('[-gam,B''; B,Pd] Zd'));
    mq2 = metadata(Q2);     xq2 = rand_funs(mq2,2);     Tq2 = rand_pts(mq2);
    ckval(act(fix_op(Q2,NAMES,VALS),xq2,Tq2),grid_act({-gam,Bpt;Bp,Pd},xq2,Tq2,NAMES,VALS),...
        tag('[-gam,B''; B,Pd]'));
    % A matrix entry into L2 (legacy Q2): the constant function.
    D21 = rand_dp(2,1,gam,h);
    R = [D21, Ap];      ckok(R,tag('[D21, A]'));
    mr = metadata(R);   xr = rand_funs(mr,2);   Tr = rand_pts(mr);
    ckval(act(fix_op(R,NAMES,VALS),xr,Tr),grid_act({D21,Ap},xr,Tr,NAMES,VALS),tag('[D21, A]'));
    % A matrix entry out of L2 (legacy Q1): the integral.
    D12 = rand_dp(1,2,gam,h);
    R = [D12; Ap];      ckok(R,tag('[D12; A]'));
    mr = metadata(R);   xr = rand_funs(mr,2);   Tr = rand_pts(mr);
    ckval(act(fix_op(R,NAMES,VALS),xr,Tr),grid_act({D12;Ap},xr,Tr,NAMES,VALS),tag('[D12; A]'));
    ckerr(@() [gam, Ap],'horzcat:dimMismatch',tag('[gam, A] rows'));
    ckerr(@() [gam; Ap],'vertcat:dimMismatch',tag('[gam; A] columns'));
    R = [Ap, dpvar(zeros(2,1))];   ck(isa(R,'copvar'),tag('[A, dpvar(zeros)] is copvar'));
    ckok(R,tag('[A, dpvar(zeros)]'));
    ckerr(@() [Ap, zeros(2,1)],'copvar:horzcatBadOperand',tag('numeric entry refused, as before'));

    % Block level: sopvar / sdopvar operands.
    Xb = X.C{2,2};      Xdb = Xd.C{2,2};
    mb = struct('vars',{v},'dom',repmat(dom1,nd,1),'space_out',true(1,nd),...
        'space_in',true(1,nd),'dim_out',2,'dim_in',2);
    xb = rand_funs(mb,2);   Tb = rand_pts(mb);  Yb = act(Xb,xb,Tb);  Ib = ident_act(mb,xb,Tb,1);
    R = gam*Xb;     ck(isa(R,'sdopvar'),tag('gam*Xb is sdopvar'));
    ckval(act(fix_op(R,NAMES,VALS),xb,Tb),lin([gv 0],Yb,Yb),tag('gam*Xb'));
    R = Xb*gam;     ckval(act(fix_op(R,NAMES,VALS),xb,Tb),lin([gv 0],Yb,Yb),tag('Xb*gam'));
    R = gam + Xb;   ck(isa(R,'sdopvar'),tag('gam+Xb is sdopvar'));
    ckval(act(fix_op(R,NAMES,VALS),xb,Tb),lin([1 gv],Yb,Ib),tag('gam+Xb'));
    R = Xb - gam;   ckval(act(fix_op(R,NAMES,VALS),xb,Tb),lin([1 -gv],Yb,Ib),tag('Xb-gam'));
    R = dpvar(0.25) + Xb;   ck(isa(R,'sopvar'),tag('dpvar(0.25)+Xb is sopvar'));
    ckval(act(R,xb,Tb),lin([1 0.25],Yb,Ib),tag('dpvar(0.25)+Xb'));
    Db = rand_dp(3,2,gam,h);    Dbv = dp_val(Db,NAMES,VALS);
    R = Db*Xb;      ck(isa(R,'sdopvar'),tag('Db*Xb is sdopvar'));
    ckval(act(fix_op(R,NAMES,VALS),xb,Tb),{Yb{1}*Dbv.'},tag('Db*Xb'));
    ckerr(@() gam*Xdb,'cdopvar:decisionTimesDecision',tag('gam*Xdb refused'));
    D21b = rand_dp(2,1,gam,h);
    R = [D21b, Xb];     ck(isa(R,'cdopvar') && isequal(size(R),[1 2]),tag('[D21, Xb] is a 1x2 container'));
    ckok(R,tag('[D21, Xb]'));
    mr = metadata(R);   xr = rand_funs(mr,2);   Tr = rand_pts(mr);
    ckval(act(fix_op(R,NAMES,VALS),xr,Tr),grid_act({D21b,Xb},xr,Tr,NAMES,VALS),tag('[D21, Xb]'));
    ckerr(@() [gam, Xb],'horzcat:dimMismatch',tag('[gam, Xb] rows'));
    ckerr(@() [gam; Xdb],'vertcat:dimMismatch',tag('[gam; Xdb] columns'));
end

% % % ======================================================= legacy, 1-D
fprintf('--- against legacy opvar/dopvar (1-D)\n');
s = polynomial({'s'});  th = polynomial({'s_dum'});    % opvar2sopvar needs the '_dum' dummy
I0 = [0 1];
NAMES = {'gam','h1','h2'};  VALS = VALS(1:3);
Xo  = rand_opvar([1 1;2 2],2,s,th,I0);
cmp = @(Lg,Rc,msg) ck(eq(opvar2copvar(fix_dopvar(Lg,NAMES,VALS)),fix_op(Rc,NAMES,VALS),1e-10),msg);
Xc  = opvar2copvar(Xo);
cmp(gam*Xo,gam*Xc,'legacy gam*X');
cmp(Xo*gam,Xc*gam,'legacy X*gam');
cmp(-gam*Xo,-gam*Xc,'legacy -gam*X');
cmp((2*gam-0.5)*Xo,(2*gam-0.5)*Xc,'legacy (2gam-0.5)*X');
cmp(gam+Xo,gam+Xc,'legacy gam+X');
cmp(Xo-gam,Xc-gam,'legacy X-gam');
cmp(gam-Xo,gam-Xc,'legacy gam-X');
Dm = rand_dp(3,3,gam,h);    Dm(1,2:3) = 0;
cmp(Dm+Xo,Dm+Xc,'legacy Dm+X');
cmp(Xo-Dm,Xc-Dm,'legacy X-Dm');
Xout = rand_opvar([0 1;2 2],2,s,th,I0);     % one output type: L2^2
D32 = rand_dp(3,2,gam,h);
cmp(D32*Xout,D32*opvar2copvar(Xout),'legacy D*X (output L2)');
Xin = rand_opvar([1 0;2 2],2,s,th,I0);      % one input type: L2^2
D23 = rand_dp(2,3,gam,h);
cmp(Xin*D23,opvar2copvar(Xin)*D23,'legacy X*D (input L2)');
Aop = rand_opvar([0 0;2 2],2,s,th,I0);
D21 = rand_dp(2,1,gam,h);   D12 = rand_dp(1,2,gam,h);
cmp([D21, Aop],[D21, opvar2copvar(Aop)],'legacy [D, A] (Q2)');
cmp([D12; Aop],[D12; opvar2copvar(Aop)],'legacy [D; A] (Q1)');
Do = rand_opvar([1 1;0 0],2,s,th,I0);   Bo = rand_opvar([0 1;2 0],2,s,th,I0);
Co = rand_opvar([1 0;0 2],2,s,th,I0);
Lq = [-gam, -Do', -Bo'; -Do, -gam, Co; -Bo, Co', Aop];
Dc = opvar2copvar(Do);  Bc = opvar2copvar(Bo);  Cc = opvar2copvar(Co);  Ac = opvar2copvar(Aop);
Rq = [-gam, -Dc', -Bc'; -Dc, -gam, Cc; -Bc, Cc', Ac];
cmp(Lq,Rq,'legacy DEMO5 shape');
DEMO6 = [-gam*eye(2), randn(2,1); randn(1,2), -gam];
cmp(mat2opvar(DEMO6,[3 3;0 0],[s,th],I0),mat2copvar_sop(DEMO6,3,{{}},[]),...
    'legacy mat2opvar(dpvar) vs mat2copvar_sop, -gam*eye(2) block matrix');

fprintf('test_dpvar_ops_sop: %d passed, %d failed\n',npass,nfail);
if nfail>0
    error('test_dpvar_ops_sop:failed','%d checks failed.',nfail);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decision values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = rand_dp(m,n,gam,h)
% A random m x n dpvar affine in gam, h1, h2, every entry nonzero.
D = randn(m,n) + gam*randn(m,n) + h(1)*randn(m,n) + h(2)*randn(m,n);
end


function V = dp_val(D,names,vals)
% Value of a spatially constant dpvar at the decision values: rows of D.C
% are (matrix row outer, decision variable inner, constant first).
assert(size(D.degmat,1)==1 && ~any(D.degmat(:)),'dp_val: D must be constant in space');
[tf,loc] = ismember(D.dvarname,names);
assert(all(tf),'dp_val: unknown decision variable');
dv = [1; reshape(vals(loc),[],1)];
m = size(D,1);
V = full(kron(speye(m),dv.')*D.C);
end


function P = fix_op(P,names,vals)
% Substitute decision values: C = unvec(A + B'*d) per gamma cell, for every
% sdopvar block, by the kernel definition in sdopvar.m. Fixed objects pass.
if isa(P,'sdopvar')
    [tf,loc] = ismember(P.Zd,names);
    assert(all(tf),'fix_op: a decision variable has no value');
    d = reshape(vals(loc),[],1);
    NL = prod([cellfun(@numel,P.ZL),1]);    NR = prod([cellfun(@numel,P.ZR),1]);
    prm = cell(size(P.params.A));
    for k = 1:numel(prm)
        a = P.params.A{k};  b = P.params.B{k};
        if isempty(a),  prm{k} = sparse(P.dims(1)*NL,P.dims(2)*NR);    continue,   end
        prm{k} = reshape(a + b.'*d,P.dims(1)*NL,P.dims(2)*NR);
    end
    P = sopvar(prm,P.vars,P.ZL,P.ZR,P.dom,P.dims);
elseif isa(P,'cdopvar')
    C = P.C;
    for k = 1:numel(C)
        if ~isempty(C{k}),  C{k} = fix_op(C{k},names,vals);  end
    end
    m = metadata(P);
    P = copvar(C,rmfield(m,'Zd'));
end
end


function Pf = fix_dopvar(P,names,vals)
% A legacy dopvar at the decision values, field by field (pi_dpvar_eval).
% Either class: an 'opvar' can hold dpvar fields ([D, A] with a dpvar D).
Pf = opvar();
Pf.I = P.I;     Pf.var1 = P.var1;   Pf.var2 = P.var2;
f = {'P','Q1','Q2'};
for k = 1:3,    Pf.(f{k}) = ev_dp(P.(f{k}),names,vals);       end
f = {'R0','R1','R2'};
for k = 1:3,    Pf.R.(f{k}) = ev_dp(P.R.(f{k}),names,vals);   end
Pf.dim = Pf.dim;
end


function y = ev_dp(x,names,vals)
if isa(x,'dpvar'),  y = pi_dpvar_eval(x,names,vals);
else,               y = x;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test functions, points and the quadrature evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xf = rand_funs(meta,deg)
% One random polynomial test function per input space: degree <= deg in
% each of the space's variables, dim_in(j) components.
N = numel(meta.dim_in);
xf = struct('vars',cell(1,N),'exps',[],'coef',[],'mix',[]);
for j = 1:N
    vj = meta.vars(meta.space_in(j,:));
    g = cell(1,numel(vj));  [g{:}] = ndgrid(0:deg);
    if isempty(vj),     E = zeros(1,0);
    else,               E = reshape(cat(numel(vj)+1,g{:}),[],numel(vj));
    end
    xf(j).vars = vj;    xf(j).exps = E;
    xf(j).coef = randn(size(E,1),meta.dim_in(j));
end
end


function F = fx(f,S)
% Values of test function f at the points S (columns: f.vars), G x p; with
% f.mix set, the values of f.mix * f.
G = size(S,1);  E = f.exps;
Mn = ones(G,size(E,1));
for e = 1:size(E,1)
    for k = 1:size(E,2),    Mn(:,e) = Mn(:,e).*S(:,k).^E(e,k);    end
end
F = Mn*f.coef;
if ~isempty(f.mix),     F = F*f.mix.';  end
end


function T = rand_pts(meta)
% Five random output points per output space, in its domain; a single
% (empty) point for R^q.
M = numel(meta.dim_out);    T = cell(1,M);
for i = 1:M
    k = find(meta.space_out(i,:));
    if isempty(k),  T{i} = zeros(1,0);  continue,   end
    d = meta.dom(k,:);
    T{i} = d(:,1).' + rand(5,numel(k)).*(d(:,2)-d(:,1)).';
end
end


function Y = act(P,xf,T)
% Output of the FIXED operator P (copvar or sopvar) on the test functions
% xf at the points T, by quadrature of the kernel definition.
if isa(P,'sopvar'),     Y = {act_block(P,xf(1),T{1},sort(P.vars.out))};     return,     end
M = size(P.C,1);    N = size(P.C,2);    Y = cell(1,M);
for i = 1:M
    to = P.vars(P.space_out(i,:));
    Y{i} = zeros(size(T{i},1),P.dim_out(i));
    for j = 1:N
        if isempty(P.C{i,j}),   continue,   end
        Y{i} = Y{i} + act_block(P.C{i,j},xf(j),T{i},to);
    end
end
end


function y = act_block(b,f,Tp,Tvars)
% (b*f)(t) = sum_gamma int K_gamma(t,s') f(s') ds', K_gamma = (I_m (x)
% ZL(t))' C_gamma (I_n (x) ZR(s')); in a shared direction k gamma_k = 1 is
% s'_k = t_k, 2 is int_a^t, 3 is int_t^b; input-only directions int_a^b.
vout = b.vars.out;  vin = b.vars.in;
m = b.dims(1);      n = b.dims(2);
n3 = numel(intersect(vin,vout));
[~,io] = ismember(vout,Tvars);      tp = Tp(:,io);
[~,i3] = ismember(vin(1:n3),vout);  % position of each shared variable in vout
[~,jf] = ismember(f.vars,vin);      % columns of f in vin order
NL = prod([cellfun(@numel,b.ZL),1]);    NR = prod([cellfun(@numel,b.ZR),1]);
[xg,wg] = gauss(4);
G = size(tp,1);     y = zeros(G,m);
for g = 1:G
    t = tp(g,:);
    zl = 1;     for k = 1:numel(vout),  zl = kron(zl,t(k).^b.ZL{k}(:));     end
    for c = 1:numel(b.params)
        C = full(b.params{c});
        if isempty(C) || ~any(C(:)),    continue,   end
        gm = ones(1,n3);
        if n3>0,    s = cell(1,n3);     [s{:}] = ind2sub([3*ones(1,n3),1],c);   gm = [s{:}];    end
        % nodes and weights per input variable
        nodes = cell(1,numel(vin));     wts = cell(1,numel(vin));
        for k = 1:numel(vin)
            a = b.dom.in(k,1);  bb = b.dom.in(k,2);
            if k<=n3
                tk = t(i3(k));
                switch gm(k)
                    case 1, nodes{k} = tk;  wts{k} = 1;             continue
                    case 2, bb = tk;
                    case 3, a = tk;
                end
            end
            nodes{k} = (a+bb)/2 + (bb-a)/2*xg;  wts{k} = (bb-a)/2*wg;
        end
        [S,W] = tensor_nodes(nodes,wts);
        zr = ones(size(S,1),1);
        for k = 1:numel(vin),   zr = rowkron(zr,S(:,k).^(b.ZR{k}(:).'));   end
        X = fx(f,S(:,jf));
        vv = (W.'*rowkron(X,zr)).';         % int kron(x, ZR) ds'
        u = C*vv;
        y(g,:) = y(g,:) + (reshape(u,NL,m).'*zl).';
    end
end
end


function Y = ident_act(meta,xf,T,c)
% c*x_i at the points of output space i (square: space i in = out).
M = numel(meta.dim_out);    Y = cell(1,M);
for i = 1:M,    Y{i} = c*fx(xf(i),T{i});    end
end


function Y = const_act(Mv,meta,xf,T)
% Action of the constant matrix Mv between the spaces of meta, by the
% definition: output-only variables constant, shared ones pointwise,
% input-only ones integrated over their domain.
M = numel(meta.dim_out);    N = numel(meta.dim_in);
ro = cumsum([0;meta.dim_out(:)]);   co = cumsum([0;meta.dim_in(:)]);
[xg,wg] = gauss(4);
Y = cell(1,M);
for i = 1:M
    vo = meta.vars(meta.space_out(i,:));
    Y{i} = zeros(size(T{i},1),meta.dim_out(i));
    for j = 1:N
        Mij = Mv(ro(i)+1:ro(i+1),co(j)+1:co(j+1));
        vi = xf(j).vars;
        for g = 1:size(T{i},1)
            nodes = cell(1,numel(vi));  wts = nodes;
            for k = 1:numel(vi)
                p = find(strcmp(vo,vi{k}));
                if ~isempty(p),     nodes{k} = T{i}(g,p);   wts{k} = 1;
                else
                    d = meta.dom(strcmp(meta.vars,vi{k}),:);
                    nodes{k} = (d(1)+d(2))/2 + (d(2)-d(1))/2*xg;  wts{k} = (d(2)-d(1))/2*wg;
                end
            end
            [S,W] = tensor_nodes(nodes,wts);
            Y{i}(g,:) = Y{i}(g,:) + (Mij*(fx(xf(j),S).'*W)).';
        end
    end
end
end


function Y = grid_act(O,xf,T,names,vals)
% Row i of a grid of operands: sum_j O{i,j} applied to x_j, a dpvar entry
% by its value as a constant block between the row's and column's spaces.
M = size(O,1);  N = size(O,2);
% Row and column spaces from the operator entries.
Y = cell(1,M);
for i = 1:M
    for j = 1:N
        o = O{i,j};
        if isa(o,'dpvar')
            % spaces of row i / column j, read off an operator in them
            vo = row_vars(O,i);     mo = size(o,1);
            mm = struct('vars',{union(vo,xf(j).vars)},'dom',[],'space_out',[],'space_in',[],...
                'dim_out',mo,'dim_in',size(o,2));
            mm.vars = reshape(mm.vars,1,[]);
            mm.dom = repmat([0 1],numel(mm.vars),1);
            mm.space_out = ismember(mm.vars,vo);   mm.space_in = ismember(mm.vars,xf(j).vars);
            yij = const_act(dp_val(o,names,vals),mm,xf(j),{T{i}});
            yij = yij{1};
        else
            of = fix_op(o,names,vals);
            yij = act(of,xf(j),T(i));
            yij = yij{1};
        end
        if isempty(Y{i}),   Y{i} = yij;     else,   Y{i} = Y{i} + yij;  end
    end
end
end


function vo = row_vars(O,i)
% Output variables of row i of the operand grid, from an operator entry.
for j = 1:size(O,2)
    o = O{i,j};
    if isa(o,'copvar') || isa(o,'cdopvar'),     vo = o.vars(o.space_out(1,:));  return
    elseif isa(o,'sopvar') || isa(o,'sdopvar'), vo = sort(o.vars.out);          return
    end
end
vo = cell(1,0);
end


function [S,W] = tensor_nodes(nodes,wts)
% Tensor product of per-variable nodes and weights.
% repelem(S,nk,1) repeats each earlier row nk times and repmat tiles the new
% nodes, so row (r-1)*nk + q has weight W(r)*w(q): kron(W,w).
S = zeros(1,0);     W = 1;
for k = 1:numel(nodes)
    nk = numel(nodes{k});
    S = [repelem(S,nk,1), repmat(nodes{k}(:),size(S,1),1)];
    W = kron(W(:),wts{k}(:));
end
end


function Z = rowkron(A,B)
% Row-wise Kronecker product: Z(q,:) = kron(A(q,:),B(q,:)). The product
% permute(A,[1 3 2]).*B is G x nb x na with (q,j,i) = A(q,i)*B(q,j), so the
% column-major reshape gives column j + (i-1)*nb, kron's order.
Z = reshape(permute(A,[1 3 2]).*B,size(A,1),[]);
end


function [x,w] = gauss(n)
% Gauss-Legendre nodes and weights on [-1,1] (Golub-Welsch).
k = 1:n-1;  beta = k./sqrt(4*k.^2-1);
[V,D] = eig(diag(beta,1)+diag(beta,-1));
[x,p] = sort(diag(D));  w = 2*V(1,p).'.^2;
end
