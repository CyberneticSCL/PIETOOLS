function T = heatNd_poincare_validate(stage,vopts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = HEATND_POINCARE_VALIDATE(STAGE,VOPTS) validation runs of the Poincare
% diagnostic (HEATND_POINCARE) against the exact constants, the stock
% reference numbers below, and - in 3-D - the predictions of the tensor
% bounds from the 1-D references at the same per-direction basis.
%
% STAGE
%   'I0'   HEATND_POINCARE_OPS I0 for N = 1, 2, 3 (default and alternative
%          BC sets / boxes); no SDP; seconds;
%   '1d'   N = 1, DD and DN: L2 {none, product, faces} d = 1..3 and energy
%          {product, faces} d = 1..2; per case a bisection (HEATND_BISECT
%          logic, maxsolve 40, rtol 1e-5, retry 'loose'), the objective mode
%          (+ re-certification at (1-1e-3) lam-hat), and a soundness probe at
%          1.001 x exact;
%   '2d'   N = 2 (DD x DN), d = 1: L2 {none, product, faces} and energy faces,
%          each full / dir 1 / dir 2, and energy product full; as '1d';
%   '2d2'  N = 2, d = 2, faces: L2 full / dir 1 / dir 2: objective,
%          re-certification and soundness only (a bisection is 4-5 min per
%          case at stock's 15-41 s per solve);
%   '3d0'  N = 3: I0 and the smoke test (L2 full, d = 1, no generator, lam
%          = 0, all rows, MOSEK defaults: the reference's measured 88.6 s);
%   '3d'   N = 3, d = 1 faces, TWO-SIDED tests of the predictions, no search:
%          VOPTS.tests (cellstr of 'L2dir<i>', 'Edir<i>', 'L2full') at
%          VOPTS.frac x prediction (default [0.99 1.01]), one solve per
%          lam, all rows (HEATND_LINDEP's QR of this At exceeded 70 GB,
%          MEASURED), tight tolerances, MOSEK time capped at VOPTS.tmax3 s
%          (a capped solve is 0, never coerced). Predictions: per direction
%          the 1-D objective-mode lam-hat at the same basis and generators for
%          that direction's BC (tensor bounds B1 = B2: equality, reference
%          study). PER-DIRECTION TESTS USE THE FACES OF THE TESTED DIRECTION
%          ONLY (VOPTS.faces3 'dir'): the B2 certificate (1-D certificate (x)
%          T_ihat^2) needs only g_i, and B1 holds for any generator set, so the
%          prediction is the same as with all 2N faces, while the all-faces
%          3-D solves are priced over 20 min (MEASURED capped solves: L2 m
%          28672, 7 x 512, ~110 s/iteration; energy m 31488, 7 x 640, ~140
%          s/iteration, 19.7 GB). A certified x is saved (x_<test>_<lam>.mat);
%          'L2full' is priced the same and is not solved by default;
%   '3dsum' N = 3: the full L2 LOWER bound B2 as a certificate, no solve:
%          the saved certified per-direction points x_i at lam_i (faces of
%          direction i only, a subcone of the all-faces cone) are summed into
%          the all-faces L2 full SDP - G0 = sum_i G0_i, face block (i,e) from
%          x_i - since sum_i Q_i(lam_i) = Q_full(sum_i lam_i); the result is
%          checked by the HEATND_SOLVE gate (rel_b <= 1e-6 on all rows of the
%          full SDP, every block PSD). VOPTS.sumfiles: the three x files.
%
% VOPTS (defaults): dir (fullfile(tempdir,'heatNd','poincare') - NOT
%   durable), tests ({}), frac ([0.99 1.01]), tight3 (true), tmax3 (1080 s
%   MOSEK time cap), faces3 ('dir' | 'all'), sumfiles ({}), maxsolve (40),
%   rtol (1e-5).
%
% Stock reference (scout ref study, stock PIETOOLS on this tree, MOSEK 11,
% bisection to 1e-4 in lam/exact; 'T_d' = tensor degree d):
%   1-D L2 none: DD T1 1.0e-4, T2 0.091, T3 1.051; DN T1 1.8e-4, T2 0.109,
%       T3 1.122. L2 product: DD 0.99997 / 9.86899 / 9.86960; DN 0.99997 /
%       2.46740 / 2.46740. Energy product: DD 9.70065 / 9.86960; DN 2.46709
%       / 2.46740.
%   2-D (opvar2d) T1: L2 none full 0.00165, dir1 1.0e-4, dir2 3.3e-4;
%       product 0.00241 / 0.00132 / 0.00048; faces [3 4 5 6] 1.99994 /
%       0.99997 / 0.99997; energy faces 12.03906 / 9.71626 / 2.46725; energy
%       product full 2.96% of lambda1. T2 faces: L2 12.33508 / 9.86930 /
%       2.46732.
%   1-D faces (not run in stock): predicted by the stock 2-D per-direction
%       faces values through B1 = B2.
%
% OUTPUT T: struct array, one row per case: N, bc, form, sel, d, gen, mu
% (exact), ref, refsrc, m, Ks, nnz, lo, hi (certified bracket), stop,
% nsolve, lamhat, obj_st, lam_recert, recert_st, sound_st, tests ([lam st
% rel_b psd_relmin t_mosek iter mosek_threads] per solve), t (s), peak (MB);
% saved to VOPTS.dir.
%
% Initial coding MMP, 09/27/2026
% MMP, 09/27/2026 (final reviews): the '3d' log prints cert_viol, the
%   quantity HEATND_SOLVE now gates -1 on (it printed cert_rel).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2,    vopts = struct();   end
df = struct('dir',fullfile(tempdir,'heatNd','poincare'),'tests',{{}},'frac',[0.99 1.01], ...
            'tight3',true,'tmax3',1080,'faces3','dir','sumfiles',{{}},'maxsolve',40,'rtol',1e-5);
fn = fieldnames(df);
for i = 1:numel(fn),    if ~isfield(vopts,fn{i}),   vopts.(fn{i}) = df.(fn{i});     end,    end
if ~exist(vopts.dir,'dir'),     mkdir(vopts.dir);   end
info = heatNd_path();
fprintf('heatNd_poincare_validate(%s): maxNumCompThreads = %d (MOSEK uses its own count, per row), %s\n', ...
        stage,info.threads,char(datetime('now')));
T = struct('N',{},'bc',{},'form',{},'sel',{},'d',{},'gen',{},'mu',{},'ref',{},'refsrc',{},'m',{},'Ks',{}, ...
           'nnz',{},'lo',{},'hi',{},'stop',{},'nsolve',{},'lamhat',{},'obj_st',{},'lam_recert',{}, ...
           'recert_st',{},'sound_st',{},'tests',{},'t',{},'peak',{});
ref = stockref();
switch stage
    case 'I0'
        cs = { {1,{'DD'},[0 1]}, {1,{'DN'},[0 1]}, {1,{'ND'},[-1 2]}, {2,{'DD','DN'},[0 1;0 1]}, ...
               {2,{'ND','DD'},[0 2;-1 1]}, {3,{'DD','DN','DN'},[0 1;0 1;0 1]}, {3,{'DN','ND','DD'},[0 1;1 2;0 3]} };
        for c = 1:numel(cs)
            [N,bc,dom] = cs{c}{:};
            t0 = tic;   O = heatNd_poincare_ops(heatNd_pie(N,0,bc,dom),true);
            fprintf('\nI0 N=%d %s %s: %s (%d checks, %.1f s)\n',N,strjoin(bc,'x'),mat2str(dom), ...
                    tern(O.I0ok,'PASS','FAIL'),numel(O.I0),toc(t0));
            for k = 1:numel(O.I0),  fprintf('  %-50s %9.2e  (tol %.0e)  %d\n',O.I0(k).check,O.I0(k).err,O.I0(k).tol,O.I0(k).ok);  end
        end
        return
    case {'1d','2d','2d2'}
        switch stage
            case '1d'
                cs = {};
                for b = {'DD','DN'}
                    for d = 1:3,  cs(end+1,:) = {1,b,'L2','full',d,'none'};     end %#ok<AGROW>
                    for d = 1:3,  cs(end+1,:) = {1,b,'L2','full',d,'product'};  end %#ok<AGROW>
                    for d = 1:3,  cs(end+1,:) = {1,b,'L2','full',d,'faces'};    end %#ok<AGROW>
                    for d = 1:2,  cs(end+1,:) = {1,b,'energy','full',d,'product'};  end %#ok<AGROW>
                    for d = 1:2,  cs(end+1,:) = {1,b,'energy','full',d,'faces'};    end %#ok<AGROW>
                end
                full_search = true;
            case '2d'
                b = {'DD','DN'};    cs = {};
                for g = {'none','product','faces'}
                    for s = {'full',1,2},   cs(end+1,:) = {2,b,'L2',s{1},1,g{1}};   end %#ok<AGROW>
                end
                for s = {'full',1,2},   cs(end+1,:) = {2,b,'energy',s{1},1,'faces'};  end %#ok<AGROW>
                cs(end+1,:) = {2,b,'energy','full',1,'product'};
                full_search = true;
            case '2d2'
                b = {'DD','DN'};    cs = {};
                for s = {'full',1,2},   cs(end+1,:) = {2,b,'L2',s{1},2,'faces'};  end %#ok<AGROW>
                full_search = false;
        end
        Oc = containers.Map();
        for c = 1:size(cs,1)
            [N,bc,form,sel,d,gen] = cs{c,:};
            if ~iscell(bc),     bc = {bc};  end
            key = strjoin(bc,'x');
            if ~isKey(Oc,key),  Oc(key) = heatNd_poincare_ops(heatNd_pie(N,0,bc));    end
            T(end+1) = run_case(Oc(key),struct('form',form,'sel',sel,'d',d,'gen',gen),ref,vopts,full_search); %#ok<AGROW>
            report(T(end));
        end
        save(fullfile(vopts.dir,sprintf('pval_%s.mat',stage)),'T');
    case '3d0'
        heatNd_poincare_validate('I0',vopts);
        O = heatNd_poincare_ops(heatNd_pie(3,0));
        t0 = tic;
        R = heatNd_poincare(O,struct('form','L2','sel','full','d',1,'gen','none'), ...
                            struct('mode','test','lams',0,'solve',struct('tight',false), ...
                                   'bisect',struct('retry',{{}})));
        T(end+1) = mkrow(O,R,NaN,'smoke: lam = 0 certifies (reference 88.6 s)');
        T(end).tests = trace7(R.test.trace);    T(end).t = toc(t0);
        report(T(end));
        save(fullfile(vopts.dir,sprintf('pval_3d0_%s.mat',char(datetime('now','Format','HHmmss')))),'T');
    case '3d'
        O = heatNd_poincare_ops(heatNd_pie(3,0));
        % 1-D references at the same per-direction basis (seconds).
        P = struct();
        for b = unique(O.bc)
            O1 = heatNd_poincare_ops(heatNd_pie(1,0,b));
            for f = {'L2','energy'}
                Ro = heatNd_poincare(O1,struct('form',f{1},'d',1,'gen','faces'), ...
                                     struct('mode','objective','verbose',false));
                P.([f{1} '_' b{1}]) = Ro.obj.lamhat;
                fprintf('1-D reference %s %s d=1 faces: lam-hat %.7f (st %+d, recert %+d)\n',f{1},b{1}, ...
                        Ro.obj.lamhat,Ro.obj.v.st,Ro.obj.recert.st);
            end
        end
        for k = 1:numel(vopts.tests)
            nm = vopts.tests{k};
            fc = [];                                % all 2N faces
            if strcmp(nm,'L2full')
                sp = struct('form','L2','sel','full');
                pr = sum(cellfun(@(b) P.(['L2_' b]),O.bc));
            else
                tok = regexp(nm,'^(L2|E)dir(\d)$','tokens','once');
                if isempty(tok),    error('heatNd_poincare_validate:test','Unknown test ''%s''.',nm),  end
                i = str2double(tok{2});
                if strcmp(tok{1},'E'),  sp = struct('form','energy','sel',i);  else,  sp = struct('form','L2','sel',i);  end
                pr = P.([tern(strcmp(tok{1},'E'),'energy','L2') '_' O.bc{i}]);
                if strcmp(vopts.faces3,'dir'),  fc = false(O.N,2);  fc(i,:) = true;     end
            end
            sp.d = 1;   sp.gen = 'faces';   sp.faces = fc;
            t0 = tic;   lams = vopts.frac*pr;
            R = heatNd_poincare(O,sp,struct('mode','sdp','lams',lams));
            so = struct('tight',vopts.tight3,'keepx',true,'param',struct('MSK_DPAR_OPTIMIZER_MAX_TIME',vopts.tmax3));
            T(end+1) = mkrow(O,R,pr,sprintf('%s: 1-D-predicted (tensor bounds), faces %s',nm,vopts.faces3));  %#ok<AGROW>
            for j = 1:numel(lams)
                v = heatNd_solve(R.Ds{j},so);
                T(end).tests(end+1,:) = [lams(j) v.st v.rel_b v.psd_relmin v.t_mosek v.iter v.mosek_threads];
                fprintf('  %s lam %.7f (%.3f x pred %.7f): st %+d  %s  rel_b %.2e psd %+.1e cert_viol %.1e  MOSEK %.1f s  %d it  %d thr  peak %.0f MB\n', ...
                        nm,lams(j),lams(j)/pr,pr,v.st,v.why,v.rel_b,v.psd_relmin,v.cert_viol,v.t_mosek,v.iter, ...
                        v.mosek_threads,v.mem_peak);
                if v.st==1                          % for the '3dsum' certificate
                    x = v.x;    bscl = R.Ds{j}.bscl;    K = R.Ds{j}.K;  lam = lams(j); %#ok<NASGU>
                    save(fullfile(vopts.dir,sprintf('x_%s_%.6f.mat',nm,lams(j))),'x','bscl','K','lam','sp','nm');
                end
                clear v
            end
            T(end).t = toc(t0);     T(end).peak = peakMB();
            report(T(end));
            save(fullfile(vopts.dir,sprintf('pval_3d_%s_%s.mat',nm,char(datetime('now','Format','HHmmss')))),'T','P');
        end
    case '3dsum'
        T = sum_certificate(vopts);
    otherwise
        error('heatNd_poincare_validate:stage','STAGE is I0, 1d, 2d, 2d2, 3d0, 3d or 3dsum.')
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = sum_certificate(vopts)
% B2 as a certificate: certified per-direction points (faces of direction i
% only, Gram blocks [G0, face (i,a), face (i,b)]) summed into the all-faces
% L2 full SDP (blocks [G0, (1,a), (1,b), ..., (N,a), (N,b)], the order
% HEATND_POINCARE declares them), in original units (x * bscl).
O = heatNd_poincare_ops(heatNd_pie(3,0));   N = O.N;
t0 = tic;
R = heatNd_poincare(O,struct('form','L2','sel','full','d',1,'gen','faces'),struct('keepF',true));
F = R.F;    n = F.K.s(1);   n2 = n^2;
if F.K.f~=0 || numel(F.K.s)~=1+2*N || any(F.K.s~=n)
    error('heatNd_poincare_validate:sum','Unexpected block structure %s.',mat2str(F.K.s))
end
xf = zeros(F.nx,1);     lam = 0;    src = {};
for i = 1:N
    S = load(vopts.sumfiles{i});
    ok = strcmp(S.sp.form,'L2') && isequal(S.sp.sel,i) && isequal(S.K.s,[n n n]) && S.K.f==0 && ...
         isequal(S.sp.faces,[false(i-1,2); true(1,2); false(N-i,2)]);
    if ~ok,     error('heatNd_poincare_validate:sum','%s is not an L2 dir-%d faces-in-dir-%d point.',vopts.sumfiles{i},i,i),  end
    x = S.x*S.bscl;     lam = lam + S.lam;  src{end+1} = sprintf('dir %d at %.7f',i,S.lam);     %#ok<AGROW>
    xf(1:n2) = xf(1:n2) + x(1:n2);                              % G0
    for e = 1:2                                                 % face (i,e)
        o = (1 + 2*(i-1) + e - 1)*n2;
        xf(o+(1:n2)) = x(e*n2+(1:n2));
    end
end
bw = F.b1 + (lam-F.lam1)*F.db;
rel = norm(F.At'*xf - bw)/norm(bw);
em = inf;   eM = -inf;
for j = 1:numel(F.K.s)
    X = reshape(xf((j-1)*n2+(1:n2)),n,n);   ev = eig((X+X')/2);
    em = min(em,ev(1));     eM = max(eM,ev(end));
end
st = double(rel<=1e-6 && em/eM>=-1e-8);
fprintf(['3dsum: full L2 d=1 all faces (m %d, Ks %s) at lam = sum = %.7f (%s): rel_b %.2e, psd_relmin %+.1e ' ...
         '-> st %+d (HEATND_SOLVE gate; no solve), %.1f s\n'],F.m,mat2str(F.K.s),lam,strjoin(src,', '),rel,em/eM,st,toc(t0));
T = struct('lam',lam,'st',st,'rel_b',rel,'psd_relmin',em/eM,'m',F.m,'Ks',F.K.s,'src',{src});
save(fullfile(vopts.dir,sprintf('pval_3dsum_%s.mat',char(datetime('now','Format','HHmmss')))),'T');
end


function mb = peakMB()
mb = NaN;
try
    p = System.Diagnostics.Process.GetCurrentProcess();     p.Refresh();
    mb = double(p.PeakPagedMemorySize64)/2^20;
catch
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = run_case(O,sp,ref,vopts,full_search)
% Bisection (optional), objective + re-certification, soundness probe.
t0 = tic;
key = sprintf('N%d|%s|%s|%s|%d|%s',O.N,strjoin(O.bc,'x'),sp.form,num2str(sp.sel),sp.d,sp.gen);
rv = NaN;   rs = '';
if isKey(ref,key),  rv = ref(key);  rs = 'stock';   end
Ro = heatNd_poincare(O,sp,struct('mode','objective','verbose',false));
s = mkrow(O,Ro,rv,rs);
s.lamhat = Ro.obj.lamhat;   s.obj_st = Ro.obj.v.st;     s.lam_recert = Ro.obj.lam_recert;
if ~isempty(Ro.obj.recert),     s.recert_st = Ro.obj.recert.st;     end
Rt = heatNd_poincare(O,sp,struct('mode','test','lams',1.001*Ro.mu,'verbose',false, ...
                     'bisect',struct('retry',{{}},'verbose',false)));
s.sound_st = Rt.test.trace(end,2);
s.tests = trace7(Rt.test.trace);
if full_search
    Rb = heatNd_poincare(O,sp,struct('mode','bisect','verbose',false,'bisect', ...
                         struct('maxsolve',vopts.maxsolve,'rtol',vopts.rtol,'retry',{{'loose'}},'verbose',false)));
    s.lo = Rb.B.lo;     s.hi = Rb.B.hi;     s.stop = Rb.B.stop;     s.nsolve = Rb.B.nsolve;
    s.tests = [s.tests; trace7(Rb.B.trace)];
end
s.t = toc(t0);
end


function s = mkrow(O,R,rv,rs)
s = struct('N',O.N,'bc',{O.bc},'form',R.spec.form,'sel',R.spec.sel,'d',R.spec.d,'gen',R.spec.gen, ...
           'mu',R.mu,'ref',rv,'refsrc',rs,'m',R.shape.m,'Ks',R.shape.Ks,'nnz',R.shape.nnz, ...
           'lo',NaN,'hi',NaN,'stop','','nsolve',0,'lamhat',NaN,'obj_st',NaN,'lam_recert',NaN, ...
           'recert_st',NaN,'sound_st',NaN,'tests',zeros(0,7),'t',0,'peak',R.mem.peak_private);
end


function M = trace7(tr)
% [lam st rel_b psd_relmin t_mosek iter mosek_threads] of each attempt.
M = tr(tr(:,12)>0,[1 2 3 4 5 6 8]);
end


function report(s)
fprintf(['%s N=%d %-9s %-6s sel=%-4s d=%d %-7s m=%5d | exact %.6f ref %-9.6g (%s) | bracket [%.7f, %.7f) %s | ' ...
         'obj %.7f (st %+d) recert %.7f st %+d | 1.001 exact: st %+d | %.1f s, peak %.0f MB, MOSEK threads %s\n'], ...
        datestr(now,'HH:MM:SS'),s.N,strjoin(s.bc,'x'),s.form,num2str(s.sel),s.d,s.gen,s.m,s.mu,s.ref,s.refsrc, ...
        s.lo,s.hi,s.stop,s.lamhat,s.obj_st,s.lam_recert,s.recert_st,s.sound_st,s.t,s.peak, ...
        mat2str(unique(s.tests(:,7))'));
for k = 1:size(s.tests,1)
    if s.N==3 || size(s.tests,1)<=3
        fprintf('     lam %.7f  st %+d  rel_b %.2e  psd %+.1e  MOSEK %.1f s  %d it  %d thr\n',s.tests(k,:));
    end
end
end


function ref = stockref()
% Stock reference numbers (header), keyed N|bc|form|sel|d|gen.
ref = containers.Map();
L = {'N1|DD|L2|full|1|none',1.0e-4; 'N1|DD|L2|full|2|none',0.091; 'N1|DD|L2|full|3|none',1.051;
     'N1|DN|L2|full|1|none',1.8e-4; 'N1|DN|L2|full|2|none',0.109; 'N1|DN|L2|full|3|none',1.122;
     'N1|DD|L2|full|1|product',0.99997; 'N1|DD|L2|full|2|product',9.86899; 'N1|DD|L2|full|3|product',9.86960;
     'N1|DN|L2|full|1|product',0.99997; 'N1|DN|L2|full|2|product',2.46740; 'N1|DN|L2|full|3|product',2.46740;
     'N1|DD|energy|full|1|product',9.70065; 'N1|DD|energy|full|2|product',9.86960;
     'N1|DN|energy|full|1|product',2.46709; 'N1|DN|energy|full|2|product',2.46740;
     'N2|DDxDN|L2|full|1|none',0.00165; 'N2|DDxDN|L2|1|1|none',1.0e-4; 'N2|DDxDN|L2|2|1|none',3.3e-4;
     'N2|DDxDN|L2|full|1|product',0.00241; 'N2|DDxDN|L2|1|1|product',0.00132; 'N2|DDxDN|L2|2|1|product',0.00048;
     'N2|DDxDN|L2|full|1|faces',1.99994; 'N2|DDxDN|L2|1|1|faces',0.99997; 'N2|DDxDN|L2|2|1|faces',0.99997;
     'N2|DDxDN|energy|full|1|faces',12.03906; 'N2|DDxDN|energy|1|1|faces',9.71626; 'N2|DDxDN|energy|2|1|faces',2.46725;
     'N2|DDxDN|energy|full|1|product',0.0296*12.337006;
     'N2|DDxDN|L2|full|2|faces',12.33508; 'N2|DDxDN|L2|1|2|faces',9.86930; 'N2|DDxDN|L2|2|2|faces',2.46732};
for k = 1:size(L,1),    ref(L{k,1}) = L{k,2};   end
end


function x = tern(c,a,b),   if c, x = a; else, x = b; end,  end
