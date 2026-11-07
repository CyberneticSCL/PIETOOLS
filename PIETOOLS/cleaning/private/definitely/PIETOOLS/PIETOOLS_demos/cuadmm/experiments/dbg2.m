% dbg2.m -- two diagnostics the baseline needs before the GPU arm runs.
%
% (A) psd_relmin is NEGATIVE (-0.09 to -0.53) on exactly the rows with K.f>0,
%     while the K.f=0 rows read +1e-11. Two candidate causes and they demand
%     different fixes, so guessing is not an option:
%       H1 SLICING. solinfo.RRx is in decvartable order and solinfo.x is
%          SeDuMi's cone order; they differ once there are free variables.
%          Slicing the wrong one by cone blocks would misread every Gram.
%          Testable directly: is RR the identity, and do x and RRx/bscl agree?
%       H2 SCALE. For the feasibility cases ||b|| is eppos2 ~ 1e-6, so the whole
%          solution sits at ~1e-5 and a relative PSD margin is measured against
%          the noise floor. Testable by raising eppos2 and re-reading the margin;
%          if H2 holds the margin improves with no code change.
%     H1 and H2 make OPPOSITE predictions for the K.f=0 rows, which are clean
%     under H1's explanation but should be equally noisy under H2's if their
%     ||b|| is equally small -- and stabpde_rd1 has ||b||=8e-6 and reads +9e-12.
%     That already leans away from H2, so measure rather than assume.
%
% (B) the four 2-D input/output cases die with "Too many input arguments" and
%     the message alone does not say which call is at fault.

cuadmm_path;
fprintf('=== A: slicing vs scale ===\n');
for eppos2 = [1e-6 1e-4 1e-2 1]
    clear stateNameGenerator
    pvar s t
    x = pde_var(1,s,[0,1]);
    PIE = initialize(convert([diff(x,t,1)==diff(x,s,2)+0.5*pi^2*x; ...
                              subs(x,s,0)==0; subs(x,s,1)==0]));
    st = lpisettings('light');
    st.sos_opts.solver='mosek'; st.sos_opts.simplify=false;
    st.eppos2 = eppos2;
    evalc('[sol,P] = PIETOOLS_PIE2PDEstability(PIE,st);');
    S = cuadmm_private('sdpshape',sol);
    RR = mkRR(sol);
    xr = sol.solinfo.RRx(:);
    xc = sol.solinfo.x(:);
    nvar = S.Kf + sum(double(S.Ks).^2);
    isI  = (size(RR,1)==size(RR,2)) && (nnz(RR-speye(size(RR,1)))==0);
    % compare the two vectors after undoing the b-normalisation
    bf=[]; for q=1:sol.expr.num, bf=[bf;sol.expr.b{q}]; end %#ok<AGROW>
    bscl = norm(full(bf));
    agree = norm(xr - xc*bscl)/max(norm(xr),eps);
    [pmR,prR] = psdmargin(xr,S);
    [pmC,prC] = psdmargin(xc,S);
    fprintf('A eppos2=%-6g normb=%-10.3e RR_is_I=%d nvar=%d len_x=%d len_RRx=%d agree=%.2e | relmin_RRx=%+.4e relmin_x=%+.4e | absmin_x=%+.3e\n', ...
        eppos2,bscl,isI,nvar,numel(xc),numel(xr),agree,prR,prC,pmC);
end

fprintf('=== B: 2-D io stack ===\n');
try
    [sol,M] = cuadmm_private('bl_b_io2','Hinf_gain_2D','light');   %#ok<ASGLU>
    fprintf('B ok m=%d\n',cuadmm_private('sdpshape',sol).m);
catch ME
    fprintf('B ERR %s\n',ME.message);
    for k=1:numel(ME.stack)
        fprintf('B   at %s line %d\n',ME.stack(k).name,ME.stack(k).line);
    end
end
fprintf('DBG2DONE\n');

function [pmin,prel] = psdmargin(x,S)
off = S.Kf;  pmin = inf;  pmax = -inf;
for k = 1:numel(S.Ks)
    N = double(S.Ks(k));
    if off+N^2 > numel(x), pmin=NaN; prel=NaN; return; end
    Xk = reshape(x(off+(1:N^2)),N,N);  ev = eig((Xk+Xk')/2);
    pmin = min(pmin,min(ev));  pmax = max(pmax,max(ev));  off = off + N^2;
end
prel = pmin/max(pmax,eps);
end

function RR = mkRR(prog)
RR = speye(prog.var.idx{1}-1);
for i=1:prog.var.num
    sz = prog.var.idx{i+1}-prog.var.idx{i};
    switch prog.var.type{i}
        case 'poly', RR = spantiblkdiag(RR,speye(sz));
        case 'sos',  RR = spblkdiag(RR,speye(sz));
    end
end
for i=1:prog.extravar.num
    RR = spblkdiag(RR,speye(prog.extravar.idx{i+1}-prog.extravar.idx{i}));
end
end
