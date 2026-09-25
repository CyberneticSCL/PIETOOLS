% poinc1d.m -- how tight can the NEGATIVITY certificate get, with the Lyapunov
% operator FIXED to the identity?
%
% Maintainer's framing: for the homogeneous heat equation P=I works, so the
% levers are all in the negativity test, and the question is really how tightly
% the LPI can certify the Poincare inequality.
%
% Concretely, for  x_t = x_ss + lam*x  on [0,1] with Dirichlet BCs and
% V = int x^2:
%     Vdot = -2 int x_s^2 + 2 lam int x^2
% so stability up to lam* is EXACTLY the Poincare inequality
%     int x_s^2 >= lam* int x^2,     lam* = pi^2 = 9.8696.
% Fixing P=I removes the Lyapunov search entirely, so the reach measured here
% is a property of the negativity certificate ALONE.
%
% Reported as the fraction lam/lam* reached, which is the quantity that
% matters and is comparable across settings and across dimension.

cuadmm_path;
LAMSTAR = pi^2;
FR = [0.50 0.75 0.90 0.99 0.999 0.9999];

% CORRECTION: override2 = 0 means the psatz term IS added (Deop = De1op+De2op)
% and that is the LIGHT DEFAULT, so an earlier 'psatzD(ovr2=0)' variant simply
% re-tested the default.  override2 = 1 is what turns psatz OFF.
V = { 'Dup1 psatzON(default)' , @(st) st
      'Dup1 psatzOFF'         , @(st) setf(st,'override2',1)
      'Dup2 psatzON'          , @(st) dup(st,2)
      'Dup2 psatzOFF'         , @(st) setf(dup(st,2),'override2',1)
      'Dup3 psatzOFF'         , @(st) setf(dup(st,3),'override2',1) };

fprintf('P1 variant|m|nblk|Ns|frac|lam|cert|rel_b|psd_relmin|feasratio|t\n');
for v = 1:size(V,1)
    lab = V{v,1};  best = 0;  Sm = [];
    for f = FR
        lam = f*LAMSTAR;
        try
            st = lpisettings('light');
            st.eppos2 = 1e-2;  st.eppos = 1e-2;
            st.sos_opts.solver = 'sedumi';   % 1-D: SeDuMi is fine and fast
            st = V{v,2}(st);
            [prog,Dop,Deop] = build_fixedP(lam,st);
            evalc('sol = lpisolve(prog,st.sos_opts);');
            S = cuadmm_private('sdpshape',sol);  if isempty(Sm), Sm = S; end
            I = sol.solinfo.info;  xv = sol.solinfo.RRx(:);
            Atf=[];bf=[];
            for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
            rb = norm(full(Atf'*xv-bf))/norm(full(bf));
            [pmin,pmax] = gramspec(xv,S);
            triv = (abs(rb-1)<=1e-6) || norm(xv)<=1e-12;
            cert = (I.numerr==0) && (I.feasratio>0.5) && (rb<1e-6) && ~triv;
            fprintf('P1 %s|%d|%d|[%s]|%.2f|%.4f|%d|%.3e|%+.3e|%+.4f|%.2f\n', ...
                lab,S.m,S.nblk,strtrim(num2str(S.Ks)),f,lam,cert,rb, ...
                pmin/max(pmax,eps),I.feasratio,0);
            if cert, best = f; else, break; end
        catch ME
            fprintf('P1 %s|ERR|%.2f|%s\n',lab,f,strrep(ME.message,newline,' '));
            break
        end
    end
    fprintf('P1REACH %s|reach=%.2f|m=%d\n',lab,best,Sm.m);
end
fprintf('P1DONE\n');


function [prog,Dop,Deop] = build_fixedP(lam,st)
% The stability LPI with the Lyapunov operator PINNED to the identity.
pvar s t
x   = pde_var('state',1,s,[0,1]);
PIE = convert([diff(x,t,1)==diff(x,s,2)+lam*x; subs(x,s,0)==0; subs(x,s,1)==0]);
PIE = initialize(PIE);
Top = PIE.T;  Aop = PIE.A;

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);

% P = I : no decision variables at all on the Lyapunov side.
opvar Iop;
Iop.I = PIE.dom;  Iop.var1 = PIE.vars(1,1);  Iop.var2 = PIE.vars(1,2);
Iop.dim = Top.dim;
n0 = Top.dim(1,1);  n1 = Top.dim(2,1);
if n0>0, Iop.P = eye(n0); end
if n1>0, Iop.R.R0 = eye(n1); end

Dop = Aop'*(Iop*Top) + (Iop*Top)'*Aop;
Dop = clean_opvar(Dop,1e-12);

[prog,De1op] = poslpivar(prog,Dop.dim,st.dd2,st.options2);
if st.override2~=1
    [prog,De2op] = poslpivar(prog,Dop.dim,st.dd3,st.options3);
    Deop = De1op+De2op;
else
    Deop = De1op;
end
prog = lpi_eq(prog,Dop+Deop,'symmetric');
end

function [pmin,pmax] = gramspec(x,S)
pmin = inf; pmax = -inf; off = S.Kf;
for k = 1:numel(S.Ks)
    N = S.Ks(k);  Xk = reshape(x(off+(1:N^2)),N,N);  Xk = (Xk+Xk')/2;
    ev = eig(Xk);  pmin = min(pmin,min(ev));  pmax = max(pmax,max(ev));
    off = off + N^2;
end
end

function st = dup(st,Dup)
n1=1; n2=1; n3=1; n4=n2+n3;
st.dd2 = {n1+Dup,   [n2+Dup-1, n3+Dup,   n4+Dup  ], [n2+Dup-1, n3+Dup,   n4+Dup  ]};
st.dd3 = {n1+Dup-1, [n2+Dup-2, n3+Dup-1, n4+Dup-1], [n2+Dup-2, n3+Dup-1, n4+Dup-1]};
end

function st = setf(st,f,v), st.(f) = v; end
