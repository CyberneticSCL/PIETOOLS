% hinf_bisect.m -- BISECTION on gamma instead of minimising it.
%
% Maintainer's observation: with an objective, IPM residuals often fail to
% converge while gamma still converges to the right value -- residuals mix
% feasibility and optimality -- and bisection is the more reliable method.
%
% This matters for the solver evaluation too: fixing gamma turns the H-infinity
% LPI into a PURE FEASIBILITY program (c = 0, and the gamma>=0 lpi_ineq
% disappears with it), which is exactly the class where cuADMM is strong and
% where its dual-infeasibility stall does not arise.  PIETOOLS_Hinf_gain
% documents this path itself ("a specific gain test ... results in a
% feasibility test instead of an optimization problem").
%
% Reference from the direct optimisation, same plant and settings:
%   light gamma = 0.05222589, heavy gamma = 0.05211894 (Mosek).
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
ROOT = fullfile(cuadmm_outdir(),'cu_bis'); if ~exist(ROOT,'dir'), mkdir(ROOT); end

for pre = {'light','heavy'}
    st = lpisettings(pre{1});  st.sos_opts.solver = 'mosek';
    % ---- sanity: bracket the known answer
    for g = [0.04 0.0522 0.06 0.10]
        [ok,rb,ne,fr,S] = feas(pre{1},st,g);
        fprintf('HB %-5s gam=%.4f feas=%d rel_b=%.3e numerr=%d feasratio=%+.3f m=%d Kf=%d Ns=[%s]\n', ...
                pre{1},g,ok,rb,ne,fr,S.m,S.Kf,strtrim(num2str(S.Ks)));
    end
    % ---- bisect
    lo = 0.04; hi = 0.10; nit = 0;
    while (hi-lo) > 1e-6
        mid = 0.5*(lo+hi);
        ok = feas(pre{1},st,mid);
        if ok, hi = mid; else, lo = mid; end
        nit = nit + 1;
    end
    fprintf('HB %-5s BISECTION gamma = %.8f  (%d solves)\n', pre{1}, hi, nit);
end
fprintf('HBDONE\n');

function [ok,rb,ne,fr,S,prog] = feas(pre,st,gamval)
pvar s t
x = pde_var('state',1,s,[0,1]);
w = pde_var('input',1);
z = pde_var('output',1);
PIE = initialize(convert([diff(x,t,1)==diff(x,s,2)+2*x+s*w;
      z == int(x,s,[0,1]); subs(x,s,0)==0; subs(x,s,1)==0]));
Top=PIE.T; Aop=PIE.A; Bwop=PIE.Bw; Czop=PIE.Cz; Dzwop=PIE.Dzw;
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
% gam is a NUMBER here -- no lpidecvar, no lpi_ineq, no lpisetobj.
gam = gamval;
[prog,R1op] = poslpivar(prog,Top.dim,st.dd1,st.options1);
if st.override1~=1
    [prog,P2op] = poslpivar(prog,Top.dim,st.dd12,st.options12);
    Rop = R1op+P2op;
else
    Rop = R1op;
end
Qdeg = get_lpivar_degs(Rop,Top);
[prog,Qop] = lpivar(prog,Top.dim,Qdeg);
prog = lpi_eq(prog,Top'*Qop-Rop);
Iw = mat2opvar(eye(size(Bwop,2)),Bwop.dim(:,2),PIE.vars,PIE.dom);
Iz = mat2opvar(eye(size(Czop,1)),Czop.dim(:,1),PIE.vars,PIE.dom);
Dop = [-gam*Iw,    Dzwop',   Bwop'*Qop;
        Dzwop,     -gam*Iz,  Czop;
        Qop'*Bwop, Czop',    Aop'*Qop+Qop'*Aop];
[prog,De1op] = poslpivar(prog,Dop.dim,st.dd2,st.options2);
if st.override2~=1
    [prog,De2op] = poslpivar(prog,Dop.dim,st.dd3,st.options3);
    Deop = De1op+De2op;
else
    Deop = De1op;
end
prog = lpi_eq(prog,Deop+Dop,'symmetric');
S = cuadmm_private('sdpshape',prog);
evalc('sol = lpisolve(prog,st.sos_opts);');
I = sol.solinfo.info; xv = sol.solinfo.RRx(:);
Atf=[];bf=[];
for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
rb = norm(full(Atf'*xv-bf))/norm(full(bf));
ne = I.numerr; fr = I.feasratio;
triv = (abs(rb-1)<=1e-6)||norm(xv)<=1e-12;
ok = (I.numerr==0)&&(I.feasratio>0.5)&&(rb<1e-6)&&~triv;
prog = sol;
end
