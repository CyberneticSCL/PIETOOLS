% odegen.m -- generator comparison on the ONLY plant found with a nonzero
% POINTWISE MULTIPLIER cell: a 2-D PDE coupled to a finite-dimensional ODE.
%
%   x_t = lap(x) + lam*x + X          on [0,1]^2, Dirichlet all four edges
%   X_t = -2X + int int x
%
% Measured: Top.dim(:,1) = [1 0 0 1], and Q.R00 is NONZERO -- the first plant
% in this whole investigation where the multiplier cell exists.  Two things
% are being tested at once and must be reported separately:
%  (a) GENERALITY: does linear4 still beat the product and no-psatz?
%  (b) STRUCTURE: this is the n1>0 regime where the degree-balance rule is
%      known to be over-determined (single-factor blocks need offset e, not
%      e/2, from the same caps).  The generator CHOICE may still be fine.
% lam* is not analytic for the coupled system, so this is a WITHIN-PLANT
% comparison only: largest lam certified by each generator set.
cuadmm_path;
LAM = [1 2 5 10 15 19];
SETS = {'none',[] ; 'product',1 ; 'linear4',[3 4 5 6]};
fprintf('OG set|m|nblk|Ns|lam|cert|rel_b|feasratio|numerr|t\n');
for is = 1:size(SETS,1)
    nm = SETS{is,1}; codes = SETS{is,2}; best = -inf; mm=-1; nb=-1; ks=[];
    for lam = LAM
        try
            t0=tic; prog = b9(lam,codes);
            so.solver='mosek'; evalc('sol = lpisolve(prog,so);'); tw=toc(t0);
            S=cuadmm_private('sdpshape',sol); if mm<0, mm=S.m; nb=S.nblk; ks=S.Ks; end
            I=sol.solinfo.info; xv=sol.solinfo.RRx(:);
            Atf=[];bf=[];
            for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
            rb=norm(full(Atf'*xv-bf))/norm(full(bf));
            triv=(abs(rb-1)<=1e-6)||norm(xv)<=1e-12;
            cert=(I.numerr==0)&&(I.feasratio>0.5)&&(rb<1e-6)&&~triv;
            fprintf('OG %s|%d|%d|[%s]|%g|%d|%.3e|%+.4f|%d|%.1f\n', ...
                nm,S.m,S.nblk,strtrim(num2str(S.Ks)),lam,cert,rb,I.feasratio,I.numerr,tw);
            if cert, best=lam; else, break; end
        catch ME
            fprintf('OG %s|ERR|%g|%s\n',nm,lam,strrep(ME.message,newline,' ')); break
        end
    end
    if isinf(best), bs='none'; else, bs=sprintf('%g',best); end
    fprintf('OGREACH %s|lam_max=%s|m=%d|nblk=%d|Ns=[%s]\n',nm,bs,mm,nb,strtrim(num2str(ks)));
end
fprintf('OGDONE\n');

function prog = b9(lam,codes)
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],[0 1;0 1]);
X = pde_var('state',1,[],[]);
sys = [diff(x,t,1) == diff(x,s1,2)+diff(x,s2,2)+lam*x + X;
       diff(X,t,1) == -2*X + int(int(x,s1,[0,1]),s2,[0,1]);
       subs(x,s1,0)==0; subs(x,s1,1)==0;
       subs(x,s2,0)==0; subs(x,s2,1)==0];
PIE = initialize(convert(sys));
Top=PIE.T; Aop=PIE.A;
st=lpisettings('heavy'); st.settings_2d.eppos=1e-2*[1;1;1;1];
s2d=st.settings_2d;
dx=s2d.LF_deg.dx; dy=s2d.LF_deg.dy; d2=s2d.LF_deg.d2;
eqd.dx={1+dx{1};1+dx{2};1+dx{3}};
eqd.dy={1+dy{1},1+dy{2},1+dy{3}};
eqd.d2=cellfun(@(c) 1+c, d2,'UniformOutput',false);
prog=lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
np=Top.dim(:,1);
Iop=opvar2d(eye(sum(np)),[np np],PIE.dom,PIE.vars);
Qop=clean_opvar((Aop'*(Iop*Top))'+Aop'*(Iop*Top),1e-12);
eq_opts=get_eq_opts_2D(Qop,s2d.eq_opts,1e-12);
[prog,Qeop]=poslpivar_2d(prog,Qop.dim,eqd,eq_opts);
for p=codes
    o=eq_opts; o.psatz=p;
    [prog,Qe]=poslpivar_2d(prog,Qop.dim,eqd,o);
    Qeop=Qeop+Qe;
end
prog=lpi_eq_2d(prog,Qop+Qeop,'symmetric');
end
