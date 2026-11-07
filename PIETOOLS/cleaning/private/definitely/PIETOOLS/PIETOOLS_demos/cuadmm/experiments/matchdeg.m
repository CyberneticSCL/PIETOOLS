% matchdeg.m -- THE CONFOUND CONTROL.
%
% In my generator table every Sigma_k block was built at the FULL eq_deg and
% THEN multiplied by g.  Since Qop + Qeop == 0 pins deg(Qop), the product
% (deg 4) overshoots the constraint by 4 while each linear generator (deg 1)
% overshoots by 1.  So "which generator" was confounded with "how far the
% block overshoots".  If the product certifies once its block degree is
% lowered to land on deg(Qop), then "the stock product certifies nothing" is
% an artifact of the degree convention, not a property of the generator.
%
% Counter-pressure, already measured: ANY reduction collapses a 2-D block from
% 264 to 104 monomials because reduce_joint_degs drags subset caps down.  So
% two effects fight -- less overshoot versus lost span.  Measure, do not argue.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
cuadmm_shadow('poslpivar');
LAMSTAR = 2*pi^2;
FR = [0.10 0.25 0.50 0.90 0.99];

V = { 'product  @eq_deg'   , 1,          0
      'product  @eq_deg-1' , 1,         -1
      'product  @eq_deg-2' , 1,         -2
      'linear4  @eq_deg'   , [3 4 5 6],  0
      'linear4  @eq_deg-1' , [3 4 5 6], -1 };

fprintf('MD variant|m|nblk|Ns|frac|cert|rel_b|feasratio|numerr\n');
for v = 1:size(V,1)
    lab=V{v,1}; codes=V{v,2}; off=V{v,3}; best=0; mm=-1; nb=-1; ks=[];
    for f = FR
        try
            prog = build7(f*LAMSTAR,codes,off);
            so.solver='mosek'; evalc('sol = lpisolve(prog,so);');
            S = cuadmm_private('sdpshape',sol); if mm<0, mm=S.m; nb=S.nblk; ks=S.Ks; end
            I = sol.solinfo.info; xv = sol.solinfo.RRx(:);
            Atf=[];bf=[];
            for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
            rb = norm(full(Atf'*xv-bf))/norm(full(bf));
            triv = (abs(rb-1)<=1e-6)||norm(xv)<=1e-12;
            cert = (I.numerr==0)&&(I.feasratio>0.5)&&(rb<1e-6)&&~triv;
            fprintf('MD %s|%d|%d|[%s]|%.2f|%d|%.3e|%+.4f|%d\n', ...
                lab,S.m,S.nblk,strtrim(num2str(S.Ks)),f,cert,rb,I.feasratio,I.numerr);
            if cert, best=f; else, break; end
        catch ME
            fprintf('MD %s|ERR|%.2f|%s\n',lab,f,strrep(ME.message,newline,' ')); break
        end
    end
    fprintf('MDREACH %s|reach=%.2f|m=%d|nblk=%d|Ns=[%s]\n',lab,best,mm,nb,strtrim(num2str(ks)));
end
fprintf('MDDONE\n');

cuadmm_shadow('off');   % remove the shadow: it stays active for the whole session otherwise

function prog = build7(lam,codes,off)
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],[0 1;0 1]);
PIE = initialize(convert([diff(x,t,1)==diff(x,s1,2)+diff(x,s2,2)+lam*x;
      subs(x,s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0]));
Top=PIE.T; Aop=PIE.A;
st = lpisettings('heavy'); st.settings_2d.eppos = 1e-2*[1;1;1;1];
s2d = st.settings_2d;
dx=s2d.LF_deg.dx; dy=s2d.LF_deg.dy; d2=s2d.LF_deg.d2;
eqd.dx = {1+dx{1}; 1+dx{2}; 1+dx{3}};
eqd.dy = {1+dy{1}, 1+dy{2}, 1+dy{3}};
eqd.d2 = cellfun(@(c) 1+c, d2, 'UniformOutput', false);
pd = eqd; if off~=0, pd = bump(eqd,off); end
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
np = Top.dim(:,1);
Iop = opvar2d(eye(sum(np)),[np np],PIE.dom,PIE.vars);
Qop = clean_opvar((Aop'*(Iop*Top))' + Aop'*(Iop*Top),1e-12);
eq_opts = get_eq_opts_2D(Qop,s2d.eq_opts,1e-12);
[prog,Qeop] = poslpivar_2d(prog,Qop.dim,eqd,eq_opts);
for p = codes
    o = eq_opts; o.psatz = p;
    [prog,Qe] = poslpivar_2d(prog,Qop.dim,pd,o);
    Qeop = Qeop + Qe;
end
prog = lpi_eq_2d(prog,Qop+Qeop,'symmetric');
end
function d = bump(d,k)
d.dx = cellfun(@(c) max(c+k,0), d.dx, 'UniformOutput', false);
d.dy = cellfun(@(c) max(c+k,0), d.dy, 'UniformOutput', false);
d.d2 = cellfun(@(c) max(c+k,0), d.d2, 'UniformOutput', false);
end
