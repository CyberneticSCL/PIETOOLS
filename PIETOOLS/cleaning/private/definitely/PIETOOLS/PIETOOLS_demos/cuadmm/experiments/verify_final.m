% verify_final.m -- regression check of the UPDATED repo patch (normalised
% generators, codes 3-6).  NO shadow on the path.  Must reproduce:
%   [0,1]^2  reach 0.99  (normalisation is a no-op on a unit side)
%   [0,1]x[0,2] reach 0.99 (raw was 0.10 -- this is the fix)
cuadmm_path;
fprintf('VF poslpivar_2d resolves to %s\n', which('poslpivar_2d'));
FR = [0.10 0.25 0.50 0.90 0.99];
C = { [0 1;0 1], 'unit 1x1' , 2*pi^2
      [0 1;0 2], 'elong 1x2', pi^2*(1+1/4) };
for ic = 1:size(C,1)
    dom=C{ic,1}; lab=C{ic,2}; ls=C{ic,3}; best=0; mm=-1;
    for f = FR
        try
            prog = b8(dom,f*ls,[3 4 5 6]);
            so.solver='mosek'; evalc('sol = lpisolve(prog,so);');
            S=cuadmm_private('sdpshape',sol); if mm<0, mm=S.m; end
            I=sol.solinfo.info; xv=sol.solinfo.RRx(:);
            Atf=[];bf=[];
            for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
            rb=norm(full(Atf'*xv-bf))/norm(full(bf));
            triv=(abs(rb-1)<=1e-6)||norm(xv)<=1e-12;
            cert=(I.numerr==0)&&(I.feasratio>0.5)&&(rb<1e-6)&&~triv;
            fprintf('VF %s|%.2f|%d|%.3e|%+.4f|%d\n',lab,f,cert,rb,I.feasratio,I.numerr);
            if cert, best=f; else, break; end
        catch ME
            fprintf('VF %s|ERR|%.2f|%s\n',lab,f,strrep(ME.message,newline,' ')); break
        end
    end
    fprintf('VFREACH %s|reach=%.2f|m=%d\n',lab,best,mm);
end
fprintf('VFDONE\n');

function prog = b8(dom,lam,codes)
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],dom);
PIE = initialize(convert([diff(x,t,1)==diff(x,s1,2)+diff(x,s2,2)+lam*x;
      subs(x,s1,dom(1,1))==0; subs(x,s1,dom(1,2))==0;
      subs(x,s2,dom(2,1))==0; subs(x,s2,dom(2,2))==0]));
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
