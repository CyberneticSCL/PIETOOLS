% normgen.m -- the elongated-domain failure: scaling, or geometry?
%
% MEASURED: linear4 reaches 0.99 on [0,1]^2 and 0.90 on the OFFSET [0.3,1.3]^2
% (still unit length), but 0.00 on [0,1]x[0,2].  The offset case is the
% diagnostic: translation is fine, LENGTH is not.  On [0,2] the generators
% (s2-0) and (2-s2) reach magnitude 2 while the s1 generators reach 1, so the
% Gram blocks are mismatched in scale.  PREDICTION: normalising each generator
% to (s-a)/(b-a), which is in [0,1] on any interval, should restore it.
% If normalisation fixes it, the failure is CONDITIONING and the patch should
% ship normalised generators.  If it does not, the obstruction is geometric
% and the result is genuinely limited to near-square domains.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
cuadmm_shadow('poslpivar');   % function-handle multipliers

FR = [0.10 0.25 0.50 0.90 0.99];
% dom, label, lamstar
C = { [0 1;0 1], 'unit  1x1', 2*pi^2
      [0 1;0 2], 'elong 1x2', pi^2*(1+1/4)
      [0 1;0 4], 'elong 1x4', pi^2*(1+1/16) };

fprintf('NG dom|gens|m|nblk|frac|lam|cert|rel_b|feasratio|numerr\n');
for ic = 1:size(C,1)
    dom = C{ic,1}; lab = C{ic,2}; lamstar = C{ic,3};
    a1=dom(1,1); b1=dom(1,2); a2=dom(2,1); b2=dom(2,2);
    RAW  = { @(x,y) x-a1, @(x,y) b1-x, @(x,y) y-a2, @(x,y) b2-y };
    NORM = { @(x,y) (x-a1)/(b1-a1), @(x,y) (b1-x)/(b1-a1), ...
             @(x,y) (y-a2)/(b2-a2), @(x,y) (b2-y)/(b2-a2) };
    for g = {{'raw',RAW},{'norm',NORM}}
        gn = g{1}{1}; gl = g{1}{2}; best=0; mm=-1; nb=-1;
        for f = FR
            try
                prog = build6(dom,f*lamstar,gl);
                so.solver='mosek'; evalc('sol = lpisolve(prog,so);');
                S = cuadmm_private('sdpshape',sol); if mm<0, mm=S.m; nb=S.nblk; end
                I = sol.solinfo.info; xv = sol.solinfo.RRx(:);
                Atf=[];bf=[];
                for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
                rb = norm(full(Atf'*xv-bf))/norm(full(bf));
                triv = (abs(rb-1)<=1e-6)||norm(xv)<=1e-12;
                cert = (I.numerr==0)&&(I.feasratio>0.5)&&(rb<1e-6)&&~triv;
                fprintf('NG %s|%s|%d|%d|%.2f|%.4f|%d|%.3e|%+.4f|%d\n', ...
                    lab,gn,S.m,S.nblk,f,f*lamstar,cert,rb,I.feasratio,I.numerr);
                if cert, best=f; else, break; end
            catch ME
                fprintf('NG %s|%s|ERR|%.2f|%s\n',lab,gn,f,strrep(ME.message,newline,' ')); break
            end
        end
        fprintf('NGREACH %s|%s|reach=%.2f|m=%d|nblk=%d\n',lab,gn,best,mm,nb);
    end
end

% --- and the well-posed Neumann variant: Dirichlet at one end, Neumann at the
%     other in s2 (pure NN leaves the constant mode in the kernel, so the
%     boundary operator is genuinely not invertible -- that build failure was
%     the physics, not the harness).
fprintf('NG --- mixed D/N in s2, lam* = pi^2*(1 + 1/4) for the quarter-wave mode\n');
dom = [0 1;0 1]; lamstar = pi^2*(1+0.25);
RAW = { @(x,y) x, @(x,y) 1-x, @(x,y) y, @(x,y) 1-y };
for gn = {'none','linear4'}
    gl = {}; if strcmp(gn{1},'linear4'), gl = RAW; end
    best=0; mm=-1;
    for f = FR
        try
            prog = build6(dom,f*lamstar,gl,'DN');
            so.solver='mosek'; evalc('sol = lpisolve(prog,so);');
            S = cuadmm_private('sdpshape',sol); if mm<0, mm=S.m; end
            I = sol.solinfo.info; xv = sol.solinfo.RRx(:);
            Atf=[];bf=[];
            for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
            rb = norm(full(Atf'*xv-bf))/norm(full(bf));
            triv = (abs(rb-1)<=1e-6)||norm(xv)<=1e-12;
            cert = (I.numerr==0)&&(I.feasratio>0.5)&&(rb<1e-6)&&~triv;
            fprintf('NG DN-s2|%s|%d|%.2f|%d|%.3e|%+.4f|%d\n',gn{1},S.m,f,cert,rb,I.feasratio,I.numerr);
            if cert, best=f; else, break; end
        catch ME
            fprintf('NG DN-s2|%s|ERR|%.2f|%s\n',gn{1},f,strrep(ME.message,newline,' ')); break
        end
    end
    fprintf('NGREACH DN-s2|%s|reach=%.2f|m=%d\n',gn{1},best,mm);
end
fprintf('NGDONE\n');

cuadmm_shadow('off');   % remove the shadow: it stays active for the whole session otherwise

function prog = build6(dom,lam,glist,bc2)
if nargin<4, bc2 = 'DD'; end
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],dom);
sys = diff(x,t,1) == diff(x,s1,2)+diff(x,s2,2)+lam*x;
sys = [sys; subs(x,s1,dom(1,1))==0; subs(x,s1,dom(1,2))==0];
switch bc2
    case 'DD', sys = [sys; subs(x,s2,dom(2,1))==0; subs(x,s2,dom(2,2))==0];
    case 'DN', sys = [sys; subs(x,s2,dom(2,1))==0; subs(diff(x,s2,1),s2,dom(2,2))==0];
end
PIE = initialize(convert(sys));
Top = PIE.T; Aop = PIE.A;
st = lpisettings('heavy'); st.settings_2d.eppos = 1e-2*[1;1;1;1];
s2d = st.settings_2d;
dx=s2d.LF_deg.dx; dy=s2d.LF_deg.dy; d2=s2d.LF_deg.d2;
eqd.dx = {1+dx{1}; 1+dx{2}; 1+dx{3}};
eqd.dy = {1+dy{1}, 1+dy{2}, 1+dy{3}};
eqd.d2 = cellfun(@(c) 1+c, d2, 'UniformOutput', false);
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
np = Top.dim(:,1);
Iop = opvar2d(eye(sum(np)),[np np],PIE.dom,PIE.vars);
Qop = clean_opvar((Aop'*(Iop*Top))' + Aop'*(Iop*Top),1e-12);
eq_opts = get_eq_opts_2D(Qop,s2d.eq_opts,1e-12);
[prog,Qeop] = poslpivar_2d(prog,Qop.dim,eqd,eq_opts);
for i = 1:numel(glist)
    o = eq_opts; o.psatz = glist{i};
    [prog,Qe] = poslpivar_2d(prog,Qop.dim,eqd,o);
    Qeop = Qeop + Qe;
end
prog = lpi_eq_2d(prog,Qop+Qeop,'symmetric');
end
