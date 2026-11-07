% psatzgrid.m -- the generator GRID, per the maintainer's correction.
%
% On [0,1] the 1-D choices are  {s(1-s)},  {s, 1-s},  or  {s, 1-s, s(1-s)}.
% In 2-D the certificate needs COMBINATIONS: one factor from each direction,
% so with three choices per direction there are up to NINE multiplier terms.
% This is NOT simply Putinar's Positivstellensatz on the four linear
% generators, which is why the grid is enumerated explicitly.
%
% P is pinned to the identity throughout, so the negativity certificate is the
% only unknown and reach is a property of the generator set alone.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
cuadmm_shadow('poslpivar');

LAMSTAR = 2*pi^2;
FR = [0.10 0.50 0.90 0.99 0.999 0.9999];

% 1-D factors per direction
A = @(x) x;  B = @(x) 1-x;  C = @(x) x.*(1-x);

SETS = {
 'CONTROL g1,g2 (handles)' , { @(x,y) C(x) , @(x,y) C(y) }
 'linear4 {s1,1-s1,s2,1-s2}', { @(x,y) A(x), @(x,y) B(x), @(x,y) A(y), @(x,y) B(y) }
 'BxB 4 corner products'   , { @(x,y) A(x).*A(y), @(x,y) A(x).*B(y), ...
                               @(x,y) B(x).*A(y), @(x,y) B(x).*B(y) }
 'CxC 9 products'          , { @(x,y) A(x).*A(y), @(x,y) A(x).*B(y), @(x,y) A(x).*C(y), ...
                               @(x,y) B(x).*A(y), @(x,y) B(x).*B(y), @(x,y) B(x).*C(y), ...
                               @(x,y) C(x).*A(y), @(x,y) C(x).*B(y), @(x,y) C(x).*C(y) }
};

fprintf('PG set|nterm|m|nblk|frac|lam|cert|rel_b|feasratio|numerr|t\n');
for k = 1:size(SETS,1)
    lab = SETS{k,1};  gl = SETS{k,2};  best = 0;  Sm = [];
    for f = FR
        lam = f*LAMSTAR;
        try
            t0 = tic;  prog = build2(lam,gl);
            so.solver='mosek';  evalc('sol = lpisolve(prog,so);');  tw = toc(t0);
            S = cuadmm_private('sdpshape',sol); if isempty(Sm), Sm=S; end
            I = sol.solinfo.info; xv = sol.solinfo.RRx(:);
            Atf=[];bf=[];
            for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
            rb = norm(full(Atf'*xv-bf))/norm(full(bf));
            triv = (abs(rb-1)<=1e-6) || norm(xv)<=1e-12;
            cert = (I.numerr==0) && (I.feasratio>0.5) && (rb<1e-6) && ~triv;
            fprintf('PG %s|%d|%d|%d|%.4f|%.4f|%d|%.3e|%+.4f|%d|%.1f\n', ...
                lab,numel(gl),S.m,S.nblk,f,lam,cert,rb,I.feasratio,I.numerr,tw);
            if cert, best=f; else, break; end
        catch ME
            fprintf('PG %s|ERR|%.4f|%s\n',lab,f,strrep(ME.message,newline,' '));
            break
        end
    end
    if isempty(Sm), Sm.m=-1; Sm.nblk=-1; end
    fprintf('PGREACH %s|nterm=%d|reach=%.4f|m=%d|nblk=%d\n',lab,numel(gl),best,Sm.m,Sm.nblk);
end
fprintf('PGDONE\n');

cuadmm_shadow('off');   % remove the shadow: it stays active for the whole session otherwise

function prog = build2(lam,glist)
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],[0,1;0,1]);
PIE = convert([diff(x,t,1)==diff(x,s1,2)+diff(x,s2,2)+lam*x;
               subs(x,s1,0)==0; subs(x,s1,1)==0;
               subs(x,s2,0)==0; subs(x,s2,1)==0]);
PIE = initialize(PIE);
Top = PIE.T; Aop = PIE.A;
st = lpisettings('heavy');  st.settings_2d.eppos = 1e-2*[1;1;1;1];
s2d = st.settings_2d;
dx=s2d.LF_deg.dx; dy=s2d.LF_deg.dy; d2=s2d.LF_deg.d2;
s2d.eq_deg.dx = {1+dx{1}; 1+dx{2}; 1+dx{3}};
s2d.eq_deg.dy = {1+dy{1}, 1+dy{2}, 1+dy{3}};
s2d.eq_deg.d2 = cellfun(@(c) 1+c, d2, 'UniformOutput', false);
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
np  = Top.dim(:,1);
Iop = opvar2d(eye(sum(np)),[np np],PIE.dom,PIE.vars);
Qop = (Aop'*(Iop*Top))' + Aop'*(Iop*Top);
Qop = clean_opvar(Qop,1e-12);
eq_opts = get_eq_opts_2D(Qop,s2d.eq_opts,1e-12);
[prog,Qeop] = poslpivar_2d(prog,Qop.dim,s2d.eq_deg,eq_opts);
for i = 1:numel(glist)
    o = eq_opts;  o.psatz = glist{i};
    [prog,Qe] = poslpivar_2d(prog,Qop.dim,s2d.eq_deg,o);
    Qeop = Qeop + Qe;
end
prog = lpi_eq_2d(prog,Qop+Qeop,'symmetric');
end
