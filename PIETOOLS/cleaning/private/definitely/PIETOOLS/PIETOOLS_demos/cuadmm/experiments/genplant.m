% genplant.m -- does the 2-D linear-generator psatz result survive beyond the
% symmetric Dirichlet Laplacian?
%
% Everything measured so far used ONE plant: scalar, isotropic, self-adjoint,
% unit square, Dirichlet on all four edges, P = I.  That plant is maximally
% favourable to edge generators, because the extremal mode sin(pi s1)sin(pi s2)
% vanishes on every face.  Each plant below removes ONE assumption, so a
% failure can be attributed.
%
% The claim under test is a WITHIN-PLANT comparison -- does {s1,1-s1,s2,1-s2}
% still beat the stock product and no-psatz? -- so exact lam* is not required;
% it is reported where analytic to make "reach" comparable across plants.
%
% Also recorded per plant: whether Qop.R22{1,1} (the POINTWISE MULTIPLIER cell)
% is nonzero.  It is identically zero for the baseline plant, and a
% variable-coefficient plant should make it nonzero -- a structurally
% different regime for every degree argument.
cuadmm_path;

P = {};
%        name                    a1    a2    b1  bc1   bc2   dom            lamstar
P{end+1} = {'P1 baseline'         ,1   ,1   ,0 ,'DD','DD',[0 1;0 1]  , 2*pi^2};
P{end+1} = {'P2 aniso a2=0.25'    ,1   ,0.25,0 ,'DD','DD',[0 1;0 1]  , pi^2*(1+0.25)};
P{end+1} = {'P3 elongated 1x2'    ,1   ,1   ,0 ,'DD','DD',[0 1;0 2]  , pi^2*(1+1/4)};
P{end+1} = {'P4 offset [.3,1.3]^2',1   ,1   ,0 ,'DD','DD',[0.3 1.3;0.3 1.3], 2*pi^2};
P{end+1} = {'P5 Neumann in s2'    ,1   ,1   ,0 ,'DD','NN',[0 1;0 1]  , pi^2};
P{end+1} = {'P6 advection b1=3'   ,1   ,1   ,3 ,'DD','DD',[0 1;0 1]  , 2*pi^2+9/4};

FR = [0.25 0.50 0.90 0.99];
SETS = {'none',[] ; 'product',1 ; 'linear4',[3 4 5 6]};

fprintf('GP plant|mult_nonzero|set|m|nblk|frac|lam|cert|rel_b|feasratio|numerr|t\n');
for ip = 1:numel(P)
    pl = P{ip};  lamstar = pl{8};
    % structural probe: is the pointwise multiplier cell nonzero?
    try
        [~,Qop0] = buildQ(pl, 0.1*lamstar);
        mc = Qop0.R22{1,1};
        if isempty(mc), mz = 0;
        elseif isa(mc,'double'), mz = double(any(mc(:)~=0));
        else, mz = double(~isempty(mc.coefficient) && any(any(mc.coefficient~=0)));
        end
    catch ME
        fprintf('GP %s|BUILD_FAILED|%s\n', pl{1}, strrep(ME.message,newline,' '));
        continue
    end
    for is = 1:size(SETS,1)
        nm = SETS{is,1}; codes = SETS{is,2}; best = 0; mm=-1; nb=-1;
        for f = FR
            lam = f*lamstar;
            try
                t0=tic; [prog,~] = buildQ(pl,lam,codes);
                so.solver='mosek'; evalc('sol = lpisolve(prog,so);'); tw=toc(t0);
                S = cuadmm_private('sdpshape',sol); if mm<0, mm=S.m; nb=S.nblk; end
                I = sol.solinfo.info; xv = sol.solinfo.RRx(:);
                Atf=[];bf=[];
                for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
                rb = norm(full(Atf'*xv-bf))/norm(full(bf));
                triv = (abs(rb-1)<=1e-6)||norm(xv)<=1e-12;
                cert = (I.numerr==0)&&(I.feasratio>0.5)&&(rb<1e-6)&&~triv;
                fprintf('GP %s|%d|%s|%d|%d|%.2f|%.4f|%d|%.3e|%+.4f|%d|%.1f\n', ...
                    pl{1},mz,nm,S.m,S.nblk,f,lam,cert,rb,I.feasratio,I.numerr,tw);
                if cert, best=f; else, break; end
            catch ME
                fprintf('GP %s|%d|%s|ERR|%.2f|%s\n',pl{1},mz,nm,f,strrep(ME.message,newline,' '));
                break
            end
        end
        fprintf('GPREACH %s|mult=%d|%s|reach=%.2f|m=%d|nblk=%d\n',pl{1},mz,nm,best,mm,nb);
    end
end
fprintf('GPDONE\n');


function [prog,Qop] = buildQ(pl,lam,codes)
if nargin<3, codes = []; end
[~,a1,a2,b1,bc1,bc2,dom,~] = deal(pl{:});
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],dom);
rhs = a1*diff(x,s1,2) + a2*diff(x,s2,2) + lam*x;
if b1~=0, rhs = rhs + b1*diff(x,s1,1); end
sys = diff(x,t,1) == rhs;
sys = [sys; bcpair(x,s1,dom(1,:),bc1); bcpair(x,s2,dom(2,:),bc2)];
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
for p = codes
    o = eq_opts; o.psatz = p;
    [prog,Qe] = poslpivar_2d(prog,Qop.dim,eqd,o);
    Qeop = Qeop + Qe;
end
prog = lpi_eq_2d(prog,Qop+Qeop,'symmetric');
end

function bc = bcpair(x,v,iv,kind)
switch kind
    case 'DD', bc = [subs(x,v,iv(1))==0;             subs(x,v,iv(2))==0];
    case 'NN', bc = [subs(diff(x,v,1),v,iv(1))==0;   subs(diff(x,v,1),v,iv(2))==0];
    case 'DN', bc = [subs(x,v,iv(1))==0;             subs(diff(x,v,1),v,iv(2))==0];
end
end
