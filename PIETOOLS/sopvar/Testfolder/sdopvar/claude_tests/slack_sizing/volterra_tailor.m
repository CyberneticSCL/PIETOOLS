function T = volterra_tailor(dws)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = VOLTERRA_TAILOR(DWS) the norm bound of the Volterra operator,
% ||T|| = 2/pi (VOLTERRA_NORM_SOP), from  min gam  s.t.  gam - T'*T >= 0,
% with the inequality posed by 'lpi_ineq_sop' at the reader's degrees and
% the weight raised by dw in DWS (default 0:3), against the stock-style
% slack of VOLTERRA_NORM_SOP (poscopvar at degree 1, plain + product at the
% SAME degree) and the legacy DEMO2 LPI (lpi_ineq with psatz 1). Accuracy
% is the bound minus 2/pi; complexity the SDP shape.
% OUTPUT: struct array T (mode, dw, bound, excess, ndv, m, Ksmax, D, w, t).
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1 || isempty(dws),    dws = 0:3;  end
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
exact = 2/pi;
a = 0;  b = 1;
opvar Top;
Top.R.R1 = 1;   Top.I = [a,b];
Top.var1 = polynomial({'s'});   Top.var2 = polynomial({'s_dum'});
sopts = struct('solver','mosek','simplify',false);
T = struct('mode',{},'dw',{},'bound',{},'excess',{},'ndv',{},'m',{},'Ksmax',{},'D',{},'w',{},'terms',{},'t',{});

% % % Legacy DEMO2
t0 = tic;
prob = lpiprogram(Top.vars,Top.I);
[prob,gam] = lpidecvar(prob,'gam');
prob = lpi_ineq(prob,gam-Top'*Top,struct('psatz',1));
prob = lpisetobj(prob,gam);
prob = lpisolve(prob,sopts);
S = cx_shape(prob);
T(end+1) = row('legacy lpi_ineq psatz 1',NaN,sqrt(double(lpigetsol(prob,gam))),exact,S,NaN,NaN,'',toc(t0));

% % % Container, stock-style slack (volterra_norm_sop): degree 1, plain + product at the same degree
Tc = opvar2copvar(Top);
[sp,dm] = copvar_space_list(Tc,'out');
dom = struct('vars',{Tc.vars},'dom',Tc.dom);
for ds = {'stock-style deg 1', 'stock-style deg 2'}
    t0 = tic;
    d1 = 1 + strcmp(ds{1},'stock-style deg 2');
    prog = lpiprogram_sop(Tc);
    [prog,gam] = lpidecvar(prog,'gam');
    [prog,W1] = poscopvar(prog,dm,sp,dom,d1);
    [prog,W2] = poscopvar(prog,dm,sp,dom,d1,struct('psatz',1));
    prog = lpi_eq_sop(prog,gam-Tc'*Tc-(W1+W2),'symmetric');
    prog = lpisetobj(prog,gam);
    prog = lpisolve(prog,sopts);
    S = cx_shape(prog);
    T(end+1) = row(ds{1},NaN,sqrt(double(lpigetsol_sop(prog,gam))),exact,S,d1,d1,'[0 1]@[0 0]',toc(t0)); %#ok<AGROW>
end

% % % Container, lpi_ineq_sop at the reader's degrees, weight raised by dw
for dw = reshape(dws,1,[])
    for ps = {'product','faces'}
        t0 = tic;
        prog = lpiprogram_sop(Tc);
        [prog,gam] = lpidecvar(prog,'gam');
        [prog,~,inf] = lpi_ineq_sop(prog,gam-Tc'*Tc,struct('dw',dw,'psatz',ps{1}));
        prog = lpisetobj(prog,gam);
        prog = lpisolve(prog,sopts);
        S = cx_shape(prog);
        lbl = sprintf('lpi_ineq_sop %s',ps{1});
        T(end+1) = row(lbl,dw,sqrt(double(lpigetsol_sop(prog,gam))),exact,S,inf.degrees.D,inf.degrees.w,...
                       sprintf('%s@%s',mat2str(inf.terms.codes),mat2str(inf.terms.offsets)),toc(t0)); %#ok<AGROW>
    end
end
fprintf('\nVolterra norm, exact 2/pi = %.6f\n',exact);
for i = 1:numel(T)
    fprintf('  %-26s dw %s  bound %.6f  excess %.2e  | ndv %5d m %4d Ks<=%3d | D %s w %s %s | %.2f s\n',...
            T(i).mode,nstr(T(i).dw),T(i).bound,T(i).excess,T(i).ndv,T(i).m,T(i).Ksmax,nstr(T(i).D),nstr(T(i).w),T(i).terms,T(i).t);
end
end


function r = row(mode,dw,bound,exact,S,D,w,terms,t)
r = struct('mode',mode,'dw',dw,'bound',bound,'excess',bound-exact,'ndv',S.ndv,'m',S.m,...
           'Ksmax',max([S.Ks,0]),'D',D,'w',w,'terms',terms,'t',t);
end


function s = nstr(x)
if isnan(x),    s = '-';    else,   s = mat2str(x);    end
end
