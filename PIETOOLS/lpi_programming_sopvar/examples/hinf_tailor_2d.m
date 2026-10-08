function T = hinf_tailor_2d(modes,level,form)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = HINF_TAILOR_2D(MODES,LEVEL,FORM) the 2-D H-infinity gain LPI (primal
% KYP, gamma a decision variable) on the container path for the plant io2
% of cx_plant (2-D reaction-diffusion at r = 15, closed-form L2 gain
% 0.1711), with the slack N sized by the stock 2-D settings
% (cx_hinf_slack2d: eq_deg/eq_opts of settings_2d plus its psatz terms) or
% by 'lpi_ineq_sop' (get_lift_degs per direction: faces at w by default).
% FORM 'Q' (default): the non-coercive form of PIETOOLS_Hinf_gain_2D_non_
% coercive, R (cx_hinf_lf2d, no eppos), Q (cx_hinf_qdeg2d), T'Q - R = 0
% and K33 = A'Q + Q'A; 'P': the coercive form of PIETOOLS_Hinf_gain_2D,
% P = R + eppos I and K33 = A'PT + T'PA. The storage is the same in every
% mode. The stock 2-D settings are a reference point, not a truth (ad hoc,
% little tested).
% MODES (default {'stock2d','new','new dw0','new dD0','new product'}; also
% 'deg D w [faces|product|none] [capJ]': lift D and weight w in both
% directions on the L2 spaces, R^q spaces at weight 0, faces by default,
% J a scalar cap on the total degree over the four axes [int, mult]),
% LEVEL lpisettings name (default 'light').
% OUTPUT: struct array T (mode, gam, excess over the closed-form gain, the
%   SDP shape, the slack's decision count and degrees, times).
%
% MEASURED (10/08/2026, light storage, MOSEK). Form 'Q': no tailored slack
% admits a finite gain. (1,1), (2,1), (1,2) with four faces, (2,2) with the
% product pair (cap 6 and uncapped) and (2,2) with four faces (cap 6: nx
% 8.9e5; uncapped: nx 9.8e5, gamma 12.4) all end with status UNKNOWN,
% residuals ~1e-9 and the objective growing to 1e1-3e3; the stock light
% slack (R's degrees + 3 in every slot) is a 1.5M-variable program, stopped
% after 9 min. Form 'P': 'new like' (the reader on P: (2,1) per direction,
% four faces) solves cleanly in 5 s to gamma 0.2665 (+56%) at nx 1.1e5;
% (3,1) faces 0.1974 (+15%) at 3.4e5; (2,2) faces 0.171718 (+0.36%) at
% 5.6e5, m 5420, 57 s (cap 6: 4.9e5, 36 s); (3,3) product 0.171777 (+0.40%)
% at 1.43e6, 491 s; (2,2), (3,2), (2,3) product fail (UNKNOWN); the stock
% light slack (7.6e5 variables, 13978 rows) was stopped by the budget at
% 459 s, objective 0.34 and falling, and the baseline suite (09/24) ran the
% stock coercive executive on this plant to 1.272 with status UNKNOWN.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1 || isempty(modes),  modes = {'stock2d','new','new dw0','new dD0','new product'};   end
if nargin<2 || isempty(level),  level = 'light';    end
if nargin<3 || isempty(form),   form = 'Q';         end
if ~any(strcmp(form,{'Q','P'})),    error('hinf_tailor_2d:form','FORM should be ''Q'' or ''P''.');    end
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
% closed-form gain of io2 (cx_cases_hinf)
nu = 1;     rr = 15;    Mm = 10;    Nn = 10;
mu_mn   = nu*pi^2*((2*(1:Mm)'-1).^2 + (2*(1:Nn)-1).^2);
mn_fact = (2*(1:Mm)'-1).*(2*(1:Nn)-1);
gex = (8/pi^2)*sqrt(sum((1./((mu_mn-rr).*mn_fact)).^2,'all'));
PIE = cx_plant('io2');
PIE = initialize(PIE);
st = cx_settings(level,'mosek');
S = cx_hinf_set2d(st,[1e-4;1e-6;1e-6;1e-6]);
Tm = cx_hinf_op2d(PIE.T);   Am = cx_hinf_op2d(PIE.A);
Bw = cx_hinf_op2d(PIE.B1);  Cz = cx_hinf_op2d(PIE.C1);    Dzw = cx_hinf_op2d(PIE.D11);
Iw = cx_hinf_op2d(mat2opvar(eye(size(PIE.B1,2)),PIE.B1.dim(:,2),PIE.vars,PIE.dom));
Iz = cx_hinf_op2d(mat2opvar(eye(size(PIE.C1,1)),PIE.C1.dim(:,1),PIE.vars,PIE.dom));
dom = struct();     dom.vars = reshape(PIE.vars(:,1).varname,1,[]);    dom.dom = PIE.dom;
fprintf('hinf_tailor_2d: io2, closed-form gain %.6f, level %s, form %s\n',gex,level,form);
T = struct('form',{},'mode',{},'gam',{},'excess',{},'numerr',{},'pinf',{},'rel_b',{},'ndv',{},'m',{},...
           'Ksmax',{},'sumKs2',{},'nN',{},'D',{},'w',{},'terms',{},'t_base',{},'t_slack',{},'t_solve',{});
for im = 1:numel(modes)
    mode = modes{im};
    t0 = tic;
    prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
    [prog,gam] = lpidecvar(prog,'gam');
    prog = lpi_ineq(prog,gam);
    prog = lpisetobj(prog,gam);
    if strcmp(form,'Q')             % non-coercive: R = T'Q >= 0 (no eppos), K33 = A'Q + Q'A
        [prog,Rm] = cx_hinf_lf2d(prog,Tm,PIE,S,false);
        [sp,dm] = cx_space_list(Tm,'out');
        [prog,Qm] = lpivar_cdopvar(prog,dm,sp,dom,cx_hinf_qdeg2d(Rm));
        prog = lpi_eq_cdopvar(prog,Tm'*Qm - Rm);
        Km = [-(gam*Iw),   Dzw',        Bw'*Qm;
               Dzw,        -(gam*Iz),   Cz;
               Qm'*Bw,     Cz',         Am'*Qm + Qm'*Am];
    else                            % coercive: P = R + eppos I, K33 = A'PT + T'PA
        [prog,Rm] = cx_hinf_lf2d(prog,Tm,PIE,S,true);
        PB = Rm*Bw;     PA = Rm*Am;
        Km = [-(gam*Iw),   Dzw',        PB'*Tm;
               Dzw,        -(gam*Iz),   Cz;
               Tm'*PB,     Cz',         PA'*Tm + Tm'*PA];
    end
    t_base = toc(t0);
    t0 = tic;
    D = NaN;    w = NaN;    terms = '';    nN = NaN;
    n0 = numel(prog.decvartable);
    switch mode
        case 'stock2d'
            prog = cx_hinf_slack2d(prog,Km,PIE,S);
            nN = numel(prog.decvartable)-n0;
        otherwise
            o = struct();
            tok = strsplit(strtrim(mode));
            if strcmp(tok{1},'deg')
                % 'deg D w [faces|product|none]': D, w in both directions on
                % the L2 space, R^q spaces at weight 0
                Dd = str2double(tok{2});    wd = str2double(tok{3});
                ps = 'faces';   J = [];
                for t = 4:numel(tok)    % 'capJ': scalar cap J on the total degree over [int, mult]
                    if strncmp(tok{t},'cap',3),     J = str2double(tok{t}(4:end));
                    else,                           ps = tok{t};
                    end
                end
                nv = numel(Km.vars);    Mk = size(Km.space_out,1);   deg = cell(1,Mk);
                for k = 1:Mk            % as get_lift_degs lays the specification out
                    own = find(Km.space_out(k,:));
                    if isempty(own),    deg{k} = struct('int',zeros(1,nv));
                    else
                        mult = zeros(1,nv);     mult(own) = Dd;
                        deg{k} = struct('int',wd*ones(1,nv),'mult',mult);
                        if ~isempty(J),     deg{k}.joint = J;  end
                    end
                end
                o.deg = deg;
                switch ps
                    case 'faces',   o.psatz = [0 3 4 5 6];  o.psatz_offset = zeros(1,5);
                    case 'product', o.psatz = [0 1];        o.psatz_offset = [0 1];
                    case 'none',    o.psatz = 0;            o.psatz_offset = 0;
                end
            else
                switch mode
                    case 'new'
                    case 'new dw0',     o.dw = 0;
                    case 'new dD0',     o.dD = 0;
                    case 'new product', o.psatz = 'product';
                    case 'new like',    o.like = Rm;      % the reader on R, laid out over K's spaces
                    otherwise,  error('hinf_tailor_2d:mode','Unknown mode ''%s''.',mode)
                end
            end
            [prog,Nm,inf] = lpi_ineq_sop(prog,-Km,o);
            nN = numel(Nm.Zd);
            if isfield(inf.degrees,'D'),    D = inf.degrees.D;  w = inf.degrees.w;
            else
                for k = 1:numel(inf.deg)
                    if isfield(inf.deg{k},'mult'),  D = inf.deg{k}.mult;    w = inf.deg{k}.int;     end
                end
            end
            terms = sprintf('%s@%s',mat2str(inf.terms.codes),mat2str(inf.terms.offsets));
    end
    t_slack = toc(t0);
    t0 = tic;   prog = lpisolve(prog,st.sos_opts);  t_solve = toc(t0);
    info = prog.solinfo.info;
    [rb,~,~] = cx_resid(prog);
    g = double(lpigetsol_sop(prog,gam));
    Sh = cx_shape(prog);
    T(end+1) = struct('form',form,'mode',mode,'gam',g,'excess',g-gex,'numerr',info.numerr,'pinf',info.pinf,'rel_b',rb,...
                      'ndv',Sh.ndv,'m',Sh.m,'Ksmax',max([Sh.Ks,0]),'sumKs2',sum(Sh.Ks.^2),'nN',nN,...
                      'D',D,'w',w,'terms',terms,'t_base',t_base,'t_slack',t_slack,'t_solve',t_solve); %#ok<AGROW>
    fprintf('  %-12s gamma %.6f (excess %+.2e) numerr %d pinf %d rel_b %.1e | ndv %7d m %6d Ks<=%4d sumKs2 %8d nN %7d | D %s w %s %s | base %.1fs slack %.1fs solve %.1fs\n',...
            mode,g,g-gex,info.numerr,info.pinf,rb,Sh.ndv,Sh.m,max([Sh.Ks,0]),sum(Sh.Ks.^2),nN,mat2str(D),mat2str(w),terms,t_base,t_slack,t_solve);
end
end
