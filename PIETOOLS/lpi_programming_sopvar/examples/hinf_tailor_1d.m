function T = hinf_tailor_1d(plants,levels,modes,form)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = HINF_TAILOR_1D(PLANTS,LEVELS,MODES,FORM) measures the slack sizing of
% the H-infinity gain LPI (primal KYP, gamma a decision variable) on the
% container path: the stock settings ('poslpivar_settings_sop' 'slack':
% dd2/options2 + dd3/options3 of lpisettings) against 'lpi_ineq_sop' sizing
% N from -K by 'get_lift_degs' (defaults and variants). FORM selects the
% storage: 'Q' (default), the non-coercive form of PIETOOLS_Hinf_gain, R =
% T'Q >= 0 with Q an lpivar and K33 = A'Q + Q'A; 'P', the coercive form of
% PIETOOLS_Hinf_gain_coercive, V = <Tx,PTx> with P = R >= 0 and K33 = A'PT
% + T'PA; 'Qd', the dual non-coercive form of PIETOOLS_Hinf_gain_dual, TQ =
% R and K33 = Q'A' + AQ; 'Pd', the dual coercive form of
% PIETOOLS_Hinf_gain_dual_coercive, K33 = TPA' + APT'; the dual forms swap
% the spaces of z and w. R (and Q) are the same in every mode, so only the
% negativity constraint K + N = 0 changes.
% For every (plant, level, mode): gamma, the
% stock executive's gamma on the same settings, the SDP shape (decision
% variables, equality rows, Gram sizes), the slack's decision count and
% degrees, and the build and solve times.
%
% INPUT (defaults): PLANTS cell of cx_plant argument cells
%   ({{'io1'},{'io1',2}}); LEVELS lpisettings names ({'light','heavy'});
%   FORM 'Q', 'Qd', 'P' or 'Pd' (plants with Tw ~= 0 are refused except in 'Q');
%   MODES cellstr among 'stock', 'new' (lpi_ineq_sop defaults: dD = dw = 1,
%   plain + product at w-1), 'new dw0', 'new dw2', 'new dD0', 'new faces',
%   'new none', 'new like' (the reader on R, opts.like, laid out over the
%   spaces of K), and 'deg D w [product|faces|none] [capJ]' (explicit lift
%   D and weight w on the L2 space, R^q spaces at weight 0, the product
%   pair by default, optional joint cap J), e.g. 'deg 2 2', 'deg 2 1 cap3'.
% OUTPUT: struct array T, one row per (plant, level, mode).
% Needs cx_plant, cx_settings, cx_shape, cx_resid (claude_tests/cx_exec).
%
% MEASURED (10/08/2026, io1, closed-form gain 0.182479080, MOSEK; details
% in lpi_programming_sopvar/README.md under lpi_ineq_sop). 'Q': every
% slack with D >= 2, w >= 2 gives 0.1824816 (+2.5e-6, a floor of this form);
% stock light (1,2) 0.182626; the reader on K overshoots to (5,3), the
% reader on R ('new like') gives the minimum (2,2). 'Qd': as 'Q' but
% reaching the closed form to 1e-7. 'P': stock light 0.1928 (+5.7%,
% UNKNOWN); lift-1 slacks diverge; the reader on K, (4,3) or (4,2) at dw 0,
% reaches the floor 1e-6 at nx 1786 (stock heavy 2332); 'new like' (2,2)
% +1.1e-5 at 810, (2,3) at the floor. 'Pd': stock 8252; every tailored
% slack primal infeasible.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1 || isempty(plants),     plants = {{'io1'},{'io1',2}};    end
if nargin<2 || isempty(levels),     levels = {'light','heavy'};     end
if nargin<3 || isempty(modes)
    modes = {'stock','new','new dw0','new dw2','new dD0','new faces','new none'};
end
if nargin<4 || isempty(form),   form = 'Q';     end
if ~any(strcmp(form,{'Q','Qd','P','Pd'})),  error('hinf_tailor_1d:form','FORM should be ''Q'', ''Qd'', ''P'' or ''Pd''.');    end
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
T = struct('plant',{},'level',{},'form',{},'mode',{},'gam_stock',{},'gam',{},'relgap',{},'numerr',{},...
           'pinf',{},'rel_b',{},'psd_relmin',{},'ndv',{},'m',{},'Ksmax',{},'sumKs2',{},'nN',{},...
           'D',{},'w',{},'terms',{},'t_slack',{},'t_solve',{});
for ip = 1:numel(plants)
    PIE = cx_plant(plants{ip}{:});
    plbl = strjoin(cellfun(@num2str,plants{ip},'UniformOutput',false),',');
    if ~strcmp(form,'Q') && ~(PIE.Tw==0)
        error('hinf_tailor_1d:Tw','Forms ''P'' and ''Pd'' are written for Tw == 0 (no boundary disturbance).');
    end
    for il = 1:numel(levels)
        st = cx_settings(levels{il},'mosek');
        switch form
            case 'Q',   [~,~,gam_s] = PIETOOLS_Hinf_gain(PIE,st);
            case 'Qd',  [~,~,gam_s] = PIETOOLS_Hinf_gain_dual(PIE,st);
            case 'P',   [~,~,gam_s] = PIETOOLS_Hinf_gain_coercive(PIE,st);
            case 'Pd',  [~,~,gam_s] = PIETOOLS_Hinf_gain_dual_coercive(PIE,st);
        end
        fprintf('\n%s / %s / form %s: stock executive gamma %.6f; stock slack dd2 %s dd3 %s (override2 %d, options2 %s, options3 %s)\n',...
                plbl,levels{il},form,gam_s,degstr(st.dd2),degstr(st.dd3),st.override2,...
                jsonencode(st.options2),jsonencode(st.options3));
        Tm = opvar2copvar(PIE.T);   Am = opvar2copvar(PIE.A);
        Bw = opvar2copvar(PIE.Bw);  Cz = opvar2copvar(PIE.Cz);  Dzw = opvar2copvar(PIE.Dzw);
        [spw,dmw] = copvar_space_list(Bw,'in');     [spz,dmz] = copvar_space_list(Cz,'out');
        Iw = eye_copvar_sop(dmw,spw,PIE.dom);       Iz = eye_copvar_sop(dmz,spz,PIE.dom);
        for im = 1:numel(modes)
            mode = modes{im};
            prog = lpiprogram_sop(Tm);
            [prog,gam] = lpidecvar(prog,'gam');
            prog = lpi_ineq(prog,gam);
            prog = lpisetobj(prog,gam);
            [prog,Rm] = poslpivar_settings_sop(prog,Tm,st,'lf','out',PIE.dom);
            if strcmp(form,'Q')         % non-coercive: R = T'Q >= 0, K33 = A'Q + Q'A
                [sp,dm] = copvar_space_list(Tm,'out');
                [prog,Qm] = lpivar_cdopvar(prog,dm,sp,PIE.dom,get_lpivar_degs_sop(Rm));
                prog = lpi_eq_sop(prog,Tm'*Qm-Rm);
                Km = [-(gam*Iw),   Dzw',        Bw'*Qm;
                       Dzw,        -(gam*Iz),   Cz;
                       Qm'*Bw,     Cz',         Am'*Qm + Qm'*Am];
            elseif strcmp(form,'Qd')    % dual non-coercive (PIETOOLS_Hinf_gain_dual): TQ = R, K33 = Q'A' + AQ
                [sp,dm] = copvar_space_list(Tm,'out');
                [prog,Qm] = lpivar_cdopvar(prog,dm,sp,PIE.dom,get_lpivar_degs_sop(Rm));
                prog = lpi_eq_sop(prog,Tm*Qm-Rm);
                Km = [-(gam*Iz),   Dzw,         Cz*Qm;
                       Dzw',       -(gam*Iw),   Bw';
                       Qm'*Cz',    Bw,          Qm'*Am' + Am*Qm];
            elseif strcmp(form,'P')     % coercive: V = <Tx,PTx>, P = R >= 0, K33 = A'PT + T'PA
                PB = Rm*Bw;     PA = Rm*Am;
                Km = [-(gam*Iw),   Dzw',        PB'*Tm;
                       Dzw,        -(gam*Iz),   Cz;
                       Tm'*PB,     Cz',         PA'*Tm + Tm'*PA];
            else                        % dual coercive (PIETOOLS_Hinf_gain_dual_coercive): K33 = TPA' + APT'
                Km = [-(gam*Iz),   Dzw,         Cz*Rm*Tm';
                       Dzw',       -(gam*Iw),   Bw';
                       Tm*Rm*Cz',  Bw,          Tm*Rm*Am' + Am*Rm*Tm'];
            end
            t0 = tic;
            D = NaN;    w = NaN;    terms = '';
            switch mode
                case 'stock'
                    [prog,Nm] = poslpivar_settings_sop(prog,Km,st,'slack','out',PIE.dom);
                    prog = lpi_eq_sop(prog,Nm+Km,'symmetric');
                otherwise
                    if strcmp(mode,'new like'),     o = struct('like',Rm);  % the reader on R, laid out over K's spaces
                    else,   o = mode_opts(mode,copvar_space_list(Km,'out'));
                    end
                    [prog,Nm,inf] = lpi_ineq_sop(prog,-Km,o);
                    if isfield(inf.degrees,'D'),    D = inf.degrees.D;  w = inf.degrees.w;
                    else    % explicit degrees: read them back from the spec of the L2 space
                        for k = 1:numel(inf.deg)
                            if isfield(inf.deg{k},'mult'),  D = inf.deg{k}.mult;    w = inf.deg{k}.int;     end
                        end
                    end
                    terms = sprintf('%s@%s',mat2str(inf.terms.codes),mat2str(inf.terms.offsets));
            end
            t_slack = toc(t0);
            t0 = tic;   prog = lpisolve(prog,st.sos_opts);  t_solve = toc(t0);
            info = prog.solinfo.info;
            [rb,~,prel] = cx_resid(prog);
            try     g = double(lpigetsol_sop(prog,gam));
            catch,  g = NaN;        % sosgetsol refuses a program the solver found infeasible
            end
            S = cx_shape(prog);
            r = struct('plant',plbl,'level',levels{il},'form',form,'mode',mode,'gam_stock',gam_s,'gam',g,...
                       'relgap',g/gam_s-1,'numerr',info.numerr,'pinf',info.pinf,'rel_b',rb,'psd_relmin',prel,...
                       'ndv',S.ndv,'m',S.m,'Ksmax',max([S.Ks,0]),'sumKs2',sum(S.Ks.^2),'nN',numel(Nm.Zd),...
                       'D',D,'w',w,'terms',terms,'t_slack',t_slack,'t_solve',t_solve);
            T(end+1) = r; %#ok<AGROW>
            fprintf('  %-10s gamma %.6f (%+.2e vs stock exec) numerr %d pinf %d rel_b %.1e psd %.1e | ndv %6d m %5d Ks<=%3d sumKs2 %7d nN %6d | D %s w %s %s | slack %.2fs solve %.2fs\n',...
                    mode,g,r.relgap,info.numerr,info.pinf,rb,prel,S.ndv,S.m,r.Ksmax,r.sumKs2,r.nN,...
                    mat2str(D),mat2str(w),terms,t_slack,t_solve);
        end
    end
end
end


function o = mode_opts(mode,sp)
% 'new ...' variants of the reader, or 'deg D w [product|faces|none] [capJ]':
% explicit lift D and weight w on every L2 space (R^q spaces at weight 0),
% product pair by default, optional joint cap J.
o = struct();
tok = strsplit(strtrim(mode));
if strcmp(tok{1},'deg')
    D = str2double(tok{2});     w = str2double(tok{3});
    ps = 'product';     J = [];
    for t = 4:numel(tok)
        if strncmp(tok{t},'cap',3),     J = str2double(tok{t}(4:end));
        else,                           ps = tok{t};
        end
    end
    deg = cell(1,numel(sp));
    for k = 1:numel(sp)
        if isempty(sp{k}),  deg{k} = struct('int',0);
        else
            deg{k} = struct('int',w,'mult',D);
            if ~isempty(J),     deg{k}.joint = J;  end
        end
    end
    o.deg = deg;
    switch ps
        case 'product',     o.psatz = [0 1];    o.psatz_offset = [0 1];
        case 'faces',       o.psatz = [0 3 4];  o.psatz_offset = [0 0 0];
        case 'quad',        o.psatz = [0 5];    o.psatz_offset = [0 1];     % the single-direction quadratic code; = 'product' in 1-D
        case 'none',        o.psatz = 0;        o.psatz_offset = 0;
        otherwise,          error('hinf_tailor_1d:mode','Unknown psatz ''%s''.',ps)
    end
    return
end
switch mode
    case 'new'
    case 'new dw0',     o.dw = 0;
    case 'new dw2',     o.dw = 2;
    case 'new dD0',     o.dD = 0;
    case 'new dD0 dw0', o.dD = 0;   o.dw = 0;
    case 'new faces',   o.psatz = 'faces';
    case 'new none',    o.psatz = 'none';
    otherwise,          error('hinf_tailor_1d:mode','Unknown mode ''%s''.',mode)
end
end


function s = degstr(dd)
if iscell(dd),  s = ['{' strjoin(cellfun(@mat2str,dd,'UniformOutput',false),',') '}'];
else,           s = mat2str(dd);
end
end
