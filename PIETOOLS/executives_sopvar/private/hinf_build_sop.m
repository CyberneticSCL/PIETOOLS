function [prog,Rout,gam,info] = hinf_build_sop(PIE,st,form,label)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,R,GAM,INFO] = HINF_BUILD_SOP(PIE,ST,FORM,LABEL) the H-infinity gain
% executives on the container path, 1-D and 2-D, one LPI per FORM:
%
%   'Q'   PIETOOLS_Hinf_gain, _2D_non_coercive:   R = T'Q >= 0 (Q an lpivar),
%         K = [-gam I_w, D', B'Q; D, -gam I_z, C; Q'B, C', A'Q + Q'A] <= 0;
%         with a boundary disturbance (Tw ~= 0) the coercive form, as stock;
%   'Qd'  PIETOOLS_Hinf_gain_dual (1-D only):     TQ = R,
%         K = [-gam I_z, D, CQ; D', -gam I_w, B'; Q'C', B, Q'A' + AQ] <= 0;
%   'P'   PIETOOLS_Hinf_gain_coercive, _2D:       P = R >= 0,
%         K = [-gam I_w, D', (PB)'T; D, -gam I_z, C; T'PB, C', (PA)'T + T'PA]
%         <= 0, with the Tw terms of the stock when Tw ~= 0;
%   'Pd'  PIETOOLS_Hinf_gain_dual_coercive, _dual_2D:  K33 = TPA' + APT'.
%
% The storage is the stock one (EXEC_TOOLS_SOP: the lf degrees of ST in
% 1-D, LF_deg / LF_opts of settings_2d in 2-D; the 2-D coercive forms add
% eppos I as the stock 2-D executives do, the 1-D ones add nothing). The
% slack is LPI_INEQ_SOP with 'like' the storage and the weight raised by
% ST.sop.dw_Q (Q forms) or dw_P (coercive forms), or the stock slack when
% ST.sop.slack is 'stock'; the program is rebuilt and re-solved by
% LPI_SOLVE_LOOP_SOP while the solver returns UNKNOWN or INFEASIBLE. With
% ST.sop.gam_fixed set, the LPI is the feasibility test at that gamma (the
% stock executives' 'gain' argument) and GAM is that value when certified,
% Inf otherwise. INFO carries the certified gain, in 1-D the numerical
% witness (PIE_WITNESS_SOP 'gain': the frequency response's maximum, a
% lower bound), their relative gap, the dual kernel of the negativity
% constraint and its alignment with the worst input, the slack degrees and
% terms, the SDP shape and the loop history.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4,    label = ['Hinf ' form];     end
st = exec_sop_settings(st);
if ~isa(PIE,'pie_struct'),  error('hinf_build_sop:pie','The PIE should be a ''pie_struct''.');  end
PIE = initialize(PIE);
tl = exec_tools_sop(PIE,st);    nv = tl.nv;
Tw0 = (PIE.Tw==0);
if nv==2 && ~Tw0 && ~(PIE.B1==0)
    error('hinf_build_sop:TwB','The PIE takes both the input w and its derivative dw/dt; not supported (as the stock 2-D executives).');
end
switch form
    case 'Q'
        if ~Tw0,    form = 'P';     label = [label ' -> P (Tw ~= 0)'];   end
    case {'Qd','Pd'}
        if ~Tw0
            error('hinf_build_sop:Tw','Hinf-dual LPI cannot be solved for systems with disturbances at the boundary.');
        end
        if strcmp(form,'Qd') && nv==2
            error('hinf_build_sop:Qd2D','The dual Q form has no stock 2-D counterpart; PIETOOLS_Hinf_gain_dual_2D_sop is the coercive dual.');
        end
    case 'P'
    otherwise
        error('hinf_build_sop:form','FORM should be ''Q'', ''Qd'', ''P'' or ''Pd''.');
end
coercive = any(strcmp(form,{'P','Pd'}));
dwbase = st.sop.dw_Q;   if coercive,    dwbase = st.sop.dw_P;   end
gfix = [];  if isfield(st.sop,'gam_fixed') && ~isempty(st.sop.gam_fixed) && st.sop.gam_fixed~=0,    gfix = st.sop.gam_fixed;    end
Tm = op2copvar_sop(PIE.T);  Am = op2copvar_sop(PIE.A);
Bw = op2copvar_sop(PIE.B1); Cz = op2copvar_sop(PIE.C1);  Dzw = op2copvar_sop(PIE.D11);
Iw = tl.ident(PIE.B1,2);    Iz = tl.ident(PIE.C1,1);
Twm = [];   if ~Tw0,    Twm = op2copvar_sop(PIE.Tw);    end
% strictness on the storage: the stock 2-D coercive executives add eppos I
% (the dual with [1e-4;0;0;1e-6] when the settings carry no eppos); the 1-D
% gain executives add nothing
e = [];
if coercive && nv==2
    e = tl.eppos;
    if strcmp(form,'Pd') && tl.eppos_default,   e = [1e-4;0;0;1e-6];    end
end
[prog,aux] = lpi_solve_loop_sop(@assemble,st,label);
cl = aux.classify;
certified = any(strcmp(cl.status,{'optimal','inaccurate'}));
if isempty(gfix),   gam = cl.obj;
elseif certified,   gam = gfix;
else,               gam = Inf;
end
Rsol = [];  Rout = [];
try
    Rsol = getsol_lpivar_sop(prog,aux.Rm);
    Rout = tl.toop(Rsol);
catch
    Rout = Rsol;
end
info = exec_info_sop(prog,aux,st,label);
info.form = form;   info.gam = gam;  info.gam_fixed = gfix;   info.nv = nv;
info.storage = Rsol;    info.K = aux.Km;    info.N = aux.Nm;
info.witness = [];  info.gain_lb = NaN;     info.omega = NaN;   info.gap = NaN;     info.dual_alignment = NaN;
if st.sop.witness && nv==1
    if any(strcmp(form,{'Q','P'})),     Wt = pie_witness_sop(PIE,'gain',st.sop,info.dual,aux.wblk);
    else,                               Wt = pie_witness_sop(PIE,'gain',st.sop);
    end
    info.witness = Wt;
    if Wt.ok
        info.gain_lb = Wt.gain;     info.omega = Wt.omega;  info.dual_alignment = Wt.dual_alignment;
        if ~isnan(gam) && Wt.gain>0,    info.gap = gam/Wt.gain - 1;     end
    end
end
if st.sop.verbose
    fprintf('  [%s] %s: gamma %.8g (certified upper bound: %d); numerical lower bound %.8g at omega %.4g; gap %.2e; dual alignment %.3f\n', ...
            label,cl.status,gam,info.certified,info.gain_lb,info.omega,info.gap,info.dual_alignment);
end

    function [prog,aux] = assemble(dw)
        prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
        if st.sop.keepdual,     prog.sopeq = {};    end     % ask lpi_eq_cdopvar to record its rows
        if isempty(gfix)
            [prog,gamv] = lpidecvar(prog,'gam');
            prog = lpi_ineq(prog,gamv);
            prog = lpisetobj(prog,gamv);
        else
            gamv = gfix;
        end
        [prog,Rm] = tl.storage(prog,Tm,'out');
        if ~isempty(e) && ~all(e==0),   Rm = Rm + tl.eye(1,e);     end
        switch form
            case 'Q'
                [sp,dm] = copvar_space_list(Tm,'out');
                [prog,Qm] = lpivar_cdopvar(prog,dm,sp,tl.dom,get_lpivar_degs_sop(Rm));
                prog = lpi_eq_sop(prog,Tm'*Qm-Rm);
                Km = [-(gamv*Iw),   Dzw',       Bw'*Qm;
                       Dzw,         -(gamv*Iz), Cz;
                       Qm'*Bw,      Cz',        Am'*Qm + Qm'*Am];
                wblk = 1;
            case 'Qd'
                [sp,dm] = copvar_space_list(Tm,'out');
                [prog,Qm] = lpivar_cdopvar(prog,dm,sp,tl.dom,get_lpivar_degs_sop(Rm));
                prog = lpi_eq_sop(prog,Tm*Qm-Rm);
                Km = [-(gamv*Iz),   Dzw,        Cz*Qm;
                       Dzw',        -(gamv*Iw), Bw';
                       Qm'*Cz',     Bw,         Qm'*Am' + Am*Qm];
                wblk = 2;
            case 'P'
                PB = Rm*Bw;     PA = Rm*Am;
                K11 = -(gamv*Iw);   K13 = PB'*Tm;   K31 = Tm'*PB;
                if ~Tw0
                    X = Twm'*PB;
                    K11 = on_registry_sop(K11,X) + X + PB'*Twm;
                    K13 = K13 + Twm'*PA;    K31 = K31 + PA'*Twm;
                end
                Km = [K11,      Dzw',       K13;
                      Dzw,      -(gamv*Iz), Cz;
                      K31,      Cz',        PA'*Tm + Tm'*PA];
                wblk = 1;
            case 'Pd'
                CP = Cz*Rm;     AP = Am*Rm;
                Km = [-(gamv*Iz),   Dzw,        CP*Tm';
                       Dzw',        -(gamv*Iw), Bw';
                       Tm*CP',      Bw,         Tm*AP' + AP*Tm'];
                wblk = 2;
        end
        if strcmp(st.sop.slack,'stock')
            [prog,Nm] = tl.slack(prog,Km);
            [prog,tag] = lpi_eq_sop(prog,Nm+Km,'symmetric');
            li = struct('tag',tag,'degrees',[],'terms',[]);
        else
            o = struct('like',Rm,'dw',dwbase+dw,'margin',st.sop.margin);
            [prog,Nm,li] = lpi_ineq_sop(prog,-Km,o);
        end
        obj = gamv;     if ~isempty(gfix),   obj = [];  end
        aux = struct('obj',obj,'Rm',Rm,'Km',Km,'Nm',Nm,'li',li,'wblk',wblk,'form',form);
    end
end
