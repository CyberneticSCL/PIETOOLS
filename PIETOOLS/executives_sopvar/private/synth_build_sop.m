function [prog,Gop,gam,info] = synth_build_sop(PIE,st,form,label)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,G,GAM,INFO] = SYNTH_BUILD_SOP(PIE,ST,FORM,LABEL) the H-infinity
% synthesis executives on the container path, the LPIs of the stock files as
% the cx_exec transcriptions wrote them (10/06/2026), one per FORM:
%
%   'ctrl'  PIETOOLS_Hinf_control (1-D):   P >= 0 (+ eppos), Z free,
%           CZ = CP + Dzu Z, AZ = AP + Bu Z,
%           K = [-gam I_z, Dzw, CZ T'; Dzw', -gam I_w, B'; T CZ', B, AZ T' + T AZ']
%           <= 0, the state feedback K = Z P^{-1} (getController);
%   'est'   PIETOOLS_Hinf_estimator, _2D:  P >= 0 (+ eppos), Z free,
%           PB = PB1 + Z Dyw, PA = PA + Z Cy,
%           K = [-gam I_w, -Dzw', -PB' T; -Dzw, -gam I_z, Cz; -T' PB, Cz', PA' T + T' PA]
%           <= 0, with the Tw terms of the stock when Tw ~= 0; the observer
%           gain L = P^{-1} Z (getObserver, getObserver_2D).
%
% The storage is the stock one (EXEC_TOOLS_SOP) plus eppos I; Z is an
% lpivar_cdopvar with ST.ddZ (1-D) or the largest entry of settings_2d.Zop_deg
% as its degree cap (2-D). The slack is LPI_INEQ_SOP with 'like' the
% storage and the weight raised by ST.sop.dw_P, or the stock slack when
% ST.sop.slack is 'stock'; the program is rebuilt while the solver returns
% UNKNOWN or INFEASIBLE (LPI_SOLVE_LOOP_SOP). With ST.sop.gam_fixed set the
% LPI is the feasibility test at that gamma. INFO carries P and Z as
% operators, the dual kernel, the slack degrees and terms, the shape and the
% history.
%
% Initial coding MMP, 10/08/2026
% MMP, 10/09/2026: In 1-D the gain is reconstructed on the container path,
%                  getController_sop / getObserver_sop (COPVAR/INV: Gohberg-
%                  Krein on the L2 block, Schur complement for R^k), from the
%                  solved containers; the opvar for closedLoopPIE / PIESIM is
%                  its conversion. The stock getController truncates K at
%                  1e-4 and its inverse was measured at residual 1e-3..1e-2
%                  on simple non-separable operators (opvar/
%                  inverse_dependency_map_2026_10_09.md). INFO gains G (the
%                  container gain) and inv (the inverse diagnostics). 2-D
%                  keeps getObserver_2D.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4,    label = ['Hinf ' form];     end
st = exec_sop_settings(st);
if ~isa(PIE,'pie_struct'),  error('synth_build_sop:pie','The PIE should be a ''pie_struct''.');  end
PIE = initialize(PIE);
tl = exec_tools_sop(PIE,st);    nv = tl.nv;
Tw0 = (PIE.Tw==0);
switch form
    case 'ctrl'
        if nv==2,   error('synth_build_sop:dim2','Optimal control of 2D PIEs is currently not supported.');   end
        if ~Tw0 || ~(PIE.Tu==0)
            error('synth_build_sop:Tw','Hinf-dual LPI cannot currently be solved for systems with disturbances or inputs at the boundary.');
        end
    case 'est'
        if nv==2 && ~Tw0 && ~(PIE.B1==0)
            error('synth_build_sop:TwB','The PIE takes both the input w and its derivative dw/dt; LPI based estimator synthesis is currently not supported.');
        end
    otherwise
        error('synth_build_sop:form','FORM should be ''ctrl'' or ''est''.');
end
gfix = [];  if isfield(st.sop,'gam_fixed') && ~isempty(st.sop.gam_fixed) && st.sop.gam_fixed~=0,    gfix = st.sop.gam_fixed;    end
eppos = tl.eppos;
Tm = op2copvar_sop(PIE.T);  Am = op2copvar_sop(PIE.A);
Bw = op2copvar_sop(PIE.B1); Cz = op2copvar_sop(PIE.C1);  Dzw = op2copvar_sop(PIE.D11);
Iw = tl.ident(PIE.B1,2);    Iz = tl.ident(PIE.C1,1);
Bu = [];    Dzu = [];   Cy = [];    Dyw = [];   Twm = [];
if strcmp(form,'ctrl'),     Bu = op2copvar_sop(PIE.B2);     Dzu = op2copvar_sop(PIE.D12);   end
if strcmp(form,'est'),      Cy = op2copvar_sop(PIE.C2);     Dyw = op2copvar_sop(PIE.D21);
    if ~Tw0,    Twm = op2copvar_sop(PIE.Tw);    end
end
[prog,aux] = lpi_solve_loop_sop(@assemble,st,label);
cl = aux.classify;
certified = any(strcmp(cl.status,{'optimal','inaccurate'}));
if isempty(gfix),   gam = cl.obj;
elseif certified,   gam = gfix;
else,               gam = Inf;
end
Gop = [];   Pop = [];   Zop = [];   Psol = [];   Zsol = [];
Gc = [];    iinv = [];                                                      % MMP, 10/09/2026
try
    Psol = getsol_lpivar_sop(prog,aux.Pm);  Zsol = getsol_lpivar_sop(prog,aux.Zm);
    Pop = tl.toop(Psol);    Zop = tl.toop(Zsol);
%   if strcmp(form,'ctrl'),     Gop = getController(Pop,Zop);               % MMP, 10/09/2026 (was)
%   elseif nv==1,               Gop = getObserver(Pop,Zop);                 % MMP, 10/09/2026 (was)
%   else,                       Gop = getObserver_2D(Pop,Zop);              % MMP, 10/09/2026 (was)
%   end                                                                     % MMP, 10/09/2026 (was)
    if nv==1 && strcmp(form,'ctrl'),    [Gc,Gop,iinv] = getController_sop(Psol,Zsol);  % MMP, 10/09/2026
    elseif nv==1,                       [Gc,Gop,iinv] = getObserver_sop(Psol,Zsol);    % MMP, 10/09/2026
    else,                               Gop = getObserver_2D(Pop,Zop);      % MMP, 10/09/2026
    end                                                                     % MMP, 10/09/2026
catch ME
    if st.sop.verbose,  fprintf('  [%s] gain extraction: %s\n',label,ME.message);  end
end
info = exec_info_sop(prog,aux,st,label);
info.form = form;   info.gam = gam;  info.gam_fixed = gfix;   info.nv = nv;
info.P = Pop;   info.Z = Zop;   info.Psol = Psol;   info.Zsol = Zsol;   info.K = aux.Km;
info.G = Gc;    info.inv = iinv;                                            % MMP, 10/09/2026
if st.sop.verbose
    fprintf('  [%s] %s: gamma %.8g (certified: %d); gain class %s\n',label,cl.status,gam,info.certified,class(Gop));
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
        [prog,Pm] = tl.storage(prog,Tm,'out');
        if ~all(eppos==0),  Pm = Pm + tl.eye(1,eppos);  end
        switch form
            case 'ctrl'
                [spu,dmu] = copvar_space_list(Bu,'in');     % Z: state -> R^nu
                [spx,dmx] = copvar_space_list(Bu,'out');
                [prog,Zm] = lpivar_cdopvar(prog,struct('out',dmu,'in',dmx),struct('out',{spu},'in',{spx}),tl.dom,tl.zdeg);
                Dzu_r = on_registry_sop(Dzu,Zm);
                CZ = Cz*Pm + Dzu_r*Zm;      AZ = Am*Pm + Bu*Zm;
                Km = [-(gamv*Iz),   Dzw,         CZ*Tm';
                       Dzw',        -(gamv*Iw),  Bw';
                       Tm*CZ',      Bw,          AZ*Tm' + Tm*AZ'];
            case 'est'
                [spx,dmx] = copvar_space_list(Cy,'in');     % Z: R^ny -> state
                [spy,dmy] = copvar_space_list(Cy,'out');
                [prog,Zm] = lpivar_cdopvar(prog,struct('out',dmx,'in',dmy),struct('out',{spx},'in',{spy}),tl.dom,tl.zdeg);
                Dyw_r = on_registry_sop(Dyw,Zm);
                PB = Pm*Bw + Zm*Dyw_r;      PA = Pm*Am + Zm*Cy;
                K11 = -(gamv*Iw);
                K13 = -PB'*Tm;              K31 = -Tm'*PB;
                if ~Tw0
                    X   = Twm'*PB;
                    K11 = on_registry_sop(K11,X) + X + PB'*Twm;
                    K13 = K13 - Twm'*PA;    K31 = K31 - PA'*Twm;
                end
                Km = [K11,      -Dzw',       K13;
                      -Dzw,     -(gamv*Iz),  Cz;
                      K31,      Cz',         PA'*Tm + Tm'*PA];
        end
        if strcmp(st.sop.slack,'stock')
            [prog,Nm] = tl.slack(prog,Km);
            [prog,tag] = lpi_eq_sop(prog,Nm+Km,'symmetric');
            li = struct('tag',tag,'degrees',[],'terms',[]);
        else
            o = struct('like',Pm,'dw',st.sop.dw_P+dw,'margin',st.sop.margin);
            [prog,Nm,li] = lpi_ineq_sop(prog,-Km,o);
        end
        obj = gamv;     if ~isempty(gfix),   obj = [];  end
        aux = struct('obj',obj,'Pm',Pm,'Zm',Zm,'Km',Km,'Nm',Nm,'li',li,'form',form);
    end
end
