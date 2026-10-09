function [prog,out,gam,info] = h2_build_sop(PIE,st,form,label)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,OUT,GAM,INFO] = H2_BUILD_SOP(PIE,ST,FORM,LABEL) the H2 executives on
% the container path, the LPIs of the stock files as the cx_exec
% transcriptions wrote them (10/06/2026), one per FORM:
%
%   'c'    PIETOOLS_H2_norm_c, _2D_c_non_coercive:
%          R >= 0, TQ = R, W >= 0 on R^nz,
%          Dneg = [-gam I_w, B'; B, Q'A' + AQ] <= 0,
%          Dpos = [W, CQ; Q'C', R] >= 0,  gam >= trace(W);
%   'o'    PIETOOLS_H2_norm_o, _2D_o_non_coercive:
%          R >= 0, T'Q = R, W >= 0 on R^nw,
%          Dneg = [-gam I_z, C; C', Q'A + A'Q] <= 0,
%          Dpos = [W, B'Q; Q'B, R] >= 0,  gam >= trace(W);
%   'cco'  PIETOOLS_H2_norm_c_coercive, _2D_c:
%          W >= 0 (2-D: + eppos on the R^n and L2[x,y] spaces, as stock),
%          (AW)T' + T(WA') + BB' <= 0, gam >= trace(C W C'), the norm bound
%          is sqrt(gam);
%   'oco'  PIETOOLS_H2_norm_o_coercive, _2D_o:
%          W >= 0 (2-D: + eppos), (A'W)T + T'(WA) + C'C <= 0,
%          gam >= trace(B' W B), the norm bound is sqrt(gam);
%   'ctrl' PIETOOLS_H2_control (1-D):   P >= 0 (+ eppos), Z, W free,
%          Dneg = [-gam I_w, B'; B, (TP)A' + (AP)T' + (B2 Z)T' + (T Z')B2'] <= 0,
%          Dpos = [W, CP + D12 Z; (.)', P] >= 0, gam >= trace(W), K = Z P^{-1};
%   'est'  PIETOOLS_H2_estimator (1-D): P >= 0 (+ eppos), Z, W free,
%          Dneg = [-gam I_z, C; C', (T'P)A + (A'P)T + (T'Z)C2 + (C2'Z')T] <= 0,
%          Dpos = [W, -(B'P + D21' Z'); (.)', P] >= 0, gam >= trace(W), L = P^{-1}Z.
%
% Each negativity / positivity constraint is LPI_INEQ_SOP with 'like' the
% storage (weight raised by ST.sop.dw_Q for the Q forms 'c', 'o', by dw_P
% otherwise), or the stock slack when ST.sop.slack is 'stock'; the trace
% inequality is the stock lpi_ineq on the 'dpvar' TRACE_RN_SOP returns. GAM
% is the value the stock executive returns (gam, or sqrt(gam) in the
% coercive forms). OUT is W (norms), K (ctrl) or L (est), as an opvar
% (opvar2d in 2-D). INFO carries, in 1-D with ST.sop.witness set, the
% numerical H2 norm of the open loop (PIE_WITNESS_SOP 'h2') for the four
% norm executives. The H2 branches of the stock are the least tested (the
% maintainer dropped H2 from the 09/26 regime); the status is reported, not
% presumed.
%
% Initial coding MMP, 10/08/2026
% MMP, 10/09/2026: The 'ctrl' / 'est' gains are reconstructed on the
%                  container path, getController_sop / getObserver_sop
%                  (COPVAR/INV), from the solved containers, as in
%                  synth_build_sop; OUT stays the opvar for closedLoopPIE /
%                  PIESIM, INFO gains G (container) and inv (diagnostics).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4,    label = ['H2 ' form];   end
st = exec_sop_settings(st);
if ~isa(PIE,'pie_struct'),  error('h2_build_sop:pie','The PIE should be a ''pie_struct''.');    end
PIE = initialize(PIE);
tl = exec_tools_sop(PIE,st);    nv = tl.nv;
if ~(PIE.Tw==0),    error('h2_build_sop:Tw','H2 norm LPI cannot be solved with disturbances at the boundary.');  end
switch form
    case 'ctrl'
        if nv==2,           error('h2_build_sop:dim2','H2 control of 2-D PIEs is not supported (as stock).');    end
        if ~(PIE.Tu==0),    error('h2_build_sop:Tu','Inputs at the boundary are not supported.');   end
        if ~(PIE.D11==0),   error('h2_build_sop:D11','Feedthrough D11 is not supported.');          end
    case 'est'
        if nv==2,           error('h2_build_sop:dim2','H2 estimation of 2-D PIEs is not supported (as stock).'); end
        if ~(PIE.D11==0),   error('h2_build_sop:D11','Feedthrough D11 is not supported.');          end
    case {'c','o','cco','oco'}
    otherwise
        error('h2_build_sop:form','FORM should be ''c'', ''o'', ''cco'', ''oco'', ''ctrl'' or ''est''.');
end
qform = any(strcmp(form,{'c','o'}));
dwbase = st.sop.dw_P;   if qform,    dwbase = st.sop.dw_Q;   end
eppos = tl.eppos;
Tm = op2copvar_sop(PIE.T);  Am = op2copvar_sop(PIE.A);
Bw = op2copvar_sop(PIE.B1); Cz = op2copvar_sop(PIE.C1);
Iw = tl.ident(PIE.B1,2);    Iz = tl.ident(PIE.C1,1);
B2m = [];   D12m = [];  C2m = [];   D21m = [];
if strcmp(form,'ctrl'),     B2m = op2copvar_sop(PIE.B2);    D12m = op2copvar_sop(PIE.D12);  end
if strcmp(form,'est'),      C2m = op2copvar_sop(PIE.C2);    D21m = op2copvar_sop(PIE.D21);  end
[prog,aux] = lpi_solve_loop_sop(@assemble,st,label);
cl = aux.classify;
gam_lpi = cl.obj;
gam = gam_lpi;
if any(strcmp(form,{'cco','oco'})),     gam = sqrt(gam_lpi);    end
% outputs
out = [];   Psol = [];  Zsol = [];  Wsol = [];  Rsol = [];  Qsol = [];
Pop = [];   Zop = [];   Wop = [];   Rop = [];   Qop = [];
Gc = [];    iinv = [];                                                      % MMP, 10/09/2026
try
    switch form
        case {'c','o'}
            Wsol = getsol_lpivar_sop(prog,aux.Wm);  Rsol = getsol_lpivar_sop(prog,aux.Rm);  Qsol = getsol_lpivar_sop(prog,aux.Qm);
            Wop = tl.toop(Wsol);    Rop = tl.toop(Rsol);    Qop = tl.toop(Qsol);    out = Wop;
        case {'cco','oco'}
            Wsol = getsol_lpivar_sop(prog,aux.Wm);  Wop = tl.toop(Wsol);    out = Wop;
        case {'ctrl','est'}
            Psol = getsol_lpivar_sop(prog,aux.Pm);  Zsol = getsol_lpivar_sop(prog,aux.Zm);  Wsol = getsol_lpivar_sop(prog,aux.Wm);
            Pop = tl.toop(Psol);    Zop = tl.toop(Zsol);    Wop = tl.toop(Wsol);
%           if strcmp(form,'ctrl'),     out = getController(Pop,Zop);       % MMP, 10/09/2026 (was)
%           else,                       out = getObserver(Pop,Zop);         % MMP, 10/09/2026 (was)
%           end                                                             % MMP, 10/09/2026 (was)
            if strcmp(form,'ctrl'),     [Gc,out,iinv] = getController_sop(Psol,Zsol);  % MMP, 10/09/2026
            else,                       [Gc,out,iinv] = getObserver_sop(Psol,Zsol);    % MMP, 10/09/2026
            end                                                             % MMP, 10/09/2026
    end
catch ME
    if st.sop.verbose,  fprintf('  [%s] output extraction: %s\n',label,ME.message);  end
end
info = exec_info_sop(prog,aux,st,label);
info.form = form;   info.gam = gam;  info.gam_lpi = gam_lpi;    info.nv = nv;
info.P = Pop;   info.Z = Zop;   info.W = Wop;   info.R = Rop;   info.Q = Qop;
info.G = Gc;    info.inv = iinv;                                            % MMP, 10/09/2026
info.witness = [];  info.h2_num = NaN;  info.gap = NaN;
if st.sop.witness && nv==1 && any(strcmp(form,{'c','o','cco','oco'}))
    Wt = pie_witness_sop(PIE,'h2',st.sop);
    info.witness = Wt;
    if Wt.ok
        info.h2_num = Wt.h2;
        if ~isnan(gam) && Wt.h2>0,  info.gap = gam/Wt.h2 - 1;   end
    end
end
if st.sop.verbose
    fprintf('  [%s] %s: value %.8g (LPI objective %.8g); numerical H2 norm %.8g; gap %.2e\n', ...
            label,cl.status,gam,gam_lpi,info.h2_num,info.gap);
end

    function [prog,aux] = assemble(dw)
        prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
        if st.sop.keepdual,     prog.sopeq = {};    end     % ask lpi_eq_cdopvar to record its rows
        [prog,gamv] = lpidecvar(prog,'gam');
        if ~strcmp(form,'o'),   prog = lpi_ineq(prog,gamv);     end       % as the stock 1-D files (H2_norm_o declares none)
        prog = lpisetobj(prog,gamv);
        Rm = [];    Qm = [];    Pm = [];    Zm = [];    Wm = [];
        switch form
            case 'c'
                [sp,dm] = copvar_space_list(Tm,'out');
                [prog,Rm] = tl.storage(prog,Tm,'out');
                [prog,Qm] = lpivar_cdopvar(prog,dm,sp,tl.dom,get_lpivar_degs_sop(Rm));
                prog = lpi_eq_sop(prog,Tm*Qm - Rm);
                [prog,Wm] = wdecl(prog,Cz,'out');
                Dneg = [-(gamv*Iw),  Bw';
                        Bw,          Qm'*Am' + Am*Qm];
                Dpos = [Wm,         Cz*Qm;
                        Qm'*Cz',    Rm];
                like = Rm;  tr = trace_rn_sop(Wm);
            case 'o'
                [sp,dm] = copvar_space_list(Tm,'out');
                [prog,Rm] = tl.storage(prog,Tm,'out');
                [prog,Qm] = lpivar_cdopvar(prog,dm,sp,tl.dom,get_lpivar_degs_sop(Rm));
                prog = lpi_eq_sop(prog,Tm'*Qm - Rm);
                [prog,Wm] = wdecl(prog,Bw,'in');
                Dneg = [-(gamv*Iz),  Cz;
                        Cz',         Qm'*Am + Am'*Qm];
                Dpos = [Wm,         Bw'*Qm;
                        Qm'*Bw,     Rm];
                like = Rm;  tr = trace_rn_sop(Wm);
            case 'cco'
                [prog,Wm] = tl.storage(prog,Tm,'out');
                Wm = add_eppos2d(Wm);
                Dneg = (Am*Wm)*Tm' + Tm*(Wm*Am') + Bw*Bw';
                Dpos = [];  like = Wm;  tr = trace_rn_sop((Cz*Wm)*Cz');
            case 'oco'
                [prog,Wm] = tl.storage(prog,Tm,'out');
                Wm = add_eppos2d(Wm);
                Dneg = (Am'*Wm)*Tm + Tm'*(Wm*Am) + Cz'*Cz;
                Dpos = [];  like = Wm;  tr = trace_rn_sop((Bw'*Wm)*Bw);
            case 'ctrl'
                [spx,dmx] = copvar_space_list(Am,'out');
                [prog,Pm] = tl.storage(prog,Am,'out');
                Pm = Pm + tl.eye(1,eppos);
                [spu,dmu] = copvar_space_list(B2m,'in');
                [prog,Zm] = lpivar_cdopvar(prog,struct('out',dmu,'in',dmx),struct('out',{spu},'in',{spx}),tl.dom,tl.zdeg);
                [spz,dmz] = copvar_space_list(Cz,'out');
                [prog,Wm] = lpivar_cdopvar(prog,dmz,spz,tl.dom,tl.zdeg);
                Dneg = [-(gamv*Iw),  Bw';
                        Bw,          (Tm*Pm)*Am' + (Am*Pm)*Tm' + (B2m*Zm)*Tm' + (Tm*Zm')*B2m'];
                Dp12 = Cz*Pm + on_registry_sop(D12m,Zm)*Zm;
                Dpos = [Wm,     Dp12;
                        Dp12',  Pm];
                like = Pm;  tr = trace_rn_sop(Wm);
            case 'est'
                [spx,dmx] = copvar_space_list(Tm,'out');
                [prog,Pm] = tl.storage(prog,Tm,'out');
                Pm = Pm + tl.eye(1,eppos);
                [spy,dmy] = copvar_space_list(C2m,'out');
                [prog,Zm] = lpivar_cdopvar(prog,struct('out',dmx,'in',dmy),struct('out',{spx},'in',{spy}),tl.dom,tl.zdeg);
                [spw,dmw] = copvar_space_list(Bw,'in');
                [prog,Wm] = lpivar_cdopvar(prog,dmw,spw,tl.dom,tl.zdeg);
                Dneg = [-(gamv*Iz),  Cz;
                        Cz',         (Tm'*Pm)*Am + (Am'*Pm)*Tm + (Tm'*Zm)*C2m + (C2m'*Zm')*Tm];
                D12 = Bw'*Pm + on_registry_sop(D21m',Zm')*Zm';
                Dpos = [Wm,       -D12;
                        -D12',    Pm];
                like = Pm;  tr = trace_rn_sop(Wm);
        end
        if strcmp(st.sop.slack,'stock')
            [prog,Nneg] = tl.slack(prog,Dneg);
            [prog,tag] = lpi_eq_sop(prog,Nneg + Dneg,'symmetric');
            li = struct('tag',tag,'degrees',[],'terms',[]);
            if ~isempty(Dpos)
                [prog,Npos] = tl.slack(prog,Dpos);
                prog = lpi_eq_sop(prog,Npos - Dpos,'symmetric');
            end
        else
            o = struct('like',like,'dw',dwbase+dw,'margin',st.sop.margin);
            [prog,~,li] = lpi_ineq_sop(prog,-Dneg,o);
            if ~isempty(Dpos)
                [prog,~,~] = lpi_ineq_sop(prog,Dpos,o);
            end
        end
        prog = lpi_ineq(prog,gamv - tr);
        aux = struct('obj',gamv,'Rm',Rm,'Qm',Qm,'Pm',Pm,'Zm',Zm,'Wm',Wm,'li',li,'form',form);
    end

    function [prog,Wm] = wdecl(prog,X,side)
        % W >= 0 on the SIDE space of X: poslpivar's default degrees in 1-D
        % (the stock poslpivar(prog,dim)); in 2-D the stock poslpivar_2d
        % defaults have no container form, so an R^n space takes the identity
        % basis and a space with variables the LF storage degrees (the trace
        % runs over the R^n block, so the value is the same whenever that
        % block exists; a distributed w gives the solver floor in both)
        if nv==1
            [prog,Wm] = poslpivar_sop(prog,X,[],struct(),side,PIE.dom);
        elseif any(any(X.(['space_' side])))
            [prog,Wm] = tl.storage(prog,X,side);
        else
            [prog,Wm] = poslpivar_sop(prog,X,[],struct(),side);
        end
    end

    function Wm = add_eppos2d(Wm)
        % the stock 2-D coercive H2 executives add eppos on R^n and L2[x,y]
        % (blkdiag(eppos(1) I, 0, 0, eppos(4) I)); the 1-D ones add nothing
        if nv==2
            e = eppos;  e(2:3) = 0;
            if ~all(e==0),  Wm = Wm + tl.eye(1,e);  end
        end
    end
end
